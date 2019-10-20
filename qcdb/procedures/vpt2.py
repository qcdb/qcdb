#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2019 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

"""Module with functions for Psi4/Cfour interface. Portions that require
calls to Boost Python psi4 module are here, otherwise in qcdb module.
Also calls to qcdb module are here and not elsewhere in driver.
Organizationally, this module isolates qcdb code from psi4 code.

"""
import os
import re
import copy
import glob
import uuid
import shelve
import shutil
import difflib
import datetime
import subprocess
import pprint
pp = pprint.PrettyPrinter(width=120)
from pathlib import Path

#from psi4.driver.p4util.exceptions import *

import qcelemental as qcel
from qcelemental.models import ResultInput
import qcengine as qcng

from ..driver import pe
from ..util import provenance_stamp
from ..keywords import Keywords, register_kwds
from ..driver.driver_util import kwargs_lower, get_package
from ..driver.driver_helpers import get_active_molecule
from ..driver.proc_table import procedures
from ..programs.cfour.harvester import backtransform, format_fjobarc, harvest_zmat, jajo2mol
from ..programs.cfour.jajo import getrec


def run_cfour_module(xmod):
    # Find environment by merging PSIPATH and PATH environment variables
    lenv = {
        'PATH': ':'.join([os.path.abspath(x) for x in os.environ.get('PSIPATH', '').split(':') if x != '']) + \
                ':' + os.environ.get('PATH') + \
                ':' + core.get_datadir() + '/basis' + \
                ':' + core.psi_top_srcdir() + '/share/basis',
        'CFOUR_NUM_CORES': os.environ.get('CFOUR_NUM_CORES'),
        'LD_LIBRARY_PATH': os.environ.get('LD_LIBRARY_PATH')
        }
    #   Filter out None values as subprocess will fault on them
    lenv = {k: v for k, v in lenv.items() if v is not None}

    # Call executable xcfour, directing cfour output to the psi4 output file
    try:
        retcode = subprocess.Popen([xmod], bufsize=0, stdout=subprocess.PIPE, env=lenv)
    except OSError as e:
        sys.stderr.write('Program %s not found in path or execution failed: %s\n' % (cfour_executable, e.strerror))
        #p4out.write('Program %s not found in path or execution failed: %s\n' % (cfour_executable, e.strerror))
        message = ('Program %s not found in path or execution failed: %s\n' % (cfour_executable, e.strerror))
        raise ValidationError(message)        

    c4out = ''
    while True:
        data = retcode.stdout.readline()
        if not data:
            break
        #if core.outfile_name() == 'stdout':
        #    sys.stdout.write(data)
        #else:
        #    p4out.write(data)
        #    p4out.flush()
        c4out += data
    #internal_p4c4_info['output'] = c4out
    return c4out


@register_kwds(pe.nu_options)
def vpt2(name, **kwargs):
    """Perform vibrational second-order perturbation computation through
    Cfour to get anharmonic frequencies. This version uses c4 for the disp
    and pt2 but gets gradients from p4.

    :type c4full: :ref:`boolean <op_py_boolean>`
    :param c4full: ``'on'`` || |dl| ``'off'`` |dr|

        Indicates whether when *name* indicates a Cfour method and *mode*
        indicates a sow/reap approach, sown files are direct ZMAT files
        and FJOBARC files are expected to reap, so that Cfour only, not
        Cfour-through-Psi4, is needed for distributed jobs.

    .. caution:: Some features are not yet implemented. Buy a developer a coffee.

       - Presently uses all gradients. Could mix in analytic 2nd-derivs.

       - Collect resutls.

       - Manage scratch / subdir better.

       - Allow CFOUR_BASIS

       - Consider forcing some tighter convcrit, c4 and p4

       - mixed ang/bohr signals

       - error by converting to ang in psi?

       - Expand CURRENT DIPOLE XYZ beyond SCF

       - Remember additional FJOBARC record TOTENER2 if EXCITE .ne. NONE

       - switch C --> S/R with recovery using shelf

    """
    from ..driver import endorsed_plugins
    kwargs = kwargs_lower(kwargs)

#    if 'options' in kwargs:
#        driver_helpers.set_options(kwargs.pop('options'))
#
#    # Bounce if name is function
#    if hasattr(name, '__call__'):
#        return name(energy, kwargs.pop('label', 'custom function'), ptype='energy', **kwargs)
#
#    # Allow specification of methods to arbitrary order
    lowername = name.lower()
    package = get_package(lowername, kwargs)
#    lowername, level = driver_helpers._parse_arbitrary_order(lowername)
#    if level:
#        kwargs['level'] = level

    # Make sure the molecule the user provided is the active one
    molecule = kwargs.pop('molecule', get_active_molecule())
    molecule.update_geometry()

#    if len(pe.nu_options.scroll) == 0:
#        #print('EMPTY OPT')
#        pe.load_nu_options()

# -----
    verbose = kwargs.pop('verbose', 0)
    scratch_messy = kwargs.pop('scratch_messy', True)  # TODO

    kwgs = {'accession': kwargs['accession'], 'verbose': verbose}

#    optstash = p4util.OptionsState(
#        ['BASIS'])

    # Option mode of operation- whether vpt2 run in one job or files farmed out
    if not('vpt2_mode' in kwargs):
        if ('mode' in kwargs):
            kwargs['vpt2_mode'] = kwargs['mode']
            del kwargs['mode']
        else:
            kwargs['vpt2_mode'] = 'continuous'

    # Switches for route through code- S/R or continuous & Psi4 or Cfour gradients
    isSowReap = True if kwargs['vpt2_mode'].lower() == 'sowreap' else False
#!BR#    isC4notP4 = bool(re.match('cfour', lowername)) or bool(re.match('c4-', lowername))
    isC4notP4 = False  # TODO until intf_psi4 hooked up to qcng
    isC4fully = True if ('c4full' in kwargs and yes.match(str(kwargs['c4full'])) and isC4notP4 and isSowReap) else False
    print('isSowReap=', isSowReap, 'isC4notP4=', isC4notP4, 'isC4fully=', isC4fully)

    cfourharness = qcng.get_program('qcdb-cfour')
    config = qcng.config.get_config(local_options={"ncores": 2})

    # Save submission directory and basis set
    current_directory = os.getcwd()
#    user_basis = core.get_global_option('BASIS')

    # Open data persistence shelf- vital for sowreap, checkpoint for continuouw
#    shelf = shelve.open(current_directory + '/' + os.path.splitext(core.outfile_name())[0] + '.shelf', writeback=True)
    shelf = shelve.open(current_directory + '/vpt2scratch.shelf', writeback=True)

    # Cfour keywords to request vpt2 analysis through findif gradients
    c000_opts = Keywords()
    pe.load_options(c000_opts)
    c000_opts.require('CFOUR', 'VIBRATION', 'FINDIF', **kwgs)
    c000_opts.require('CFOUR', 'FREQ_ALGORITHM', 'PARALLEL', **kwgs)
    c000_opts.require('CFOUR', 'ANH_ALGORITHM', 'PARALLEL', **kwgs)
    c000_opts.require('CFOUR', 'ANHARMONIC', 'VPT2', **kwgs)
    c000_opts.require('CFOUR', 'FD_PROJECT', 'OFF', **kwgs)


    # When a Psi4 method is requested for vpt2, a skeleton of
    #   computations in Cfour is still required to hang the gradients
    #   upon. The skeleton is as cheap as possible (integrals only
    #   & sto-3g) and set up here.
#!BR#    if isC4notP4:
#!BR#        skelname = lowername
#!BR#    else:
    if True:
        skelname = 'c4-scf'
#        core.set_global_option('BASIS', 'STO-3G')
    #    P4  'c4-scf'/'cfour'CALC_LEVEL      lowername  # temporary
    #    C4  lowername                       cfour{}  # temporary

    if 'status' not in shelf:
        shelf['status'] = 'initialized'
        shelf['linkage'] = os.getpid()
        shelf['zmat'] = {}  # Cfour-generated ZMAT files with finite difference geometries
        shelf['fjobarc'] = {}  # Cfour- or Psi4-generated ascii files with packaged gradient results
        shelf['results'] = {}  # models.Result
        shelf.sync()
    else:
        pass
        # how decide whether to use. keep precedent of intco.dat in mind

#    # Construct and move into directory job scratch / cfour scratch / harm
#    psioh = core.IOManager.shared_object()
#    psio = core.IO.shared_object()
#    os.chdir(psioh.get_default_path())  # psi_scratch
#    cfour_tmpdir = kwargs['path'] if 'path' in kwargs else \
#        'psi.' + str(os.getpid()) + '.' + psio.get_default_namespace() + \
#        '.cfour.' + str(uuid.uuid4())[:8]
#    if not os.path.exists(cfour_tmpdir):
#        os.mkdir(cfour_tmpdir)
#    os.chdir(cfour_tmpdir)  # psi_scratch/cfour
#    if not os.path.exists('harm'):
#        os.mkdir('harm')
#    os.chdir('harm')  # psi_scratch/cfour/harm

#    psioh.set_specific_retention(32, True)  # temporary, to track p4 scratch
    #shelf['status'] = 'anharm_jobs_sown'  # temporary to force backtrack
    print('STAT', shelf['status'])  # temporary

    resi = ResultInput(
        **{
            'driver': 'energy',  # to prevent qcdb imposition of analytic hessian
            'extras': {
                'qcdb:options': copy.deepcopy(c000_opts), #pe.nu_options),
            },
            'model': {
                'method': 'c4-scf', #'hf',
                #'basis': '6-31g',
                'basis': 'sto-3g',
            },
            'molecule': molecule.to_schema(dtype=2),
            'provenance': provenance_stamp(__name__),
        })

    # Generate the ZMAT input file in scratch
    cfourrec = cfourharness.qcdb_build_input(resi, config)
    shelf['genbas'] = cfourrec['infiles']['GENBAS']
    shelf['zmat']['000-000'] = cfourrec['infiles']['ZMAT']
    shelf.sync()

#    with open('ZMAT', 'w') as handle:
#        cfour_infile = write_zmat(skelname, 1)
#        handle.write(cfour_infile)
#    print('\n====== Begin ZMAT input for CFOUR ======')
#    print(open('ZMAT', 'r').read())
#    print('======= End ZMAT input for CFOUR =======\n')

    # Check existing shelf consistent with generated ZMAT, store
#    if ('000-000' in shelf['zmat']) and (shelf['zmat']['000-000'] != cfour_infile):
#        diff = difflib.Differ().compare(shelf['zmat']['000-000'].splitlines(), cfour_infile.splitlines())
#        raise ValidationError("""Input file translated to Cfour ZMAT does not match ZMAT stored in shelf.\n\n""" +
#            '\n'.join(list(diff)))

    # Reset basis after Cfour skeleton seeded
#    core.set_global_option('BASIS', user_basis)

    if shelf['status'] == 'initialized':
        print('{:_^45}'.format('  VPT2 Setup: Harmonic  '))

        # Generate the displacements that will form the harmonic freq
        scrkwgs = {'scratch_directory': config.scratch_directory, 'scratch_messy': True, 'scratch_suffix': '_000'}
        success, dexe = qcng.util.execute(['xjoda'], cfourrec['infiles'], [], **scrkwgs)
        partial = dexe['stdout']

        scrkwgs.update({'scratch_name': Path(dexe['scratch_directory']).name, 'scratch_exist_ok': True})
        scrkwgs.update({'scratch_messy': scratch_messy})
        success, dexe = qcng.util.execute(['xsymcor'], {}, ['zmat*'], **scrkwgs)
        partial += dexe['stdout']
    
        print(partial)  # partial.out

        # Read the displacements that will form the harmonic freq
        zmats0N = ['000-' + item[-3:] for item in dexe['outfiles']['zmat*']]
        for zm_2 in zmats0N:
            _, zm2 = zm_2.split('-')
            shelf['zmat'][zm_2] = dexe['outfiles']['zmat*']['zmat' + zm2]
            shelf.sync()
            print(f'  CFOUR scratch file zmat{zm2} for {zm_2} has been read\n')
            #print('%s\n' % shelf['zmat'][zm_2])

        # S/R: Write distributed input files for harmonic freq
        if isSowReap:
            os.chdir(current_directory)
            inputSansMol = p4util.format_currentstate_for_input(gradient, lowername, allButMol=True, **kwargs)
            for zm12 in zmats0N:
                zm1, zm2 = zm12.split('-')

                ifile = vpt2_sow_files(zm12, shelf['linkage'], isC4notP4, isC4fully,
                    shelf['zmat'][zm12], inputSansMol, shelf['genbas'])

                with open('VPT2-' + zm12 + '.in', 'w') as handle:
                    handle.write(ifile)

            msg = vpt2_instructions('harmonic', current_directory, zmats0N)
            core.print_out(msg)
            print(msg)

        shelf['status'] = 'harm_jobs_sown'

        # S/R: Pause for distributed calculations
        if isSowReap:
            shelf.close()
            return 0.0

    if shelf['status'] == 'harm_jobs_sown':
        zmats0N = [item for item in sorted(shelf['zmat'].keys()) if (item[:3] == '000' and item[-3:] != '000')]

        # S/R: Check that distributed calcs all completed correctly
        if isSowReap:
            msg = vpt2_instructions('harmonic', current_directory, zmats0N)
            core.print_out(msg)
            isOk, msg = sown_jobs_status(current_directory, 'VPT2', zmats0N, reap_job_validate,
                shelf['linkage'], ['CURRENT ENERGY', 'CURRENT DIPOLE', 'CURRENT GRADIENT'])
            core.print_out(msg)
            print(msg)
            if not isOk:
                shelf.close()
                return 0.0

        # Collect all results from gradients forming the harmonic freq
        for zm12 in zmats0N:
            zm1, zm2 = zm12.split('-')
            if zm12 not in shelf['fjobarc']:
                print('{:_^45}'.format(f'  VPT2 Computation: {zm12}  '))

                fjobarc = vpt2_reaprun_files(zm12, shelf['linkage'], isSowReap, isC4notP4, isC4fully,
                    shelf['zmat'][zm12], #current_directory, psioh.get_default_path(), cfour_tmpdir,
                    lowername, kwargs, shelf['genbas'], config, package, scratch_messy=scratch_messy)
                shelf['fjobarc'][zm12] = fjobarc
                shelf.sync()
        shelf['status'] = 'harm_jobs_reaped'

    if shelf['status'] == 'harm_jobs_reaped':
        zmats0N = [item for item in sorted(shelf['zmat'].keys()) if (item[:3] == '000' and item[-3:] != '000')]

        print('{:_^45}'.format('  VPT2 Results: Harmonic  '))
        for k, v in shelf.items():
            print('   {:_^20}'.format(k))
            #pp.pprint(v)

        #scrkwgs = {'scratch_directory': config.scratch_directory, 'scratch_messy': True, 'scratch_suffix': '_000'}
        #scrkwgs.update({'scratch_name': Path(dexe['scratch_directory']).name, 'scratch_exist_ok': True})
        #scrkwgs.update({'scratch_messy': scratch_messy})

        # Process the gradients into harmonic freq
        scrkwgs = {'scratch_directory': config.scratch_directory, 'scratch_messy': True, 'scratch_suffix': '_000med'}
        success, dexe = qcng.util.execute(['xjoda'], {'ZMAT': shelf['zmat']['000-000'], 'GENBAS': shelf['genbas']}, [], **scrkwgs)
        harmout = dexe['stdout']

        scrkwgs.update({'scratch_name': Path(dexe['scratch_directory']).name, 'scratch_exist_ok': True})
        success, dexe = qcng.util.execute(['xsymcor'], {}, [], **scrkwgs)
        print('xsymcor', success)
        harmout += dexe['stdout']

        for zm12 in zmats0N:
            zm1, zm2 = zm12.split('-')
            success, dexe = qcng.util.execute(['xja2fja'], {'FJOBARC': shelf['fjobarc'][zm12]}, [], **scrkwgs)
            print(zm12, 'xja2fja', success)
            harmout += dexe['stdout']

            success, dexe = qcng.util.execute(['xsymcor'], {}, [], **scrkwgs)
            print(zm12, 'xsymcor', success)
            harmout += dexe['stdout']

        success, dexe = qcng.util.execute(['xjoda'], {}, [], **scrkwgs)
        print('xjoda', success)
        harmout += dexe['stdout']

        for zm in Path(dexe['scratch_directory']).glob('zmat*'):
            print('Removing', zm)
            os.remove(zm)

        scrkwgs.update({'scratch_messy': scratch_messy})
        success, dexe = qcng.util.execute(['xcubic'], {}, ['zmat*'], **scrkwgs)
        print('xcubic', success)
        harmout += dexe['stdout']
        #print('HARMOUT')
        #print(harmout)

        pp.pprint(shelf['zmat'].keys())    
        pp.pprint(shelf['fjobarc'].keys())    



#        os.chdir(psioh.get_default_path() + cfour_tmpdir + '/harm')  # psi_scratch/cfour/harm
#        harmout = run_cfour_module('xjoda')
#        harmout += run_cfour_module('xsymcor')
#        for zm12 in zmats0N:
#            zm1, zm2 = zm12.split('-')
#            with open('FJOBARC', 'w') as handle:
#                handle.write(shelf['fjobarc'][zm12])
#            harmout += run_cfour_module('xja2fja')
#            harmout += run_cfour_module('xsymcor')
#            shutil.move('FJOBARC', 'fja.' + zm12)
#            try:
#                os.remove('zmat' + zm2)
#            except OSError:
#                pass
#        harmout += run_cfour_module('xjoda')
#        harmout += run_cfour_module('xcubic')
#        core.print_out(harmout)
#        with open('harm.out', 'w') as handle:
#            handle.write(harmout)

        # Generate displacements along harmonic normal modes
        #zmatsN0 = [item[-3:] for item in sorted(shelf['zmat'].keys()) if (item[:3] == '000' and item[-3:] != '000')]
        for fl, contents in dexe['outfiles']['zmat*'].items():
            zmN_ = fl[-3:] + '-000'
            shelf['zmat'][zmN_] = contents
            shelf.sync()
            print(f'  CFOUR scratch file {fl} for {zmN_} has been read\n')
            #print('%s\n' % shelf['zmat'][zm12])

        zmatsN0 = [item for item in sorted(shelf['zmat'].keys()) if (item[:3] != '000' and item[-3:] == '000')]
        for zmN0 in zmatsN0:
            zm1, _ = zmN0.split('-')

            # Collect displacements along the normal coordinates generated by the harmonic freq.
            #   Further harmonic freqs are to be run at each of these to produce quartic force field.
            #   To carry these out, generate displacements for findif by gradient at each displacement.

            scrkwgs = {'scratch_directory': config.scratch_directory, 'scratch_messy': True, 'scratch_suffix': f'_{zmN0}'}
            success, dexe = qcng.util.execute(['xjoda'], {'ZMAT': shelf['zmat'][zmN0], 'GENBAS': shelf['genbas']}, [], **scrkwgs)

            scrkwgs.update({'scratch_name': Path(dexe['scratch_directory']).name, 'scratch_exist_ok': True})
            scrkwgs.update({'scratch_messy': scratch_messy})
            success, dexe = qcng.util.execute(['xsymcor'], {}, ['zmat*'], **scrkwgs)

            for fl, contents in dexe['outfiles']['zmat*'].items():
                zm12 = zm1 + '-' + fl[-3:]
                shelf['zmat'][zm12] = contents
                shelf.sync()
                print('  CFOUR scratch file %s for %s has been read\n' % (fl, zm12))
                #print('%s\n' % shelf['zmat'][zm12])

#        zmatsN0 = [item[-3:] for item in sorted(glob.glob('zmat*'))]
#        os.chdir('..')  # psi_scratch/cfour
#        for zm1 in zmatsN0:
#            zm12 = zm1 + '-000'
#            with open(psioh.get_default_path() + cfour_tmpdir + '/harm/zmat' + zm1, 'r') as handle:
#                shelf['zmat'][zm12] = handle.read()
#                shelf.sync()
#                core.print_out('  CFOUR scratch file %s for %s has been read\n' % ('zmat' + zm1, zm12))
#                core.print_out('%s\n' % shelf['zmat'][zm12])
#
#            # Collect displacements along the normal coordinates generated by the harmonic freq.
#            #   Further harmonic freqs are to be run at each of these to produce quartic force field.
#            #   To carry these out, generate displacements for findif by gradient at each displacement.
#            if os.path.exists(zm1):
#                shutil.rmtree(zm1)
#            os.mkdir(zm1)
#            os.chdir(zm1)  # psi_scratch/cfour/004
#            with open('ZMAT', 'w') as handle:
#                handle.write(shelf['zmat'][zm12])
#            shutil.copy2('../harm/GENBAS', 'GENBAS')  # ln -s $ecpdir/ECPDATA $j/ECPDATA
#            with open('partial.out', 'w') as handle:
#                handle.write(run_cfour_module('xjoda'))
#                handle.write(run_cfour_module('xsymcor'))
#
#            # Read the displacements that will form the anharmonic freq
#            zmatsNN = [item[-3:] for item in sorted(glob.glob('zmat*'))]
#            for zm2 in zmatsNN:
#                zm12 = zm1 + '-' + zm2
#                with open(psioh.get_default_path() + cfour_tmpdir + '/' + zm1 + '/zmat' + zm2, 'r') as handle:
#                    shelf['zmat'][zm12] = handle.read()
#                    shelf.sync()
#                    core.print_out('  CFOUR scratch file %s for %s has been read\n' % ('zmat' + zm2, zm12))
#                    core.print_out('%s\n' % shelf['zmat'][zm12])
#            os.chdir('..')  # psi_scratch/cfour

        zmatsNN = [item for item in sorted(shelf['zmat'].keys()) if (item[:3] != '000' and item[-3:] != '000')]

        # S/R: Write distributed input files for anharmonic freq
        if isSowReap:
            os.chdir(current_directory)
            inputSansMol = p4util.format_currentstate_for_input(gradient, lowername, allButMol=True, **kwargs)
            for zm12 in zmatsNN:
                zm1, zm2 = zm12.split('-')

                ifile = vpt2_sow_files(zm12, shelf['linkage'], isC4notP4, isC4fully,
                    shelf['zmat'][zm12], inputSansMol, shelf['genbas'])
                # GENBAS needed here

                with open('VPT2-' + zm12 + '.in', 'w') as handle:
                    handle.write(ifile)

            msg = vpt2_instructions('anharmonic', current_directory, zmatsNN)
            core.print_out(msg)
            print(msg)

        shelf['status'] = 'anharm_jobs_sown'

        # S/R: Pause for distributed calculations
        if isSowReap:
            shelf.close()
            return 0.0

    if shelf['status'] == 'anharm_jobs_sown':
        zmatsNN = [item for item in sorted(shelf['zmat'].keys()) if (item[:3] != '000' and item[-3:] != '000')]

        # S/R: Check that distributed calcs all completed correctly
        if isSowReap:
            msg = vpt2_instructions('anharmonic', current_directory, zmatsNN)
            core.print_out(msg)
            isOk, msg = sown_jobs_status(current_directory, 'VPT2', zmatsNN, 
                reap_job_validate, shelf['linkage'], 
                ['CURRENT ENERGY', 'CURRENT DIPOLE', 'CURRENT GRADIENT'])
            core.print_out(msg)
            print(msg)
            if not isOk:
                shelf.close()
                return 0.0

        # Collect all results from gradients forming the anharmonic freq
        for zmNN in zmatsNN:
            zm1, zm2 = zmNN.split('-')
            if zmNN not in shelf['fjobarc']:
                print('{:_^45}'.format(f'  VPT2 Computation: {zmNN}'))

                fjobarc = vpt2_reaprun_files(zmNN, shelf['linkage'], isSowReap, isC4notP4, isC4fully,
                    shelf['zmat'][zmNN],
                    lowername, kwargs, shelf['genbas'], config, package, scratch_messy=scratch_messy)
                shelf['fjobarc'][zmNN] = fjobarc
                shelf.sync()
        shelf['status'] = 'anharm_jobs_reaped'

    if shelf['status'] == 'anharm_jobs_reaped':
        zmats0N = [item for item in sorted(shelf['zmat'].keys()) if (item[:3] == '000' and item[-3:] != '000')]
        zmatsN0 = [item for item in sorted(shelf['zmat'].keys()) if (item[:3] != '000' and item[-3:] == '000')]
        zmatsNN = [item for item in sorted(shelf['zmat'].keys()) if (item[:3] != '000' and item[-3:] != '000')]

        print('{:_^45}'.format('  VPT2 Results: Harmonic  '))

        # Process the gradients into harmonic freq
        scrkwgs = {'scratch_directory': config.scratch_directory, 'scratch_messy': True, 'scratch_suffix': '_000final'}
        success, dexe = qcng.util.execute(['xjoda'], {'ZMAT': shelf['zmat']['000-000'], 'GENBAS': shelf['genbas']}, [], **scrkwgs)
        anharmout = dexe['stdout']

        scrkwgs.update({'scratch_name': Path(dexe['scratch_directory']).name, 'scratch_exist_ok': True})
        success, dexe = qcng.util.execute(['xsymcor'], {}, [], **scrkwgs)
        anharmout += dexe['stdout']

        for zm12 in zmats0N:
            zm1, zm2 = zm12.split('-')
            success, dexe = qcng.util.execute(['xja2fja'], {'FJOBARC': shelf['fjobarc'][zm12]}, [], **scrkwgs)
            anharmout += dexe['stdout']
            success, dexe = qcng.util.execute(['xsymcor'], {}, [], **scrkwgs)
            anharmout += dexe['stdout']

        success, dexe = qcng.util.execute(['xjoda'], {}, [], **scrkwgs)
        anharmout += dexe['stdout']
        scrkwgs.update({'scratch_messy': scratch_messy})
        success, dexe = qcng.util.execute(['xcubic'], {}, ['zmat*', 'JOBARC', 'JAINDX'], as_binary=['JOBARC', 'JAINDX'], **scrkwgs)
        anharmout += dexe['stdout']

        jobarc0 = dexe['outfiles']['JOBARC']
        jaindx0 = dexe['outfiles']['JAINDX']

        # Process the gradients into harmonic freq at each normco displaced point
        os.chdir('..')  # psi_scratch/cfour
        for zm1_ in zmatsN0:
            zm1, _ = zm1_.split('-')

            scrkwgs = {'scratch_directory': config.scratch_directory, 'scratch_messy': True, 'scratch_suffix': f'_{zm1}final'}
            success, dexe = qcng.util.execute(['xjoda'], {'ZMAT': shelf['zmat'][zm1_], 'GENBAS': shelf['genbas']}, [], **scrkwgs)
            anharmout = dexe['stdout']

            scrkwgs.update({'scratch_name': Path(dexe['scratch_directory']).name, 'scratch_exist_ok': True})
            success, dexe = qcng.util.execute(['xsymcor'], {}, [], **scrkwgs)
            anharmout += dexe['stdout']

            for zm12 in [item for item in zmatsNN if (item[:3] == zm1 and item[-3:] != '000')]:
                _, zm2 = zm12.split('-')
                print(zm12, shelf['fjobarc'][zm12])

                success, dexe = qcng.util.execute(['xja2fja'], {'FJOBARC': shelf['fjobarc'][zm12]}, [], **scrkwgs)
                anharmout += dexe['stdout']

                success, dexe = qcng.util.execute(['xsymcor'], {}, [], **scrkwgs)
                anharmout += dexe['stdout']

                os.remove(Path(dexe['scratch_directory']) / 'FJOBARC')

            success, dexe = qcng.util.execute(['xjoda'], {}, [], **scrkwgs)
            anharmout += dexe['stdout']

            scrkwgs.update({'scratch_messy': scratch_messy})
            success, dexe = qcng.util.execute(['xja2fja'], {}, ['FJOBARC'], **scrkwgs)
            anharmout += dexe['stdout']
            shelf['fjobarc'][zm1_] = dexe['outfiles']['FJOBARC']
            shelf.sync()

            print('PARTIAL', zm1_, '\n', anharmout)


        # Process the harmonic freqs at normco displacements into anharmonic freq
        print('{:_^45}'.format('  VPT2 Results: Anharmonic  '))

        pprint.pprint(shelf['zmat'].keys())
        pprint.pprint(shelf['fjobarc'].keys())

        scrkwgs = {'scratch_directory': config.scratch_directory, 'scratch_messy': True, 'scratch_suffix': '_000anharm'}
        success, dexe = qcng.util.execute(['ls'], {'JOBARC': jobarc0, 'JAINDX': jaindx0}, [], as_binary=['JOBARC', 'JAINDX'], **scrkwgs)
        anharmout = ''
        scrkwgs.update({'scratch_name': Path(dexe['scratch_directory']).name, 'scratch_exist_ok': True})

        for zm1_ in zmatsN0:
            zm1, _ = zm1_.split('-')
            success, dexe = qcng.util.execute(['xja2fja'], {'FJOBARC': shelf['fjobarc'][zm1_]}, [], **scrkwgs)
            print(zm1_, 'xja2fja', success)
            anharmout += dexe['stdout']

            success, dexe = qcng.util.execute(['xcubic'], {}, [], **scrkwgs)
            print(zm1_, 'xcubic', success)
            anharmout += dexe['stdout']

        print(anharmout)  # anharm.out

        shelf['status'] = 'vpt2_completed'

    # Finish up
    shelf.close()


def vpt2_sow_files(item, linkage, isC4notP4, isC4fully, zmat, inputSansMol, inputGenbas):
    """Provided with the particular displacement number *item* and the
    associated *zmat* file contents and *linkage*, and common contents
    *inputSansMol*, returns contents of input file to be sown.

    """
    inputReapOrders = r"""
print_variables()

print_out('VPT2 RESULT: linkage {0} for item {1} yields CURRENT ENERGY being %r\n' % (variable('CURRENT ENERGY')))
print_out('VPT2 RESULT: linkage {0} for item {1} yields CURRENT GRADIENT being %r\n' % (p4util.mat2arr(core.get_gradient())))
print_out('VPT2 RESULT: linkage {0} for item {1} yields CURRENT DIPOLE being [%r, %r, %r]\n' % (variable('CURRENT DIPOLE X'), variable('CURRENT DIPOLE Y'), variable('CURRENT DIPOLE Z')))
""".format(linkage, item)

    # Direct Cfour for gradients
    if isC4fully:
        inputString = zmat
        with open('VPT2-GENBAS', 'w') as handle:
            handle.write(inputGenbas)

#!BR#    # Cfour for gradients
#!BR#    elif isC4notP4:
#!BR#        # GENBAS needed here
#!BR#        inputString = 'extracted_genbas = """\n' + inputGenbas.replace('\n\n', '\nblankline\n') + '\n"""\n\n'
#!BR#        inputString += """cfour {\n%s\n}\n\nenergy('cfour', genbas=extracted_genbas)\n\n""" % (zmat)
#!BR#        inputString += inputReapOrders
#!BR#        inputString += r"""
#!BR#print_out('VPT2 RESULT: linkage {0} for item {1} yields CURRENT MOLECULE being %r\n' % (get_active_molecule().create_psi4_string_from_molecule()))
#!BR#""".format(linkage, item)

    # Psi4 for gradients
    else:
        inputString = p4util.format_molecule_for_input(
            qcdb.cfour.harvest_zmat(zmat).create_psi4_string_from_molecule(),
            name='disp' + item[:3] + item[-3:])
        inputString += inputSansMol
        inputString += inputReapOrders

    return inputString


#def vpt2_reaprun_files(item, linkage, isSowReap, isC4notP4, isC4fully, zmat, outdir, scrdir, c4scrdir, lowername, kwargs):
def vpt2_reaprun_files(item, linkage, isSowReap, isC4notP4, isC4fully, zmat, lowername, kwargs, genbas, config, package, scratch_messy):
    """Provided with the particular displacement number *item* and the
    associated *zmat* file with geometry and *linkage*, returns the
    FJOBARC contents. Depending on the mode settings of *isC4notP4*,
    *isSowReap*, and *isC4fully*, either runs (using *lowername* and
    *kwargs*) or reaps contents. *outdir* is where psi4 was invoked,
    *scrdir* is the psi4 scratch directory, and *c4scrdir* is Cfour
    scratch directory within.

    """
    # Extract qcel.models.Molecule at findif orientation
    zmmol = harvest_zmat(zmat)

    # Cfour S/R Direct for gradients
    if isC4fully:
        with open('VPT2-' + item + '.fja', 'r') as handle:
            fjobarc = handle.read()

#!BR#    # Cfour for gradients
#!BR#    elif isC4notP4:
#!BR#
#!BR#        # S/R: Reap results from output file
#!BR#        if isSowReap:
#!BR#            isOk, msg, results = reap_job_validate(outdir, 'VPT2', item, linkage,
#!BR#                ['CURRENT ENERGY', 'CURRENT DIPOLE', 'CURRENT GRADIENT', 'CURRENT MOLECULE'])
#!BR#            if not isOk:
#!BR#                raise ValidationError(msg)
#!BR#
#!BR#            fje = results['CURRENT ENERGY']
#!BR#            fjgrd = results['CURRENT GRADIENT']
#!BR#            fjdip = [item / constants.dipmom_au2debye for item in results['CURRENT DIPOLE']]
#!BR#            c4mol = qcdb.Molecule(results['CURRENT MOLECULE'])
#!BR#            c4mol.update_geometry()
#!BR#
#!BR#        # C: Run the job and collect results
#!BR#        else:
#!BR#            # Prepare Cfour skeleton calc directory
#!BR#            {'ZMAT': zmat, 'GENBAS': shelf['genbas']}, 
#!BR#            os.chdir(scrdir + c4scrdir)  # psi_scratch/cfour
#!BR#            if os.path.exists('scr.' + item):
#!BR#                shutil.rmtree('scr.' + item)
#!BR#            os.mkdir('scr.' + item)
#!BR#            os.chdir('scr.' + item)  # psi_scratch/cfour/scr.000-004
#!BR#            with open('ZMAT', 'w') as handle:
#!BR#                handle.write(zmat)
#!BR#            shutil.copy2('../harm/GENBAS', 'GENBAS')
#!BR#
#!BR#            #os.chdir(scrdir + '/scr.' + item)
#!BR#            #run_cfour_module('xja2fja')
#!BR#            #with open('FJOBARC', 'r') as handle:
#!BR#            #    fjobarc = handle.read()
#!BR#
#!BR#            # Run Cfour calc using ZMAT & GENBAS in scratch, outdir redirects to outfile
#!BR#            os.chdir(outdir)  # current_directory
#!BR#            core.get_active_molecule().set_name('blank_molecule_psi4_yo')
#!BR#            energy('cfour', path=c4scrdir + '/scr.' + item)
#!BR##            os.chdir(scrdir + '/scr.' + item)
#!BR#
#!BR#            fje = core.variable('CURRENT ENERGY')
#!BR#            fjgrd = p4util.mat2arr(core.get_gradient())
#!BR#            fjdip = [core.variable('CURRENT DIPOLE X') / constants.dipmom_au2debye,
#!BR#                     core.variable('CURRENT DIPOLE Y') / constants.dipmom_au2debye,
#!BR#                     core.variable('CURRENT DIPOLE Z') / constants.dipmom_au2debye]
#!BR#            # TODO switch DIPOLE to vector
#!BR#            c4mol = qcdb.Molecule(core.get_active_molecule().create_psi4_string_from_molecule())
#!BR#            c4mol.update_geometry()
#!BR#
#!BR#        # Get map btwn ZMAT and C4 orientation, then use it, grad and dipole to forge FJOBARC file
#!BR#        fjobarc = format_fjobarc(fje,
#!BR#            *backtransform(chgeMol=zmmol, permMol=c4mol), gradient=fjgrd, dipole=fjdip)

    # Psi4 for gradients
    else:

        # Run Cfour skeleton calc and extract qcdb.Molecule at needed C4 orientation
        scrkwgs = {'scratch_directory': config.scratch_directory, 'scratch_messy': True, 'scratch_suffix': f'_{item}skel'}
        success, dexe = qcng.util.execute(['xjoda'], {'ZMAT': zmat, 'GENBAS': genbas},  [], **scrkwgs)
        skel = dexe['stdout']

        scrkwgs.update({'scratch_name': Path(dexe['scratch_directory']).name, 'scratch_exist_ok': True})
        success, dexe = qcng.util.execute(['xvmol'], {}, [], **scrkwgs)
        skel += dexe['stdout']

        scrkwgs.update({'scratch_messy': scratch_messy})
        success, dexe = qcng.util.execute(['xvmol2ja'], {}, ['JOBARC', 'JAINDX'], as_binary=['JOBARC', 'JAINDX'], **scrkwgs)
        skel += dexe['stdout']

        print(skel)  # partial.out

        print(f"""  CFOUR scratch file 'JOBARC (binary)' for {item} has been read""")
        c4mol = jajo2mol(getrec([b'COORD   ', b'ATOMCHRG', b'MAP2ZMAT', b'IFLAGS  '], dexe['outfiles']['JOBARC'], dexe['outfiles']['JAINDX']))

        # S/R: Reap results from output file
        if isSowReap:
            isOk, msg, results = reap_job_validate(outdir, 'VPT2', item, linkage,
                ['CURRENT ENERGY', 'CURRENT DIPOLE', 'CURRENT GRADIENT'])
            if not isOk:
                raise ValidationError(msg)

            fje = results['CURRENT ENERGY']
            fjgrd = results['CURRENT GRADIENT']
            fjdip = [item / constants.dipmom_au2debye for item in results['CURRENT DIPOLE']]

        # C: Run the job and collect results
        else:

            resi = ResultInput(
                **{
                    'driver': 'gradient',
                    'extras': {
                        'qcdb:options': copy.deepcopy(pe.nu_options),
                    },
                    'model': {
                        'method': lowername,
                        'basis': '(auto)',
                    },
                    'molecule': zmmol,
                    'provenance': provenance_stamp(__name__),
                })
            res = qcng.compute(resi, 'qcdb-' + package, raise_error=True)
            #pp.pprint(res.dict())

            fje = res.extras['qcdb:qcvars']['CURRENT ENERGY']['data']
            fjgrd = res.extras['qcdb:qcvars']['CURRENT GRADIENT']['data']
            #au2debye = qcel.constants.get('dipmom_au2debye', return_tuple=True).data
            # TODO switch DIPOLE to vector
            fjdip = [float(res.extras['qcdb:qcvars']['CURRENT DIPOLE X']['data']) / qcel.constants.dipmom_au2debye,
                     float(res.extras['qcdb:qcvars']['CURRENT DIPOLE Y']['data']) / qcel.constants.dipmom_au2debye,
                     float(res.extras['qcdb:qcvars']['CURRENT DIPOLE Z']['data']) / qcel.constants.dipmom_au2debye]

        # Transform results into C4 orientation (defined by c4mol) & forge FJOBARC file
        fjobarc = format_fjobarc(fje, *backtransform(chgeMol=zmmol, permMol=c4mol, chgeGrad=fjgrd, chgeDip=fjdip))

    return fjobarc


def vpt2_instructions(stage, dir, zmats):
    """Stores all the instructions to the user for running
    :py:func:`~wrappers_cfour.vpt2` in sowreap mode. Depending on the
    *stage*, Pieces together instruction strings for the appropriate
    *stage* individualized by working directory *dir* and sown inputs
    *zmats* information.

    """
    stepFiles = ''
    for zm12 in sorted(zmats):
        stepFiles += """             psi4 %-27s %-27s\n""" % ('VPT2-' + zm12 + '.in', 'VPT2-' + zm12 + '.out')

    step0 = """
    The vpt2 sow/reap procedure has been selected through mode='sowreap'. This
    output file, the corresponding input file, and the data persistence file
    must not be edited by the user over the course of the sow/reap procedure.
    Throughout, psi4 can be invoked to move to the next stage of the procedure
    or to tally up the 'sown' jobs. This output file is overwritten each time
    psi4 is invoked, but all results and instructions accumulate.

    This procedure involves two stages of distributed calculations, harmonic and
    anharmonic, and a mimimum of three invokations of psi4 on the original input
    file (including the one that initially generated this text). From the input
    geometry (0), displacements are generated for which gradients are required.
    Input files for these are 'sown' in the current directory (1). Upon
    completion, their output files are 'reaped' into a harmonic force field (2).
    At displacements along the normal coordinates, further displacements are
    generated for which gradients are required. Input files for these are again
    'sown' in the current directory (3). Upon completion, their output files are
    'reaped' into an anharmonic force field (4), terminating the vpt2 procedure.
    Follow the instructions below to continue.

    (0)  Read Only
    --------------
       %s
       %s
       %s

""" % (dir + '/' + os.path.splitext(core.outfile_name())[0] + '.in',
       dir + '/' + core.outfile_name(),
       dir + '/' + os.path.splitext(core.outfile_name())[0] + '.shelf')
    step1 = """
    (1)  Sow
    --------
       Run all of the VPT2-000-*.in input files on any variety of computer
       architecture. The output file names must be as given below (default).

"""
    step2 = """
    (2)  Reap
    ---------
       Gather all the resulting output files in this directory along with the
       three read-only files from (0). Invoke psi4 again. The job will be
       trivial in length (unless sto-3g integrals on the molecule are costly)
       and give results for the harmonic frequency stage in this output file. It
       will also supply the next set of instructions.

             psi4 %-27s %-27s

""" % (os.path.splitext(core.outfile_name())[0] + '.in', core.outfile_name())
    step3 = """
    (3)  Sow
    --------
       Run all of the VPT2-*-*.in input files on any variety of computer
       architecture. The output file names must be as given below (default).

"""
    step4 = """
    (4)  Reap
    ---------
       Gather all the resulting output files in this directory along with the
       three read-only files from (0). Invoke psi4 again. The job will be
       trivial in length (unless sto-3g integrals on the molecule are costly)
       and give results for the harmonic and anharmonic frequency stages in this
       output file.

             psi4 %-27s %-27s

""" % (os.path.splitext(core.outfile_name())[0] + '.in', core.outfile_name())

    if stage == 'harmonic':
        instructions = step0 + step1 + stepFiles + step2
    elif stage == 'anharmonic':
        instructions = step0 + step3 + stepFiles + step4

    return instructions


def sown_jobs_status(dir, prefix, zmats, validate_func=None, linkage=None, keys=None):
    """Evaluate the output file status of jobs in *zmats* which should
    exist at *dir* + '/' + prefix + '-' + job + '.out'. Returns string with
    formatted summary of job status and boolean of whether all complete.
    Return boolean *isOk* signals whether all *zmats* have completed and,
    if *validate_func* present, are validated.

    """
    isOk = True
    msgError = ''
    instructions = '\n'
    instructions += '{:_^45}'.format(prefix + ' Status: ' + datetime.datetime.now().strftime("%Y-%m-%d %H:%M"), strNotOutfile=True)
    instructions += '\n'
    for job in sorted(zmats):
        outfile = dir + '/' + prefix + '-' + job + '.out'
        fjafile = dir + '/' + prefix + '-' + job + '.fja'
        formatArgs = [prefix + '-' + job, '', '', '', '']
        if os.path.isfile(outfile):
            with open(outfile, 'r') as handle:
                for line in handle:
                    if line.find('Buy a developer a beer!') > -1:
                        formatArgs[3] = 'Completed'
                        if reap_job_validate is not None:
                            isOkJob, msg, temp = reap_job_validate(dir, prefix, job, linkage, keys)
                            if isOkJob:
                                formatArgs[4] = '& Validated'
                            else:
                                isOk = False
                                msgError += msg
                                formatArgs[4] = 'INVALID'
                        break
                else:
                    isOk = False
                    formatArgs[2] = 'Running'
        elif os.path.isfile(fjafile):
            formatArgs[3] = 'Completed'
        else:
            isOk = False
            formatArgs[1] = 'Waiting'
        instructions += """             {0:<27} {1:^10} {2:^10} {3:^10} {4:^10}\n""".format(*formatArgs)
    instructions += '\n' + msgError + '\n\n'

    return isOk, instructions


def reap_job_validate(dir, prefix, item, linkage, keys):
    """For a given output file whose path is constructed with
    *dir* + '/' + *prefix* + '-' + *item* + '.out', tests that the file
    exists and has *prefix* RESULTS lines for each piece of information
    requested in list *keys* and that those lines correspond to the
    appropriate *linkage* and *item*. Returns *keys* along with their
    scanned values in dict *reapings*, along with error and success
    messages in *instructions* and a boolean *isOk* indicating whether
    all *keys* reaped sucessfully.

    """
    isOk = True
    instructions = ''
    reapings = {}
    outfile = dir + '/' + prefix + '-' + item + '.out'

    try:
        with open(outfile, 'r') as handle:
            for line in handle:
                if line.find(prefix + ' RESULT:') == 0:
                    sline = line.split()
                    if sline[2:7] == ['linkage', str(linkage), 'for', 'item', item]:
                        yieldsAt = line.find('yields')
                        beingAt = line.find('being')
                        if beingAt > yieldsAt > -1:
                            key = line[yieldsAt + 6:beingAt].strip()
                            val = line[beingAt + 5:].strip()
                            if key in keys:
                                reapings[key] = eval(val)
                                #core.print_out('  CFOUR scratch file %s for %s has been read\n' % ('JOBARC', zm12))
                        else:
                            isOk = False
                            instructions += """Outfile file %s
    has corrupted sowreap result line:\n%s\n\n""" % (outfile, line)
                    else:
                        isOk = False
                        instructions += """Outfile file %s
    has sowreap result of either incompatible linkage (observed: %s, expected: %s)
    or incompatible job affiliation (observed: %s, expected: %s).\n\n""" % \
                            (outfile, sline[3], linkage, sline[6], item)
            else:
                if len(reapings) != len(keys):
                    isOk = False
                    instructions += """Output file %s
    has missing results (observed: %s, expected: %s).\n\n""" % \
                        (outfile, reapings.keys(), keys)
    except IOError:
        isOk = False
        instructions += """Output file %s
    that was judged present and complete at the beginning of this
    job is now missing. Replace it and invoke psi4 again.\n\n""" % (outfile)

    # return file contents in instructions
    return isOk, instructions, reapings
