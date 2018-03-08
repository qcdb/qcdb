#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2017 The Psi4 Developers.
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

"""Module with a *procedures* dictionary specifying available quantum
chemical methods and functions driving the main quantum chemical
functionality, namely single-point energies, geometry optimizations,
properties, and vibrational frequency calculations.

"""
from __future__ import print_function
from __future__ import absolute_import
#   import os
#   import re
#   import sys
#   import shutil

import copy
import pprint
pp = pprint.PrettyPrinter(width=120)

#   import numpy as np
#   
#   import qcdb

from . import pe
from . import driver_util
from . import driver_helpers
#from . import cbs_driver
##   from psi4.driver import driver_nbody
##   from psi4.driver import p4util
from .proc_table import procedures
##   from psi4.driver.procrouting import *
##   from psi4.driver.p4util.exceptions import *
##   # never import wrappers or aliases into this file



#   def _energy_is_invariant(gradient, stationary_criterion=1.e-2):
#       """Polls options and probes `gradient` to return whether current method
#       and system expected to be invariant to translations and rotations of
#       the coordinate system.
#   
#       """
#       stationary_point = gradient.rms() < stationary_criterion  # 1.e-2 pulled out of a hat
#   
#       mol = core.get_active_molecule()
#       efp_present = hasattr(mol, 'EFP')
#   
#       translations_projection_sound = (not core.get_option('SCF', 'EXTERN') and
#                                        not core.get_option('SCF', 'PERTURB_H') and
#                                        not efp_present)
#       rotations_projection_sound = (translations_projection_sound and 
#                                     stationary_point)
#   
#       return translations_projection_sound, rotations_projection_sound
   
   
#   MOVED   #def energy(name, **kwargs):
#   MOVED   #    r"""Function to compute the single-point electronic energy."""
#   MOVED   #
#   MOVED   #    from . import endorsed_plugins
#   MOVED   #    kwargs = driver_util.kwargs_lower(kwargs)
#   MOVED   #
#   MOVED   #    # Bounce if name is function
#   MOVED   #    if hasattr(name, '__call__'):
#   MOVED   #        return name(energy, kwargs.pop('label', 'custom function'), ptype='energy', **kwargs)
#   MOVED   #
#   MOVED   #    # Allow specification of methods to arbitrary order
#   MOVED   #    lowername = name.lower()
#   MOVED   #    lowername, level = driver_helpers._parse_arbitrary_order(lowername)
#   MOVED   #    if level:
#   MOVED   #        kwargs['level'] = level
#   MOVED   #
#   MOVED   #    # Make sure the molecule the user provided is the active one
#   MOVED   #    molecule = kwargs.pop('molecule', driver_helpers.get_active_molecule())
#   MOVED   #    molecule.update_geometry()
#   MOVED   #
#   MOVED   #    if len(pe.nu_options.scroll) == 0:
#   MOVED   #        print('EMPTY OPT')
#   MOVED   #        pe.load_nu_options()
#   MOVED   #
#   MOVED   ##    # Bounce to CP if bsse kwarg
#   MOVED   ##    if kwargs.get('bsse_type', None) is not None:
#   MOVED   ##        return driver_nbody.nbody_gufunc(energy, name, ptype='energy', **kwargs)
#   MOVED   #
#   MOVED   #    # Bounce to CBS if "method/basis" name
#   MOVED   #    if '/' in lowername:
#   MOVED   #        return cbs_driver._cbs_gufunc(energy, name, ptype='energy', molecule=molecule, **kwargs)
#   MOVED   #
#   MOVED   #    # Commit to procedures['energy'] call hereafter
#   MOVED   #    return_wfn = kwargs.pop('return_wfn', False)
#   MOVED   #    pe.active_qcvars = {}
#   MOVED   #
#   MOVED   ##    #for precallback in hooks['energy']['pre']:
#   MOVED   ##    #    precallback(lowername, **kwargs)
#   MOVED   ##
#   MOVED   ##    optstash = driver_util._set_convergence_criterion('energy', lowername, 6, 8, 6, 8, 6)
#   MOVED   ##
#   MOVED   ##    # Before invoking the procedure, we rename any file that should be read.
#   MOVED   ##    # This is a workaround to do restarts with the current PSI4 capabilities
#   MOVED   ##    # before actual, clean restarts are put in there
#   MOVED   ##    # Restartfile is always converted to a single-element list if
#   MOVED   ##    # it contains a single string
#   MOVED   ##    if 'restart_file' in kwargs:
#   MOVED   ##        restartfile = kwargs['restart_file']  # Option still available for procedure-specific action
#   MOVED   ##        if restartfile != list(restartfile):
#   MOVED   ##            restartfile = [restartfile]
#   MOVED   ##        # Rename the files to be read to be consistent with psi4's file system
#   MOVED   ##        for item in restartfile:
#   MOVED   ##            name_split = re.split(r'\.', item)
#   MOVED   ##            filenum = name_split[len(name_split) - 1]
#   MOVED   ##            try:
#   MOVED   ##                filenum = int(filenum)
#   MOVED   ##            except ValueError:
#   MOVED   ##                filenum = 32  # Default file number is the checkpoint one
#   MOVED   ##            psioh = core.IOManager.shared_object()
#   MOVED   ##            psio = core.IO.shared_object()
#   MOVED   ##            filepath = psioh.get_file_path(filenum)
#   MOVED   ##            namespace = psio.get_default_namespace()
#   MOVED   ##            pid = str(os.getpid())
#   MOVED   ##            prefix = 'psi'
#   MOVED   ##            targetfile = filepath + prefix + '.' + pid + '.' + namespace + '.' + str(filenum)
#   MOVED   ##            shutil.copy(item, targetfile)
#   MOVED   #
#   MOVED   #    print('QWER', pe.nu_options.print_changed())
#   MOVED   #    package = driver_util.get_package(lowername, kwargs)
#   MOVED   #    #for k, v in pkgprefix.items():
#   MOVED   #    #    if lowername.startswith(k):
#   MOVED   #    #        package = v
#   MOVED   #    #        break
#   MOVED   #    #else:
#   MOVED   #    #    package = kwargs.get('package', 'psi4')
#   MOVED   #    #print('\nENE calling', 'procedures', package, lowername, 'with', lowername, molecule, pe.nu_options, kwargs)
#   MOVED   #    #jobrec = procedures['energy'][package][lowername](lowername, molecule=molecule, options=pe.active_options, **kwargs)
#   MOVED   #    jobrec = procedures['energy'][package][lowername](lowername, molecule=molecule, options=pe.nu_options, ptype='energy', **kwargs)
#   MOVED   #
#   MOVED   ##    for postcallback in hooks['energy']['post']:
#   MOVED   ##        postcallback(lowername, wfn=wfn, **kwargs)
#   MOVED   ##
#   MOVED   ##    optstash.restore()
#   MOVED   #    #jobrec.pop('raw_output')  # just to moderate printint to screen
#   MOVED   #    pp.pprint(jobrec)
#   MOVED   #    pe.active_qcvars = copy.deepcopy(jobrec['qcvars'])
#   MOVED   #
#   MOVED   #    if return_wfn:  # TODO current energy safer than wfn.energy() for now, but should be revisited
#   MOVED   #
#   MOVED   ##        # TODO place this with the associated call, very awkward to call this in other areas at the moment
#   MOVED   ##        if lowername in ['efp', 'mrcc', 'dmrg', 'psimrcc']:
#   MOVED   ##            core.print_out("\n\nWarning! %s does not have an associated derived wavefunction." % name)
#   MOVED   ##            core.print_out("The returned wavefunction is the incoming reference wavefunction.\n\n")
#   MOVED   ##        elif 'sapt' in lowername:
#   MOVED   ##            core.print_out("\n\nWarning! %s does not have an associated derived wavefunction." % name)
#   MOVED   ##            core.print_out("The returned wavefunction is the dimer SCF wavefunction.\n\n")
#   MOVED   #
#   MOVED   #        return (float(jobrec['qcvars']['CURRENT ENERGY'].data), jobrec)
#   MOVED   #    else:
#   MOVED   #        return float(jobrec['qcvars']['CURRENT ENERGY'].data)
#   MOVED   #        # float() is for decimal.Decimal


def hessian(name, **kwargs):
#    r"""Function complementary to :py:func:`~frequency`. Computes force
#    constants, deciding analytic, finite difference of gradients, or
#    finite difference of energies.
#
#    :returns: :py:class:`~psi4.core.Matrix` |w--w| Total non-mass-weighted electronic Hessian in Hartrees/Bohr/Bohr.
#
#    :returns: (:py:class:`~psi4.core.Matrix`, :py:class:`~psi4.core.Wavefunction`) |w--w| Hessian and wavefunction when **return_wfn** specified.
#
#    :examples:
#
#    >>> # [1] Frequency calculation without thermochemical analysis
#    >>> hessian('mp3')
#
#    >>> # [2] Frequency calc w/o thermo analysis getting the Hessian
#    >>> #     in file, core.Matrix, and np.array forms
#    >>> set hessian_write on
#    >>> H, wfn = hessian('ccsd', return_wfn=True)
#    >>> wfn.hessian().print_out()
#    >>> np.array(H)
#
#    """
    from . import endorsed_plugins
    kwargs = driver_util.kwargs_lower(kwargs)
    text = ''

#    # Bounce to CP if bsse kwarg (someday)
#    if kwargs.get('bsse_type', None) is not None:
#        raise ValidationError("Hessian: Cannot specify bsse_type for hessian yet.")
#
#    # Figure out what kind of gradient this is
#    if hasattr(name, '__call__'):
#        if name.__name__ in ['cbs', 'complete_basis_set']:
#            gradient_type = 'cbs_wrapper'
#        else:
#            # Bounce to name if name is non-CBS function
#            gradient_type = 'custom_function'
#    elif '/' in name:
#        gradient_type = 'cbs_gufunc'
#    else:
#        gradient_type = 'conventional'
#
#    if gradient_type != 'conventional':
#        raise ValidationError("Hessian: Does not yet support more advanced input or custom functions.")

    lowername = name.lower()
    package = driver_util.get_package(lowername, kwargs)

#    # Check if this is a CBS extrapolation
#    if "/" in lowername:
#        return driver_cbs._cbs_gufunc('hessian', lowername, **kwargs)

    return_wfn = kwargs.pop('return_wfn', False)
#    core.clean_variables()
    dertype = 2

#    # Prevent methods that do not have associated energies
#    if lowername in energy_only_methods:
#        raise ValidationError("hessian('%s') does not have an associated hessian" % name)

#    optstash = p4util.OptionsState(
#        ['FINDIF', 'HESSIAN_WRITE'],
#        ['FINDIF', 'FD_PROJECT'],
#        )

#    # Allow specification of methods to arbitrary order
#    lowername, level = driver_util._parse_arbitrary_order(lowername)
#    if level:
#        kwargs['level'] = level

    dertype = driver_util.find_derivative_type('hessian', package, lowername, kwargs.pop('freq_dertype', kwargs.pop('dertype', None)))

    # Make sure the molecule the user provided is the active one
    molecule = kwargs.pop('molecule', driver_helpers.get_active_molecule())
    molecule.update_geometry()

    if len(pe.nu_options.scroll) == 0:
        print('EMPTY OPT')
        pe.load_nu_options()

#    # S/R: Mode of operation- whether finite difference freq run in one job or files farmed out
#    freq_mode = kwargs.pop('mode', 'continuous').lower()
#    if freq_mode == 'continuous':
#        pass
#    elif freq_mode == 'sow':
#        if dertype == 2:
#            raise ValidationError("""Frequency execution mode 'sow' not valid for analytic Hessian calculation.""")
#    elif freq_mode == 'reap':
#        freq_linkage = kwargs.get('linkage', None)
#        if freq_linkage is None:
#            raise ValidationError("""Frequency execution mode 'reap' requires a linkage option.""")
#    else:
#        raise ValidationError("""Frequency execution mode '%s' not valid.""" % (freq_mode))

    # Set method-dependent scf convergence criteria (test on procedures['energy'] since that's guaranteed)
#    optstash_conv = driver_util._set_convergence_criterion('energy', lowername, 8, 10, 8, 10, 8)

#    # Select certain irreps
#    irrep = kwargs.get('irrep', -1)
#    if irrep == -1:
#        pass  # do all irreps
#    else:
#        irrep = driver_util.parse_cotton_irreps(irrep, molecule.schoenflies_symbol())
#        irrep -= 1  # A1 irrep is externally 1, internally 0
#        if dertype == 2:
#            core.print_out("""hessian() switching to finite difference by gradients for partial Hessian calculation.\n""")
#            dertype = 1

#    # At stationary point?
#    G0 = gradient(lowername, molecule=molecule, **kwargs)
#    translations_projection_sound, rotations_projection_sound = _energy_is_invariant(G0)
#    core.print_out('\n  Based on options and gradient (rms={:.2E}), recommend {}projecting translations and {}projecting rotations.\n'.
#                   format(G0.rms(), '' if translations_projection_sound else 'not ',
#                   '' if rotations_projection_sound else 'not '))
#    if not core.has_option_changed('FINDIF', 'FD_PROJECT'):
#        core.set_local_option('FINDIF', 'FD_PROJECT', rotations_projection_sound)

    # Does an analytic procedure exist for the requested method?
    if dertype == 2:
        text += """qcdb.hessian() will perform analytic frequency computation.\n"""

        # We have the desired method. Do it.
        jobrec = procedures['hessian'][package][lowername](lowername, molecule=molecule, options=pe.nu_options, ptype='hessian', **kwargs)
#        wfn.set_gradient(G0)
#        optstash.restore()
#        optstash_conv.restore()
#
#        # TODO: check that current energy's being set to the right figure when this code is actually used
#        core.set_variable('CURRENT ENERGY', wfn.energy())
        # TODO who's setting CURRENT?
        pp.pprint(jobrec)
        pe.active_qcvars = copy.deepcopy(jobrec['qcvars'])

        if return_wfn:
            return (jobrec['qcvars']['CURRENT HESSIAN'].data, jobrec)
        else:
            return jobrec['qcvars']['CURRENT HESSIAN'].data

    elif dertype == 1:
        raise FeatureNotImplemented("""hessian(dertype=1)""")

#        core.print_out("""hessian() will perform frequency computation by finite difference of analytic gradients.\n""")
#
#        # Shifting the geometry so need to copy the active molecule
#        moleculeclone = molecule.clone()
#
#        # Obtain list of displacements
#        displacements = core.fd_geoms_freq_1(moleculeclone, irrep)
#        moleculeclone.reinterpret_coordentry(False)
#        moleculeclone.fix_orientation(True)
#
#        # Record undisplaced symmetry for projection of displaced point groups
#        core.set_parent_symmetry(molecule.schoenflies_symbol())
#
#        ndisp = len(displacements)
#        print(""" %d displacements needed.""" % ndisp)
#        gradients = []
#        energies = []
#
#        # S/R: Write instructions for sow/reap procedure to output file and reap input file
#        if freq_mode == 'sow':
#            instructionsO = """\n#    The frequency sow/reap procedure has been selected through mode='sow'. In addition\n"""
#            instructionsO += """#    to this output file (which contains no quantum chemical calculations), this job\n"""
#            instructionsO += """#    has produced a number of input files (FREQ-*.in) for individual components\n"""
#            instructionsO += """#    and a single input file (FREQ-master.in) with a frequency(mode='reap') command.\n"""
#            instructionsO += """#    These files may look very peculiar since they contain processed and pickled python\n"""
#            instructionsO += """#    rather than normal input. Follow the instructions below (repeated in FREQ-master.in)\n"""
#            instructionsO += """#    to continue.\n#\n"""
#            instructionsO += """#    Alternatively, a single-job execution of the hessian may be accessed through\n"""
#            instructionsO += """#    the frequency wrapper option mode='continuous'.\n#\n"""
#            core.print_out(instructionsO)
#
#            instructionsM = """\n#    Follow the instructions below to carry out this frequency computation.\n#\n"""
#            instructionsM += """#    (1)  Run all of the FREQ-*.in input files on any variety of computer architecture.\n"""
#            instructionsM += """#       The output file names must be as given below (these are the defaults when executed\n"""
#            instructionsM += """#       as `psi4 FREQ-1.in`, etc.).\n#\n"""
#            for rgt in range(ndisp):
#                pre = 'FREQ-' + str(rgt + 1)
#                instructionsM += """#             psi4 -i %-27s -o %-27s\n""" % (pre + '.in', pre + '.out')
#            instructionsM += """#\n#    (2)  Gather all the resulting output files in a directory. Place input file\n"""
#            instructionsM += """#         FREQ-master.in into that directory and run it. The job will be minimal in\n"""
#            instructionsM += """#         length and give summary results for the frequency computation in its output file.\n#\n"""
#            instructionsM += """#             psi4 -i %-27s -o %-27s\n#\n\n""" % ('FREQ-master.in', 'FREQ-master.out')
#
#            with open('FREQ-master.in', 'wb') as fmaster:
#                fmaster.write('# This is a psi4 input file auto-generated from the hessian() wrapper.\n\n'.encode('utf-8'))
#                fmaster.write(p4util.format_molecule_for_input(moleculeclone).encode('utf-8'))
#                fmaster.write(p4util.format_options_for_input(moleculeclone, **kwargs))
#                p4util.format_kwargs_for_input(fmaster, lmode=2, return_wfn=True, freq_dertype=1, **kwargs)
#                fmaster.write(("""retE, retwfn = %s('%s', **kwargs)\n\n""" % (frequency.__name__, lowername)).encode('utf-8'))
#                fmaster.write(instructionsM.encode('utf-8'))
#            core.print_out(instructionsM)
#
#        for n, displacement in enumerate(displacements):
#            rfile = 'FREQ-%s' % (n + 1)
#
#            # Build string of title banner
#            banners = ''
#            banners += """core.print_out('\\n')\n"""
#            banners += """p4util.banner(' Hessian Computation: Gradient Displacement %d ')\n""" % (n + 1)
#            banners += """core.print_out('\\n')\n\n"""
#
#            if freq_mode == 'continuous':
#
#                # print progress to file and screen
#                core.print_out('\n')
#                p4util.banner('Loading displacement %d of %d' % (n + 1, ndisp))
#                print(""" %d""" % (n + 1), end=('\n' if (n + 1 == ndisp) else ''))
#                sys.stdout.flush()
#
#                # Load in displacement into the active molecule (xyz coordinates only)
#                moleculeclone.set_geometry(displacement)
#
#                # Perform the gradient calculation
#                G, wfn = gradient(lowername, molecule=moleculeclone, return_wfn=True, **kwargs)
#                gradients.append(wfn.gradient())
#                energies.append(core.get_variable('CURRENT ENERGY'))
#
#                # clean may be necessary when changing irreps of displacements
#                core.clean()
#
#            # S/R: Write each displaced geometry to an input file
#            elif freq_mode == 'sow':
#                moleculeclone.set_geometry(displacement)
#
#                # S/R: Prepare molecule, options, kwargs, function call and energy save
#                #      forcexyz in molecule writer S/R enforcement of !reinterpret_coordentry above
#                with open('%s.in' % (rfile), 'wb') as freagent:
#                    freagent.write('# This is a psi4 input file auto-generated from the hessian() wrapper.\n\n')
#                    freagent.write(p4util.format_molecule_for_input(moleculeclone, forcexyz=True).encode('utf-8'))
#                    freagent.write(p4util.format_options_for_input(moleculeclone, **kwargs).encode('utf-8'))
#                    kwargs['return_wfn'] = True
#                    p4util.format_kwargs_for_input(freagent, **kwargs)
#                    freagent.write("""G, wfn = %s('%s', **kwargs)\n\n""" % (gradient.__name__, lowername))
#                    freagent.write("""core.print_out('\\nHESSIAN RESULT: computation %d for item %d """ % (os.getpid(), n + 1))
#                    freagent.write("""yields electronic gradient %r\\n' % (p4util.mat2arr(wfn.gradient())))\n\n""")
#                    freagent.write("""core.print_out('\\nHESSIAN RESULT: computation %d for item %d """ % (os.getpid(), n + 1))
#                    freagent.write("""yields electronic energy %20.12f\\n' % (get_variable('CURRENT ENERGY')))\n\n""")
#
#            # S/R: Read energy from each displaced geometry output file and save in energies array
#            elif freq_mode == 'reap':
#                exec(banners)
#                core.set_variable('NUCLEAR REPULSION ENERGY', moleculeclone.nuclear_repulsion_energy())
#                pygrad = p4util.extract_sowreap_from_output(rfile, 'HESSIAN', n, freq_linkage, True, label='electronic gradient')
#                p4mat = core.Matrix.from_list(pygrad)
#                p4mat.print_out()
#                gradients.append(p4mat)
#                energies.append(p4util.extract_sowreap_from_output(rfile, 'HESSIAN', n, freq_linkage, True))
#
#        # S/R: Quit sow after writing files. Initialize skeleton wfn to receive grad for reap
#        if freq_mode == 'sow':
#            optstash.restore()
#            optstash_conv.restore()
#            if return_wfn:
#                return (None, None)
#            else:
#                return None
#        elif freq_mode == 'reap':
#            wfn = core.Wavefunction.build(molecule, core.get_global_option('BASIS'))
#
#        # Assemble Hessian from gradients
#        #   Final disp is undisp, so wfn has mol, G, H general to freq calc
#        H = core.fd_freq_1(molecule, gradients, irrep)  # TODO or moleculeclone?
#        wfn.set_hessian(H)
#        wfn.set_gradient(G0)
#        wfn.set_frequencies(core.get_frequencies())
#
#        # The last item in the list is the reference energy, return it
#        core.set_variable('CURRENT ENERGY', energies[-1])
#
#        core.set_parent_symmetry('')
#        optstash.restore()
#        optstash_conv.restore()
#
#        if return_wfn:
#            return (wfn.hessian(), wfn)
#        else:
#            return wfn.hessian()

    else:
        raise FeatureNotImplemented("""hessian(dertype=0)""")
#        core.print_out("""hessian() will perform frequency computation by finite difference of analytic energies.\n""")
#
#        # Set method-dependent scf convergence criteria (test on procedures['energy'] since that's guaranteed)
#        optstash.restore()
#        optstash_conv.restore()
#        optstash_conv = driver_util._set_convergence_criterion('energy', lowername, 10, 11, 10, 11, 10)
#
#        # Shifting the geometry so need to copy the active molecule
#        moleculeclone = molecule.clone()
#
#        # Obtain list of displacements
#        displacements = core.fd_geoms_freq_0(moleculeclone, irrep)
#        moleculeclone.fix_orientation(True)
#        moleculeclone.reinterpret_coordentry(False)
#
#        # Record undisplaced symmetry for projection of diplaced point groups
#        core.set_parent_symmetry(molecule.schoenflies_symbol())
#
#        ndisp = len(displacements)
#
#        # This version is pretty dependent on the reference geometry being last (as it is now)
#        print(' %d displacements needed.' % ndisp)
#        energies = []
#
#        # S/R: Write instructions for sow/reap procedure to output file and reap input file
#        if freq_mode == 'sow':
#            instructionsO = """\n#    The frequency sow/reap procedure has been selected through mode='sow'. In addition\n"""
#            instructionsO += """#    to this output file (which contains no quantum chemical calculations), this job\n"""
#            instructionsO += """#    has produced a number of input files (FREQ-*.in) for individual components\n"""
#            instructionsO += """#    and a single input file (FREQ-master.in) with a frequency(mode='reap') command.\n"""
#            instructionsO += """#    These files may look very peculiar since they contain processed and pickled python\n"""
#            instructionsO += """#    rather than normal input. Follow the instructions below (repeated in FREQ-master.in)\n"""
#            instructionsO += """#    to continue.\n#\n"""
#            instructionsO += """#    Alternatively, a single-job execution of the hessian may be accessed through\n"""
#            instructionsO += """#    the frequency wrapper option mode='continuous'.\n#\n"""
#            core.print_out(instructionsO)
#
#            instructionsM = """\n#    Follow the instructions below to carry out this frequency computation.\n#\n"""
#            instructionsM += """#    (1)  Run all of the FREQ-*.in input files on any variety of computer architecture.\n"""
#            instructionsM += """#       The output file names must be as given below (these are the defaults when executed\n"""
#            instructionsM += """#       as `psi4 FREQ-1.in`, etc.).\n#\n"""
#            for rgt in range(ndisp):
#                pre = 'FREQ-' + str(rgt + 1)
#                instructionsM += """#             psi4 -i %-27s -o %-27s\n""" % (pre + '.in', pre + '.out')
#            instructionsM += """#\n#    (2)  Gather all the resulting output files in a directory. Place input file\n"""
#            instructionsM += """#         FREQ-master.in into that directory and run it. The job will be minimal in\n"""
#            instructionsM += """#         length and give summary results for the frequency computation in its output file.\n#\n"""
#            instructionsM += """#             psi4 -i %-27s -o %-27s\n#\n\n""" % ('FREQ-master.in', 'FREQ-master.out')
#
#            with open('FREQ-master.in', 'wb') as fmaster:
#                fmaster.write('# This is a psi4 input file auto-generated from the hessian() wrapper.\n\n'.encode('utf-8'))
#                fmaster.write(p4util.format_molecule_for_input(moleculeclone).encode('utf-8'))
#                fmaster.write(p4util.format_options_for_input(moleculeclone, **kwargs))
#                p4util.format_kwargs_for_input(fmaster, lmode=2, return_wfn=True, freq_dertype=0, **kwargs)
#                fmaster.write(("""retE, retwfn = %s('%s', **kwargs)\n\n""" % (frequency.__name__, lowername)).encode('utf-8'))
#                fmaster.write(instructionsM.encode('utf-8'))
#            core.print_out(instructionsM)
#
#        for n, displacement in enumerate(displacements):
#            rfile = 'FREQ-%s' % (n + 1)
#
#            # Build string of title banner
#            banners = ''
#            banners += """core.print_out('\\n')\n"""
#            banners += """p4util.banner(' Hessian Computation: Energy Displacement %d ')\n""" % (n + 1)
#            banners += """core.print_out('\\n')\n\n"""
#
#            if freq_mode == 'continuous':
#
#                # print progress to file and screen
#                core.print_out('\n')
#                p4util.banner('Loading displacement %d of %d' % (n + 1, ndisp))
#                print(""" %d""" % (n + 1), end=('\n' if (n + 1 == ndisp) else ''))
#                sys.stdout.flush()
#
#                # Load in displacement into the active molecule
#                moleculeclone.set_geometry(displacement)
#
#                # Perform the energy calculation
#                E, wfn = energy(lowername, return_wfn=True, molecule=moleculeclone, **kwargs)
#                energies.append(core.get_variable('CURRENT ENERGY'))
#
#                # clean may be necessary when changing irreps of displacements
#                core.clean()
#
#            # S/R: Write each displaced geometry to an input file
#            elif freq_mode == 'sow':
#                moleculeclone.set_geometry(displacement)
#
#                # S/R: Prepare molecule, options, kwargs, function call and energy save
#                with open('%s.in' % (rfile), 'wb') as freagent:
#                    freagent.write('# This is a psi4 input file auto-generated from the gradient() wrapper.\n\n')
#                    freagent.write(p4util.format_molecule_for_input(moleculeclone, forcexyz=True).encode('utf-8'))
#                    freagent.write(p4util.format_options_for_input(moleculeclone, **kwargs).encode('utf-8'))
#                    p4util.format_kwargs_for_input(freagent, **kwargs)
#                    freagent.write("""electronic_energy = %s('%s', **kwargs)\n\n""" % (energy.__name__, lowername))
#                    freagent.write("""core.print_out('\\nHESSIAN RESULT: computation %d for item %d """ % (os.getpid(), n + 1))
#                    freagent.write("""yields electronic energy %20.12f\\n' % (electronic_energy))\n\n""")
#
#            # S/R: Read energy from each displaced geometry output file and save in energies array
#            elif freq_mode == 'reap':
#                exec(banners)
#                core.set_variable('NUCLEAR REPULSION ENERGY', moleculeclone.nuclear_repulsion_energy())
#                energies.append(p4util.extract_sowreap_from_output(rfile, 'HESSIAN', n, freq_linkage, True))
#
#        # S/R: Quit sow after writing files. Initialize skeleton wfn to receive grad for reap
#        if freq_mode == 'sow':
#            optstash.restore()
#            optstash_conv.restore()
#            if return_wfn:
#                return (None, None)
#            else:
#                return None
#        elif freq_mode == 'reap':
#        #    core.set_variable('CURRENT ENERGY', energies[-1])
#            wfn = core.Wavefunction.build(molecule, core.get_global_option('BASIS'))
#
#        # Assemble Hessian from energies
#        H = core.fd_freq_0(molecule, energies, irrep)
#        wfn.set_hessian(H)
#        wfn.set_gradient(G0)
#        wfn.set_frequencies(core.get_frequencies())
#
#        # The last item in the list is the reference energy, return it
#        core.set_variable('CURRENT ENERGY', energies[-1])
#
#        core.set_parent_symmetry('')
#        optstash.restore()
#        optstash_conv.restore()
#
#        if return_wfn:
#            return (wfn.hessian(), wfn)
#        else:
#            return wfn.hessian()


#def frequency(name, **kwargs):
#    r"""Function to compute harmonic vibrational frequencies.
#
#    :aliases: frequencies(), freq()
#
#    :returns: *float* |w--w| Total electronic energy in Hartrees.
#
#    :returns: (*float*, :py:class:`~psi4.core.Wavefunction`) |w--w| energy and wavefunction when **return_wfn** specified.
#
#    :type name: string
#    :param name: ``'scf'`` || ``'mp2'`` || ``'ci5'`` || etc.
#
#        First argument, usually unlabeled. Indicates the computational method
#        to be applied to the system.
#
#    :type molecule: :ref:`molecule <op_py_molecule>`
#    :param molecule: ``h2o`` || etc.
#
#        The target molecule, if not the last molecule defined.
#
#    :type return_wfn: :ref:`boolean <op_py_boolean>`
#    :param return_wfn: ``'on'`` || |dl| ``'off'`` |dr|
#
#        Indicate to additionally return the :py:class:`~psi4.core.Wavefunction`
#        calculation result as the second element (after *float* energy) of a tuple.
#        Arrays of frequencies and the Hessian can be accessed through the wavefunction.
#
#    :type func: :ref:`function <op_py_function>`
#    :param func: |dl| ``gradient`` |dr| || ``energy`` || ``cbs``
#
#        Indicates the type of calculation to be performed on the molecule.
#        The default dertype accesses ``'gradient'`` or ``'energy'``, while
#        ``'cbs'`` performs a multistage finite difference calculation.
#        If a nested series of python functions is intended (see :ref:`sec:intercalls`),
#        use keyword ``freq_func`` instead of ``func``.
#
#    :type mode: string
#    :param mode: |dl| ``'continuous'`` |dr| || ``'sow'`` || ``'reap'``
#
#        For a finite difference of energies or gradients frequency, indicates
#        whether the calculations required to complete the frequency are to be run
#        in one file (``'continuous'``) or are to be farmed out in an
#        embarrassingly parallel fashion (``'sow'``/``'reap'``)/ For the latter,
#        run an initial job with ``'sow'`` and follow instructions in its output file.
#        For maximum flexibility, ``return_wfn`` is always on in ``'reap'`` mode.
#
#    :type dertype: :ref:`dertype <op_py_dertype>`
#    :param dertype: |dl| ``'hessian'`` |dr| || ``'gradient'`` || ``'energy'``
#
#        Indicates whether analytic (if available- they're not), finite
#        difference of gradients (if available) or finite difference of
#        energies is to be performed.
#
#    :type irrep: int or string
#    :param irrep: |dl| ``-1`` |dr| || ``1`` || ``'b2'`` || ``'App'`` || etc.
#
#        Indicates which symmetry block (:ref:`Cotton <table:irrepOrdering>` ordering) of vibrational
#        frequencies to be computed. ``1``, ``'1'``, or ``'a1'`` represents
#        :math:`a_1`, requesting only the totally symmetric modes.
#        ``-1`` indicates a full frequency calculation.
#
#    .. note:: Analytic hessians are only available for RHF. For all other methods, Frequencies will
#        proceed through finite differences according to availability of gradients or energies.
#
#    .. _`table:freq_gen`:
#
#    +-------------------------+---------------------------------------------------------------------------------------------------------------+
#    | name                    | calls method                                                                                                  |
#    +=========================+===============================================================================================================+
#    | scf                     | Hartree--Fock (HF) :ref:`[manual] <sec:scf>`                                                                  |
#    +-------------------------+---------------------------------------------------------------------------------------------------------------+
#
#    :examples:
#
#    >>> # [1] Frequency calculation for all modes through highest available derivatives
#    >>> frequency('ccsd')
#
#    >>> # [2] Frequency calculation for b2 modes through finite difference of gradients
#    >>> #     printing lowest mode frequency to screen and Hessian to output
#    >>> E, wfn = frequencies('scf', dertype=1, irrep=4, return_wfn=True)
#    >>> print wfn.frequencies().get(0, 0)
#    >>> wfn.hessian().print_out()
#
#    >>> # [3] Frequency calculation at default conditions and Hessian reuse at STP
#    >>> E, wfn = freq('mp2', return_wfn=True)
#    >>> set t 273.15
#    >>> set p 100000
#    >>> thermo(wfn, wfn.frequencies())
#
#    """
#    kwargs = p4util.kwargs_lower(kwargs)
#
#    # Bounce (someday) if name is function
#    if hasattr(name, '__call__'):
#        raise ValidationError("Frequency: Cannot use custom function")
#
#    lowername = name.lower()
#
#    old_global_basis = None
#    if "/" in lowername:
#        if ("+" in lowername) or ("[" in lowername) or (lowername.count('/') > 1):
#            raise ValidationError("Frequency: Cannot extrapolate or delta correct frequencies yet.")
#        else:
#            old_global_basis = core.get_global_option("BASIS")
#            lowername, new_basis = lowername.split('/')
#            core.set_global_option('BASIS', new_basis)
#
#    if kwargs.get('bsse_type', None) is not None:
#        raise ValdiationError("Frequency: Does not currently support 'bsse_type' arguements")
#
#    return_wfn = kwargs.pop('return_wfn', False)
#
#    # are we in sow/reap mode?
#    freq_mode = kwargs.get('mode', 'continuous').lower()
#    if freq_mode not in ['continuous', 'sow', 'reap']:
#        raise ValidationError("""Frequency execution mode '%s' not valid.""" % (freq_mode))
#
#    # Make sure the molecule the user provided is the active one
#    molecule = kwargs.pop('molecule', core.get_active_molecule())
#    molecule.update_geometry()
#
#    # Compute the hessian
#    H, wfn = hessian(lowername, return_wfn=True, molecule=molecule, **kwargs)
#
#    # S/R: Quit after getting new displacements
#    if freq_mode == 'sow':
#        return 0.0
#
#    # Project final frequencies?
#    translations_projection_sound, rotations_projection_sound = _energy_is_invariant(wfn.gradient())
#    project_trans = kwargs.get('project_trans', translations_projection_sound)
#    project_rot = kwargs.get('project_rot', rotations_projection_sound)
#
#    irrep = kwargs.get('irrep', None)
#    vibinfo = vibanal_wfn(wfn, irrep=irrep, project_trans=project_trans, project_rot=project_rot)
#    vibonly = qcdb.vib.filter_nonvib(vibinfo)
#    wfn.set_frequencies(core.Vector.from_array(qcdb.vib.filter_omega_to_real(vibonly['omega'].data)))
#    wfn.frequency_analysis = vibinfo
#
#    for postcallback in hooks['frequency']['post']:
#        postcallback(lowername, wfn=wfn, **kwargs)
#
#    # Reset old global basis if needed
#    if not old_global_basis is None:
#        core.set_global_option("BASIS", old_global_basis)
#
#    if return_wfn:
#        return (core.get_variable('CURRENT ENERGY'), wfn)
#    else:
#        return core.get_variable('CURRENT ENERGY')
#
#
#def vibanal_wfn(wfn, hess=None, irrep=None, molecule=None, project_trans=True, project_rot=True):
#    # TODO should go back to private
#
#    if hess is None:
#        nmwhess = np.asarray(wfn.hessian())
#    else:
#        nmwhess = hess
#
#    mol = wfn.molecule()
#    geom = np.asarray(mol.geometry())
#    symbols = [mol.symbol(at) for at in range(mol.natom())]
#
#    if molecule is not None:
#        molecule.update_geometry()
#        if mol.natom() != molecule.natom():
#            raise ValidationError('Impostor molecule trying to be analyzed! natom {} != {}'.format(mol.natom(), molecule.natom()))
#        if abs(mol.nuclear_repulsion_energy() - molecule.nuclear_repulsion_energy()) > 1.e-6:
#            raise ValidationError('Impostor molecule trying to be analyzed! NRE {} != {}'.format(mol.nuclear_repulsion_energy(), molecule.nuclear_repulsion_energy()))
#        if not np.allclose(np.asarray(mol.geometry()), np.asarray(molecule.geometry()), atol=1.e-6):
#            core.print_out('Warning: geometry center/orientation mismatch. Normal modes may not be in expected coordinate system.')
#        #    raise ValidationError('Impostor molecule trying to be analyzed! geometry\n{}\n   !=\n{}'.format(
#        #        np.asarray(mol.geometry()), np.asarray(molecule.geometry())))
#        mol = molecule
#            
#    m = np.asarray([mol.mass(at) for at in range(mol.natom())])
#    irrep_labels = mol.irrep_labels()
#
#    vibinfo, vibtext = qcdb.vib.harmonic_analysis(nmwhess, geom, m, wfn.basisset(), irrep_labels,
#                                                  project_trans=project_trans, project_rot=project_rot)
#
#    core.print_out(vibtext)
#    core.print_out(qcdb.vib.print_vibs(vibinfo, shortlong=True, normco='x', atom_lbl=symbols))
#
#    if core.has_option_changed('THERMO', 'ROTATIONAL_SYMMETRY_NUMBER'):
#        rsn = core.get_option('THERMO', 'ROTATIONAL_SYMMETRY_NUMBER')
#    else:
#        rsn = mol.rotational_symmetry_number()
#
#    if irrep is None:
#        therminfo, thermtext = qcdb.vib.thermo(vibinfo,
#                                      T=core.get_option("THERMO", "T"),  # 298.15
#                                      P=core.get_option("THERMO", "P"),  # 101325.
#                                      multiplicity=mol.multiplicity(),
#                                      molecular_mass=np.sum(m),
#                                      sigma=rsn,
#                                      rotor_type=mol.rotor_type(),
#                                      rot_const=np.asarray(mol.rotational_constants()),
#                                      E0=core.get_variable('CURRENT ENERGY'))  # someday, wfn.energy()
#
#        core.set_variable("ZPVE", therminfo['ZPE_corr'].data)
#        core.set_variable("THERMAL ENERGY CORRECTION", therminfo['E_corr'].data)
#        core.set_variable("ENTHALPY CORRECTION", therminfo['H_corr'].data)
#        core.set_variable("GIBBS FREE ENERGY CORRECTION", therminfo['G_corr'].data)
#
#        core.set_variable("ZERO K ENTHALPHY", therminfo['ZPE_tot'].data)
#        core.set_variable("THERMAL ENERGY", therminfo['E_tot'].data)
#        core.set_variable("ENTHALPY", therminfo['H_tot'].data)
#        core.set_variable("GIBBS FREE ENERGY", therminfo['G_tot'].data)
#
#        core.print_out(thermtext)
#    else:
#        core.print_out('  Thermochemical analysis skipped for partial frequency calculation.\n')
#
#    return vibinfo


## Aliases
#opt = optimize
#freq = frequency
#frequencies = frequency
#prop = properties
