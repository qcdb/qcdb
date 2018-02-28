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
"""Module with functions that interface with Grimme's DFTD3 code."""
from __future__ import absolute_import
from __future__ import print_function
import sys
import collections
#import os
#import re
#
#try:
#    from psi4.driver.p4util.exceptions import *
#    from psi4 import core
#    isP4regime = True
#except ImportError:
#    from .exceptions import *
#    isP4regime = False
#from .dashparam import *
#from .molecule import Molecule

import numpy as np

from ..util import update_with_error, der0th, der1st
from ..exceptions import *
#from .dashparam import dash_server, dashcoeff
from .. import dashparam
from .. import molparse
from .. import __version__

import pprint


"""Compute dispersion correction using Grimme's DFTD3 executable.

Function to call Grimme's dftd3 program to compute the -D correction
of level `dashlvl` using parameters for the functional `func`.
`dashparam` can supply a full set of dispersion parameters in the
absence of `func` or individual overrides in the presence of `func`.

The DFTD3 executable must be independently compiled and found in
:envvar:`PATH` or :envvar:`PSIPATH`.

Parameters
----------
mol : qcdb.Molecule or psi4.core.Molecule or str
    Molecule on which to run dispersion calculation. Both qcdb and
    psi4.core Molecule classes have been extended by this method, so
    either allowed. Alternately, a string that can be instantiated
    into a qcdb.Molecule.
func : str or None
    Density functional (Psi4, not Turbomole, names) for which to
    load parameters from dashcoeff[dashlvl][func]. This is not
    passed to DFTD3 and thus may be a dummy or `None`. Any or all
    parameters initialized can be overwritten via `dashparam`.
dashlvl : {'d2p4', 'd2gr', 'd3zero', 'd3bj', 'd3mzero', d3mbj', 'd', 'd2', 'd3', 'd3m'}
    Flavor of a posteriori dispersion correction for which to load
    parameters and call procedure in DFTD3. Must be a keys in
    dashcoeff dict (or a key in dashalias that resolves to one).
dashparam : dict, optional
    Dictionary of the same keys as dashcoeff[dashlvl] used to
    override any or all values initialized by
    dashcoeff[dashlvl][func].
dertype : {None, 0, 'none', 'energy', 1, 'first', 'gradient'}, optional
    Maximum derivative level at which to run DFTD3. For large
    `mol`, energy-only calculations can be significantly more
    efficient. Also controls return values, see below.
verbose : bool, optional
    When `True`, additionally include DFTD3 output in output.

Returns
-------
energy : float, optional
    When `dertype` is 0, energy [Eh].
gradient : list of lists of floats or psi4.core.Matrix, optional
    When `dertype` is 1, (nat, 3) gradient [Eh/a0].
(energy, gradient) : float and list of lists of floats or psi4.core.Matrix, optional
    When `dertype` is unspecified, both energy [Eh] and (nat, 3) gradient [Eh/a0].

Notes
-----
research site: https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3
Psi4 mode: When `psi4` the python module is importable at `import qcdb`
           time, Psi4 mode is activated, with the following alterations:
           * output goes to output file
           * gradient returned as psi4.core.Matrix, not list o'lists
           * scratch is written to randomly named subdirectory of psi scratch
           * psivar "DISPERSION CORRECTION ENERGY" is set
           * `verbose` triggered when PRINT keywork of SCF module >=3

"""


def run_dftd3(molrec,
              functional=None,
              dashlevel=None,
              dashparams=None,
              dertype=None,
              verbose=1):
    """

    Required Input Fields
    ---------------------
    dashlevel
    dashparams
    functional
    molecule
    do_gradient

    molecule
    provenance
        creator
        version
        routine

    Optional Input Fields
    ---------------------

    Output Fields
    -------------
    error
    success
    raw_output


    """
    jobrec = from_arrays_qc(
        functional=functional,
        dashlevel=dashlevel,
        dashparams=dashparams,
        dertype=dertype)
    jobrec['molecule'] = molrec

    jobrec['error'] = ''
    jobrec['success'] = False
    jobrec['raw_output'] = None

    prov = {}
    prov['creator'] = 'QCDB'
    prov['version'] = __version__
    prov['routine'] = sys._getframe().f_code.co_name
    jobrec['provenance'] = [prov]

    try:
        dftd3_driver(jobrec)
        jobrec['success'] = True
    except Exception as err:
        jobrec['success'] = False
        #json_data["error"] += repr(error)

    return jobrec


def from_arrays_qc(functional=None,
                   dashlevel=None,
                   dashparams=None,
                   dertype=None,
                   verbose=1):

    jobrec = {}

    processed = _validate_and_fill_dertype(dertype=dertype)
    update_with_error(jobrec, processed)

    processed = _validate_and_fill_dashparam(
        functional=functional, dashlevel=dashlevel, dashparams=dashparams)
    update_with_error(jobrec, processed)

    return jobrec


def dftd3_driver(jobrec):
    import json

    try:
        jobrec['dashlevel']
        jobrec['dashparams']
        jobrec['functional']
        jobrec['molecule']
        jobrec['do_gradient']
    except KeyError as err:
        raise KeyError(
            'Required fields missing from ({})'.format(jobrec.keys())) from err

    # this is what the dftd3 program needs, not what the job needs
    # * form dftd3_parameters string that governs dispersion calc
    # * form dftd3_geometry string that supplies geometry to dispersion calc
    # * form command and arguments
    dftd3rec = {}

    dftd3par = dashparam.dftd3_coeff_formatter(jobrec['dashlevel'],
                                               jobrec['dashparams'])
    dftd3rec['dftd3par'] = dftd3par

    dftd3_geometry = molparse.to_string(
        jobrec['molecule'], dtype='xyz', units='Angstrom', ghost_format='')
    dftd3rec['dftd3_geometry'] = dftd3_geometry

    command = ['dftd3', './dftd3_geometry.xyz']
    if jobrec['do_gradient'] is True:
        command.append('-grad')
    dftd3rec['command'] = command

    print('IN')
    #    pprint.pprint(dftd3rec)

    # test json roundtrip
    jdftd3rec = json.dumps(dftd3rec)
    dftd3rec = json.loads(jdftd3rec)

    subprocess_dftd3(dftd3rec)  # updates dftd3rec

    print('OUT')
    #    pprint.pprint(dftd3rec)

    #if maxder == 1:
    #require 'dftd3grad'
    #    derivfile = './dftd3_gradient'

    #for reqd in ['stdout'
    #try:
    jobrec = dftd3_harvest(jobrec, dftd3rec)

    print('CALC')
    pprint.pprint(jobrec)

    ## Print program output to file if verbose
    #if not verbose and isP4regime:
    #    verbose = True if core.get_option('SCF', 'PRINT') >= 3 else False
    #if verbose:

    #    text = '\n  ==> DFTD3 Output <==\n'
    #    text += out.decode('utf-8')
    #    if dertype != 0:
    #        with open(derivfile, 'r') as handle:
    #            text += handle.read().replace('D', 'E')
    #        text += '\n'
    #    if isP4regime:
    #        core.print_out(text)
    #    else:
    #        print(text)

    ## return -D & d(-D)/dx
    #if dertype == -1:
    #    return dashd, dashdderiv
    #elif dertype == 0:
    #    return dashd
    #elif dertype == 1:
    #    return psi_dashdderiv

    #for reqd in ['molecule'

    # Find out if running from Psi4 for scratch details and such
    # try:
    #     import psi4
    # except ImportError as err:
    #     isP4regime = False
    # else:
    #     isP4regime = True

    #if isP4regime:
    #    psioh = core.IOManager.shared_object()
    #    psio = core.IO.shared_object()
    #    os.chdir(psioh.get_default_path())
    #    dftd3_tmpdir = 'psi.' + str(os.getpid()) + '.' + psio.get_default_namespace() + \
    #        '.dftd3.' + str(uuid.uuid4())[:8]


#    # We seem to have a problem with one atom, force the correct result
#    if numAtoms == 1:
#        dashd = 0.0
#        dashdderiv = core.Matrix(1, 3)

#    if dertype == -1:
#        return dashd, dashdderiv
#    elif dertype == 0:
#        return dashd
#    elif dertype == 1:
#        return dashdderiv

#if args["json"]:

#def cli_dftd3(
#    with open(args["input"], 'r') as f:
#        json_data = json.load(f)
#
#    psi4.extras._success_flag_ = True
#    psi4.extras.exit_printing()
#    psi4.json_wrapper.run_json(json_data)
#
#    with open(args["input"], 'w') as f:
#        json.dump(json_data, f)
#
#    if args["output"] != "stdout":
#        os.unlink(args["output"])
#
#    sys.exit()


def subprocess_dftd3(dftd3rec):
    """Minimal localized DFTD3 call and harvest from and into `dftd3rec`.

    Required Input Fields
    ---------------------
    command : list
        Command and arguments to execute. 
    dftd3par : str
        Parameter file contents defining job run.
    dftd3_geometry
        XYZ of real atoms file contents defining job geometry.

    Optional Input Fields
    ---------------------
    scratch_location : str, optional
        Override the default scratch location.
        Note that this is PARENT dir.
    scratch_messy : bool, optional
        If present and `True`, scratch is left behind.

    Output Fields
    -------------
    stdout : str
        Main output file that gets written to stdout.
    dftd3_gradient : str, optional
        If `-grad` present in `command`, contents of file with real atom gradient.

    """
    import os
    import uuid
    import shutil
    import socket
    import subprocess

    try:
        dftd3rec['command']
        dftd3rec['dftd3par']
        dftd3rec['dftd3_geometry']
    except KeyError as err:
        raise KeyError('Required fields missing from ({})'.format(
            dftd3rec.keys())) from err

    current_directory = os.getcwd()

    # move ~/.dftd3par.<hostname> out of the way so it won't interfere
    defaultfile = os.path.expanduser(
        '~') + '/.dftd3par.' + socket.gethostname()
    defmoved = False
    if os.path.isfile(defaultfile):
        os.rename(defaultfile, defaultfile + '_hide')
        defmoved = True

    # find environment by merging PSIPATH and PATH environment variables
    # * filter out None values as subprocess will fault on them
    lenv = {
        'PATH': ':'.join([os.path.abspath(x) for x in os.environ.get('PSIPATH', '').split(':') if x != '']) + \
                ':' + os.environ.get('PATH'),
        'LD_LIBRARY_PATH': os.environ.get('LD_LIBRARY_PATH')
        }
    lenv = {k: v for k, v in lenv.items() if v is not None}

    # set up unique scratch directory and move in
    if 'scratch_location' in dftd3rec:
        basedir = dftd3rec['scratch_location']
    else:
        basedir = os.environ['HOME'] + os.sep
    dftd3_tmpdir = basedir + 'dftd3_' + str(uuid.uuid4())[:8]
    if not os.path.exists(dftd3_tmpdir):
        os.mkdir(dftd3_tmpdir)
    os.chdir(dftd3_tmpdir)

    # write governing inputs
    paramfileold = 'dftd3_parameters'  # older patched name
    with open(paramfileold, 'w') as handle:
        handle.write(dftd3rec['dftd3par'])
    paramfile = '.dftd3par.local'  # new mainline name
    with open(paramfile, 'w') as handle:
        handle.write(dftd3rec['dftd3par'])
    geomfile = './dftd3_geometry.xyz'
    with open(geomfile, 'w') as handle:
        handle.write(dftd3rec['dftd3_geometry'])

    # call `dftd3` program
    try:
        spcall = subprocess.Popen(
            dftd3rec['command'], stdout=subprocess.PIPE, env=lenv)
    except OSError as err:
        raise OSError('Command (`{}`) failed with PATH ({})'.format(
            ' '.join(command), lenv['PATH'])) from err

    # recover output data
    out, err = spcall.communicate()
    dftd3rec['stdout'] = out.decode('utf-8')

    if '-grad' in dftd3rec['command']:
        derivfile = './dftd3_gradient'
        with open(derivfile, 'r') as handle:
            dftd3rec['dftd3_gradient'] = handle.read()

    # clean up files and remove scratch directory
    if 'scratch_messy' not in dftd3rec or dftd3rec['scratch_messy'] is False:
        os.unlink(paramfileold)
        os.unlink(paramfile)
        os.unlink(geomfile)
        if '-grad' in dftd3rec['command']:
            os.unlink(derivfile)

        os.chdir('..')
        try:
            shutil.rmtree(dftd3_tmpdir)
        except OSError as err:
            raise OSError('Unable to remove dftd3 temporary directory: {}'.
                          format(dftd3_tmpdir)) from err

    if defmoved is True:
        os.rename(defaultfile + '_hide', defaultfile)

    os.chdir(current_directory)

    return dftd3rec


def dftd3_harvest(jobrec, dftd3rec):
    """Processes raw results from read-only `dftd3rec` into QCAspect fields in returned `jobrec`."""

    try:
        jobrec['molecule']['real']
        jobrec['do_gradient']
    except KeyError as err:
        raise KeyError(
            'Required fields missing from ({})'.format(jobrec.keys())) from err

    try:
        dftd3rec['stdout']
        if jobrec['do_gradient'] is True:
            dftd3rec['dftd3_gradient']
    except KeyError as err:
        raise KeyError('Required fields missing from ({})'.format(
            dftd3rec.keys())) from err

    # parse energy output (could go further and break into E6, E8, E10 and Cn coeff)
    for ln in dftd3rec['stdout'].splitlines():
        if re.search('DFTD3 V', ln):
            version = ln.replace('DFTD3', '').replace('|', '').strip().lower()
        if re.match(' Edisp /kcal,au', ln):
            ene = float(ln.split()[3])
        if re.match(' normal termination of dftd3', ln):
            break
    else:
        raise Dftd3Error('Unsuccessful run. Possibly -D variant not available in dftd3 version.')

    # parse gradient output
    if 'dftd3_gradient' in dftd3rec:
        real = jobrec['molecule']['real']
        fnat = real.shape[0]
        rnat = np.sum(real)
        srealgrad = dftd3rec['dftd3_gradient'].replace('D', 'E')
        realgrad = np.fromstring(
            srealgrad, count=3 * rnat, sep=' ').reshape((-1, 3))

        ireal = np.argwhere(real).reshape((-1))
        fullgrad = np.zeros((fnat, 3))
        fullgrad[ireal, :] = realgrad

        # TODO if dftd3rec['do_gradient']: raise Dftd3Error

    #QCAspect = collections.namedtuple('QCAspect', 'lbl unit data comment')
    calcinfo = []
    calcinfo.append(QCAspect('DISPERSION CORRECTION ENERGY', '[Eh]', ene, ''))
    calcinfo.append(
        QCAspect('DISPERSION CORRECTION GRADIENT', '[Eh/a0]', fullgrad, ''))
    calcinfo = {info.lbl: info for info in calcinfo}

    # amalgamate output
    text = dftd3rec['stdout']
    text += '\n  <<<  DFTD3 {} {} Results  >>>'.format('', '') #name.lower(), calledby.capitalize()))  # banner()
    text += qcdb.print_variables(calcinfo)
    jobrec['raw_output'] = text

    jobrec['qcvars'] = calcinfo
    prov = {}
    prov['creator'] = 'DFTD3'
    prov['routine'] = sys._getframe().f_code.co_name,
    prov['version'] = version,
    jobrec['provenance'].append(prov)

    return jobrec

    # Prepare results for Psi4
    if isP4regime and dertype != 0:
        print('P4-2?', isP4regime)
        core.set_variable('DISPERSION CORRECTION ENERGY', dashd)
        psi_dashdderiv = core.Matrix.from_list(dashdderiv)


def _validate_and_fill_dertype(dertype=None):
    """Issue `do_gradient=True` unless `dertype` is energy or nonsense."""

    if dertype is None:
        do_gradient = True
    elif der0th.match(str(dertype)):
        do_gradient = False
    elif der1st.match(str(dertype)):
        do_gradient = True
    else:
        raise ValidationError(
            """Requested derivative level ({}) not valid for run_dftd3.""".
            format(dertype))

    return ({'do_gradient': do_gradient})


def _validate_and_fill_dashparam(dashlevel, functional=None, dashparams=None):
    """Take the three paths of empirical dispersion parameter information
    (DFT functional, dispersion correction level, and particular parameters)
    and populate the parameter array.

    """
    if functional is None and dashparams is None:
        raise ValidationError(
            """Can't guess -D parameters without functional ({}) or dashparams ({})""".
            format(functional, dashparams))

    if functional is None:
        dftd3_params = {}
    else:
        functionaleff, dashleveleff, dftd3_params = dashparam.dash_server(
            functional, dashlevel, return_triplet=True)

    if dashparams is None:
        dashparams = {}
    else:
        dashleveleff, allowed_params = dashparam.dash_server(
            func=None, dashlvl=dashlevel, return_levelkeys=True)
        if not set(dashparams.keys()).issubset(allowed_params):
            raise ValidationError(
                'Requested keys ({}) not among allowed ({}) for dispersion level ({})'.
                format(dashparams.keys(), allowed_params, dashlevel))

    dftd3_params.update(dashparams)

    for func, params in dashparam.dashcoeff[dashleveleff].items():
        if dftd3_params == params:
            functionaleff = func
            break
    else:
        functionaleff = ''

    return {
        'dashlevel': dashleveleff,
        'dashparams': dftd3_params,
        'functional': functionaleff
    }


#fullgrad = np.insert(realgrad, ghosts, 0., axis=0)

#    # Create (if necessary) and update qcdb.Molecule
#    try:
#        mol.update_geometry()
#    except AttributeError as err:
#        mol = Molecule(mol)
#        mol.update_geometry()

#    # Create (if necessary) and update qcdb.Molecule
#    if isinstance(mol, (Molecule, core.Molecule)):
#        # 1st: called on a qcdb.Molecule
#        # 2nd: called on a python export of a psi4.Molecule (py-side through Psi4's driver)
#        pass
#    elif isinstance(mol, basestring):
#        # called on a string representation of a psi4.Molecule (c-side through psi4.Dispersion)
#        mol = Molecule(mol)
#    else:
#        raise ValidationError("""Argument mol must be psi4string or qcdb.Molecule""")
#    mol.update_geometry()

#try:
#    # Attach method to libmints psi4.Molecule class
#    core.Molecule.run_dftd3 = run_dftd3
#except (NameError, AttributeError):
#    # But don't worry if that doesn't work b/c
#    #   it'll get attached to qcdb.Molecule class
#    pass
