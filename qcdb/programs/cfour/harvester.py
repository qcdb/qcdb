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

import uuid
import struct
from collections import defaultdict

import numpy as np
import qcelemental as qcel
from qcelemental.models import Molecule

from ...exceptions import ValidationError
#from ..molecule import Molecule
from ...util import conv_float2negexp

#def nu_muster_memory(mem, ropts, verbose=1):
#    """Transform input *mem* in MB into psi4-type options for cfour.
#
#    """
#    text = ''
#    accession = 4567
#
#    # prepare memory keywords to be set as c-side keywords
#    #options['CFOUR']['CFOUR_MEMORY_SIZE']['value'] = int(mem)
#    #options['CFOUR']['CFOUR_MEM_UNIT']['value'] = 'MB'
#    ropts.require('CFOUR', 'MEMORY_SIZE', int(mem), accession=accession, verbose=verbose)
#    ropts.require('CFOUR', 'MEM_UNIT', 'MB', accession=accession, verbose=verbose)
#
#    #for item in options['CFOUR']:
#    #    options['CFOUR'][item]['clobber'] = True
#    #return text, options
#
#    return ''
#
#def muster_memory(mem):
#    """Transform input *mem* in MB into psi4-type options for cfour.
#
#    """
#    text = ''
#
#    # prepare memory keywords to be set as c-side keywords
#    options = defaultdict(lambda: defaultdict(dict))
#    options['CFOUR']['CFOUR_MEMORY_SIZE']['value'] = int(mem)
#    options['CFOUR']['CFOUR_MEM_UNIT']['value'] = 'MB'
#
#    for item in options['CFOUR']:
#        options['CFOUR'][item]['clobber'] = True
#    return text, options

#   Ways of modifying a computation
#   global:     set global c-side option
#   local:      set local c-side option
#   kwarg:      set kwarg
#   i-local:    set global=local c-side option to an interface module
#   ro-def:     code uses default entirely specified by read_options
#   module-def: code uses default that is complex mixture of read_options settings
#   i-def:      interfaced code uses defaults not entirely expressed in read_options
#   driver-def: driver code sets complex defaults
#
#   Pure psi4 operation
#   kwarg ~= local > global > driver-def > module-def > ro-def
#
#   Interfaced psi4 operation
#   kwarg ~= i-local > local > global > driver-def > i-def

#   P4 infrastructure replacing interfaced infrastructure (mol, basis, mem) where unavoidable overlap in how things are specified (mult in mol{} vs keyword) is treated as a clobber & complain if conflict VS P4 infrastructure as an aliased/convenient leak into interfaced infrastructure (psi) and is strictly no clobber or complain.


def muster_inherited_options(ropts, verbose=1):
    import sys
    accession = sys._getframe().f_code.co_name + '_' + str(uuid.uuid4())[:8]
    kwgs = {'accession': accession, 'verbose': verbose}
    do_translate = ropts.scroll['QCDB']['TRANSLATE_QCDB'].value

    # qcdb/memory [B] --> cfour/memory_size [MB]
    qopt = ropts.scroll['QCDB']['MEMORY']
    if do_translate or qopt.is_required():
        mem = int(0.000001 * qopt.value)
        print('\n\nMEMORY', mem, '\n\n')
        ropts.suggest('CFOUR', 'MEMORY_SIZE', mem, **kwgs)
        ropts.suggest('CFOUR', 'MEM_UNIT', 'MB', **kwgs)

    # qcdb/puream --> cfour/spherical
    ropts.suggest('CFOUR', 'SPHERICAL', ropts.scroll['QCDB']['PUREAM'].value, **kwgs)

    # qcdb/reference --> cfour/reference
    # TODO ref or scf__ref?
    qref = ropts.scroll['QCDB']['SCF__REFERENCE'].value
    if qref in ['RHF', 'UHF', 'ROHF']:
    #ref = {'RHF': 'RHF',
    #       'UHF': 'UHF',
    #       'ROHF': 'ROHF'}[ropts.scroll['QCDB']['REFERENCE'].value]
        ropts.suggest('CFOUR', 'REFERENCE', qref, **kwgs)

    # qcdb/scf__d_convergence --> cfour/scf_conv
    qopt = ropts.scroll['QCDB']['SCF__D_CONVERGENCE']
    if qopt.disputed():
        conv = conv_float2negexp(qopt.value)
    #conv = conv_float2negexp(ropts.scroll['QCDB']['SCF__D_CONVERGENCE'].value)
        ropts.suggest('CFOUR', 'SCF_CONV', conv, **kwgs)

    # qcdb/scf__maxiter --> cfour/scf_maxcyc
    ropts.suggest('CFOUR', 'SCF_MAXCYC', ropts.scroll['QCDB']['SCF__MAXITER'].value, **kwgs)

    # qcdb/scf__damping_percentage --> cfour/scf_damping
    damp = int(10 * ropts.scroll['QCDB']['SCF__DAMPING_PERCENTAGE'].value)
    ropts.suggest('CFOUR', 'SCF_DAMPING', damp, **kwgs)


def muster_psi4options(opt):
    """Translate psi4 keywords *opt* that have been explicitly set into
    their Cfour counterparts. Since explicitly set Cfour module keyword
    values will always be used preferentially to these inferred from
    psi4, the 'clobber' property is set to False.

    """
    text = ''
    options = defaultdict(lambda: defaultdict(dict))

    if 'GLOBALS' in opt:
        if 'PUREAM' in opt['GLOBALS']:
            options['CFOUR']['CFOUR_SPHERICAL']['value'] = \
                opt['MINTS']['PUREAM']['value']

    if 'SCF' in opt:
        if 'REFERENCE' in opt['SCF']:
            options['CFOUR']['CFOUR_REFERENCE']['value'] = \
                {'RHF': 'RHF',
                 'UHF': 'UHF',
                 'ROHF': 'ROHF'}[opt['SCF']['REFERENCE']['value']]

        if 'D_CONVERGENCE' in opt['SCF']:
            options['CFOUR']['CFOUR_SCF_CONV']['value'] = \
                conv_float2negexp(opt['SCF']['D_CONVERGENCE']['value'])

        if 'MAXITER' in opt['SCF']:
            options['CFOUR']['CFOUR_SCF_MAXCYC']['value'] = \
                opt['SCF']['MAXITER']['value']

        if 'DAMPING_PERCENTAGE' in opt['SCF']:
            options['CFOUR']['CFOUR_SCF_DAMPING']['value'] = \
                int(10 * opt['SCF']['DAMPING_PERCENTAGE']['value'])

    for item in options['CFOUR']:
        options['CFOUR'][item]['clobber'] = False
    return text, options

# Philosophy break:
#   Specification options
#   Massaging options

#   * No program's defaults should be tampered with w/o provocation

#   want all defaults applied to all programs, so p4 scf_conv is 5 and c4 scf_conv is 5
#   want separate regimes, so conv 6 covers all the p4 parts and cfour_conv = 8 covers the c4 parts
#   want mixture, so basis gets applied to c4 but others don't
#   first case, when options specified explicitly

#   [scf][d_convergence]    [cfour][cfour_scf_conv]     what happens?
#   8 from opt()            7 by default
#   6 from set {...}        7 by default                6 (guideline that psi4 format converts when clear)
#   8 from opt()            5 from set {...}            5 (local trumps)
#   6 from set {...}        5 from set {...}            5 (local trumps)
#
#   energy(name)            [cfour][cfour_calc_level]
#   c4-scf                  SCF by default
#   c4-scf                  CCSD from set {...}


def nu_muster_modelchem(name, dertype, ropts, verbose=1):
    lowername = name.lower()
    accession=2345

    if dertype == 0:
        if lowername == 'c4-cfour':
            pass  # permit clean operation of sandwich mode
        else:
            ropts.require('CFOUR', 'DERIV_LEVEL', 'ZERO', accession=accession, verbose=verbose)
    elif dertype == 1:
        ropts.require('CFOUR', 'DERIV_LEVEL', 'FIRST', accession=accession, verbose=verbose)
    elif dertype == 2:
        #ropts.require('CFOUR', 'DERIV_LEVEL', 'SECOND', accession=accession, verbose=verbose)
        ropts.require('CFOUR', 'VIBRATION', 'EXACT', accession=accession, verbose=verbose)
    else:
        raise ValidationError(f"""Requested Cfour dertype '{dertype}' is not available.""")

    if lowername == 'c4-cfour':
        pass
    elif lowername in ['c4-scf', 'c4-hf']:
        ropts.require('CFOUR', 'CALC_LEVEL', 'SCF', accession=accession, verbose=verbose)

    elif lowername == 'c4-mp2':
        ropts.require('CFOUR', 'CALC_LEVEL', 'MP2', accession=accession, verbose=verbose)

    elif lowername == 'c4-mp3':
        ropts.require('CFOUR', 'CALC_LEVEL', 'MP3', accession=accession, verbose=verbose)

    elif lowername == 'c4-mp4(sdq)':
        ropts.require('CFOUR', 'CALC_LEVEL', 'SDQ-MP4', accession=accession, verbose=verbose)

    elif lowername == 'c4-mp4':
        ropts.require('CFOUR', 'CALC_LEVEL', 'MP4', accession=accession, verbose=verbose)

    elif lowername == 'c4-cc2':
        ropts.require('CFOUR', 'CALC_LEVEL', 'CC2', accession=accession, verbose=verbose)

    elif lowername == 'c4-ccsd':
        ropts.require('CFOUR', 'CALC_LEVEL', 'CCSD', accession=accession, verbose=verbose)
        ropts.suggest('CFOUR', 'CC_PROGRAM', 'ECC', accession=accession, verbose=verbose)

    elif lowername == 'c4-ccsd-dboc':
        ropts.require('CFOUR', 'CALC_LEVEL', 'CCSD', accession=accession, verbose=verbose)
        ropts.require('CFOUR', 'DERIV_LEVEL', 'FIRST', accession=accession, verbose=verbose)
        ropts.require('CFOUR', 'DBOC', 'ON', accession=accession, verbose=verbose)
        ropts.suggest('CFOUR', 'CC_PROGRAM', 'ECC', accession=accession, verbose=verbose)

    elif lowername == 'c4-cc3':
        ropts.require('CFOUR', 'CALC_LEVEL', 'CC3', accession=accession, verbose=verbose)

    elif lowername == 'c4-ccsd(t)':
        # Can't use (T) b/c bug in xsymcor lops it off
        ropts.require('CFOUR', 'CALC_LEVEL', 'CCSD[T]', accession=accession, verbose=verbose)
        ropts.suggest('CFOUR', 'CC_PROGRAM', 'ECC', accession=accession, verbose=verbose)

    elif lowername == 'c4-ccsdt':
        # TODO, CC_PROG needs defaulting on a per-reference basis
        ropts.require('CFOUR', 'CALC_LEVEL', 'CCSDT', accession=accession, verbose=verbose)
        ropts.suggest('CFOUR', 'CC_PROGRAM', 'ECC', accession=accession, verbose=verbose)

    elif lowername == 'c4-ccsdt(q)':
        ropts.require('CFOUR', 'CALC_LEVEL', 'CCSDT(Q)', accession=accession, verbose=verbose)
        ropts.suggest('CFOUR', 'CC_PROGRAM', 'NCC', accession=accession, verbose=verbose)

    elif lowername == 'c4-ccsdtq':
        ropts.require('CFOUR', 'CALC_LEVEL', 'CCSDTQ', accession=accession, verbose=verbose)
        ropts.suggest('CFOUR', 'CC_PROGRAM', 'NCC', accession=accession, verbose=verbose)

    else:
        raise ValidationError(f"""Requested Cfour computational method '{lowername}' is not available.""")

#    # Set clobbering
#    if 'CFOUR_DERIV_LEVEL' in options['CFOUR']:
#        options['CFOUR']['CFOUR_DERIV_LEVEL']['clobber'] = True
#        options['CFOUR']['CFOUR_DERIV_LEVEL']['superclobber'] = True
#    if 'CFOUR_CALC_LEVEL' in options['CFOUR']:
#        options['CFOUR']['CFOUR_CALC_LEVEL']['clobber'] = True
#        options['CFOUR']['CFOUR_CALC_LEVEL']['superclobber'] = True
#    if 'CFOUR_CC_PROGRAM' in options['CFOUR']:
#        options['CFOUR']['CFOUR_CC_PROGRAM']['clobber'] = False

#    return ''


def muster_modelchem(name, dertype):
    """Transform calculation method *name* and derivative level *dertype*
    into options for cfour. While deliberately requested pieces,
    generally |cfour__cfour_deriv_level| and |cfour__cfour_calc_level|,
    are set to complain if contradicted ('clobber' set to True), other
    'recommended' settings, like |cfour__cfour_cc_program|, can be
    countermanded by keywords in input file ('clobber' set to False).
    Occasionally, want these pieces to actually overcome keywords in
    input file ('superclobber' set to True).

    """
    text = ''
    lowername = name.lower()
    options = defaultdict(lambda: defaultdict(dict))

    if dertype == 0:
        if lowername == 'c4-cfour':
            pass  # permit clean operation of sandwich mode
        else:
            options['CFOUR']['CFOUR_DERIV_LEVEL']['value'] = 'ZERO'
    elif dertype == 1:
        options['CFOUR']['CFOUR_DERIV_LEVEL']['value'] = 'FIRST'
    elif dertype == 2:
        options['CFOUR']['CFOUR_DERIV_LEVEL']['value'] = 'SECOND'
    else:
        raise ValidationError("""Requested Cfour dertype %d is not available.""" % (dertype))

    if lowername == 'c4-cfour':
        pass
    elif lowername in ['c4-scf', 'c4-hf']:
        options['CFOUR']['CFOUR_CALC_LEVEL']['value'] = 'SCF'

    elif lowername == 'c4-mp2':
        options['CFOUR']['CFOUR_CALC_LEVEL']['value'] = 'MP2'

    elif lowername == 'c4-mp3':
        options['CFOUR']['CFOUR_CALC_LEVEL']['value'] = 'MP3'

    elif lowername == 'c4-mp4(sdq)':
        options['CFOUR']['CFOUR_CALC_LEVEL']['value'] = 'SDQ-MP4'

    elif lowername == 'c4-mp4':
        options['CFOUR']['CFOUR_CALC_LEVEL']['value'] = 'MP4'

    elif lowername == 'c4-cc2':
        options['CFOUR']['CFOUR_CALC_LEVEL']['value'] = 'CC2'

    elif lowername == 'c4-ccsd':
        options['CFOUR']['CFOUR_CALC_LEVEL']['value'] = 'CCSD'
        options['CFOUR']['CFOUR_CC_PROGRAM']['value'] = 'ECC'

    elif lowername == 'c4-cc3':
        options['CFOUR']['CFOUR_CALC_LEVEL']['value'] = 'CC3'

    elif lowername == 'c4-ccsd(t)':
        # Can't use (T) b/c bug in xsymcor lops it off
        #options['CFOUR']['CFOUR_CALC_LEVEL']['value'] = 'CCSD(T)'
        options['CFOUR']['CFOUR_CALC_LEVEL']['value'] = 'CCSD[T]'
        options['CFOUR']['CFOUR_CC_PROGRAM']['value'] = 'ECC'

    elif lowername == 'c4-ccsdt':
        options['CFOUR']['CFOUR_CALC_LEVEL']['value'] = 'CCSDT'
        options['CFOUR']['CFOUR_CC_PROGRAM']['value'] = 'ECC'

    elif lowername == 'c4-ccsdt(q)':
        options['CFOUR']['CFOUR_CALC_LEVEL']['value'] = 'CCSDT(Q)'
        options['CFOUR']['CFOUR_CC_PROGRAM']['value'] = 'NCC'

    elif lowername == 'c4-ccsdtq':
        options['CFOUR']['CFOUR_CALC_LEVEL']['value'] = 'CCSDTQ'
        options['CFOUR']['CFOUR_CC_PROGRAM']['value'] = 'NCC'

    else:
        raise ValidationError("""Requested Cfour computational methods %d is not available.""" % (lowername))

    # Set clobbering
    if 'CFOUR_DERIV_LEVEL' in options['CFOUR']:
        options['CFOUR']['CFOUR_DERIV_LEVEL']['clobber'] = True
        options['CFOUR']['CFOUR_DERIV_LEVEL']['superclobber'] = True
    if 'CFOUR_CALC_LEVEL' in options['CFOUR']:
        options['CFOUR']['CFOUR_CALC_LEVEL']['clobber'] = True
        options['CFOUR']['CFOUR_CALC_LEVEL']['superclobber'] = True
    if 'CFOUR_CC_PROGRAM' in options['CFOUR']:
        options['CFOUR']['CFOUR_CC_PROGRAM']['clobber'] = False

    return text, options


def cfour_list():
    """Return an array of Cfour methods with energies. Appended
    to procedures['energy'].

    """
    val = []
    val.append('c4-cfour')
    val.append('c4-scf')
    val.append('c4-hf')
    val.append('c4-mp2')
    val.append('c4-mp3')
    val.append('c4-mp4(sdq)')
    val.append('c4-mp4')
    val.append('c4-cc2')
    val.append('c4-ccsd')
    val.append('c4-ccsd-dboc')
    val.append('c4-cc3')
    val.append('c4-ccsd(t)')
    val.append('c4-ccsdt')
    val.append('c4-ccsdt(q)')
    val.append('c4-ccsdtq')
    return val


def cfour_gradient_list():
    """Return an array of Cfour methods with analytical gradients.
    Appended to procedures['gradient'].

    """
    val = []
    val.append('cfour')
    val.append('c4-scf')
    val.append('c4-hf')
    val.append('c4-mp2')
    val.append('c4-mp3')
    val.append('c4-mp4(sdq)')
    val.append('c4-mp4')
    val.append('c4-cc2')
    val.append('c4-ccsd')
    val.append('c4-cc3')
    val.append('c4-ccsd(t)')
    val.append('c4-ccsdt')
    return val


def cfour_hessian_list():
    """Return an array of Cfour methods with analytical Hessians.
    Appended to procedures['hessian'].

    """
    val = []
    val.append('cfour')
    val.append('c4-scf')
    val.append('c4-hf')
    val.append('c4-mp2')
    val.append('c4-ccsd')
    val.append('c4-ccsd(t)')
    return val


def cfour_psivar_list():
    """Return a dict with keys of most Cfour methods and values of dicts
    with the PSI Variables returned by those methods. Used by cbs()
    wrapper to avoid unnecessary computations in compound methods.
    Result is appended to ``VARH``.

    """
    VARH = {}
    VARH['c4-scf'] = {
                         'c4-scf': 'SCF TOTAL ENERGY'}
    VARH['c4-hf'] = {
                          'c4-hf': 'HF TOTAL ENERGY'}
    VARH['c4-mp2'] = {
                          'c4-hf': 'HF TOTAL ENERGY',
                         'c4-mp2': 'MP2 TOTAL ENERGY'}
    VARH['c4-mp3'] = {
                          'c4-hf': 'HF TOTAL ENERGY',
                         'c4-mp2': 'MP2 TOTAL ENERGY',
                       'c4-mp2.5': 'MP2.5 TOTAL ENERGY',
                         'c4-mp3': 'MP3 TOTAL ENERGY'}
    VARH['c4-mp4(sdq)'] = {
                          'c4-hf': 'HF TOTAL ENERGY',
                         'c4-mp2': 'MP2 TOTAL ENERGY',
                       'c4-mp2.5': 'MP2.5 TOTAL ENERGY',
                         'c4-mp3': 'MP3 TOTAL ENERGY',
                    'c4-mp4(sdq)': 'MP4(SDQ) TOTAL ENERGY'}
    VARH['c4-mp4'] = {
                          'c4-hf': 'HF TOTAL ENERGY',
                         'c4-mp2': 'MP2 TOTAL ENERGY',
                       'c4-mp2.5': 'MP2.5 TOTAL ENERGY',
                         'c4-mp3': 'MP3 TOTAL ENERGY',
                    'c4-mp4(sdq)': 'MP4(SDQ) TOTAL ENERGY',
                         'c4-mp4': 'MP4(SDTQ) TOTAL ENERGY'}
    VARH['c4-cc2'] = {
                          'c4-hf': 'HF TOTAL ENERGY',
                         'c4-mp2': 'MP2 TOTAL ENERGY',
                         'c4-cc2': 'CC2 TOTAL ENERGY'}
    VARH['c4-ccsd'] = {
                          'c4-hf': 'HF TOTAL ENERGY',
                         'c4-mp2': 'MP2 TOTAL ENERGY',
                        'c4-ccsd': 'CCSD TOTAL ENERGY'}
    VARH['c4-cc3'] = {
                          'c4-hf': 'HF TOTAL ENERGY',
                         'c4-mp2': 'MP2 TOTAL ENERGY',
                         'c4-cc3': 'CC3 TOTAL ENERGY'}
    VARH['c4-ccsd(t)'] = {
                          'c4-hf': 'HF TOTAL ENERGY',
                         'c4-mp2': 'MP2 TOTAL ENERGY',
                        'c4-ccsd': 'CCSD TOTAL ENERGY',
                     'c4-ccsd(t)': 'CCSD(T) TOTAL ENERGY'}
    VARH['c4-ccsdt'] = {
                          'c4-hf': 'HF TOTAL ENERGY',
                         'c4-mp2': 'MP2 TOTAL ENERGY',
                        'c4-ccsd': 'CCSD TOTAL ENERGY',
                       'c4-ccsdt': 'CCSDT TOTAL ENERGY'}

    return VARH


def harvest_zmat(zmat):
    """Parses the contents of the Cfour ZMAT file into array and
    coordinate information. The coordinate info is converted into a
    rather dinky Molecule (no fragment, but does read charge, mult,
    unit). Return qcdb.Molecule. Written for findif zmat* where
    geometry always Cartesian and Bohr.

    """
    zmat = zmat.splitlines()[1:]  # skip comment line
    Nat = 0
    readCoord = True
    isBohr = ''
    charge = 0
    mult = 1
    molxyz = ''
    for line in zmat:
        if line.strip() == '':
            readCoord = False
        elif readCoord:
            molxyz += line + '\n'
            Nat += 1
        else:
            if line.find('CHARGE') > -1:
                idx = line.find('CHARGE')
                charge = line[idx + 7:]
                idxc = charge.find(',')
                if idxc > -1:
                    charge = charge[:idxc]
                charge = int(charge)
            if line.find('MULTIPLICITY') > -1:
                idx = line.find('MULTIPLICITY')
                mult = line[idx + 13:]
                idxc = mult.find(',')
                if idxc > -1:
                    mult = mult[:idxc]
                mult = int(mult)
            if line.find('UNITS=BOHR') > -1:
                isBohr = ' bohr'

    molxyz = f'{Nat}{isBohr}\n{charge} {mult}\n' + molxyz
    mol = Molecule(validate=False, **qcel.molparse.to_schema(qcel.molparse.from_string(molxyz, dtype='xyz+', fix_com=True, fix_orientation=True)["qm"], dtype=2))

    return mol


#def backtransform(chgeMol, permMol, chgeGrad=None, chgeDip=None):
#def format_fjobarc(fje, fjelem, fjcoord, fjgrd, map, fjdip):
def format_fjobarc(energy, amap, elem, coordinates, gradient, dipole):
    """Takes the key results from a gradient computation (*energy*,
    element Z list *elem*, *coordinates*, *gradient*,
    *dipole*, and atom ordering *amap*) and writes a string *fja*
    that exactly mimics the contents of a Cfour FJOBARC file.

    """
    fja = 'TOTENERG\n'
    fja += '%15d%15d\n' % (struct.unpack("ii", struct.pack("d", energy)))
    fja += 'COORD\n'
    Nat = len(coordinates)
    flatcoord = []
    for at in range(Nat):
        for xyz in range(3):
            flatcoord.append(coordinates[amap[at]][xyz])
    for idx in range(len(flatcoord)):
        if abs(flatcoord[idx]) < 1.0E-14:  # TODO
            flatcoord[idx] = 0.0
        fja += '%15d%15d' % (struct.unpack("ii", struct.pack("d", flatcoord[idx])))
        if idx % 2 == 1:
            fja += '\n'
    if len(flatcoord) % 2 == 1:
        fja += '\n'
    fja += 'MAP2ZMAT\n'
    for idx in range(Nat):
        fja += '%15d%15d' % (struct.unpack("ii", struct.pack("l", amap[idx] + 1)))
        if idx % 2 == 1:
            fja += '\n'
    if Nat % 2 == 1:
        fja += '\n'
    fja += 'GRD FILE\n'
    fja += '%5d%20.10f\n' % (Nat, 0.0)
    for at in range(Nat):
        fja += '%20.10f%20.10f%20.10f%20.10f\n' % (elem[at], coordinates[at][0], coordinates[at][1], coordinates[at][2])
    for at in range(Nat):
        fja += '%20.10f%20.10f%20.10f%20.10f\n' % (elem[at], gradient[at][0], gradient[at][1], gradient[at][2])
    fja += 'DIPOL FILE\n'
    fja += '%20.10f%20.10f%20.10f\n' % (dipole[0], dipole[1], dipole[2])

    return fja


def backtransform(chgeMol: 'Molecule', permMol: 'Molecule', chgeGrad=None, chgeDip=None):
    """Return `chgeMol` and `chgeGrd` to the native Cfour orientation embodied by `permMol`. Currently for vpt2.

<         p4c4 = OrientMols(p4Mol, outMol)
<         oriCoord = p4c4.transform_coordinates2(outMol)
---
>         # TODO watch out - haven't seen atom_map=False yet  # May 2019 we have now
>         rmsd, mill, amol = outMol.B787(p4Mol, atoms_map=True, mols_align=True, verbose=0)
>         oriCoord = mill.align_coordinates(outMol.geometry(np_out=True))


       rmsd, mill, amol = grdMol.B787(p4Mol, atoms_map=False, mols_align=True, verbose=0)

        oriCoord = mill.align_coordinates(grdMol.geometry(np_out=True))
        oriGrad = mill.align_gradient(np.array(grdGrad))
        if dipolDip is None:
            oriDip = None
        else:
            oriDip = mill.align_vector(np.array(dipolDip))

        if fcmHess is None:
            oriHess = None
        else:
            oriHess = mill.align_hessian(np.array(fcmHess))


    """
    # Set up array reorientation object -- opposite than usual
    amol, data = chgeMol.align(permMol, atoms_map=False, mols_align=True, verbose=0)
    mill = data['mill']

    oriCoord = mill.align_coordinates(chgeMol.geometry)
    oriElem = mill.align_atoms(np.array(chgeMol.atomic_numbers))
    oriElemMap = mill.atommap

    oriGrad = None if chgeGrad is None else mill.align_gradient(np.array(chgeGrad))
    oriDip = None if chgeDip is None else mill.align_vector(np.array(chgeDip))

    if chgeGrad is not None and chgeDip is not None:
        return oriElemMap, oriElem, oriCoord, oriGrad, oriDip
    else:
        return oriElemMap, oriElem, oriCoord


#def backtransform_grad(p4Mol, c4Mol, p4Grd, p4Dip):
#    """Here, p4Mol and p4Grd need to be turned into the native Cfour
#    orientation embodied by c4Mol. Currently for vpt2.
#
#    """
#    # Set up array reorientation object
#    p4c4 = OrientMols(c4Mol, p4Mol)  # opposite than usual
#    oriCoord = p4c4.transform_coordinates2(p4Mol)
#    oriGrad = p4c4.transform_gradient(p4Grd)
#    p4Elem = []
#    for at in range(p4Mol.natom()):
#        p4Elem.append(p4Mol.Z(at))
#    oriElem = p4c4.transform_elementlist(p4Elem)
#    oriElemMap = p4c4.Catommap
#    oriDip = p4c4.transform_vector(p4Dip)
#
#    #print p4c4
#    #print '    <<<   Input C4 Mol   >>>'
#    #c4Mol.print_out()
#    #print '    <<<   Input P4 Mol   >>>'
#    #p4Mol.print_out()
#    #print '    <<<   Input P4 Grad   >>>'
#    #if p4Grd is not None:
#    #    for item in p4Grd:
#    #        print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))
#    #print '    <<<   Rotated P4 Coord   >>>'
#    #if oriCoord is not None:
#    #    for item in oriCoord:
#    #        print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))
#    #print '    <<<   Rotated P4 Elem   >>>'
#    #if oriElem is not None:
#    #    for item in oriElem :
#    #        print('       %16.8f' % (item))
#    #print '    <<<   Rotated P4 Dip  >>>'
#    #if oriDip is not None:
#    #    print('       %16.8f %16.8f %16.8f' % (oriDip[0], oriDip[1], oriDip[2]))
#    #print '    <<<   Rotated P4 Grad   >>>'
#    #if oriGrad is not None:
#    #    for item in oriGrad:
#    #        print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))
#
#    return oriElemMap, oriElem, oriCoord, oriGrad, oriDip
#    #return oriElem, oriCoord, oriGrad, oriElemMap, oriDip


def jajo2mol(jajodic):
    """Returns a Molecule from entries in dictionary *jajodic* extracted
    from JAINDX and JOBARC.

    """
    zmap = jajodic[b'MAP2ZMAT']
    elem = jajodic[b'ATOMCHRG']
    coord = jajodic[b'COORD   ']
    Nat = len(elem)

    molxyz = '%d bohr\n\n' % (Nat)
    # TODO chgmult, though not really necessary for reorientation
    for at in range(Nat):
        posn = zmap[at] - 1
        el = 'GH' if elem[posn] == 0 else qcel.periodictable.to_E(elem[posn])
        posn *= 3
        molxyz += '%s %21.15f %21.15f %21.15f\n' % (el, coord[posn], coord[posn + 1], coord[posn + 2])
    mol = Molecule(validate=False, **qcel.molparse.to_schema(qcel.molparse.from_string(molxyz, dtype='xyz+', fix_com=True, fix_orientation=True)["qm"], dtype=2))

    return mol
