import os
import re
import uuid
from typing import Dict, List, Union

import qcelemental as qcel

from ...driver import pe
from ...exceptions import ValidationError
from ...util import conv_float2negexp


def muster_molecule(molrec: Dict, ropts: 'Keywords', verbose: int = 1) -> str:
    kwgs = {'accession': uuid.uuid4(), 'verbose': verbose}

    molcmd, moldata = qcel.molparse.to_string(molrec, dtype='cfour', units='Bohr', return_data=True)

    for key, val in moldata['keywords'].items():
        ropts.require('CFOUR', key, val, **kwgs)

    return molcmd


def muster_basisset(molrec: Dict, ropts: 'Keywords', native_puream: bool, verbose: int = 1) -> str:
    """Function to print the BASIS=SPECIAL block for Cfour according
    to the active atoms in Molecule. Special short basis names
    are used by qcdb GENBAS-writer in accordance with
    Cfour constraints.

    """
    accession = uuid.uuid4()

    text = [f"""{elem.upper()}:CD_{iat + 1}""" for iat, elem in enumerate(molrec['elem'])]
    text.append('')
    text.append('')
    text = '\n'.join(text)

    ropts.require('CFOUR', 'BASIS', 'SPECIAL', accession=accession, verbose=verbose)
    #ropts.suggest('CFOUR', 'SPHERICAL', native_puream, accession=accession, verbose=verbose)
    ropts.suggest('QCDB', 'PUREAM', native_puream, accession=accession, verbose=verbose)
    # cfour or qcdb keywords here?
    # req or sugg for puream?

    #options['CFOUR']['CFOUR_SPHERICAL']['superclobber'] = True

    return text


def extract_basis_from_genbas(basis: str, elem: Union[str, List], exact: bool = True, verbose: int = 1) -> str:
    """

    Parameters
    ----------
    basis : str
        Single basis for which to search GENBAS. Must be of correct casing.
    elem : str or list
        Element(s) symbols for which to search for `basis`. Any case will do.
    exact : bool, optional
        When `True`, searches only for exact `basis`. However, for a nominal
        basis like "6-31G*", Cfour searches GENBAS for "6-31G" for hydrogen,
        for example, so when `False`, extra basis sets that start like `basis`
        are additionally returned.

    Returns
    -------
    str
        Reformed GENBAS text for requested subset of elements and basis.

    no error raised if elem missing
    """
    if isinstance(elem, str):
        uelems = [elem.upper()]
    else:
        uelems = set(el.upper() for el in elem)

    #with open('/home/psilocaluser/gits/qccddb/share/qcdb/basis/GENBAS', 'r') as handle:
    #with open('GENBAS', 'r') as handle:
    library_genbas_loc = os.sep.join([pe.data_dir, 'basis', 'GENBAS'])
    with open(library_genbas_loc, 'r') as handle:
        genbas = handle.read()

    toks = re.split('(^[A-Z]{1,2}:.*$)', genbas, flags=re.MULTILINE)
    allbas = dict(zip(toks[1::2], toks[2::2]))

    wantedbas = {}
    for basline, basblock in allbas.items():
        baskey = basline.split()[0].split(':')  # ['CO', 'qz2p']
        if exact is True:
            if baskey[1] == basis and baskey[0] in uelems:  # perfect match
                wantedbas[basline] = basblock
        else:
            if baskey[1].startswith(basis[:5]) and baskey[0] in uelems:  # loose match to accomodate composing
                wantedbas[basline] = basblock

    wanted_genbas = ''.join(f'{k}\n{v}\n' for k, v in wantedbas.items())
    if verbose >= 2:
        print('Bases plucked: {}'.format(wantedbas.keys()))

    return wanted_genbas


def muster_modelchem(name: str, dertype: int, ropts: 'Keywords', verbose: int = 1) -> None:
    """Transform calculation method *name* and derivative level *dertype*
    into options for cfour. While deliberately requested pieces,
    generally |cfour__cfour_deriv_level| and |cfour__cfour_calc_level|,
    are set to complain if contradicted ('clobber' set to True), other
    'recommended' settings, like |cfour__cfour_cc_program|, can be
    countermanded by keywords in input file ('clobber' set to False).
    Occasionally, want these pieces to actually overcome keywords in
    input file ('superclobber' set to True).

    """
    lowername = name.lower()
    accession = 2345

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

    elif lowername == 'c4-cisd':
        ropts.require('CFOUR', 'CALC_LEVEL', 'CISD', accession=accession, verbose=verbose)

    elif lowername == 'c4-qcisd':
        ropts.require('CFOUR', 'CALC_LEVEL', 'QCISD', accession=accession, verbose=verbose)

    elif lowername == 'c4-qcisd(t)':
        ropts.require('CFOUR', 'CALC_LEVEL', 'QCISD(T)', accession=accession, verbose=verbose)

    elif lowername == 'c4-lccd':
        ropts.require('CFOUR', 'CALC_LEVEL', 'LCCD', accession=accession, verbose=verbose)

    elif lowername == 'c4-lccsd':
        ropts.require('CFOUR', 'CALC_LEVEL', 'LCCSD', accession=accession, verbose=verbose)

    elif lowername == 'c4-cc2':
        ropts.require('CFOUR', 'CALC_LEVEL', 'CC2', accession=accession, verbose=verbose)

    elif lowername == 'c4-ccd':
        ropts.require('CFOUR', 'CALC_LEVEL', 'CCD', accession=accession, verbose=verbose)

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

    elif lowername == 'c4-ccsd+t(ccsd)':
        if ropts.scroll['CFOUR']['CC_PROGRAM'].value == 'ECC':
            ropts.require('CFOUR', 'CALC_LEVEL', 'CCSD[T]', accession=accession, verbose=verbose)
        else:
            ropts.require('CFOUR', 'CALC_LEVEL', 'CCSD+T(CCSD)', accession=accession, verbose=verbose)

    elif lowername == 'c4-ccsd(t)':
        # Can't use (T) b/c bug in xsymcor lops it off
        ropts.require('CFOUR', 'CALC_LEVEL', 'CCSD[T]', accession=accession, verbose=verbose)
        ropts.suggest('CFOUR', 'CC_PROGRAM', 'ECC', accession=accession, verbose=verbose)

    elif lowername == 'c4-a-ccsd(t)':
        ropts.require('CFOUR', 'CALC_LEVEL', 'CCSD(T)_L', accession=accession, verbose=verbose)

    elif lowername == 'c4-ccsdt-1a':
        ropts.require('CFOUR', 'CALC_LEVEL', 'CCSDT-1', accession=accession, verbose=verbose)

    elif lowername == 'c4-ccsdt-1b':
        # note mixed case
        ropts.require('CFOUR', 'CALC_LEVEL', 'CCSDT-1b', accession=accession, verbose=verbose)

    elif lowername == 'c4-ccsdt-2':
        ropts.require('CFOUR', 'CALC_LEVEL', 'CCSDT-2', accession=accession, verbose=verbose)

    elif lowername == 'c4-ccsdt-3':
        ropts.require('CFOUR', 'CALC_LEVEL', 'CCSDT-3', accession=accession, verbose=verbose)

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


def muster_inherited_keywords(ropts: 'Keywords', verbose: int = 1) -> None:
    """Translate psi4 keywords *opt* that have been explicitly set into
    their Cfour counterparts. Since explicitly set Cfour module keyword
    values will always be used preferentially to these inferred from
    psi4, the 'clobber' property is set to False.

    """
    import sys
    accession = sys._getframe().f_code.co_name + '_' + str(uuid.uuid4())[:8]
    kwgs = {'accession': accession, 'verbose': verbose}
    do_translate = ropts.scroll['QCDB']['TRANSLATE_QCDB'].value

    # qcdb/memory [B] --> cfour/memory_size [QW]
    qopt = ropts.scroll['QCDB']['MEMORY']
    if do_translate or qopt.is_required():
        mem = int(qopt.value / 8.0)
        print('\n\nMEMORY', qopt.value, mem, '\n\n')
        ropts.suggest('CFOUR', 'MEMORY_SIZE', mem, **kwgs)
        ropts.suggest('CFOUR', 'MEM_UNIT', "INTEGERWORDS", **kwgs)

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

    # qcdb/freeze_core --> cfour/frozen_core
    ropts.suggest("CFOUR", "FROZEN_CORE", ropts.scroll["QCDB"]["FREEZE_CORE"].value, **kwgs)


if __name__ == '__main__':
    extract_basis_from_genbas('qz2p', ['Co', 'H', 'H'])
    extract_basis_from_genbas('qz2p', 'co')
