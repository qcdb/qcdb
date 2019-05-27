import os
import re
import uuid
import collections

import qcelemental as qcel

from ..driver import pe

try:
    basestring
except NameError:
    basestring = str


def format_molecule_for_cfour(molrec, ropts, verbose=1):
    accession=3456

    units = 'Bohr'
    molcmd = qcel.molparse.to_string(molrec, dtype='cfour', units=units)

    ropts.require('CFOUR', 'CHARGE', int(molrec['molecular_charge']), accession=accession, verbose=verbose)
    ropts.require('CFOUR', 'MULTIPLICITY', molrec['molecular_multiplicity'], accession=accession, verbose=verbose)
    ropts.require('CFOUR', 'UNITS', units.upper(), accession=accession, verbose=verbose)
    ropts.require('CFOUR', 'COORDINATES', 'CARTESIAN', accession=accession, verbose=verbose)

    #ropts['CFOUR']['CFOUR_CHARGE']['clobber'] = True
    #ropts['CFOUR']['CFOUR_MULTIPLICITY']['clobber'] = True
    #ropts['CFOUR']['CFOUR_UNITS']['clobber'] = True
    #ropts['CFOUR']['CFOUR_COORDINATES']['clobber'] = True

    return molcmd

def format_basis_for_cfour(molrec, ropts, native_puream, verbose=1): #puream):
    """Function to print the BASIS=SPECIAL block for Cfour according
    to the active atoms in Molecule. Special short basis names
    are used by qcdb GENBAS-writer in accordance with
    Cfour constraints.

    """
    accession = uuid.uuid4()

    text = []
    for iat, elem in enumerate(molrec['elem']):
        text.append("""{}:CD_{}""".format(elem.upper(), iat + 1))
    text.append('')
    text.append('')
    text = '\n'.join(text)

    ropts.require('CFOUR', 'BASIS', 'SPECIAL', accession=accession, verbose=verbose)
    #ropts.suggest('CFOUR', 'SPHERICAL', native_puream, accession=accession, verbose=verbose)
    ropts.suggest('QCDB', 'PUREAM', native_puream, accession=accession, verbose=verbose)
    # cfour or qcdb keywords here?
    # req or sugg for puream?

    #options['CFOUR']['CFOUR_BASIS']['value'] = 'SPECIAL'
    #options['CFOUR']['CFOUR_SPHERICAL']['value'] = puream
    #options['CFOUR']['CFOUR_BASIS']['clobber'] = True
    #options['CFOUR']['CFOUR_SPHERICAL']['clobber'] = True
    #options['CFOUR']['CFOUR_BASIS']['superclobber'] = True
    #options['CFOUR']['CFOUR_SPHERICAL']['superclobber'] = True

    return text

def old_format_basis_for_cfour(molrec, puream):
    """Function to print the BASIS=SPECIAL block for Cfour according
    to the active atoms in Molecule. Special short basis names
    are used by qcdb GENBAS-writer in accordance with
    Cfour constraints.

    """
    text = []
    for iat, elem in enumerate(molrec['elem']):
        text.append("""{}:CD_{}""".format(elem.upper(), iat + 1))
    text.append('')
    text.append('')
    text = '\n'.join(text)

    options = collections.defaultdict(lambda: collections.defaultdict(dict))
    options['CFOUR']['CFOUR_BASIS']['value'] = 'SPECIAL'
    options['CFOUR']['CFOUR_SPHERICAL']['value'] = puream

    options['CFOUR']['CFOUR_BASIS']['clobber'] = True
    options['CFOUR']['CFOUR_SPHERICAL']['clobber'] = True

    options['CFOUR']['CFOUR_BASIS']['superclobber'] = True
    options['CFOUR']['CFOUR_SPHERICAL']['superclobber'] = True

    return text, options


def extract_basis_from_genbas(basis, elem, exact=True, verbose=1):
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
    #allbas = dict(zip(toks[1::2], toks[::2]))
    allbas = dict(zip(toks[1::2], toks[2::2]))
    #print('AAA', toks[0])
    #print('BBB', toks[1])
    #print('CCC', toks[2])
    #print('YYY', toks[-2])
    #print('ZZZ', toks[-1])

    wantedbas = {}
    for basline, basblock in allbas.items():
        baskey = basline.split()[0].split(':')  # ['CO', 'qz2p']
        if exact is True:
            if baskey[1] == basis and baskey[0] in uelems:  # perfect match
                wantedbas[basline] = basblock
        else:
            if baskey[1].startswith(basis[:5]) and baskey[0] in uelems:  # loose match to accomodate composing
                wantedbas[basline] = basblock

    wanted_genbas = ''.join('{}\n{}\n'.format(k, v) for k, v in wantedbas.items())
    if verbose >= 2:
        print('Bases plucked: {}'.format(wantedbas.keys()))

    return wanted_genbas


if __name__ == '__main__':
    extract_basis_from_genbas('qz2p', ['Co', 'H', 'H'])
    extract_basis_from_genbas('qz2p', 'co')
