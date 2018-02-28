import os
import re

from ..driver import pe

try:
    basestring
except NameError:
    basestring = str


def format_basis_for_cfour(molrec, puream):
    """Function to print the BASIS=SPECIAL block for Cfour according
    to the active atoms in Molecule. Special short basis names
    are used by Psi4 libmints GENBAS-writer in accordance with
    Cfour constraints.

    """
    text = []
    for iat, elem in enumerate(molrec['elem']):
        text.append("""{}:Q4_{}""".format(elem.upper(), iat))
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


def extract_basis_from_genbas(basis, elem):
    """

    Parameters
    ----------
    basis : str
        Single basis for which to search GENBAS. Must be of correct casing.
    elem : str or list
        Element(s) symbols for which to search for `basis`. Any case will do.

    Returns
    -------
    str
        Reformed GENBAS text for requested subset of elements and basis.

    """
    if isinstance(elem, basestring):
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
        if baskey[1] == basis and baskey[0] in uelems:
            wantedbas[basline] = basblock

    wanted_genbas = ''.join('{}\n{}\n'.format(k, v) for k, v in wantedbas.items())

    return wanted_genbas


if __name__ == '__main__':
    extract_basis_from_genbas('qz2p', ['Co', 'H', 'H'])
    extract_basis_from_genbas('qz2p', 'co')
