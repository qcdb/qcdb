import copy
import pprint

import numpy as np
import qcelemental as qcel

import qcdb

from .utils import *

_arrays_prov_stamp = {"creator": "QCElemental", "version": "1.0", "routine": "qcelemental.molparse.from_arrays"}
_string_prov_stamp = {"creator": "QCElemental", "version": "1.0", "routine": "qcelemental.molparse.from_string"}

fullans1a = {
    "geom": np.array([0.0, 0.0, 0.0, 1.0, 0.0, 0.0]),
    "elea": np.array([16, 1]),
    "elez": np.array([8, 1]),
    "elem": np.array(["O", "H"]),
    "mass": np.array([15.99491462, 1.00782503]),
    "real": np.array([True, True]),
    "elbl": np.array(["", ""]),
    "units": "Angstrom",
    "fix_com": True,
    "fix_orientation": False,
    "fragment_separators": [],
    "fragment_charges": [0.0],
    "fragment_multiplicities": [2],
    "molecular_charge": 0.0,
    "molecular_multiplicity": 2,
}

# subject6 = """
#    0 1
#    O1    0         0     0.118720
#    h2   -0.753299, 0.0, -0.474880
#
#    H3    0.753299, 0.0, -0.474880
#
#    --
#    efp h2O -2.12417561  1.22597097 -0.95332054 -2.902133 1.734999 -1.953647
# --
# efp ammoniA
#     0.98792    1.87681    2.85174
# units au
#     1.68798    1.18856    3.09517
#
#     1.45873    2.55904    2.27226
#
# """
#
# ans6 = {'units': 'Bohr',
#        'geom': [0., 0., 0.118720, -0.753299, 0.0, -0.474880, 0.753299, 0.0, -0.474880],
#        'elbl': ['O1', 'h2', 'H3'],
#        'fragment_charges': [0.],
#        'fragment_multiplicities': [1],
#        'fragment_separators': [],
#        'fragment_files': ['h2O', 'ammoniA'],
#        'geom_hints': [[-2.12417561,  1.22597097, -0.95332054, -2.902133, 1.734999, -1.953647],
#                           [0.98792,    1.87681,    2.85174, 1.68798 ,   1.18856  ,  3.09517, 1.45873  ,  2.55904  ,  2.27226]],
#        'hint_types': ['xyzabc', 'points'],
#        }
#
# fullans6 = {'qm': {'geom': np.array([0., 0., 0.118720, -0.753299, 0.0, -0.474880, 0.753299, 0.0, -0.474880]),
#                   'elea': np.array([16, 1, 1]),
#                   'elez': np.array([8, 1, 1]),
#                   'elem': np.array(['O', 'H', 'H']),
#                   'mass': np.array([ 15.99491462,   1.00782503,   1.00782503]),
#                   'real': np.array([True, True, True]),
#                   'elbl': np.array(['1', '2', '3']),
#                   'units': 'Bohr',
#                   'fix_com': True,
#                   'fix_orientation': True,
#                   'fix_symmetry': 'c1',
#                   'fragment_charges': [0.],
#                   'fragment_multiplicities': [1],
#                   'fragment_separators': [],
#                   'molecular_charge': 0.,
#                   'molecular_multiplicity': 1},
#           'efp': {'fragment_files': ['h2o', 'ammonia'],
#                   'geom_hints': [[-2.12417561,  1.22597097, -0.95332054, -2.902133, 1.734999, -1.953647],
#                                  [0.98792,    1.87681,    2.85174, 1.68798 ,   1.18856  ,  3.09517, 1.45873  ,  2.55904  ,  2.27226]],
#                   'hint_types': ['xyzabc', 'points'],
#                   'units': 'Bohr',
#                   'fix_com': True,
#                   'fix_orientation': True,
#                   'fix_symmetry': 'c1',
#        }}

# @using_pylibefp
# def test_psi4_qmefp_6d():
#    subject = subject6
#
#    fullans = copy.deepcopy(fullans6)
#    fullans['efp']['geom'] = np.array([-2.22978429,  1.19270015, -0.99721732, -1.85344873,  1.5734809 ,
#        0.69660583, -0.71881655,  1.40649303, -1.90657336,  0.98792   ,
#        1.87681   ,  2.85174   ,  2.31084386,  0.57620385,  3.31175679,
#        1.87761143,  3.16604791,  1.75667803,  0.55253064,  2.78087794,
#        4.47837555])
#    fullans['efp']['elea'] = np.array([16, 1, 1, 14, 1, 1, 1])
#    fullans['efp']['elez'] = np.array([8, 1, 1, 7, 1, 1, 1])
#    fullans['efp']['elem'] = np.array(['O', 'H', 'H', 'N', 'H', 'H', 'H'])
#    fullans['efp']['mass'] = np.array([15.99491462, 1.00782503, 1.00782503, 14.00307400478, 1.00782503, 1.00782503, 1.00782503])
#    fullans['efp']['real'] = np.array([True, True, True, True, True, True, True])
#    fullans['efp']['elbl'] = np.array(['_a01o1', '_a02h2', '_a03h3', '_a01n1', '_a02h2', '_a03h3', '_a04h4'])
#    fullans['efp']['fragment_separators'] = [3]
#    fullans['efp']['fragment_charges'] = [0., 0.]
#    fullans['efp']['fragment_multiplicities'] = [1, 1]
#    fullans['efp']['molecular_charge'] = 0.
#    fullans['efp']['molecular_multiplicity'] = 1
#    fullans['efp']['hint_types'] = ['xyzabc', 'xyzabc']
#    fullans['efp']['geom_hints'][1] = [1.093116487139866, 1.9296501432128303, 2.9104336205167156, -1.1053108079381473, 2.0333070957565544, -1.488586877218809]
#
#    final, intermed = qcdb.molparse.from_string(subject, return_processed=True)
#
#    import pylibefp
#    efpobj = pylibefp.from_dict(final['efp'])
#    efpfinal = efpobj.to_dict()
#    efpfinal = qcdb.molparse.from_arrays(speclabel=False, domain='efp', **efpfinal)
#
#    assert compare_molrecs(fullans['qm'], final['qm'], 4, sys._getframe().f_code.co_name + ': full qm')
#    assert compare_molrecs(fullans['efp'], efpfinal, 4, sys._getframe().f_code.co_name + ': full efp')
#
#    import json
#    from qcdb.util import unnp
#    print(json.dumps(unnp(final['qm']), sort_keys=True, indent=4))
#    print(json.dumps(unnp(efpfinal), sort_keys=True, indent=4))


def test_qmol_11a():
    fullans = copy.deepcopy(fullans1a)
    fullans["provenance"] = _string_prov_stamp

    asdf = qcdb.Molecule(fullans)

    assert compare_molrecs(fullans, asdf.to_dict(), tnm(), relative_geoms="align", atol=1.0e-4)
    assert compare_integers(2, asdf.natom(), tnm())


def test_qmol_11b():
    fullans = copy.deepcopy(fullans1a)
    fullans["provenance"] = _arrays_prov_stamp

    asdf = qcdb.Molecule(geom=[0.0, 0.0, 0.0, 1.0, 0.0, 0.0], elez=[8, 1], fix_com=True)

    assert compare_molrecs(fullans, asdf.to_dict(), tnm(), relative_geoms="align", atol=1.0e-4)
    assert compare_integers(2, asdf.natom(), tnm())


def test_qmol_11c():
    fullans = copy.deepcopy(fullans1a)
    fullans["provenance"] = _string_prov_stamp

    asdf = qcdb.Molecule("""nocom\n8 0 0 0\n1 1 0 0""", dtype="psi4")

    assert compare_molrecs(fullans, asdf.to_dict(), tnm(), relative_geoms="align", atol=1.0e-4)
    assert compare_integers(2, asdf.natom(), tnm())


def test_qmol_11d():
    fullans = copy.deepcopy(fullans1a)
    fullans["provenance"] = _string_prov_stamp

    asdf = qcdb.Molecule("""nocom\n8 0 0 0\n1 1 0 0""", dtype="psi4+")

    assert compare_molrecs(fullans, asdf.to_dict(), tnm(), relative_geoms="align", atol=1.0e-4)
    assert compare_integers(2, asdf.natom(), tnm())


def test_qmol_11e():
    fullans = copy.deepcopy(fullans1a)
    fullans["provenance"] = _string_prov_stamp

    asdf = qcdb.Molecule("""2\n\nO 0 0 0 \n1 1 0 0 """, dtype="xyz", fix_com=True)

    assert compare_molrecs(fullans, asdf.to_dict(), tnm(), relative_geoms="align", atol=1.0e-4)
    assert compare_integers(2, asdf.natom(), tnm())


def test_qmol_11f():
    fullans = copy.deepcopy(fullans1a)
    fullans["provenance"] = _arrays_prov_stamp

    asdf = qcdb.Molecule.from_dict(fullans)

    assert compare_molrecs(fullans, asdf.to_dict(), tnm(), relative_geoms="align", atol=1.0e-4)
    assert compare_integers(2, asdf.natom(), tnm())


def test_qmol_11g():
    fullans = copy.deepcopy(fullans1a)
    fullans["provenance"] = _arrays_prov_stamp

    asdf = qcdb.Molecule.from_arrays(geom=[0.0, 0.0, 0.0, 1.0, 0.0, 0.0], elez=[8, 1], fix_com=True)

    assert compare_molrecs(fullans, asdf.to_dict(), tnm(), relative_geoms="align", atol=1.0e-4)
    assert compare_integers(2, asdf.natom(), tnm())


def test_qmol_11h():
    fullans = copy.deepcopy(fullans1a)
    fullans["provenance"] = _string_prov_stamp

    asdf = qcdb.Molecule.from_string("""nocom\n8 0 0 0\n1 1 0 0""")

    assert compare_molrecs(fullans, asdf.to_dict(), tnm(), relative_geoms="align", atol=1.0e-4)
    assert compare_integers(2, asdf.natom(), tnm())


def test_qmol_11i():
    fullans = copy.deepcopy(fullans1a)
    fullans["provenance"] = _string_prov_stamp

    asdf = qcdb.Molecule.from_string("""nocom\n8 0 0 0\n1 1 0 0""")

    assert compare_molrecs(fullans, asdf.to_dict(), tnm(), relative_geoms="align", atol=1.0e-4)
    assert compare_integers(2, asdf.natom(), tnm())


def test_qmol_11j():
    fullans = copy.deepcopy(fullans1a)
    fullans["provenance"] = _string_prov_stamp

    asdf = qcdb.Molecule.from_string("""2\n\nO 0 0 0 \n1 1 0 0 """, fix_com=True)

    assert compare_molrecs(fullans, asdf.to_dict(), tnm(), relative_geoms="align", atol=1.0e-4)
    assert compare_integers(2, asdf.natom(), tnm())


# @using_psi4_molrec
# def test_pmol_11k():
#    import psi4
#    asdf = psi4.core.Molecule.from_dict(fullans1a)
#    assess_mol_11(asdf, '[16] psi4.core.Molecule.from_dict(dict)')
#
#    assert compare_molrecs(fullans, asdf.to_dict(), 4, tnm(), relative_geoms='align')
#    assert compare_integers(2, asdf.natom(), tnm())
#
# @using_psi4_molrec
# def test_pmol_11l():
#    import psi4
#    asdf = psi4.core.Molecule.from_arrays(geom=[ 0.,  0.,  0.,  1.,  0.,  0.], elez=[8, 1], fix_com=True)
#    assess_mol_11(asdf, '[17] psi4.core.Molecule.from_arrays(geom, elez)')
#
#    assert compare_molrecs(fullans, asdf.to_dict(), 4, tnm(), relative_geoms='align')
#    assert compare_integers(2, asdf.natom(), tnm())
#
# @using_psi4_molrec
# def test_pmol_11m():
#    import psi4
#    asdf = psi4.core.Molecule.from_string("""nocom\n8 0 0 0\n1 1 0 0""")
#    assess_mol_11(asdf, '[18] psi4.core.Molecule.from_string(str, dtype="psi4")')
#
#    assert compare_molrecs(fullans, asdf.to_dict(), 4, tnm(), relative_geoms='align')
#    assert compare_integers(2, asdf.natom(), tnm())
#
# @using_psi4_molrec
# def test_pmol_11n():
#    import psi4
#    asdf = psi4.core.Molecule.from_string("""nocom\n8 0 0 0\n1 1 0 0""")
#    assess_mol_11(asdf, '[19] psi4.core.Molecule.from_string(str, dtype="psi4+")')
#
#    assert compare_molrecs(fullans, asdf.to_dict(), 4, tnm(), relative_geoms='align')
#    assert compare_integers(2, asdf.natom(), tnm())
#
# @using_psi4_molrec
# def test_pmol_11o():
#    import psi4
#    asdf = psi4.core.Molecule.from_string("""2\n\nO 0 0 0 \n1 1 0 0 """, fix_com=True)
#    assess_mol_11(asdf, '[20] psi4.core.Molecule.from_string(str, dtype="xyz")')
#
#    assert compare_molrecs(fullans, asdf.to_dict(), 4, tnm(), relative_geoms='align')
#    assert compare_integers(2, asdf.natom(), tnm())


def test_qmol_12():
    fullans = copy.deepcopy(fullans1a)
    fullans["provenance"] = _arrays_prov_stamp

    asdf = qcdb.Molecule(geom=[0.0, 0.0, 0.0, 1.0, 0.0, 0.0], elez=[8, 1], fix_com=True)

    assert compare_molrecs(fullans, asdf.to_dict(), tnm(), relative_geoms="align", atol=1.0e-4)
    assert compare_integers(2, asdf.natom(), tnm())

    import json

    smol = json.dumps(asdf.to_dict(np_out=False))
    dmol = json.loads(smol)

    asdf2 = qcdb.Molecule(dmol)

    assert compare_molrecs(fullans, asdf2.to_dict(), tnm(), relative_geoms="align", atol=1.0e-4)
    assert compare_integers(2, asdf2.natom(), tnm())
