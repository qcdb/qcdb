import copy
import pprint

import pytest
import qcelemental as qcel

import qcdb
from qcdb.programs.dftd3 import runner as dftd3

from .addons import *
from .utils import *

eneyne = """
C   0.000000  -0.667578  -2.124659
C   0.000000   0.667578  -2.124659
H   0.923621  -1.232253  -2.126185
H  -0.923621  -1.232253  -2.126185
H  -0.923621   1.232253  -2.126185
H   0.923621   1.232253  -2.126185
--
C   0.000000   0.000000   2.900503
C   0.000000   0.000000   1.693240
H   0.000000   0.000000   0.627352
H   0.000000   0.000000   3.963929
"""

db3lypd3bj = {'dashlevel': 'd3bj', 'dashparams': {'s8': 1.9889, 's6': 1.0, 'a2': 4.4211, 'a1': 0.3981}, 'do_gradient': True, 'functional': 'b3lyp'}

dpbed3zero = {'dashlevel': 'd3zero', 'dashparams': {'s6': 1.0,  's8': 0.722, 'sr6': 1.217, 'alpha6': 14.0}, 'do_gradient': False, 'functional': 'pbe'}


def test_recon_1a():
    res = dftd3.from_arrays_qc(functional='b3lyp', dashlevel='d3bj')
    assert compare_recursive(db3lypd3bj, res, atol=1.e-4)
    assert compare_strings('B3LYP-D3BJ', compute_key(res), 'key')


def test_recon_1b():
    res = dftd3.from_arrays_qc(functional='b3LYP', dashlevel='D3bj', dertype=None)
    assert compare_recursive(db3lypd3bj, res, atol=1.e-4)
    assert compare_strings('B3LYP-D3BJ', compute_key(res), 'key')


def test_recon_1c():
    res = dftd3.from_arrays_qc(dashparams={'s8': 1.9889, 's6': 1.0, 'a2': 4.4211, 'a1': 0.3981}, dashlevel='d3bj', dertype=1)
    assert compare_recursive(db3lypd3bj, res, atol=1.e-4)
    assert compare_strings('B3LYP-D3BJ', compute_key(res), 'key')


def test_recon_1d():
    res = dftd3.from_arrays_qc(functional='b3lyp', dashlevel='d3bJ', dashparams={'a2': 4.4211}, dertype='first')
    assert compare_recursive(db3lypd3bj, res, atol=1.e-4)
    assert compare_strings('B3LYP-D3BJ', compute_key(res), 'key')


def test_recon_1e():
    ans = copy.deepcopy(db3lypd3bj)
    ans['functional'] = ''
    ans['dashparams']['a2'] = 5.4211

    res = dftd3.from_arrays_qc(functional='b3lyp', dashlevel='d3bJ', dashparams={'a2': 5.4211}, dertype='first')
    assert compare_recursive(ans, res, atol=1.e-4)
    assert compare_strings('-D3BJ', compute_key(res), 'key')


def test_recon_1f():
    with pytest.raises(qcdb.ValidationError):
        dftd3.from_arrays_qc(functional='b3lyp', dashlevel='d3bJ', dashparams={'a3': 5.4211}, dertype='first')


def test_recon_2a():
    res = dftd3.from_arrays_qc(functional='pbe', dashlevel='d3zero', dertype=0)
    assert compare_recursive(dpbed3zero, res, atol=1.e-4)
    assert compare_strings('PBE-D3ZERO', compute_key(res), 'key')


def test_recon_2b():
    res = dftd3.from_arrays_qc(functional='pbe', dashlevel='d3', dertype='Energy')
    assert compare_recursive(dpbed3zero, res, atol=1.e-4)
    assert compare_strings('PBE-D3ZERO', compute_key(res), 'key')


def test_3():
    sys = qcel.molparse.from_string(eneyne)['qm']

    res = dftd3.run_dftd3(molrec=sys, functional='b3lyp', dashlevel='d3bj')
    assert compare_strings('B3LYP-D3BJ', compute_key(res), 'key')

    #assert False


def compute_key(pjrec):
    return '-'.join([pjrec['functional'], pjrec['dashlevel']]).upper()









    #"""dftd3/energy"""
    #! Exercises the various DFT-D corrections, both through python directly and through c++

ref_d2         = [-0.00390110, -0.00165271, -0.00058118]
ref_d3zero     = [-0.00285088, -0.00084340, -0.00031923]
ref_d3bj       = [-0.00784595, -0.00394347, -0.00226683]

ref_pbe_d2     = [-0.00278650, -0.00118051, -0.00041513]
ref_pbe_d3zero = [-0.00175474, -0.00045421, -0.00016839]
ref_pbe_d3bj   = [-0.00475937, -0.00235265, -0.00131239]

seneyne = """
C   0.000000  -0.667578  -2.124659
C   0.000000   0.667578  -2.124659
H   0.923621  -1.232253  -2.126185
H  -0.923621  -1.232253  -2.126185
H  -0.923621   1.232253  -2.126185
H   0.923621   1.232253  -2.126185
--
C   0.000000   0.000000   2.900503
C   0.000000   0.000000   1.693240
H   0.000000   0.000000   0.627352
H   0.000000   0.000000   3.963929
"""

@using_dftd3
def test_10_qmol():
    eneyne = qcdb.set_molecule(seneyne)
    eneyne.update_geometry()

    E, G = eneyne.run_dftd3('b3lyp', 'd2')
    assert compare_values(ref_d2[0], E, 7, 'Q: Ethene-Ethyne -D2')


@using_dftd3
@using_psi4
def test_10_pmol():
    import psi4
    eneyne = psi4.geometry(seneyne)
    eneyne.update_geometry()

    E, G = eneyne.run_dftd3('b3lyp', 'd2')
    assert compare_values(ref_d2[0], E, 7, 'P: Ethene-Ethyne -D2')

@using_dftd3
def test_11_energy():
    eneyne = qcdb.set_molecule(seneyne)
    eneyne.update_geometry()

    E, jrec = qcdb.energy('d3-b3lyp-d2', return_wfn=True)
    assert compare_values(ref_d2[0], E, 7, 'P: Ethene-Ethyne -D2')
    assert compare_values(ref_d2[0], jrec['qcvars']['DISPERSION CORRECTION ENERGY'].data, 7, tnm())
    assert compare_values(ref_d2[0], jrec['qcvars']['B3LYP-D2 DISPERSION CORRECTION ENERGY'].data, 7, tnm())

@using_dftd3
def test_11_b():
    eneyne = qcdb.set_molecule(seneyne)
    eneyne.update_geometry()
    mA = eneyne.extract_subsets(1)
    mB = eneyne.extract_subsets(2)

    E, jrec = qcdb.energy('d3-b3lyp-d3bj', return_wfn=True, molecule=mA)
    assert compare_values(ref_d3bj[1], E, 7, tnm())
    assert compare_values(ref_d3bj[1], jrec['qcvars']['DISPERSION CORRECTION ENERGY'].data, 7, tnm())
    assert compare_values(ref_d3bj[1], jrec['qcvars']['B3LYP-D3(BJ) DISPERSION CORRECTION ENERGY'].data, 7, tnm())
