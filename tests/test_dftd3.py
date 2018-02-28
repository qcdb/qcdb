import sys
import copy

import pytest

from utils import *

import qcdb
from qcdb.iface_dftd3 import dftd3


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
    assert compare_dicts(db3lypd3bj, res, 4, sys._getframe().f_code.co_name)
    

def test_recon_1b():
    res = dftd3.from_arrays_qc(functional='b3LYP', dashlevel='D3bj', dertype=None)
    assert compare_dicts(db3lypd3bj, res, 4, sys._getframe().f_code.co_name)
    

def test_recon_1c():
    res = dftd3.from_arrays_qc(dashparams={'s8': 1.9889, 's6': 1.0, 'a2': 4.4211, 'a1': 0.3981}, dashlevel='d3bj', dertype=1)
    assert compare_dicts(db3lypd3bj, res, 4, sys._getframe().f_code.co_name)
    

def test_recon_1d():
    res = dftd3.from_arrays_qc(functional='b3lyp', dashlevel='d3bJ', dashparams={'a2': 4.4211}, dertype='first')
    assert compare_dicts(db3lypd3bj, res , 4, sys._getframe().f_code.co_name)


def test_recon_1e():
    ans = copy.deepcopy(db3lypd3bj)
    ans['functional'] = ''
    ans['dashparams']['a2'] = 5.4211

    res = dftd3.from_arrays_qc(functional='b3lyp', dashlevel='d3bJ', dashparams={'a2': 5.4211}, dertype='first')
    assert compare_dicts(ans, res, 4, sys._getframe().f_code.co_name)
    

def test_recon_1f():
    with pytest.raises(qcdb.ValidationError):
        dftd3.from_arrays_qc(functional='b3lyp', dashlevel='d3bJ', dashparams={'a3': 5.4211}, dertype='first')
    

def test_recon_2a():
    res = dftd3.from_arrays_qc(functional='pbe', dashlevel='d3zero', dertype=0)
    assert compare_dicts(dpbed3zero, res, 4, sys._getframe().f_code.co_name)
    

def test_recon_2b():
    res = dftd3.from_arrays_qc(functional='pbe', dashlevel='d3', dertype='Energy')
    assert compare_dicts(dpbed3zero, res, 4, sys._getframe().f_code.co_name)
    

def test_3():
    sys = qcdb.molparse.from_string(eneyne)['qm']

    res = dftd3.run_dftd3(molrec=sys, functional='b3lyp', dashlevel='d3bj')

    #assert False



