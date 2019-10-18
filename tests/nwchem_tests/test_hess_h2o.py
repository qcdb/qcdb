import os
import sys
from ..addons import *
from ..utils import *
import qcdb
import numpy as np

from qcengine.progrmas.util.hessparse import load_hessian

def check_hess(return_value):
    scf =   -74.888142460799
    grads   =   np.array([[0.000000,   0.000000,   0.058550],
                        [-0.140065,   0.000000,  -0.029275]
                        [ 0.140065,   0.000000,  -0.029275]])

    assert compare_values(scf, qcdb.variable('HF TOTAL ENERGY'), 5, 'scf')
    assert compare_arrays(grads, qcdb.variable('CURRENT GRADIENT'), 5, 'scf grad')

@using_nwchem
def test_hess():
    h2o = qcdb.set_molecule('''
        O      0.00000000    0.00000000    0.00000000
        H      0.00000000    1.93042809   -1.10715266
        H      0.00000000   -1.93042809   -1.10715266
        units au''')

    qcdb.set_options({
        'basis' :   'sto-3g',
        'scf__e_convergence' :  1e-6,
        #'nwchem_driver__tight': True
        })
    val = qcdb.hessian('nwc-scf')
    check_hess(val)
