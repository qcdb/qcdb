#! gradient HF/cc-pvdz on water

import os
import sys
from ..utils import *
from ..addons import *
import qcdb
import numpy as np
import pprint

h2o = qcdb.set_molecule('''
        O     0.000000000000    0.000000000000   -0.065638538099
        H     0.000000000000   -0.757480611647    0.520865616174
        H     0.000000000000    0.757480611647    0.520865616174
        ''')
print(h2o)


def check_hf(return_value):
    ref = -76.026760737428
    nre = 9.187334240165
    grads = np.array([[0.000000, -0.000000, -0.015191], [-0.000000, -0.010677, 0.007595],
        [0.000000, 0.010677, 0.007595]])
    assert compare_values(ref, qcdb.get_variable('HF TOTAL ENERGY'), 5, 'hf ref')
    assert compare_values(nre, qcdb.get_variable('NUCLEAR REPULSION ENERGY'), 5, 'nre')
    pprint.pprint(grads)
    pprint.pprint(qcdb.get_variable('CURRENT GRADIENT'))
    
    assert compare_arrays(grads, qcdb.get_variable('CURRENT GRADIENT'), 5, 'hf grad')
@using_nwchem
def test_1_hf():
    qcdb.set_options({
        'basis': 'cc-pvdz',
        'memory': '400 mb',
        'nwchem_scf__rhf': True,
        'scf__e_convergence': 1.0e-8,
        'nwchem_scf__nopen': 0,
        #'nwchem_task_hf'   : 'gradient'
    })
    print('Testing hf...')
    val = qcdb.gradient('nwc-hf')
    check_hf(val)
    pprint.pprint(val)
