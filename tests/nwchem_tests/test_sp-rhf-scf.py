#! single-point HF/cc-pvdz on water
#nopen = 0 is singlet (Default)

import os
import sys
from utils import *
from addons import *
import qcdb

h2o = qcdb.set_molecule('''
        O     0.000000000000    0.000000000000   -0.065638538099
        H     0.000000000000   -0.757480611647    0.520865616174
        H     0.000000000000    0.757480611647    0.520865616174
        ''')
print(h2o)

def check_hf(return_value, is_df):
    if is_df:
        ref         =       -76.026760737428
    else:
        ref         =       -76.026760737428
    
    assert compare_values(ref, qcdb.get_variable('HF TOTAL ENERGY'), 5, 'total scf')

@using_nwchem
def test_1_hf():
    qcdb.set_options({
        'basis'     :   'cc-pvdz',
        'memory'    :   '400 mb',
        'nwchem_scf':   'rhf',
        'nwchem_scf_thresh':    1.0e-8,
        'nwchem_scf_nopen' :    0
        })
    print('Testing hf...')
    val = qcdb.energy('nwc-hf')
    check_hf(val, is_df=True)
