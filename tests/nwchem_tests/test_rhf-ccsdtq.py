#! single-point CCSDTQ/6-31g* on water
import os
import sys
from utils import *
from addons import *
import qcdb

h2o= qcdb.set_molecule('''
        O 0.000000000000    0.000000000000   -0.065638538099
        H 0.000000000000   -0.757480611647    0.520865616174
        H 0.000000000000    0.757480611647    0.520865616174
        ''')
print(h2o)

def check_ccsdtq(return_value):
    ref         =       -76.010496307079
    nre         =         9.187334240165
    ccsdtq      =       -76.210368641377713
    ccsdtq_corl =        -0.199872334299139
        
    assert compare_values(ref, qcdb.get_variable('HF TOTAL ENERGY'), 6, 'hf ref')  #TEST
    assert compare_values(nre, qcdb.get_variable('NUCLEAR REPULSION ENERGY'), 5, 'nre')
    assert compare_values(ccsdtq, qcdb.get_variable('CCSDTQ TOTAL ENERGY'), 6, 'CCSDTQ')  #TEST
    assert compare_values(ccsdtq_corl, qcdb.get_variable('CCSDTQ CORRELATION ENERGY'), 6, 'CCSDTQ corl')  #TEST

@using_nwchem
def test_1_ccsdtq():
    qcdb.set_options({
        'basis': '6-31g*',
        #'memory': '2000 mb',
        'nwchem_memory': '[total, 2000, global, 1700,mb]',
        #'nwchem_verify'    : True,
        'nwchem_tce_dft': False,
        'nwchem_tce_module': 'CCSDTQ',
        'nwchem_tce' : True,
        'nwchem_tce_thresh': 1.0e-7
        })
    print('Testing CCSDTQ (df)...')
    val = qcdb.energy('nwc-ccsdtq')
    check_ccsdtq(val)
