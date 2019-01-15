# single-point DFT-CCSDT/sto-3g on water

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


def check_dft(return_value, is_dft):
    if is_dft:
        dft_ref = -95.217126584667
        nre = 9.187334240165
        ccsdt_tot = -95.263267683066744
        ccsdt_corl = -0.046141098399779
    else:
        dft_ref = -95.217126584667
        ccsdt_tot = -95.263267683066744
        ccsdt_corl = -0.046141098399779
        nre = 9.187334240165


#    assert compare_values(dft_ref, qcdb.get_variable('DFT TOTAL ENERGY'), 5, 'dft ref')
    assert compare_values(ccsdt_tot, qcdb.get_variable('CCSDT TOTAL ENERGY'), 5, 'ccsdt total')
    assert compare_values(ccsdt_corl, qcdb.get_variable('CCSDT CORRELATION ENERGY'), 5, 'ccsdt corl')
    assert compare_values(nre, qcdb.get_variable('NUCLEAR REPULSION ENERGY'), 'nre')


@using_nwchem
def test_1_dft():
    qcdb.set_options({
        'memory': '600 mb',
        'basis': 'sto-3g',
        #'nwchem_dft__convergence__density': 1.0e-12,
        #'nwchem_dft_xc': 'b3lyp',
        #add grid options
        'nwchem_tce_dft': True,
        #'nwchem_tce_on' : True,
        'nwchem_tce': 'CCSDT',
        'nwchem_tce_thresh': 1.0e-12,
        'nwchem_task_tce': 'energy'
    })
    #print(jrec['qcvars'])
    print('Testing DFT energy...')
    val = qcdb.energy('nwc-ccsdt')
    check_dft(val, is_dft=True)


@using_nwchem
def test_2_ccsdt():
    qcdb.set_options({
        'basis': 'sto-3g',
        'memory': '600 mb',
        'nwchem_tce_dft': True,
        #'nwchem_tce_on' : True,
        'nwchem_tce': 'CCSDT',
        'nwchem_tce_thresh': 1.0e-12,
        'nwchem_task_tce': 'energy'
    })
    # print(jrec['qcvars'])
    print('Test CCSDT energy ...')
    val = qcdb.energy('nwc-ccsdt')
    check_dft(val, is_dft=False)
