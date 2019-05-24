# single-point DFT-CCSDT/sto-3g on water

import os
import sys
from ..utils import *
from ..addons import *
import qcdb

h2o = qcdb.set_molecule('''
        O     0.000000000000    0.000000000000   -0.065638538099
        H     0.000000000000   -0.757480611647    0.520865616174
        H     0.000000000000    0.757480611647    0.520865616174
        ''')

print(h2o)


def check_dft(return_value, is_dft):
    if is_dft:
        ref = -95.217126584667
        nre = 9.187334240165
        ccsdt_tot = -95.263267683066744
        ccsdt_corl = -0.046141098399779
        assert compare_values(ref, qcdb.get_variable('DFT TOTAL ENERGY'), 5, 'dft ref')
    else: #change values for non-dft
        ref = -74.963048525888
        ccsdt_tot = -75.012605624645104
        ccsdt_corl = -0.049557089026546
        nre = 9.187333574703
        assert compare_values(ref, qcdb.get_variable('HF TOTAL ENERGY'), 5, 'dft ref')
    assert compare_values(ccsdt_tot, qcdb.get_variable('CCSDT TOTAL ENERGY'), 5, 'ccsdt total')
    assert compare_values(ccsdt_corl, qcdb.get_variable('CCSDT CORRELATION ENERGY'), 5, 'ccsdt corl')
    assert compare_values(nre, qcdb.get_variable('NUCLEAR REPULSION ENERGY'), 5,  'nre')


@using_nwchem
def test_1_dft():
    qcdb.set_options({
        'memory': '6000 mb',
        'basis': 'sto-3g',
        #'nwchem_dft__convergence__density': 1.0e-12,
        #'nwchem_dft__xc_b3lyp': True,
        #add grid options
        'nwchem_tce__dft': True,
        'qc_module' : 'tce',
        'nwchem_tce__ccsdt': True,
        'nwchem_tce__thresh': 1.0e-12,
    })
    print('Testing CCSDT-DFT energy...')
    val = qcdb.energy('nwc-ccsdt')
    check_dft(val, is_dft=True)


@using_nwchem
def test_2_scf():
    qcdb.set_options({
        'basis': 'sto-3g',
        'memory': '6000 mb',
        'nwchem_tce__scf': True,
        'qc_module' :   'tce',
        'nwchem_tce__ccsdt': True,
        'nwchem_tce__thresh': 1.0e-12,
    })
    print('Test CCSDT-SCF energy ...')
    val = qcdb.energy('nwc-ccsdt')
    check_dft(val, is_dft=False)
