import os
import sys
from ..addons import *
from ..utils import *
import qcdb
import numpy as np

h2o = qcdb.set_molecule('''
        O      0.000000000000     0.000000000000    -0.123909374404
        H      0.000000000000     1.429936611037     0.983265845431
        H      0.000000000000    -1.429936611037     0.983265845431
        ''')

print(h2o)

def check_ccsd(return_value, is_df):
    if is_df:
        hf          =   -74.963623744947 #tce dft
        ccsd_corl   =   -0.049337013315753
        ccsd_tot    =   -75.012242419560370
        ccsd_grad   =   np.array[['0.000000,   0.000017,   0.111908'],
                        [0.000000,  -0.047847,  -0.055961],
                        [0.000000,   0.047830,  -0.055947]]
    else:
        hf          =   -74.962905406245 #scf
        ccsd_corl   =   -0.049337013315753
        ccsd_tot    =   -75.012242419560370
        ccsd_grad   =   np.array[['0.000000,   0.000017,   0.111908'],
                        [0.000000,  -0.047847,  -0.055961],
                        [0.000000,   0.047830,  -0.055947]]

    assert compare_values(hf, qcdb.get_variable('HF TOTAL ENERGY'), 5, 'hf ref')
    assert compare_values(ccsd_corl, qcdb.get_variable('CCSD CORRELATION ENERGY'), 5, 'ccsd corl')
    assert compare_values(ccsd_tot, qcdb.get_variable('CCSD TOTAL ENERGY'), 5, 'ccsd tot')
    assert compare_arrays(ccsd_grad, qcdb.get_variable('CURRENT GRADIENT'), 5, 'ccsd grad')

@using_nwchem
def test_1_ccsd_df():
    qcdb.set_options({
        'basis'     :   'sto-3g',
        'memory'    :   '1500 mb',
        'qc_module' :   'tce',
        #'scf__e_convergence': 1.0e-10,
        'nwchem_dft__rdft':   True,
        #'nwchem_scf__thresh'     :   1.0e-10,
        #'nwchem_scf__tol2e'      :   1.0e-10,
        'nwchem_tce__ccsd'     : True
        })
    print('Testing ccsd(dft)...')
    val = qcdb.gradient('nwc-ccsd')
    check_ccsd(val, is_df=True)

@using_nwchem
def test_2_ccsd_df():
    qcdb.set_options({
        'basis'     :   'sto-3g',
        'memory'    :   '1500 mb',
        'qc_module' :   'tce',
        'scf__e_convergence': 1.0e-10,
        'nwchem_scf__rhf':   True,
        #'nwchem_scf__thresh'     :   1.0e-10,
        'nwchem_scf__tol2e'      :   1.0e-10,
        'nwchem_tce__ccsd'     : True
        })
    print('Testing ccsd(hf)...')
    val = qcdb.gradient('nwc-ccsd')
    check_ccsd(val, is_df=False)
