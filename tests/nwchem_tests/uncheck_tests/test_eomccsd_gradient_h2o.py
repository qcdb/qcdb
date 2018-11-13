#equation of motion couple clustered singles and doubles gradient task through the Tensor Contraction Engine

import os
import sys
from utils import *
from addons import *
import qcdb

h2o = qcdb.set_molecule('''
        O      0.000000000000     0.000000000000    -0.123909374404
        H      0.000000000000     1.429936611037     0.983265845431
        H      0.000000000000    -1.429936611037     0.983265845431
        ''')
print(h2o)

def check_eomccsd(return_value, is_df):
    if is_df:
        ref         =       -74.962905406162 #hf ref
        nre         =         9.196934380443
        ccsd_corl   =        -0.049337013314511
        ccsd_tot    =       -75.012242419476792
        eom1_ht     =         0.457480376862104
        eom1_ev     =        12.448679742772001
    else:
        ref         =       -74.963623744855 #last iteration
        nre         =         9.160284701624
        ccsd_corl   =       -0.049944017032077
        ccsd_tot    =       -75.013567761886847
        eom1_ht     =         0.454242037050082
        eom1_ev     =        12.360559995440656

    assert compare_values(ref, qcdb.get_variable('HF TOTAL ENERGY'), 5, 'ref')
    assert compare_values(nre, qcdb.get_variable('NUCLEAR REPULSION ENERGY'), 5, 'nre')
    assert compare_values(ccsd_corl, qcdb.get_variable('CCSD CORRELATION ENERGY'), 5, 'ccsd corl')
    assert compare_values(ccsd_tot, qcdb.get_variable('CCSD TOTAL ENERGY'), 5, 'ccsd tot')
    assert compare_values(eom1_ht, qcdb.get_variable('EOM-CCSD ROOT 0 EXCITATION ENERGY'), 5, 'eom-ccsd')

@using_nwchem
def test_1_eomccsd():
    qcdb.set_options({
        'basis'     :       'sto-3g',
        'memory'    :       '1500 mb',
        'nwchem_scf_thresh': 1.0e-10,
        'nwchem_scf_tol2e' : 1.0e-10,
        'nwchem_scf':   'RHF',
        'nwchem_tce':   'CCSD',
        'nwchem_tce_nroots': 1,
        'nwchem_task_tce':  'gradient'
        })
    print('Testing...')
    val = qcdb.gradient('nwc-eom-ccsd')
    check_eomccsd(val, is_df=True)
@using_nwchem
def test_2_eomccsd(): 
    qcdb.set_options({
        'basis'     :       'sto-3g',
        'memory'    :       '1500 mb',
        'nwchem_scf_thresh': 1.0e-10,
        'nwchem_scf_tol2e' : 1.0e-10,
        'nwchem_scf':   'RHF',
        'nwchem_tce':   'CCSD',
        'nwchem_tce_nroots': 1,
        'nwchem_task_tce':  'gradient'
        })
    print('Testing...')
    val = qcdb.gradient('nwc-eom-ccsd')
    check_eomccsd(val, is_df=False)

