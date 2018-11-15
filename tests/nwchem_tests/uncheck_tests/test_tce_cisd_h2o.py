#Tensor Contraction Engine on Configuration Interaction on singles and doubles on water

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

def check_tce_cisd(return_value, is_df):
    if is_df:
        ref         =       -74.962905406171 #hf ref
        nre         =         9.196934380443
        cisd_tot    =       -75.011657140348433
        cisd_corl   =        -0.048751734177409
    else:
        ref         =       -74.962905406171 #hf ref
        nre         =         9.196934380443
        cisd_tot    =       -75.011657140348433
        cisd_corl   =        -0.048751734177409

    assert compare_values(ref, qcdb.get_variable('HF TOTAL ENERGY'), 5, 'hf ref')
    assert compare_values(nre, qcdb.get_variable('NUCLEAR REPULSION ENERGY'), 5, 'nre')
    assert compare_values(cisd_tot, qcdb.get_variable('CISD TOTAL ENERGY'), 5, 'cisd tot')
    assert compare_values(cisd_corl, qcdb.get_variable('CISD CORRELATION ENERGY'), 5, 'cisd corl')

@using_nwchem
def test_1_hf():
    qcdb.set_options({
        'memory'        :       '1500 mb',
        'basis'         :       'sto-3g',
        'nwchem_scf'    :       'RHF',
        'nwchem_scf_thresh':    1.0e-10,
        #'nwchem_scf_tol2e':    1.0e-10
        })
    print('Testing hf...')
    val = qcdb.energy('nwc-tce')
    check_tce_cisd(val, is_df=True)

@using_nwchem
def test_2_cisd():
    qcdb.set_options({
        'memory'        :       '1500 mb',
        'basis'         :       'sto-3g',
        'nwchem_scf'    :       'RHF',
        'nwchem_scf_thresh':    1.0e-10,
        #'nwchem_scf_tol2e':    1.0e-10,
        'nwchem_tce'    :       'CISD',
        #'nwchem_task_tce':     'energy'
        })
    print('Testing CISD...')
    val = qcdb.energy('nwc-tce')
    check_tce_cisd(val, is_df=True)
