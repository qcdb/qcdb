#Tensor Contraction Engine (TCE) Moller-Plesset Perturbation Theory third order (MP3/MBPT3) for water

import os
import sys
from utils import *
from addons import *
import qcdb

h2o= qcdb.set_molecule('''
        O      0.000000000000     0.000000000000    -0.123909374404
        H      0.000000000000     1.429936611037     0.983265845431
        H      0.000000000000    -1.429936611037     0.983265845431
        ''')
print(h2o)

def check_tce_mp3(val, is_df):
    if is_df:
        ref     =       -74.962905406171
        nre     =         9.196934380443
        mp3tot  =       -75.007969402203344      
        mp3corl =        -0.009584825225284
    else:
        ref     =       -74.962905406171 #hf ref
        nre     =         9.196934380443
        mp3tot  =       -75.007969402203344      
        mp3corl =        -0.009584825225284

    assert compare_values(ref, qcdb.get_variable('HF TOTAL ENERGY'), 5, 'hf ref')
    assert compare_values(nre, qcdb.get_variable('NUCLEAR REPULSION ENERGY'), 5, 'nre')
    assert compare_values(mp3tot, qcdb.get_variable('MP3 TOTAL ENERGY'), 5, 'tce mp3/mbpt3')
    assert compare_values(mp3corl, qcdb.get_variable('MP3 CORRELATION ENERGY'), 5, 'tce mp3/mbpt3 corl')

@using_nwchem
def test_1_hf():
    qcdb.set_options({
        'basis'     :       'sto-3g',
        'memory'    :       '1500 mb',
        'nwchem_scf_thresh': 1.0e-10,
        'nwchem_scf':       'rhf',
#        'nwchem_tol2e':     1.0e-10,
        'nwchem_tce':       'MP3',
#        'nwchem_tce_on':    True,
        'nwchem_tce_dft':   False,
#        'task_tce'  :       'energy'
        })
    print('Testing hf...')
    val = qcdb.energy('nwc-mp3')
    check_tce_mp2(val, is_df = True)

@using_nwchem
def test_2_tce_mp2():
    qcdb.set_options({
        'basis'     :       'sto-3g',
        'memory'    :       '1500 mb',
        'nwchem_scf_thresh': 1.0e-10,
        'nwchem_scf':       'rhf',
#        'nwchem_tol2e':     1.0e-10,
        'nwchem_tce':       'MP3',
#        'nwchem_tce_on':    True,
        'nwchem_tce_dft':   False,
#        'task_tce'  :       'energy'
        })
    print('Testing tce-mp3...')
    val = qcdb.energy('nwc-mp3')
    check_tce_mp2(val, is_df = True)
