#Multiconfiguration Self-consistent theory on CH2

import os
import sys
from utils import *
from addons import *
import qcdb

ch2 = qcdb.set_molecule('''
          C    0   0     0
          H    0  1.87  -0.82
          ''')

print(ch2)

def check_mcscf(val, is_df):
    if is_df:
        ref     =       -38.857160053058 #hf ref
        nre     =         6.144298248886
        mcscf   =       -38.958740226767
    else:
        ref     =       -38.857160053058 #hf ref
        nre     =         6.144298248886
        mcscf   =       -38.958757123933

    assert compare_values(ref, qcdb.get_variable('HF TOTAL ENERGY'), 5, 'hf ref')
    assert compare_values(nre, qcdb.get_variable('NUCLEAR REPULSION ENERGY'), 5, 'nre')
    assert compare_values(mcscf, qcdb.get_variable('MCSCF TOTAL ENERGY'), 5, 'mcscf')

@using_nwchem
def test_1_hf():
    qcdb.set_options({
        'basis'     :       '6-31g**'
        })
    print('Testing hf...')
    val = qcdb.energy('nwc-hf')
    check_mcscf(val, is_df=True)

@using_nwchem
def test_2_mcscf_3b1():
    qcdb.set_options({
        'basis'     :       '6-31g**',
        'nwchem_mcscf_active'   : 6,
        'nwchem_mcscf_actelec'  : 6,
        'nwchem_mcscf_state'    :'3b1'
        })
    print('Testing MCSCF 3b1...')
    val = qcdb.energy('nwc-mcscf')
    check_mcscf(val, is_df=True)

@using_nwchem
def test_3_mcscf_1a1():
    qcdb.set_options({
        'basis'     :       '6-31g**',
        'nwchem_mcscf_state': '1a1'
        })
    print('Testing MCSCF at 1a1...')
    val = qcdb.energy('nwc-mcscf')
    check_mcscf(val, is_df=False)
