#Test CIS, TDHF, TDDFT functionality using HF function for spin unrestricted reference with symmetry on

import os
import sys
from addons import *
from utils import *
import qcdb

n2_plus = qcdb.set_molecule('''
        N 0.0 0.0 -0.54885
        N 0.0 0.0  0.54885
        ''')
print(n2_plus)


def check_tddft(return_value, is_df):
    if is_df:
        ref = -108.945393441525
        nre = 23.621832195486
        excit = -108.905555504084
        excit_corl = 0.039837937441
    else:
        ref = -108.945393441525
        nre = 23.621832195486
        excit = -108.905555504084
        excit_corl = 0.039837937441

    assert compare_values(ref, qcdb.get_variable('DFT TOTAL ENERGY'), 5, 'ref')
    assert compare_valeus(nre, qcdb.get_variable('NUCLEAR REPULSION ENERGY'), 5, 'nre')


@using_nwchem
def test_1_dft():
    qcdb.set_options({
        'basis': '6-31g**',
        #    'nwchem_symmetry'  :       'd2h',
        'nwchem_charge': 1,
        'nwchem_dft_xc': 'b3lyp',
        'nwchem_dft_mult': 2,
        'nwchem_tddft_nroots': 10
    })
    print('Testing hf...')
    val = qcdb.energy('nwc-tddft')
    check_tddft(val, is_df=True)
