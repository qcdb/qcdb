#! pbe/sto-3g H2O DFT energy
import os
import sys
from ..utils import *
from ..addons import *
import qcdb

h2o = qcdb.set_molecule('''
        O     0.000000000000    0.000000000000   -0.068516219310
        H     0.000000000000   -0.790689573744    0.543701060724
        H     0.000000000000    0.790689573744    0.543701060724
        ''')
print(h2o)


def check_dft(return_value, is_df):
    if is_df:
        ref = -75.234018772521
    else:
        ref = -74.964662543238

    assert compare_values(ref, qcdb.get_variable('DFT TOTAL ENERGY'), 5, 'DFT total')  #TEST
    assert compare_values(ref, qcdb.get_variable('HF TOTAL ENERGY'), 5, 'HF total')  #TEST


@using_nwchem
def test_1_dft_tot():
    qcdb.set_options({
        'basis': 'sto-3g',
        'memory': '3000 mb',
        'e_convergence': 1.0e-7,
        'scf__d_convergence': 1.0e-7,
        'nwchem_charge': 0,
        'nwchem_dft__direct': True,
        'nwchem_dft__xc_pbe0': True,
        #'nwchem_dft__convergence__energy': 1.0e-7,
        #'nwchem_dft__convergence__density': 1.0e-7,
        #'nwchem_task_dft':'energy'
    })
    print('Testing DFT total energy...')
    val = qcdb.energy('nwc-pbe')
    check_dft(val, is_df=True)


@using_nwchem
def test_2_dft(): #options okay here, pulling scf hf energy
    qcdb.set_options({
        'basis': 'sto-3g',
        'memory': '3000 mb',
        'e_convergence': 1.0e-7,
        'scf__d_convergence': 1.0e-7,
        'nwchem_charge': 0,
        #'nwchem_dft__direct': True,
        #'nwchem_dft__convergence__energy': 1.0e-7,
        #'nwchem_dft__convergence__density': 1.0e-7
    })
    print('Testing HF energy...')
    val2 = qcdb.energy('nwc-hf')
    check_dft(val2, is_df=False)
