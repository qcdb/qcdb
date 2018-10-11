#! single-point CCSDTQ/6-31g* on water
import os
import sys
from utils import *
from addons import *
import qcdb


print( '        <<< Literal input to NWChem >>>')

print('''nwchem {
memory total 2000 global 1700 mb verify
geometry nocenter noautosym
O     0.000000000000    0.000000000000   -0.065638538099
H     0.000000000000   -0.757480611647    0.520865616174
H     0.000000000000    0.757480611647    0.520865616174
symmetry c2v
end

basis
* library 6-31g*
end

scf
 rhf
 thresh 1.0e-7 #Not the same as e_convergence and d_convergence
end
tce
 scf
 ccsdtq
 thresh 1.0e-7
end

task tce energy
}

energy('nwchem')
''')
h2o= qcdb.set_molecule('''
        O 0.000000000000    0.000000000000   -0.065638538099
        H 0.000000000000   -0.757480611647    0.520865616174
        H 0.000000000000    0.757480611647    0.520865616174
        ''')
print(h2o)

def check_ccsdtq(return_value, is_df):
    if is_df:
        ref         =       -76.010496307018
        ccsdtq      =       -76.210368642101955
        ccsdtq_corl =        -0.199872335037341
    else:
        print("Does not match")
        
    assert compare_values(ref, qcdb.get_variable('SCF TOTAL ENERGY'), 6, 'SCF')  #TEST
    assert compare_values(ccsdtq, qcdb.get_variable('CCSDTQ TOTAL ENERGY'), 6, 'CCSDTQ')  #TEST
    assert compare_values(ccsdtq_corl, qcdb.get_variable('CCSDTQ CORRELATION ENERGY'), 6, 'CCSDTQ corl')  #TEST

@using_nwchem
def test_1_hf():
    qcdb.set_options({
        'basis': '6-31g*',
        'memory': '2000 mb',
        #'nwchem_total_memory': '2000 mb',
        #'nwchem_global_memory': '1700 mb',
        #'nwchem_symmetry': 'c2v',
        'nwchem_scf': 'RHF',
        'nwchem_scf_thresh': 1.0e-7
        })
    print("Testing SCF energy (df)...")
    val = qcdb.energy('nwc-hf')
    check_ccsdtq(val, is_df=True)

def test_2_ccsdtq():
    qcdb.set_options({
        'basis': '6-31g*',
        'memory': '2000 mb',
        #'nwchem_total_memory': '2000 mb',
        #'nwchem_global_memory': '1700 mb',
        #'nwchem_symmetry': 'c2v',
        'nwchem_tce_dft': False,
        'nwchem_tce': 'CCSDTQ',
        'nwchem_tce_thresh': 1.0e-7
        })
    print('Testing CCSDTQ (df)...')
    val = qcdb.energy('nwc-ccsdtq')
    check_ccstdq(val, is_df=True)

#clean()
#clean_variables()
#nwchem {}

print('        <<< Translation of NWChem input to Psi4 format to NWChem >>>')
print('''
molecule {
O
H 1 R
H 1 R 2 A

R=0.958
A=104.5
symmetry c2v
}

set {
basis 6-31g*
reference rhf
nwchem_memory [total,2000,global,1700,mb,verify]
nwchem_scf_thresh 1.0e-7
nwchem_tce_thresh 1.0e-7
nwchem_tce ccsdtq
}

energy('nwchem-tce')
''')
