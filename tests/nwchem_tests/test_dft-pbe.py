#! pbe/sto-3g H2O DFT energy  
import os
import sys
from utils import *
from addons import *
import qcdb

print('        <<< Literal nwchem.nw to NWChem >>>')


print('''nwchem {
memory 300 mb
geometry units Angstrom
O     0.000000000000    0.000000000000   -0.068516219310
H     0.000000000000   -0.790689573744    0.543701060724
H     0.000000000000    0.790689573744    0.543701060724
end
basis spherical
 * library sto-3g
end
charge 0

dft
direct
xc xpbe96 cpbe96
grid lebedev 99 14
convergence energy 1e-7 density 1e-7
end

task dft energy
}

energy ('nwchem')
''')
h2o= qcdb.set_molecule('''
        O     0.000000000000    0.000000000000   -0.068516219310
        H     0.000000000000   -0.790689573744    0.543701060724
        H     0.000000000000    0.790689573744    0.543701060724
        ''')
print(h2o)

def check_dft(return_value, is_df):
    if is_df:
        ref     =   -75.234018772521
    else:
        ref     =   -74.964662543238
        
    assert compare_values(ref, qcdb.get_variable('DFT TOTAL ENERGY'), 5, 'DFT total')  #TEST
    assert compare_values(ref, qcdb.get_variable('HF TOTAL ENERGY'), 5, 'SCF total')  #TEST

@using_nwchem
def test_1_dft_tot():
    qcdb.set_options({
        'basis'         : 'sto-3g',
        'memory'        : '300 mb',
        'nwchem_charge' : 0,
        'nwchem_dft_direct': True,
        #'nwchem_dft__convergence__energy': 1.0e-7,
        #'nwchem_dft__convergence__density': 1.0e-7,
        #'nwchem_task_dft':'energy'
        })
    print('Testing DFT total energy...')
    val = qcdb.energy('nwc-dft')
    check_dft(val, is_df=True)
    
def test_2_dft():
    qcdb.set_options({
        'basis'         : 'sto-3g',
        'memory'        : '300 mb',
        'nwchem_charge' : 0,
        'nwchem_dft_direct': True,
        #'nwchem_dft__convergence__energy': 1.0e-7,
        #'nwchem_dft__convergence__density': 1.0e-7
        })
    print('Testing DFT current energy...')
    val2 = qcdb.energy('nwc-hf')
    check_dft(val, is_df=True)

#nwchem {}
#clean ()
#clean_variables()

print('        <<< Thorough Psi4 format >>>')

print('''
memory 300 mb

molecule h2o {
0 1
O
H 1 1.0
H 1 1.0 2 104.5
}

set {
print 2
basis sto-3g

df_scf_guess false #Doing nothing 
scf_type direct
e_convergence 7
d_convergence 7
dft_spherical_points 590
dft_radial_points 99
reference rks
}

set dft_functional pbe

nwchem {}

energy ('nwchem-dft')
''')
#compare_values(-75.234018772562, get_variable('dft total energy'), 5, 'DFT')  #TEST
#compare_values(-75.234018772562, get_variable('current energy'), 5, 'DFT')  #TEST

#clean()
#clean_variables()
#nwchem {}

# these make a clean slate for the next job
#revoke_global_option_changed('NWCHEM_DFT_DIRECT')
#revoke_global_option_changed('NWCHEM_DFT')
#revoke_global_option_changed('NWCHEM_DFT_XC')
#revoke_global_option_changed('NWCHEM_DFT_GRID')
#revoke_global_option_changed('NWCHEM_DFT_CONVERGENCE')

print('        <<< Translation of nwchem.nw to Psi4 format to NWChem >>>')

print('''
memory 300 mb

molecule h2o {
0 1
O
H 1 1.0
H 1 1.0 2 104.5
}

set {
print 2
basis sto-3g

df_scf_guess false
nwchem_dft_direct true
nwchem_dft_grid [lebedev, 99, 14]
nwchem_dft_convergence [energy, 1.e-7, density, 1.e-7]
}

set nwchem_dft_xc [xpbe96, cpbe96]

nwchem {}

energy ('nwchem-dft')
''')
#compare_values(-75.234018772562, get_variable('dft total energy'), 5, 'DFT')  #TEST
#compare_values(-75.234018772562, get_variable('current energy'), 5, 'DFT')  #TEST
