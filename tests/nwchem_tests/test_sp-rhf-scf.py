#! single-point HF/cc-pvdz on water
#nopen = 0 is singlet (Default)

import os
import sys
from utils import *
from addons import *
import qcdb

h2o = qcdb.set_molecule('''
        O     0.000000000000    0.000000000000   -0.065638538099
        H     0.000000000000   -0.757480611647    0.520865616174
        H     0.000000000000    0.757480611647    0.520865616174
        ''')
print(h2o)

def check_hf(return_value, is_df):
    if is_df:
        ref         =       -76.026760737428
    else:
        ref         =       -76.026760737428
    
    assert compare_values(ref, qcdb.get_variable('HF TOTAL ENERGY'), 5, 'total scf')

@using_nwchem
def test_1_hf():
    qcdb.set_options({
        'basis'     :   'cc-pvdz',
        'memory'    :   '400 mb',
        'nwchem_scf':   'rhf',
        'nwchem_scf_thresh':    1.0e-8,
        'nwchem_scf_nopen' :    0
        })
    print('Testing hf...')
    val = qcdb.energy('nwc-hf')
    check_hf(val, is_df=True)

print('''        <<< Literal nwchem.nw to NWChem >>>'

memory 400 mb

nwchem {
geometry
O     0.000000000000    0.000000000000   -0.065638538099
H     0.000000000000   -0.757480611647    0.520865616174
H     0.000000000000    0.757480611647    0.520865616174
end

basis spherical
* library cc-pvdz
end

scf
 rhf
 thresh 1.0e-8 #Not the same as e_convergence and d_convergence
 nopen 0 #singlet
end

task scf energy
}

energy('nwchem')
''')

#clean()
#clean_variables()

print('''        <<< Translation of nwchem.nw to Psi4 format to NWChem >>>'

memory 400 mb

molecule {
O
H 1 R
H 1 R 2 A

R=0.958
A=104.5
}

set {
basis cc-pvdz
nwchem_scf_thresh 1.0e-8 
nwchem_scf rhf 
nwchem_scf_nopen 0
}

#To get the same energy from Psi4, add commands below
#set df_scf_guess false #Doing nothing for NWChem
#set scf_type direct  #Doing nothing for NWChem

nwchem {}  # clear literal block since running sequentially

energy('nwchem-scf')

''')

#clean()
#clean_variables()
#nwchem {}

# these make a clean slate for the next job
#revoke_global_option_changed('NWCHEM_SCF_THRESH')
#revoke_global_option_changed('NWCHEM_SCF')
#revoke_global_option_changed('NWCHEM_SCF_NOPEN')

print('''        <<< Thorough Psi4 format >>>'

memory 400 mb

molecule {
O
H 1 R
H 1 R 2 A

R=0.958
A=104.5
}

set basis cc-pvdz
set nwchem_scf_thresh 1.0e-8 
set reference rhf

energy('nwchem-scf')

''')

#clean()
#clean_variables()

print('''        <<< Translation of nwchem.nw to Psi4 format to NWChem >>>'

memory 400 mb

molecule {
O
H 1 R
H 1 R 2 A

R=0.958
A=104.5
}

set {
basis cc-pvdz
nwchem_scf_thresh 1.0e-8
nwchem_scf rhf #or reference rhf
nwchem_scf_nopen 0
}

nwchem {}  # clear literal block since running sequentially

energy('nwchem-scf')

''')
#clean()
#clean_variables()
#nwchem {}

print('''        <<< Translation of nwchem.nw to Psi4 format to NWChem >>>'

molecule {
O
H 1 R
H 1 R 2 A

R=0.958
A=104.5
}

set basis cc-pvdz
set nwchem_scf_thresh 1.0e-8 
set reference rhf

energy('nwchem-scf')

''')
