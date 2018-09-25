#! gradient HF/cc-pvdz on water

import os
import sys
import utils
import addons
import qcdb

print('''        <<< Literal nwchem.nw to NWChem >>>'

nwchem {
memory 400 mb
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
 thresh 1.0e-8 
 nopen 0 #singlet
end

task scf gradient

}

gradient('nwchem')
''')
h2o = qcdb.set_molecule('''
        O     0.000000000000    0.000000000000   -0.065638538099
        H     0.000000000000   -0.757480611647    0.520865616174
        H     0.000000000000    0.757480611647    0.520865616174
        ''')
print(h2o)

def check_hf(return_value, is_df):
    if is_df:
        ref     =       -76.026760737428
    else:
        ref     =       -76.026760737428

    assert compare_values(ref, qcdb.get_variable('TOTAL SCF ENERGY'), 5, 'scf')

#@using_nwchem:
def test_1_hf():
    qcdb.set_options({
        'basis'     :   'cc-pvdz',
        'memory'    :   '400 mb',
        'nwchem_scf':   'rhf',
        'nwchem_scf_thresh': 1.0e-8,
        'nwchem_scf_nopen' : 0,
        })
    print('Testing hf...')
    val = qcdb.energy('nwc-hf')
    check_hf(val, is_df=True)

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

gradient('nwchem-scf')
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
set reference rhf
set nwchem_scf_thresh 1.0e-8

gradient('nwchem-scf')

clean()
clean_variables()
''')
