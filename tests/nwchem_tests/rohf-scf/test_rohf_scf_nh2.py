#! single-point HF/cc-pVDZ (Cartesian) on NH2 
import os
import sys
from utils import *
from addons import *
import qcdb

print('        <<< Literal nwchem.nw to NWChem >>>')
print('''
memory 300 mb

nwchem {
geometry
zmat
 O
 H 1 0.96
 H 1 0.96 2 104.5
end

basis cartesian
* library cc-pVDZ
end

scf
rohf
nopen 1
thresh 1.0e-8
end

task scf energy
}

energy('nwchem')
''')
h2o= qcdb.set_molecule('''
        O     0.000000000000    0.000000000000   -0.068516219310
        H     0.000000000000   -0.790689573744    0.543701060724
        H     0.000000000000    0.790689573744    0.543701060724
        ''')
print(h2o)

def check_hf(return_value, is_df):
    if is_df:
        ref= -76.010538615956
    else:
        assert compare_values(ref, qcdb.get_variable('TOTAL SCF ENERGY'), 5, 'scf total') 
        
@using_nwchem
def test_1_rohf():
    qcdb.set_options({
        'basis': 'cc-pVDZ',
        'memory': '300 mb',
        'nwchem_scf': 'ROHF',
        'nwchem_scf_nopen': 1,
        'nwchem_scf_thresh': 1.0e-8})
    print('Testing HF energy ...')
    val = qcdb.energy('nwc-hf')
    check_hf(val, is_df=True)
#clean()
#clean_variables()
#nwchem {}

print('        <<< Translation of nwchem.nw to Psi4 format to NWChem >>>')

#banner('ROHF energy calculation')
print('''
molecule {
0 2
O
H 1 R
H 1 R 2 A

R=0.96
A=104.5
}

set {
basis cc-pVDZ
puream false
nwchem_scf rohf
nwchem_scf_thresh 1.0e-8
}
''')
#To get the same energy from Psi4, add commands below
#df_scf_guess false #Doing nothing for NWChem
#scf_type direct 
#}


#clean()
#clean_variables()
#nwchem {}

print( '        <<< Thorough Psi4 format >>>')

#banner('ROHF energy calculation')
print('''
molecule {
0 2
O
H 1 R
H 1 R 2 A

R=0.96
A=104.5
}

set {
basis cc-pvdz
puream false
reference rohf 
e_convergence 8
scf_type direct
}

energy('nwchem-scf')

''') #energy might not be needed in print here

