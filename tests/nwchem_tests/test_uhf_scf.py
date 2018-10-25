#! single-point UHF/cc-pVDZ  on NH2 
import os
import sys
from addons import *
from utils import *
import qcdb
nh2= qcdb.set_molecule('''
           N        0.08546       -0.00020       -0.05091
           H       -0.25454       -0.62639        0.67895
           H       -0.25454       -0.31918       -0.95813
           ''')
print(nh2)
def check_uhf_hf(return_value, is_df=True):
    if is_df:
        ref     =       -55.566057523877
    else:
        print("Does not match")

    assert compare_values(ref, qcdb.get_variable('HF TOTAL ENERGY'), 2, 'scf')

@using_nwchem
def test_1_hf():
    qcdb.set_options({
        'basis'     : 'cc-pvdz',
        'memory'    : '400 mb',
        'nwchem_scf': 'uhf',
        'nwchem_scf_nopen': 1,
        'nwchem_scf_thresh': 1.0e-8
        })
    print('Testing hf...')
    val = qcdb.energy('nwc-hf')
    check_uhf_hf(val, is_df=True)
print('''        <<< Literal nwchem.nw to NWChem >>>'

nwchem {
memory 400 mb 

geometry
zmat
 N
 H 1 1.008
 H 1 1.008 2 105.0
end
end

basis spherical
* library cc-pVDZ
end

scf
uhf
nopen 1
thresh 1.0e-8
end

task scf energy
}

energy('nwchem')
''')

#clean()
#clean_variables()
#nwchem {}

print('''        <<< Translation of nwchem.nw to Psi4 format to NWChem >>>'

banner('UHF-SCF energy calculation')

molecule {
0 2
N
H 1 R
H 1 R 2 A

R=1.008
A=105.0
}

set {
basis cc-pVDZ
nwchem_scf uhf
nwchem_scf_nopen 1
nwchem_scf_thresh 1.0e-8
#To get the same energy from Psi4, add commands below
#df_scf_guess false #Doing nothing for NWChem
#scf_type direct  #Doing nothing for NWChem
}

energy('nwchem-scf')
''')

#clean()
#clean_variables()
#nwchem {}

print('        <<< Translation of nwchem.nw to Psi4 format to NWChem >>>')

#banner('UHF-SCF energy calculation')
print('''
molecule {
0 2
N
H 1 R
H 1 R 2 A

R=1.008
A=105.0
}

set {
basis cc-pvdz
reference uhf 
nwchem_scf_nopen 1
nwchem_scf_thresh 1.0e-8
}

energy('nwchem-scf')

''')
