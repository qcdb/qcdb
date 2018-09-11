#! single-point HF/cc-pVDZ (Cartesian) on NH2 
import os
import sys
import utils
import addons
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
        O
        H 1 0.96
        H 1 0.96 2 104.5
        ''')
print(h2o)
#print(qcdb.get_active_options().print_changed())

qcdb.set_options({'basis': 'cc-pVDZ',
                 'memory': '300 mb',
                 'nwchem_scf': 'ROHF',
                 'nwchem_scf_nopen': 1,
                 'nwchem_scf_thresh': 1.0e-8})

print(qcdb.get_active_options().print_changed())
def check_rohf(return_value, is_df):
    if is_df:
        ref= -76.010538615956
        one= -123.058841737821
        two= 37.851104681667
        nre= 9.197198440198
    else:
        print("Does not match")

    assert compare_values(ref, qcdb.get_variable('Total SCF energy'), 5, 'scf total')
    assert compare_values(one, qcdb.get_variable('One-electron energy'), 5, 'one electron')
    assert compare_values(two, qcdb.get_variable('Two-electron energy'), 5, 'two electron')
    assert compare_values(nre, qcdb.get_variable('Nuclear repulsion energy'), 5, 'nuclear repulsion')

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

qcdb.energy('nwchem-scf')


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

