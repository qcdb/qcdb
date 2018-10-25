#! Geometry optimization HF/6-31g* on water  
import os
import sys
from utils import *
from addons import *
import qcdb

print('        <<< Literal nwchem input to NWChem >>>')

print('''
memory 400 mb 

nwchem {
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
 thresh 1.0e-8 
 direct 
end

task scf optimize
}

energy('nwchem')

clean()
clean_variables()
nwchem {}  ''') 

h2o=qcdb.set_molecule('''
        O     0.000000000000    0.000000000000   -0.065638538099
        H     0.000000000000   -0.757480611647    0.520865616174
        H     0.000000000000    0.757480611647    0.520865616174
        ''')
print(h2o)

def check_rhf(return_value, is_df):
    if is_df:
        ref     =   -76.010496306999
        nre     =     9.187334240165
    else:
        ref     =   -76.010738270124
        nre     =     9.347020370478
    assert compare_values(ref, qcdb.get_variable('HF TOTAL ENERGY'), 5, 'hf ref')
    assert compare_values(nre, qcdb.get_variable('NUCLEAR REPULSION ENERGY'), 5, 'nre')


@using_nwchem
def test_1_hf():
    qcdb.set_options({
        'basis'     :   '6-31g*',
        'memory'    :   '400 mb',
        #'nwchem_geometry_center'    : False,
        #'nwchem_geometry_autosym'   : False,
        #'symmetry'   :   'c2v',
        'nwchem_scf'        :   'RHF',
        'nwchem_scf_thresh' :   1.0e-8,
        'nwchem_scf_direct' :   True,
        'nwchem_task_scf'   :   'optimize'
        })
    print('Testing HF...')
    val = qcdb.energy('nwc-hf')
    check_rhf(val, is_df=True)

print('        <<< Translation of nwchem.nw to Psi4 format to NWChem >>>')
print('''
memory 400 mb

molecule {
O     0.000000000000    0.000000000000   -0.065638538099
H     0.000000000000   -0.757480611647    0.520865616174
H     0.000000000000    0.757480611647    0.520865616174
symmetry c2v
}

set {
basis 6-31g*
nwchem_scf_thresh 1.0e-8 
nwchem_scf rhf 
}

#To get the same energy from Psi4, add commands below
set df_scf_guess false #Doing nothing for NWChem
set scf_type direct  
#NWChem default convergence criterian 
set g_convergence GAU

nwchem {}  # clear literal block since running sequentially

optimize('nwchem-scf')

clean()
clean_variables()
nwchem {}
''')
