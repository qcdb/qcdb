#! single-point CCSDTQ/6-31g* on water
import os
import sys
import utils
import addons
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
qcdb.set_options({
    'basis': '6-31g',
    'nwchem_total_memory': '2000 mb',
    'nwchem_global_memory': '1700 mb',
    'nwchem_symmetry': 'c2v',
    'nwchem_scf': 'RHF',
    'nwchem_scf_thresh': 1.0e-7,
    #'nwchem_tce_job': 'scf ccsdtq',
    #'nwchem_tce_thresh': 1.0e-7,
    'nwchem_task_tce': 'energy'})
print(qcdb.get_active_options().print_changed()
def check_rhf_ccsdtq(return_value, is_df):
    if is_df:
        ref=
        ccsdtq=
        ccsdtq_corl=
    else:
        print("Does not match")
        
    assert compare_values(-76.010496307065, qcdb.get_variable('scf total energy'), 6, 'SCF')  #TEST
    assert compare_values(-76.210368642101955, qcdb.get_variable('ccsdtq total energy'), 6, 'CCSDTQ')  #TEST
    assert compare_values(-0.199872335037341, qcdb.get_variable('ccsdtq correlation energy'), 6, 'CCSDTQ corl')  #TEST

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

assert compare_values(-76.010496307065, qcdb.get_variable('scf total energy'), 6, 'SCF')  #TEST
assert compare_values(-76.210368642101955, qcdb.get_variable('ccsdtq total energy'), 6, 'CCSDTQ')  #TEST
assert compare_values(-0.199872335037341, qcdb.get_variable('ccsdtq correlation energy'), 6, 'CCSDTQ corl')  #TEST


