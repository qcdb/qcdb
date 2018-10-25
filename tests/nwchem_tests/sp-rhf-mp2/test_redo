#! single-point MP2/cc-pvdz on water
import os
import sys
import utils
import addons
import qcdb

print("        <<< Literal nwchem.nw to NWChem >>>")
print("""
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
 thresh 1.0e-8 
 nopen 0 
end

mp2
tight
end

task mp2 energy
}

energy('nwchem')


clean()
clean_variables()
""")
h2o= qcdb.set_molecule('''
        O     0.000000000000    0.000000000000   -0.065638538099
        H     0.000000000000   -0.757480611647    0.520865616174
        H     0.000000000000    0.757480611647    0.520865616174
        ''')
print(h2o)
qcdb.set_options({
    'basis': 'cc-pvdz',
    'memory': '400 mb',
    'nwchem_scf': 'RHF',
    'nwchem_scf_thresh': '1.0e-8',
    'nwchem_scf_nopen': 0,
    'nwchem_mp2_tight': True,
    'nwchem_task_mp2': 'energy'
    })
print(qcdb.get_active_options().print_changed())
def check_mp2(return_value, is_df):
    if is_df:
        ref=-76.026760737428
        one=-123.137568018334
        two=37.923473040741
        mp2_total=-76.230777733733
        scs_mp2_total=-76.226922314540
    else:
        ref=-76.026760737428
        one=-123.137568018334
        two=37.923473040741
        mp2_total=-76.230777733733
        scs_mp2_total=-76.226922314540

    assert compare_values(ref, qcdb.get_variable('Total SCF energy'), 5, 'scf total')
    assert compare_values(one, qcdb.get_variable('One-electron energy'), 5, 'one electron scf energy')
    assert compare_values(two, qcdb.get_variable('Two-electron energy'), 5, 'two electron scf energy')
    assert compare_values(mp2_total, qcdb.get_variable('Total MP2 energy'), 5, 'mp2 energy')
    assert compare_values(scs_mp2_total, qcdb.get_variable('Total SCS-MP2 energy'), 5, 'scs-mp2 energy')

print("        <<< Translation of nwchem.nw to Psi4 format to NWChem >>>")
print("""
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
nwchem_mp2_tight true
}

nwchem {}  # clear literal block since running sequentially

energy('nwchem-mp2')


clean()
clean_variables()
nwchem {}
""")
# these make a clean slate for the next job
#revoke_global_option_changed('NWCHEM_SCF_THRESH')
#revoke_global_option_changed('NWCHEM_SCF')
#revoke_global_option_changed('NWCHEM_SCF_NOPEN')
#revoke_global_option_changed('NWCHEM_MP2_TIGHT')

print("      <<< Thorough Psi4 format >>>")
print("""
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
set nwchem_mp2_tight true 

energy('nwchem-mp2')


clean()
clean_variables()
""")

