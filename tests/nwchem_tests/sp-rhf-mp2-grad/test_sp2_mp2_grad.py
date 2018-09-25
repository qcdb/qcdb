#! Gradient MP2/cc-pvdz on water
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
 nopen 0 
end

mp2
tight
end

task mp2 gradient 
}

gradient('nwchem')


clean()
clean_variables()''')

h2o = qcdb.set_molecule('''
        O     0.000000000000    0.000000000000   -0.065638538099
        H     0.000000000000   -0.757480611647    0.520865616174
        H     0.000000000000    0.757480611647    0.520865616174
        ''')
print(h2o)
def check_mp2(return_value, is_df, a5050=True):
    if is_df:
        ref         =       -76.026760737428
        mp2_tot     =       -76.230777733733
        mp2_corl    =        -0.204016996305
        scs_corl    =        -0.200161577112
        scs_tot     =       -76.226922314540
        ss_corl     =        -0.051529405908
        os_corl     =        -0.152487590397
        #a5050       =       
    else:
        ref         =       -76.026760737428
        mp2_tot     =       -76.226922314540
        mp2_corl    =        -0.204016996305

    assert compare_values(ref, qcdb.get_variable('TOTAL SCF ENERGY'), 5, 'scf')
    assert compare_values(mp2_tot, qcdb.get_variable('TOTAL MP2 ENERGY'), 5, 'mp2 tot')
    assert compare_values(mp2_corl, qcdb.get_variable('MP2 CORRELATION ENERGY'), 5, 'mp2 corl')
    assert compare_values(scs_tot, qcdb.get_variable('SCS-MP2 TOTAL ENERGY'), 5, 'scs mp2 tot')
    assert compare_values(scs_corl, qcdb.get_variable('SCS-MP2 CORRELATION ENERGY'), 5, 'scs mp2 corl')
    assert compare_values(a5050, qcdb.get_variable('custom SCS-MP2 TOTAL ENERGY'), 5, 'custom scs mp2')

#@using_nwchem:
def test_1_mp2():
    qcdb.set_options({
        'basis'         :       'cc-pvdz',
        'memory'        :       '400 mb',
        'nwchem_scf'    :       'rhf',
        'nwchem_scf_thresh':    1.0e-4,
        'nwchem_scf_nopen':     0,
        'nwchem_mp2_tight':     True,
        'nwchem_task_mp2':      'gradient'
        })
    print('Testing mp2 ...')
    val = qcdb.gradient('nwc-mp2')
    check_mp2_tot(val, is_df=True)

def test_2_hf():
    qcdb.set_options({
        'basis'         :       'cc-pvdz',
        'memory'        :       '400 mb',
        'nwchem_scf'    :       'rhf',
        'nwchem_scf_thresh' :   1.0e-4,
        'nwchem_scf_nopen':     0
        })
    print('Testing hf...')
    val = qcdb.energy('nwc-hf')
    check_hf(val, is_df=True)
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
nwchem_mp2_tight true
}

nwchem {}  # clear literal block since running sequentially

gradient('nwchem-mp2')


clean()
clean_variables()
nwchem {}
''')
# these make a clean slate for the next job
#revoke_global_option_changed('NWCHEM_SCF_THRESH')
#revoke_global_option_changed('NWCHEM_SCF')
#revoke_global_option_changed('NWCHEM_SCF_NOPEN')
#revoke_global_option_changed('NWCHEM_MP2_TIGHT')

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
set nwchem_mp2_tight true 

gradient('nwchem-mp2')


clean()
clean_variables()
''')
