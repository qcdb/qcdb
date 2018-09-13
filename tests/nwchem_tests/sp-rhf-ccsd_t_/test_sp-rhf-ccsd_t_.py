#! single-point CCSD(T)/cc-pvdz on water
import os
import sys
import utils
import addons
import qcdb

print('        <<< Literal input to NWChem >>>')
print('''
memory 600 mb

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
 thresh 1.0e-12 #Not the same as e_convergence and d_convergence
end
tce
 scf
 ccsd(t)
 thresh 1.0e-12
end

task tce energy
}

energy('nwchem')

clean()
clean_variables()
nwchem {}
''')

h2o= qcdb.set_molecule('''
        O     0.000000000000    0.000000000000   -0.065638538099
        H     0.000000000000   -0.757480611647    0.520865616174
        H     0.000000000000    0.757480611647    0.520865616174
        ''')

print(h2o)

qcdb.set_options({
    'basis': 'cc-pvdz',
    'memory': '600 mb',
    'nwchem_scf': 'rhf',
    'nwchem_scf_thresh': 1.0e-12,
    'nwchem_tce_dft': False,
    'nwchem_tce': 'ccsd(t)',
    'nwchem_tce_thresh': 1.0e-12,
    'nwchem_task_tce': 'energy'
    })
def check_rhf_ccsd_t_(return_value, is_df):
    if is_df:
        ref         =   -76.026744421192
        ccsdcorl    =   -0.213350416141872
        ccsdtot     =   -76.240094837333771
        ccsd_t_corr =   -0.003062727448805
        ccsd_t_corl =   -0.216413143590677
        ccsd_t_tot  =   -76.243157564782578
    else:
        print("Does not match- Test fail")

    assert compare_values(ref, qcdb.get_variable('SCF TOTAL ENERGY'), 5, 'rhf ref')
    assert compare_values(ccsdcorl, qcdb.get_variable('CCSD CORRELATION ENERGY'), 5, 'ccsd corl')
    assert compare_values(ccsdtot, qcdb.get_variable('CCSD TOTAL ENERGY'), 5, 'ccsd total')
    assert compare_values(ccsd_t_corr, qcdb.get_variable('CCSD(T) CORRECTION ENERGY'), 5, 'ccsd(t) corr')
    assert compare_values(ccsd_t_corl, qcdb.get_variable('CCSD(T) CORRELATION ENERGY'), 5, 'ccsd(t) corl')
    assert compare_values(ccsd_t_tot, qcdb.get_variable('CCSD(T) TOTAL ENERGY'), 5, 'ccsd(t) tot')

#@using_nwchem
def test_1_df_rhf():
    qcdb.set_options({
        'basis': 'cc-pvdz',
        'memory': '600 mb',
        'nwchem_scf': 'rhf',
        'nwchem_scf_thresh': 1.0e-12
        })
    print('     Testing rhf ...')
    val = qcdb.energy('nwchem-scf')
    check_scf(val, is_df=True) #where to find check functions #TODO
def test_2_df_ccsd():
    qcdb.set_options({
        'basis': 'cc-pvdz',
        'memory':'600 mb',
        'nwchem_tce_dft': False,
        'nwchem_tce': 'ccsd(t)',
        'nwchem_tce_thresh': 1.0e-12,
        'nwchem_task_tce': 'energy'
        })
    print('Testing ccsd...')
    val=qcdb.energy('nwchem-tce-ccsd')
    check_ccsd(val, is_df=True)
    print('Testing ccsd(t)')
    check_ccsdt(val, is_df=True)


print( '        <<< Translation of NWChem input to Psi4 format to NWChem >>>')
print('''
molecule {
O
H 1 R
H 1 R 2 A

R=0.958
A=104.5
}

set {
basis cc-pvdz
reference rhf
nwchem_scf_thresh 1.0e-12
nwchem_tce_thresh 1.0e-12
nwchem_tce ccsd(t)
}

energy('nwchem-tce')''')
