#! single-point CCSD(T)/cc-pvdz on water
# CCSD(T) calculation through "task ccsd(t)" is cheaper than "task tce" 
# and able to harvest MP2 energy as well. 
import os
import sys
import utils
import addons
import qcdb

print( '        <<< Literal input to NWChem >>>' )
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

task ccsd(t) energy
}
''')

h2o=qcdb.set_molecule('''
        O     0.000000000000    0.000000000000   -0.065638538099
        H     0.000000000000   -0.757480611647    0.520865616174
        H     0.000000000000    0.757480611647    0.520865616174
        ''')
print(h2o)
def check_ccsd_t_2(return_value, is_df, is_5050=False):
    if is_df:
        ref         =   -76.026760737428
        ccsd_corl   =    -0.213341272556805
        ccsd_tot    =   -76.240102009984767
        scscorl     =    -0.264498694126312
        scstot      =   -76.291259431554266
        mp2corl     =    -0.204016996303923
        mp2tot      =   -76.230777733731884
        ss          =    -0.046033728720216
        os          =    -0.167307543836588
        ccsd_t_corl =    -0.003059541651480
        ccsd_t_tot  =   -76.243161551636248
        a5050corl   =   0.5*(ss + os)
        a5050tot    =   a5050corl + ref
    else:
        ref         =   -76.026760737427963

    assert compare_values(ref, qcdb.get_variable('SCF TOTAL ENERGY'), 5, 'mp2 ref')
    assert compare_values(ss, qcdb.get_variable('CCSD SAME-SPIN CORRELATION ENERGY'), 5, 'ccsd ss')
    assert compare_values(os, qcdb.get_variable('CCSD OPPOSITE-SPIN CORRELATION ENERGY'), 5, 'ccsd os')
    assert compare_values(mp2corl, qcdb.get_variable('MP2 CORRELATION ENERGY'), 5, 'mp2 corl')
    assert compare_values(mp2tot, qcdb.get_variable('MP2 TOTAL ENERGY'), 5, 'mp2 tot')
    assert compare_values(scscorl, qcdb.get_variable('SCS-CCSD CORRELATION ENERGY'), 5, 'scs ccsd corl')
    assert compare_values(scstot, qcdb.get_variable('SCS-CCSD TOTAL ENERGY'), 5, 'scs ccsd tot')
    assert compare_values(ccsd_t_corl, qcdb.get_variable('CCSD(T) CORRELATION ENERGY'), 5, 'ccsd(t) corl')
    assert compare_values(ccsd_t_tot, qcdb.get_variable('CCSD(T) TOTAL ENERGY'), 5, 'ccsd(t) tot')
    if is_5050:
        assert compare_values(a5050corl, qcdb.get_variable('custom SCS-MP2 CORRELATION ENERGY'), 5, 'mp2 scscorl')
        assert compare_values(a5050tot, qcdb.get_variable('custom SCS-MP2 TOTAL ENERGY'), 5, 'mp2 scstot')
    assert compare_values(mp2tot, return_value, 5, 'mp2 return')
    assert compare_values(ref, qcdb.get_variable('CURRENT REFERENCE ENERGY'), 5, 'mp2 ref')
   
#@using_nwchem:


#energy('nwchem')

#clean()
#clean_variables()
#nwchem {}


print('        <<< Translation of NWChem input to Psi4 format to NWChem >>>'

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

energy('nwchem-ccsd(t)')