#! single-point CCSD(T)/cc-pvdz on water
# CCSD(T) calculation through "task ccsd(t)" is cheaper than "task tce" 
# and able to harvest MP2 energy as well. 
import os
import sys
from utils import *
from addons import *
import qcdb

h2o=qcdb.set_molecule('''
        O     0.000000000000    0.000000000000   -0.065638538099
        H     0.000000000000   -0.757480611647    0.520865616174
        H     0.000000000000    0.757480611647    0.520865616174
        ''')

print(h2o)
def check_ccsd_t_2(return_value, is_df, is_5050):
    if is_df:
        ref         =   -76.026760737428
        nre         =     9.187334240165
        ccsd_corr   =    -0.213341272556805
        ccsd_tot    =   -76.240102009984767
        scscorl     =    -0.264498694126312
        scstot      =   -76.291259431554266
        mp2corl     =    -0.204016996303923
        mp2tot      =   -76.230777733731884
        ss          =    -0.046033728720216
        os          =    -0.167307543836588
        ccsd_t_corl =    -0.003059541651480
        ccsd_t_tot  =   -76.243161551636248
        t_corr      =    -0.003145008014227
        a5050corl   =   0.5*(ss + os)
        a5050tot    =   a5050corl + ref
    else:
        ref         =   -76.026760737427963
        nre         =     9.187334240165
        ccsd_corr   =    -0.213341272556805
        ccsd_tot    =   -76.240102009984767
        scscorl     =    -0.264498694126312
        scstot      =   -76.291259431554266
        mp2corl     =    -0.204016996303923
        mp2tot      =   -76.230777733731884
        ss          =    -0.046033728720216
        os          =    -0.167307543836588
        t_corr      =    -0.003059541651480
        ccsd_t_tot  =   -76.243247017998996
        ccsd_t_corl =    -0.003059541651480
        a5050corl   =   0.5*(ss + os)
        a5050tot    =   a5050corl + ref

    assert compare_values(ref, qcdb.get_variable('HF TOTAL ENERGY'), 5, 'hf ref')
    assert compare_values(nre, qcdb.get_variable('NUCLEAR REPULSION ENERGY'), 5, 'nre')
    assert compare_values(ss, qcdb.get_variable('CCSD SAME-SPIN CORRELATION ENERGY'), 5, 'ccsd ss')
    assert compare_values(os, qcdb.get_variable('CCSD OPPOSITE-SPIN CORRELATION ENERGY'), 5, 'ccsd os')
    assert compare_values(mp2corl, qcdb.get_variable('MP2 CORRELATION ENERGY'), 5, 'mp2 corl')
    assert compare_values(mp2tot, qcdb.get_variable('MP2 TOTAL ENERGY'), 5, 'mp2 tot')
    assert compare_values(scscorl, qcdb.get_variable('SCS-MP2 CORRELATION ENERGY'), 5, 'scs ccsd corl')
    assert compare_values(scstot, qcdb.get_variable('SCS-MP2 TOTAL ENERGY'), 5, 'scs ccsd tot')
    assert compare_values(t_corr, qcdb.get_variable('CCSD(T) CORRELATION ENERGY'), 5, 'ccsd(t) corl')
    assert compare_values(ccsd_t_tot, qcdb.get_variable('CCSD(T) TOTAL ENERGY'), 5, 'ccsd(t) tot')
    if is_5050:
        assert compare_values(a5050corl, qcdb.get_variable('CUSTOM SCS-MP2 CORRELATION ENERGY'), 5, 'mp2 scscorl')
        assert compare_values(a5050tot, qcdb.get_variable('CUSTOM SCS-MP2 TOTAL ENERGY'), 5, 'mp2 scstot')
    assert compare_values(mp2tot, return_value, 5, 'mp2 return')
    assert compare_values(t_corr, qcdb.get_variable('(T) CORRECTION ENERGY'), 5, 'ccsd(t) correction') 
@using_nwchem
def test_1_hf():
    qcdb.set_options({
        'basis'     :   'cc-pvdz',
        'memory'    :   '600 mb',
        'nwchem_scf':   'rhf',
        'nwchem_scf_thresh': 1.0e-12
        })
    print("Testing hf...")
    val = qcdb.energy('nwc-hf')
    check_ccsd_t_2(val, is_df=True, is_5050=False)

@using_nwchem
def test_2_mp2():
    qcdb.set_options({
        'basis'     :   'cc-pvdz',
        'memory'    :   '600 mb'
        })
    print("Testing MP2...")
    val = qcdb.energy('nwc-mp2')
    check_ccsd_t_2(val, is_df=True, is_5050=False)

@using_nwchem
def test_3_mp2_5050():
    qcdb.set_options({
        'basis'     :   'cc-pvdz',
        'memory'    :   '600 mb'
        })
    print("testing mp2...")
    val = qcdb.energy('nwc-mp2')
    check_ccsd_t_2(val, is_df=False, is_5050=True)

@using_nwchem
def test_4_ccsd():
    qcdb.set_options({
        'basis'     :   'cc-pvdz',
        'memory'    :   '600 mb',
        #'nwchem_tce_dft'    :   False,
        #'nwchem_tce':   'ccsd(t)',
        #'nwchem_tce_on': False,
        #'nwchem_tce_thresh' :   1.0e-12
        })
    print("Testing CCSD energy ...")
    val = qcdb.energy('nwc-ccsd')
    check_ccsd_t_2(val, is_df=True, is_5050=False)

@using_nwchem
def test_5_ccsd_t():
    qcdb.set_options({
        'basis'     :   'cc-pvdz',
        'memory'    :   '600 mb',
        'nwchem_tce_dft'    :   False,
        'nwchem_tce':   'ccsd(t)',
        #'nwchem_tce_on': False,
        'nwchem_tce_thresh' :   1.0e-12
        })
    print("Testing CCSD(T) energy ...")
    val = qcdb.energy('nwc-ccsd(t)')
    check_ccsd_t_2(val, is_df=True, is_5050=False)
@using_nwchem
def test_6_custom_ccsd():
    qcdb.set_options({
        'basis'     :   'cc-pvdz',
        'memory'    :   '600 mb',
        'nwchem_tce_dft'    :   False,
        'nwchem_tce':   'ccsd(t)',
        #'nwchem_tce_on': True,
        'nwchem_tce_thresh' :   1.0e-12
        })
    print("Testing custom SCS MP2 energy ...")
    val = qcdb.energy('nwc-ccsd(t)')
    check_ccsd_t_2(val, is_df=False, is_5050=True)
