#non-TCE CCSD(T) energy for N2
import os
import sys
from ..addons import *
from ..utils import *
import qcdb

n2 = qcdb.set_molecule('''
        n 0 0 -0.5
        n 0 0  0.5 ''')

print(n2)

def check_ccsd_t_(return_value, is5050):
    hf          =   -108.929838333552
    nre         =     25.929685200985
    mp2_tot     =   -109.213451339803015
    mp2_corl    =    -0.283613006250685
    ccsd_corr   =    -0.292081644010763
    ccsd_tot    =   -109.221919977563090
    scs_ccsd_corl   =   -0.361897318565201
    scs_ccsd_tot    =   -109.291735652117524
    os          =   -0.227464720235995 * 1.27
    ss          =   -0.064616923774767 * 1.13
    a5050_corl  =   0.5*(ss + os)  
    a5050_tot   =   a5050_corl + scs_ccsd_tot
    t_corr      =   -0.009217759273267
    ccsd_t_corl =   ccsd_corr + t_corr
    ccsd_t_tot  =   -109.231137736836359

    assert compare_values(hf, qcdb.get_variable('HF TOTAL ENERGY'), 5, 'hf ref')
    assert compare_values(nre, qcdb.get_variable('NUCLEAR REPULSION ENERGY'), 5, 'nre')
    assert compare_values(mp2_tot, qcdb.get_variable('MP2 TOTAL ENERGY'), 5, 'mp2 tot')
    assert compare_values(mp2_corl, qcdb.get_variable('MP2 CORRELATION ENERGY'), 5, 'mp2 corl')
    assert compare_values(ccsd_tot, qcdb.get_variable('CCSD TOTAL ENERGY'), 5, 'ccsd tot')
    assert compare_values(t_corr, qcdb.get_variable('(T) CORRECTION ENERGY'), 5, 't corr')
    assert compare_values(ccsd_t_tot, qcdb.get_variable('CCSD(T) TOTAL ENERGY'), 5, 'ccsd(t) tot')
    assert compare_values(ccsd_t_corl, qcdb.get_variable('CCSD(T) CORRELATION ENERGY'), 5, 'ccsd(t) corl')
    if is5050:
        assert compare_values(a5050_corl, qcdb.get_variable('CUSTOM SCS-CCSD CORRELATION ENERGY'), 5, 'custom scs corl')
        assert compare_values(a5050_tot, qcdb.get_variable('CUSTOM SCS-CCSD TOTAL ENERGY'), 5, 'custom scs tot')


@using_nwchem
def test_1_a5050():
    qcdb.set_options({
        'basis' : 'cc-pvdz',
        'nwchem_ccsd_thresh'    :   1.0e-8
        })
    print('testing ccsd(t)...')
    val = qcdb.energy('nwc-ccsd(t)')
    check_ccsd_t_(val, is5050=True)

@using_nwchem
def test_2_a5050_no():
    qcdb.set_options({
        'basis' : 'cc-pvdz',
        'nwchem_ccsd_thresh'    : 1.0e-8
        })
    val = qcdb.energy('nwc-ccsd(t)')
    check_ccsd_t_(val, is5050=False)
