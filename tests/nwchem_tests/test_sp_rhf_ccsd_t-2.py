#! single-point CCSD(T)/cc-pvdz on water
# CCSD(T) calculation through "task ccsd(t)" is cheaper than "task tce"
# and able to harvest MP2 energy as well.
import os
import sys
from ..utils import *
from ..addons import *
import qcdb


def check_ccsd_t_2(return_value):
        ref = -76.026760737428
        nre = 9.187334240165
        ccsd_corr = -0.213341272556805
        ccsd_tot = -76.240102009984767
        scscorl = -0.264498694126312
        scstot = -76.291259431554266
        mp2corl = -0.204016996303923
        mp2tot = -76.230777733731884
        ss = -0.046033728720216 * 1.130000000000000
        os = -0.167307543836588 * 1.270000000000000
        ccsd_t_tot = -76.243161551636248
        t_corr = -0.003059541972533
        ccsd_t_corl= ccsd_corr + t_corr
        a5050corl = 0.5 * (ss + os)
        a5050tot = a5050corl + scstot
        
        assert compare_values(ref, qcdb.variable('HF TOTAL ENERGY'), 5, 'hf ref')
        assert compare_values(nre, qcdb.variable('NUCLEAR REPULSION ENERGY'), 5, 'nre')
        assert compare_values(mp2corl, qcdb.variable('MP2 CORRELATION ENERGY'), 5, 'mp2 corl')
        assert compare_values(mp2tot, qcdb.variable('MP2 TOTAL ENERGY'), 5, 'mp2 tot')
        #assert compare_values(scscorl, qcdb.variable('SCS-CCSD CORRELATION ENERGY'), 5, 'scs ccsd corl')
        #assert compare_values(scstot, qcdb.variable('SCS-CCSD TOTAL ENERGY'), 5, 'scs ccsd tot')
        assert compare_values(ccsd_t_corl, qcdb.variable('CCSD(T) CORRELATION ENERGY'), 5, 'ccsd(t) corl')
        assert compare_values(ccsd_t_tot, qcdb.variable('CCSD(T) TOTAL ENERGY'), 5, 'ccsd(t) tot')
#        if is_5050:
#            assert compare_values(a5050corl, qcdb.variable('CUSTOM SCS-CCSD CORRELATION ENERGY'), 5, 'mp2 scscorl')
#            assert compare_values(a5050tot, qcdb.variable('CUSTOM SCS-CCSD TOTAL ENERGY'), 5, 'mp2 scstot')
       # assert compare_values(mp2tot, return_value, 5, 'mp2 return')
        assert compare_values(t_corr, qcdb.variable('(T) CORRECTION ENERGY'), 5, 'ccsd(t) correction')


@using_nwchem
def test_1_a5050_no():
    h2o = qcdb.set_molecule('''
        O     0.000000000000    0.000000000000   -0.065638538099
        H     0.000000000000   -0.757480611647    0.520865616174
        H     0.000000000000    0.757480611647    0.520865616174
        ''')

    qcdb.set_options({
        'basis': 'cc-pvdz',
        'memory': '600 mb',
        'nwchem_scf__rhf': True,
        'nwchem_scf__thresh': 1.0e-12,
    })
    print("Testing ccsd(t)...")
    val = qcdb.energy('nwc-ccsd(t)')
    check_ccsd_t_2(val)


#@using_nwchem
#def test_2_a5050():
#    qcdb.set_options({
#        'basis': 'cc-pvdz',
#        'memory': '600 mb',
#        'nwchem_scf__rhf': True,
#        'nwchem_scf__thresh': 1.0e-12,
#        #'nwchem_tce__dft': False,
#        #'nwchem_tce__module': 'ccsd(t)',
#        #'nwchem_tce': False,
#        #'nwchem_tce__thresh': 1.0e-12
#    })
#    print("testing mp2...")
#    val = qcdb.energy('nwc-ccsd(t)')
#    check_ccsd_t_2(val, is_5050=True)

