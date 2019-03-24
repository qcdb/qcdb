#! single-point MP2/cc-pvdz on water
import os
import sys
from ..utils import *
from ..addons import *
import qcdb

h2o = qcdb.set_molecule('''
        O     0.000000000000    0.000000000000   -0.065638538099
        H     0.000000000000   -0.757480611647    0.520865616174
        H     0.000000000000    0.757480611647    0.520865616174
        ''')
print(h2o)


def check_mp2(return_value, is5050):
    ref = -76.026760737428
    nre = 9.187334240165
    mp2_tot = -76.230777733733
    scs_mp2_tot = -76.226922314540
    scs_corl = -0.200161577112
    os = -0.152487590397 * 1.200000000000
    ss = -0.051529405908 * 0.333333333333
    a5050corl = 0.5 * (os + ss)
    a5050tot = a5050corl + ref
    assert compare_values(ref, qcdb.get_variable('HF TOTAL ENERGY'), 5, 'scf total')
    assert compare_values(mp2_tot, qcdb.get_variable('MP2 TOTAL ENERGY'), 5, 'mp2 energy')
    assert compare_values(scs_mp2_tot, qcdb.get_variable('SCS-MP2 TOTAL ENERGY'), 5, 'scs-mp2')
    assert compare_values(scs_corl, qcdb.get_variable('SCS-MP2 CORRELATION ENERGY'), 5, 'scs-mp2 corl')
    if is5050:
        assert compare_values(a5050corl, qcdb.get_variable('CUSTOM SCS-MP2 CORRELATION ENERGY'), 5, 'custom scs corl')
        assert compare_values(a5050tot, qcdb.get_variable('CUSTOM SCS-MP2 TOTAL ENERGY'), 5, 'custom scs-mp2 tot')


@using_nwchem
def test_1_a5050_no():
    qcdb.set_options({
        'basis': 'cc-pvdz', 
        'memory': '400 mb', 
        'nwchem_mp2__tight': True, 
        #'nwchem_task_mp2': 'energy'
        })
    print('Testing mp2(df)...')
    # print(jrec['qcvars'])
    val = qcdb.energy('nwc-mp2')  # return_wfn=True)
    check_mp2(val, is5050=False)

@using_nwchem
def test_2_a5050():
    qcdb.set_options({
        'basis': 'cc-pvdz', 
        'memory': '400 mb', 
        'nwchem_mp2__tight': True, 
        #'nwchem_task_mp2': 'energy'  # not the way to specify anymore
        })
    print('Testing scs-mp2(df)...')
    val = qcdb.energy('nwc-mp2')
    check_mp2(val, is5050=True)
