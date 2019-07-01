#! Gradient MP2/cc-pvdz on water
import os
import sys
from ..utils import *
from ..addons import *
import qcdb
import numpy as np


def check_mp2(val, is_df, is5050):
    if is_df:
        ref = -76.026760737428
        mp2_tot = -76.230777733733
        mp2_corl = -0.204016996305
        scs_corl = -0.200161577112
        scs_tot = -76.226922314540
        ss_corl = -0.051529405908 * 0.333333333333
        os_corl = -0.152487590397 * 1.200000000000
        a5050corl = 0.5 * (ss_corl + os_corl)
        a5050tot = a5050corl + ref
        ref_grad = np.array([[0.000000, 0.000000, 0.012114], [-0.001793, 0.000000, -0.006057],
                             [0.001793, 0.000000, -0.006057]])
    else:
        ref = -76.026760737428
        mp2_tot = -76.226922314540
        mp2_corl = -0.204016996305
        scs_corl = -0.200161577112
        scs_tot = -76.226922314540
        ss_corl = -0.051529405908 * 0.333333333333
        os_corl = -0.152487590397 * 1.200000000000
        a5050corl = 0.5 * (ss_corl + os_corl)
        a5050tot = a5050corl + ref
        ref_grad = np.array([[0.000000, 0.000000, 0.012114], [-0.001793, 0.000000, -0.006057],
                             [0.001793, 0.000000, -0.006057]])

    assert compare_values(ref, qcdb.get_variable('HF TOTAL ENERGY'), 5, 'scf')
    assert compare_values(mp2_tot, qcdb.get_variable('MP2 TOTAL ENERGY'), 5, 'mp2 tot')
    assert compare_values(mp2_corl, qcdb.get_variable('MP2 CORRELATION ENERGY'), 5, 'mp2 corl')
#    assert compare_values(scs_tot, qcdb.get_variable('SCS-MP2 TOTAL ENERGY'), 5, 'scs mp2 tot')
#    assert compare_values(scs_corl, qcdb.get_variable('SCS-MP2 CORRELATION ENERGY'), 5, 'scs mp2 corl')
    assert compare_arrays(ref_grad, qcdb.get_variable('CURRENT GRADIENT'), 5, 'mp2 grad')
 #   if is5050:
  #      assert compare_values(a5050corl, qcdb.get_variable('CUSTOM SCS-MP2 CORRELATION ENERGY'), 5,
                              #'custom scs mp2 corl')
   #     assert compare_values(a5050tot, qcdb.get_variable('CUSTOM SCS-MP2 TOTAL ENERGY'), 5, 'custom scs-mp2 tot')


@using_nwchem
def test_1_mp2():
    h2o = qcdb.set_molecule('''
        O     0.000000000000    0.000000000000   -0.065638538099
        H     0.000000000000   -0.757480611647    0.520865616174
        H     0.000000000000    0.757480611647    0.520865616174
        ''')

    qcdb.set_options({
        'basis': 'cc-pvdz',
        'memory': '400 mb',
        'scf__e_convergence': 1.0e-4,
        'nwchem_scf__rhf': True,
        #'nwchem_scf__thresh': 1.0e-4,
        'nwchem_scf__nopen': 0,
        'nwchem_mp2__tight': True,
    #    'nwchem_task_mp2': 'gradient'
    })
    print('Testing mp2 ...')
    val = qcdb.gradient('nwc-mp2')
    check_mp2(val, is_df=True, is5050=False)


@using_nwchem
def test_2_hf():
    h2o = qcdb.set_molecule('''
        O     0.000000000000    0.000000000000   -0.065638538099
        H     0.000000000000   -0.757480611647    0.520865616174
        H     0.000000000000    0.757480611647    0.520865616174
        ''')

    qcdb.set_options({
        'basis': 'cc-pvdz',
        'memory': '400 mb',
        'scf__e_convergence': 1.0e-4,
        'nwchem_scf__rhf': True,
        #'nwchem_scf__thresh': 1.0e-4,
        'nwchem_scf__nopen': 0,
        'nwchem_mp2__tight': True,
    })
    print('Testing hf...')
    val = qcdb.gradient('nwc-mp2')
    check_mp2(val, is_df=True, is5050=False)


#@using_nwchem
#def test_3_mp2_custom():
#    qcdb.set_options({
#        'basis': 'cc-pvdz',
#        'memory': '400 mb',
#        'scf__e_convergence': 1.0e-4,
#        'nwchem_scf__rhf': True,
        #'nwchem_scf__thresh': 1.0e-4,
#        'nwchem_scf__nopen': 0,
#        'nwchem_mp2__tight': True,
    #    'nwchem_task_mp2': 'gradient'
#    })
#    print('Testing mp2 ...')
#    val = qcdb.gradient('nwc-mp2')
#    check_mp2(val, is_df=False, is5050=True)

@using_nwchem
def test_4_mp2_array():
    h2o = qcdb.set_molecule('''
        O     0.000000000000    0.000000000000   -0.065638538099
        H     0.000000000000   -0.757480611647    0.520865616174
        H     0.000000000000    0.757480611647    0.520865616174
        ''')

    qcdb.set_options({
        'basis': 'cc-pvdz',
        'memory': '400 mb',
        'scf__e_convergence': 1.0e-4,
        'nwchem_scf__rhf': True,
        'nwchem_scf__thresh': 1.0e-4,
        'nwchem_scf__nopen': 0,
        'nwchem_mp2__tight': True,
    #    'nwchem_task_mp2': 'gradient'
    })
    print('Testing mp2 ...')
    val = qcdb.gradient('nwc-mp2')
    check_mp2(val, is_df=True, is5050=False)
