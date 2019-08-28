import os
import sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))
from .utils import *
from .addons import *

import qcdb

#! conventional and density-fitting mp2 test of mp2 itself and setting scs-mp2


def check_mp2(return_value, is_df, is_5050=False):
    if is_df:
        ref       =    -76.0167614256151865
        mp2ss     =     -0.0527406422061238
        mp2os     =     -0.1562926850310142
        mp2corl   =     -0.2090333272371381
        mp2tot    =    -76.2257947528523232
        scscorl   =     -0.2051314361059251
        scstot    =    -76.2218928617211162
        a5050corl =      0.5 * (mp2ss + mp2os)
        a5050tot  =      a5050corl + ref
    else:
        ref       =   -76.01678947133706
        mp2ss     =    -0.05268120425816
        mp2os     =    -0.15637564436589
        mp2corl   =    -0.20905684862405
        mp2tot    =   -76.22584631996111
        scscorl   =    -0.20521117465845
        scstot    =   -76.22200064599551
        a5050corl =     0.5 * (mp2ss + mp2os)
        a5050tot  =     a5050corl + ref

    assert compare_values(ref, qcdb.variable('SCF TOTAL ENERGY'), 5, 'mp2 ref')
    assert compare_values(mp2ss, qcdb.variable('MP2 SAME-SPIN CORRELATION ENERGY'), 5, 'mp2 ss')
    assert compare_values(mp2os, qcdb.variable('MP2 OPPOSITE-SPIN CORRELATION ENERGY'), 5, 'mp2 os')
    assert compare_values(mp2corl, qcdb.variable('MP2 CORRELATION ENERGY'), 5, 'mp2 corl')
    assert compare_values(mp2tot, qcdb.variable('MP2 TOTAL ENERGY'), 5, 'mp2 tot')
    assert compare_values(scscorl, qcdb.variable('SCS-MP2 CORRELATION ENERGY'), 5, 'mp2 scscorl')
    assert compare_values(scstot, qcdb.variable('SCS-MP2 TOTAL ENERGY'), 5, 'mp2 scstot')
    if is_5050:
        assert compare_values(a5050corl, qcdb.variable('custom SCS-MP2 CORRELATION ENERGY'), 5, 'mp2 scscorl')
        assert compare_values(a5050tot, qcdb.variable('custom SCS-MP2 TOTAL ENERGY'), 5, 'mp2 scstot')

    assert compare_values(ref, qcdb.variable('CURRENT REFERENCE ENERGY'), 5, 'mp2 ref')
    assert compare_values(mp2corl, qcdb.variable('CURRENT CORRELATION ENERGY'), 5, 'mp2 corl')
    assert compare_values(mp2tot, qcdb.variable('CURRENT ENERGY'), 5, 'mp2 tot')
    assert compare_values(mp2tot, return_value, 5, 'mp2 return')

@pytest.fixture
def h2o():
    return """
O
H 1 1.0
H 1 1.0 2 90.0
"""


@using_psi4
def test_1_df_mp2(h2o):
    qcdb.set_molecule(h2o)
    qcdb.set_options({
        'basis': 'cc-pvdz'
    })

    print('   Testing mp2 (df) ...')
    val = qcdb.energy('mp2')
    check_mp2(val, is_df=True)


@using_psi4
def test_2_conv_mp2(h2o):
    qcdb.set_molecule(h2o)
    qcdb.set_options({
        'basis': 'cc-pvdz',
        'psi4_mp2_type': 'conv'
    })

    print('   Testing mp2 (conv) ...')
    val = qcdb.energy('mp2')
    check_mp2(val, is_df=False)


@using_psi4
def test_3_df_scs_mp2(h2o):
    qcdb.set_molecule(h2o)
    qcdb.set_options({
        'basis': 'cc-pvdz',
        'psi4_mp2_os_scale': 1.2,
        'psi4_mp2_ss_scale': 0.33333333333333333,
        'psi4_mp2_type': 'df',
    })

#set mp2_type df

    print('   Testing explicit scs mp2 (df) ...')
    val = qcdb.energy('mp2')
    check_mp2(val, is_df=True)


@using_psi4
def test_4_conv_scs_mp2(h2o):
    qcdb.set_molecule(h2o)
    qcdb.set_options({
        'basis': 'cc-pvdz',
        'psi4_mp2_os_scale': 1.2,
        'psi4_mp2_ss_scale': 0.33333333333333333,
        'psi4_mp2_type': 'conv',
    })

    print('   Testing explicit scs mp2 (conv) ...')
    val = qcdb.energy('mp2')
    check_mp2(val, is_df=False)


@using_psi4
def test_5_df_custom_scs_mp2(h2o):
    qcdb.set_molecule(h2o)
    qcdb.set_options({
        'basis': 'cc-pvdz',
        'psi4_mp2_os_scale':  0.5,
        'psi4_mp2_ss_scale':  0.5,
    })

#set mp2_type df

    print('   Testing user-def scs mp2 (df) ...')
    val = qcdb.energy('mp2')
    check_mp2(val, is_df=True, is_5050=True)


@using_psi4
def test_6_conv_custom_scs_mp2(h2o):
    qcdb.set_molecule(h2o)
    qcdb.set_options({
        'basis': 'cc-pvdz',
        'psi4_mp2_os_scale':  0.5,
        'psi4_mp2_ss_scale':  0.5,
        'psi4_mp2_type': 'conv',
    })
#
#set mp2_type conv

    print('   Testing user-def scs mp2 (conv) ...')
    val = qcdb.energy('mp2')
    check_mp2(val, is_df=False, is_5050=True)

