import os
import sys

import numpy as np

sys.path.insert(1, os.path.join(sys.path[0], '..'))
from utils import *
from addons import *


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

    assert compare_values(ref, qcdb.get_variable('SCF TOTAL ENERGY'), 5, 'mp2 ref')
    assert compare_values(mp2ss, qcdb.get_variable('MP2 SAME-SPIN CORRELATION ENERGY'), 5, 'mp2 ss')
    assert compare_values(mp2os, qcdb.get_variable('MP2 OPPOSITE-SPIN CORRELATION ENERGY'), 5, 'mp2 os')
    assert compare_values(mp2corl, qcdb.get_variable('MP2 CORRELATION ENERGY'), 5, 'mp2 corl')
    assert compare_values(mp2tot, qcdb.get_variable('MP2 TOTAL ENERGY'), 5, 'mp2 tot')
    assert compare_values(scscorl, qcdb.get_variable('SCS-MP2 CORRELATION ENERGY'), 5, 'mp2 scscorl')
    assert compare_values(scstot, qcdb.get_variable('SCS-MP2 TOTAL ENERGY'), 5, 'mp2 scstot')
    if is_5050:
        assert compare_values(a5050corl, qcdb.get_variable('custom SCS-MP2 CORRELATION ENERGY'), 5, 'mp2 scscorl')
        assert compare_values(a5050tot, qcdb.get_variable('custom SCS-MP2 TOTAL ENERGY'), 5, 'mp2 scstot')

    assert compare_values(ref, qcdb.get_variable('CURRENT REFERENCE ENERGY'), 5, 'mp2 ref')
    assert compare_values(mp2corl, qcdb.get_variable('CURRENT CORRELATION ENERGY'), 5, 'mp2 corl')
    assert compare_values(mp2tot, qcdb.get_variable('CURRENT ENERGY'), 5, 'mp2 tot')
    assert compare_values(mp2tot, return_value, 5, 'mp2 return')

h2o = qcdb.set_molecule("""
O
H 1 1.0
H 1 1.0 2 90.0
""")

@pytest.fixture
def h2o_y():
    return """
H 0 0 1
O 0 0 0
H 1 0 0
no_com
no_reorient
"""

@pytest.fixture
def h2o_z():
    return """
O 0 0 0
H 0 1 0
H 1 0 0
no_com
no_reorient
"""

dip_z = np.array([   1.52054865,  1.52054865,  0.])
grad_z = np.array([[-0.04159183, -0.04159183, -0.        ],
                   [ 0.01821547,  0.02337636,  0.        ],
                   [ 0.02337636,  0.01821547,  0.        ]])
dip_y = np.array([   1.52054865, 0.,  1.52054865])
grad_y = np.array([[ 0.01821547, 0.,  0.02337636],
                   [-0.04159183, 0., -0.04159183],
                   [ 0.02337636, 0.,  0.01821547]])



@using_cfour
def test_2_conv_mp2(h2o_z):
    qcdb.set_molecule(h2o_z)
    qcdb.set_options({
        'basis': 'cc-pvdz',
        'psi4_mp2_type': 'conv'
    })

    print('   Testing mp2 (conv) ...')
    val = qcdb.energy('c4-mp2')
    check_mp2(val, is_df=False)


@using_cfour
def test_4_conv_scs_mp2(h2o_z):
    qcdb.set_molecule(h2o_z)
    qcdb.set_options({
        'basis': 'cc-pvdz',
        'mp2_os_scale': 1.2,
        'mp2_ss_scale': 0.33333333333333333,
        'psi4_mp2_type': 'conv',
    })

    print('   Testing explicit scs mp2 (conv) ...')
    val = qcdb.energy('c4-mp2')
    check_mp2(val, is_df=False)


@using_cfour
def test_scale(h2o_z):
    h2o_z = qcdb.Molecule(h2o_z)
    qcdb.set_options({
        'basis': 'cc-pvdz',
        'cfour_spin_scal': 'on',
        'cfour_reference': 'uhf',
        'cfour_calc_level': 'mp2',
        'cfour_deriv_level': 'first',
        'CFOUR_DIFF_TYPE': 'relaxed',
    })

    val, jrec = qcdb.energy('c4-cfour', return_wfn=True, molecule=h2o_z)
    dip = np.array([float(jrec['qcvars']['CURRENT DIPOLE {}'.format(i)].data) for i in 'XYZ'])

    tnm = sys._getframe().f_code.co_name
    assert compare_arrays(grad_z, jrec['qcvars']['CURRENT GRADIENT'].data, 5, tnm + ' grad')
    assert compare_arrays(dip_z, dip, 5, tnm + ' dipole')

@using_cfour
def test_scale_2(h2o_y):
    h2o_y = qcdb.Molecule(h2o_y)
    qcdb.set_options({
        'basis': 'cc-pvdz',
        'cfour_spin_scal': 'on',
        'cfour_reference': 'uhf',
        'cfour_calc_level': 'mp2',
        'cfour_deriv_level': 'first',
        'CFOUR_DIFF_TYPE': 'relaxed',
    })

    val, jrec = qcdb.energy('c4-cfour', return_wfn=True, molecule=h2o_y)
    dip = np.array([float(jrec['qcvars']['CURRENT DIPOLE {}'.format(i)].data) for i in 'XYZ'])

    tnm = sys._getframe().f_code.co_name
    assert compare_arrays(grad_y, jrec['qcvars']['CURRENT GRADIENT'].data, 5, tnm + ' grad')
    assert compare_arrays(dip_y, dip, 5, tnm + ' dipole')



@using_cfour
def test_6_conv_custom_scs_mp2(h2o_z):
    qcdb.set_molecule(h2o_z)
    qcdb.set_options({
        'basis': 'cc-pvdz',
        'mp2_os_scale':  0.5,
        'mp2_ss_scale':  0.5,
        'psi4_mp2_type': 'conv',
    })

    print('   Testing user-def scs mp2 (conv) ...')
    val = qcdb.energy('c4-mp2')
    check_mp2(val, is_df=False, is_5050=True)
#    assert False

