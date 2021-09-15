import os

import pytest

import qcdb

from .utils import *


@pytest.fixture
def h2o():
    return """
    O
    H 1 R
    H 1 R 2 A
    
    R=0.958
    A=104.5
"""

@pytest.fixture
def nh2():
    return """
0 2
N
H 1 R
H 1 R 2 A

R=1.008
A=105.0
"""

@using("cfour")
def test_sp_rhf_ccsd(h2o):
    """cfour/sp-rhf-ccsd/input.dat 
    #! single point CCSD/qz2p on water

    """
    h2o = qcdb.set_molecule(h2o)
    qcdb.set_options({
        #'cfour_CALC_level': 'CCSD',
        'cfour_BASIS': 'qz2p',
        'cfour_SCF_CONV': 12,
        'cfour_CC_CONV': 12
    })

    e, jrec = qcdb.energy('c4-ccsd', return_wfn=True, molecule=h2o)

    scftot = -76.062748460117
    mp2tot = -76.332940127333
    ccsdcorl = -0.275705491773
    ccsdtot = -76.338453951890
    assert compare_values(scftot, qcdb.variable('scf total energy'), 6, tnm() + 'SCF')
    assert compare_values(mp2tot, qcdb.variable('mp2 total energy'), 6, tnm() + 'MP2')
    assert compare_values(ccsdcorl, qcdb.variable('ccsd correlation energy'), 6, tnm() + 'CCSD corl')
    assert compare_values(ccsdtot, qcdb.variable('ccsd total energy'), 6, tnm() + 'CCSD')


@using("cfour")
def test_sp_rhf_ccsd_ao(h2o):
    """cfour/sp-rhf-ccsd-ao/input.dat 
    #! single point CCSD/qz2p on water

    """
    h2o = qcdb.set_molecule(h2o)
    qcdb.set_options({
        #'cfour_CALC_level': 'CCSD',
        'cfour_BASIS': 'qz2p',
        'cfour_abcdtype': 'aobasis',
        'cfour_SCF_CONV': 12,
        'cfour_CC_CONV': 12
    })

    e, jrec = qcdb.energy('c4-ccsd', return_wfn=True, molecule=h2o)

    scftot = -76.062748460117
    mp2tot = -76.332940127333
    ccsdcorl = -0.275705491773
    ccsdtot = -76.338453951890
    assert compare_values(scftot, qcdb.variable('scf total energy'), 6, tnm() + 'SCF')
    assert compare_values(mp2tot, qcdb.variable('mp2 total energy'), 6, tnm() + 'MP2')
    assert compare_values(ccsdcorl, qcdb.variable('ccsd correlation energy'), 6, tnm() + 'CCSD corl')
    assert compare_values(ccsdtot, qcdb.variable('ccsd total energy'), 6, tnm() + 'CCSD')


@using("cfour")
def test_sp_uhf_ccsd(nh2):
    """cfour/sp-uhf-ccsd/input.dat
    #! single-point CCSD/qz2p on NH2

    """
    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options({
        #'cfour_CALC_level': 'CCSD',
        'cfour_BASIS': 'qz2p',
        'cfour_REFerence': 'UHF',
        'cfour_occupation': [[3,1,1,0],[3,0,1,0]],
        'cfour_SCF_CONV': 12,
        'cfour_CC_CONV': 12,
    })

    qcdb.energy('c4-ccsd', molecule=nh2)

    scftot = -55.5893469688
    mp2tot = -55.784877360093
    ccsdcorl = -0.213298055172
    assert compare_values(scftot, qcdb.variable('scf total energy'), 6, tnm() + 'SCF')
    assert compare_values(mp2tot, qcdb.variable('mp2 total energy'), 6, tnm() + 'MP2')
    assert compare_values(ccsdcorl, qcdb.variable('ccsd correlation energy'), 6, tnm() + 'CCSD')


@using("cfour")
def test_sp_rohf_ccsd(nh2):

    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options({
        #cfour_CALC_level': 'CCSD',
        'cfour_BASIS': 'qz2p',
        'cfour_REFerence': 'ROHF',
        'cfour_OCCUPATION': [[3,1,1,0],[3,0,1,0]],
        'cfour_SCF_CONV': 12,
        'cfour_PRINT': 2
    })

    qcdb.energy('c4-ccsd', molecule=nh2)

    scftot = -55.5847372601
    smp2corl = -0.002983751786
    osmp2corl = -0.155770420921
    ssmp2corl = -0.041785354569
    mp2corl = -0.200539527276
    mp2tot = -55.785276787341
    osccsdcorl = -0.1745752
    ssccsdcorl = -0.0432743
    ccsdcorl = -0.217849506326
    ccsdtot = -55.802586766392
    compare_values(scftot, qcdb.variable('scf total energy'), 6, tnm() + 'SCF')
    # not printed compare_values(smp2corl, qcdb.variable('mp2 singles energy'), 6, tnm() + 'MP2 singles')
    # not printed compare_values(osmp2corl, qcdb.variable('mp2 opposite-spin correlation energy'), 6, tnm() + 'MP2 OS corl')
    # not printed compare_values(ssmp2corl, qcdb.variable('mp2 same-spin correlation energy'), 6, tnm() + 'MP2 SS corl')
    # not printed compare_values(mp2corl, qcdb.variable('mp2 correlation energy'), 6, tnm() + 'MP2 corl')
    # not printed compare_values(mp2tot, qcdb.variable('mp2 total energy'), 6, tnm() + 'MP2')

    # * Summer 2021: have a value for ccsd opposite-spin correlation below, but probably not collected from current version with ROHF
    # compare_values(osccsdcorl, qcdb.variable('ccsd opposite-spin correlation energy'), 6, tnm() + 'CCSD OS corl')
    compare_values(ssccsdcorl, qcdb.variable('ccsd same-spin correlation energy'), 6, tnm() + 'CCSD SS corl')
    compare_values(ccsdcorl, qcdb.variable('ccsd correlation energy'), 6, tnm() + 'CCSD corl')
    compare_values(ccsdtot, qcdb.variable('ccsd total energy'), 6, tnm() + 'CCSD')
    compare_values(ccsdcorl, qcdb.variable('current correlation energy'), 6, tnm() + 'Current corl')
    compare_values(ccsdtot, qcdb.variable('current energy'), 6, tnm() + 'Current')




#if __name__ == '__main__':
##    test_sp_rhf_mp2()
##    test_sp_uhf_mp2()
#    test_sp_rohf_mp2_sc()
