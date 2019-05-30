import os
import sys
#sys.path.insert(1, os.path.join(sys.path[0], '..'))
from .utils import *
from .addons import *

import qcdb

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


def check_rhf(tnm, fc=False):

    scftot = -76.062748460117
    if fc:
        mp2corl = -0.245151851737
        mp2tot = -76.307900311342
        ccsdcorl = -0.250330561135
        ccsdtot = -76.313079021252
        ccsdtcorr = -0.007096510226
        ccsdttot = -76.320175531445

    else:
        mp2corl = -0.270191667216
        mp2tot = -76.332940127333
        ccsdcorl = -0.275705491773
        ccsdtot = -76.338453951890
        ccsdtcorr = -0.007263597996
        ccsdttot = -76.345717549886

    assert compare_values(scftot, qcdb.get_variable('scf total energy'), 6, tnm + ' SCF')
    assert compare_values(mp2corl, qcdb.get_variable('mp2 correlation energy'), 6, tnm + ' MP2 corl')
    assert compare_values(mp2tot, qcdb.get_variable('mp2 total energy'), 6, tnm + ' MP2')
    assert compare_values(ccsdcorl, qcdb.get_variable('ccsd correlation energy'), 6, tnm + ' CCSD corl')
    assert compare_values(ccsdtot, qcdb.get_variable('ccsd total energy'), 6, tnm + ' CCSD')
    assert compare_values(ccsdtcorr, qcdb.get_variable('(t) correction energy'), 6, tnm + ' (T)')
    assert compare_values(ccsdttot, qcdb.get_variable('ccsd(t) total energy'), 6, tnm + ' CCSD(T)')


def check_uhf(tnm, fc):
    if not fc:
        scftot = -55.5893469688
        mp2corl = -0.195530391306
        mp2tot = -55.784877360382
        ccsdcorl = -0.2132980551718224
        ccsdtot = -55.802645024298
        ccsdtcorr = -0.005166587884
        ccsdttot = -55.807811611842

    assert compare_values(scftot, qcdb.get_variable('scf total energy'), 6, tnm + ' SCF')
    assert compare_values(mp2corl, qcdb.get_variable('mp2 correlation energy'), 6, tnm + ' MP2 corl')
    assert compare_values(mp2tot, qcdb.get_variable('mp2 total energy'), 6, tnm + ' MP2')
    assert compare_values(ccsdcorl, qcdb.get_variable('ccsd correlation energy'), 6, tnm + ' CCSD corl')
    assert compare_values(ccsdtot, qcdb.get_variable('ccsd total energy'), 6, tnm + ' CCSD')
    assert compare_values(ccsdtcorr, qcdb.get_variable('(t) correction energy'), 6, tnm + ' (T)')
    assert compare_values(ccsdttot, qcdb.get_variable('ccsd(t) total energy'), 6, tnm + ' CCSD(T)')


def check_rohf(tnm, fc, prog):

    scftot = -55.5847372601
    if fc:
        smp2corl = -0.002943388250
        osmp2corl = -0.134798097365
        ssmp2corl = -0.038374828493
        mp2corl = -0.176116314109
        mp2tot = -55.760853573995
        ccsdtot = -55.777569551624
        if prog == 'vcc':
            osccsdcorl = -0.152973207228
            ssccsdcorl = -0.039859075346
        elif prog == 'ecc':
            osccsdcorl = -0.1563275
            ssccsdcorl = -0.0365048
        ccsdttot = -55.782617433232
        ccsdcorl = -0.1928322917376748
        ccsdtcorr = -0.005047881608

    else:
        smp2corl = -0.002983751786
        osmp2corl = -0.155770420921
        ssmp2corl = -0.041785354569
        mp2corl = -0.200539527276
        mp2tot = -55.785276787341
        ccsdtot = -55.802586766392
        if prog == 'vcc':
            osccsdcorl = -0.174575174753
            ssccsdcorl = -0.043274331574
        elif prog == 'ecc':
            osccsdcorl = -0.1779580
            ssccsdcorl = -0.0398915
        ccsdttot = -55.807820706828
        ccsdcorl = -0.217849506326
        ccsdtcorr = -0.005233940436

    assert compare_values(scftot, qcdb.get_variable('scf total energy'), 6, tnm + ' SCF')
    # not printed assert compare_values(mp2corl, qcdb.get_variable('mp2 correlation energy'), 6, tnm + ' MP2 corl')
    # not printed assert compare_values(mp2tot, qcdb.get_variable('mp2 total energy'), 6, tnm + 'MP2')
    assert compare_values(ccsdcorl, qcdb.get_variable('ccsd correlation energy'), 6, tnm + ' CCSD corl')
    assert compare_values(ccsdtot, qcdb.get_variable('ccsd total energy'), 6, tnm + ' CCSD')
    assert compare_values(ccsdtcorr, qcdb.get_variable('(t) correction energy'), 6, tnm + ' (T)')
    assert compare_values(ccsdttot, qcdb.get_variable('ccsd(t) total energy'), 6, tnm + ' CCSD(T)')

    #assert compare_values(osmp2corl, qcdb.get_variable('mp2 opposite-spin correlation energy'), 6, tnm + ' MP2 OS corl')
    #assert compare_values(ssmp2corl, qcdb.get_variable('mp2 same-spin correlation energy'), 6, tnm + ' MP2 SS corl')
    assert compare_values(osccsdcorl, qcdb.get_variable('ccsd opposite-spin correlation energy'), 6, tnm + ' CCSD OS corl')
    assert compare_values(ssccsdcorl, qcdb.get_variable('ccsd same-spin correlation energy'), 6, tnm + ' CCSD SS corl')


@using_cfour
def test_sp_rhf_ccsd_t_ao_ecc(h2o):
    """cfour/sp-rhf-ccsd_t_-ao-ecc/input.dat"""

    h2o = qcdb.set_molecule(h2o)
    qcdb.set_options({
        #'cfour_CALC_level': 'CCSD(T)',
        'cfour_BASIS': 'qz2p',
        'cfour_abcdtype': 'aobasis',
        'cfour_cc_program': 'ecc',
        'cfour_SCF_CONV': 12,
        'cfour_CC_CONV': 12
    })

    e, jrec = qcdb.energy('c4-ccsd(t)', return_wfn=True, molecule=h2o)

    check_rhf(tnm=sys._getframe().f_code.co_name, fc=False)
    assert 'CC_PROGRAM           ICCPRO         ECC' in jrec['stdout']
    assert 'ABCDTYPE             IABCDT         AOBASIS' in jrec['stdout']



@using_cfour
def test_sp_rhf_ccsd_t_ao(h2o):

    h2o = qcdb.set_molecule(h2o)
    qcdb.set_options({
        'cfour_BASIS': 'qz2p',
        'cfour_abcdtype': 'aobasis',
        'cfour_SCF_CONV': 12,
        'cfour_CC_CONV': 12,
        'cfour_cc_program': 'vcc'
    })

    e, jrec = qcdb.energy('c4-ccsd(t)', return_wfn=True, molecule=h2o)

    check_rhf(tnm=sys._getframe().f_code.co_name, fc=False)
    assert 'CC_PROGRAM           ICCPRO         VCC' in jrec['stdout']
    assert 'ABCDTYPE             IABCDT         AOBASIS' in jrec['stdout']


@using_cfour
def test_sp_rhf_ccsd_t_ecc(h2o):

    h2o = qcdb.set_molecule(h2o)
    qcdb.set_options({
        'cfour_BASIS': 'qz2p',
        'cfour_SCF_CONV': 12,
        'cfour_CC_CONV': 12,
        'cfour_cc_program': 'ecc'
    })

    e, jrec = qcdb.energy('c4-ccsd(t)', return_wfn=True, molecule=h2o)

    check_rhf(tnm=sys._getframe().f_code.co_name, fc=False)
    assert 'CC_PROGRAM           ICCPRO         ECC' in jrec['stdout']
    assert 'ABCDTYPE             IABCDT         STANDARD' in jrec['stdout']


@using_cfour
def test_sp_rhf_ccsd_t_(h2o):

    h2o = qcdb.set_molecule(h2o)
    qcdb.set_options({
        'cfour_BASIS': 'qz2p',
        'cfour_SCF_CONV': 12,
        'cfour_CC_CONV': 12,
        'cfour_cc_program': 'vcc'
    })

    e, jrec = qcdb.energy('c4-ccsd(t)', return_wfn=True, molecule=h2o)

    check_rhf(tnm=sys._getframe().f_code.co_name, fc=False)
    assert 'CC_PROGRAM           ICCPRO         VCC' in jrec['stdout']
    assert 'ABCDTYPE             IABCDT         STANDARD' in jrec['stdout']


@using_cfour
def test_sp_rhf_ccsd_t_fc(h2o):

    h2o = qcdb.set_molecule(h2o)
    qcdb.set_options({
        'cfour_BASIS': 'qz2p',
        'cfour_abcdtype': 'aobasis',
        'cfour_dropmo' : [1],
        'cfour_SCF_CONV': 12,
        'cfour_CC_CONV': 12,
#        'cfour_cc_program': 'ecc'
    })

    e, jrec = qcdb.energy('c4-ccsd(t)', return_wfn=True, molecule=h2o)

    check_rhf(tnm=sys._getframe().f_code.co_name, fc=True)
    assert 'CC_PROGRAM           ICCPRO         ECC' in jrec['stdout']
    assert 'ABCDTYPE             IABCDT         AOBASIS' in jrec['stdout']


@using_cfour
def test_sp_rhf_ccsd_t_ncc(h2o):

    h2o = qcdb.set_molecule(h2o)
    qcdb.set_options({
        'cfour_BASIS': 'qz2p',
        'cfour_abcdtype': 'aobasis',
        'cfour_cc_program': 'ncc',
        'cfour_SCF_CONV': 12,
        'cfour_CC_CONV': 12
    })

    e, jrec = qcdb.energy('c4-ccsd(t)', return_wfn=True, molecule=h2o)

    check_rhf(tnm=sys._getframe().f_code.co_name, fc=False)
    assert 'CC_PROGRAM           ICCPRO         NCC' in jrec['stdout']
    assert 'ABCDTYPE             IABCDT         AOBASIS' in jrec['stdout']


@using_cfour
def test_sp_uhf_ccsd_t_ao_ecc(nh2):

    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options({
        'cfour_BASIS': 'qz2p',
        'cfour_abcdtype': 'aobasis',
        'cfour_reference': 'uhf',
        'cfour_occupation': [[3,1,1,0], [3,0,1,0]],
        'cfour_SCF_CONV': 12,
        'cfour_CC_CONV': 12,
        'cfour_cc_program': 'ecc'
    })

    e, jrec = qcdb.energy('c4-ccsd(t)', return_wfn=True, molecule=nh2)

    check_uhf(tnm=sys._getframe().f_code.co_name, fc=False)
    assert 'CC_PROGRAM           ICCPRO         ECC' in jrec['stdout']
    assert 'ABCDTYPE             IABCDT         AOBASIS' in jrec['stdout']


@using_cfour
def test_sp_uhf_ccsd_t_ao(nh2):

    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options({
        'cfour_BASIS': 'qz2p',
        'cfour_abcdtype': 'aobasis',
        'cfour_reference': 'uhf',
        'cfour_occupation': [[3,1,1,0], [3,0,1,0]],
        'cfour_SCF_CONV': 12,
        'cfour_CC_CONV': 12,
        'cfour_cc_program': 'vcc'
    })

    e, jrec = qcdb.energy('c4-ccsd(t)', return_wfn=True, molecule=nh2)

    check_uhf(tnm=sys._getframe().f_code.co_name, fc=False)
    assert 'CC_PROGRAM           ICCPRO         VCC' in jrec['stdout']
    assert 'ABCDTYPE             IABCDT         AOBASIS' in jrec['stdout']


@using_cfour
def test_sp_uhf_ccsd_t_ecc(nh2):

    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options({
        'cfour_BASIS': 'qz2p',
        'cfour_reference': 'uhf',
        'cfour_occupation': [[3,1,1,0], [3,0,1,0]],
        'cfour_SCF_CONV': 12,
        'cfour_CC_CONV': 12,
        'cfour_cc_program': 'ecc'
    })

    e, jrec = qcdb.energy('c4-ccsd(t)', return_wfn=True, molecule=nh2)

    check_uhf(tnm=sys._getframe().f_code.co_name, fc=False)
    assert 'CC_PROGRAM           ICCPRO         ECC' in jrec['stdout']
    assert 'ABCDTYPE             IABCDT         STANDARD' in jrec['stdout']


@using_cfour
def test_sp_uhf_ccsd_t_(nh2):

    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options({
        'cfour_BASIS': 'qz2p',
        'cfour_reference': 'uhf',
        'cfour_occupation': [[3,1,1,0], [3,0,1,0]],
        'cfour_SCF_CONV': 12,
        'cfour_CC_CONV': 12,
        'cfour_cc_program': 'vcc'
    })

    e, jrec = qcdb.energy('c4-ccsd(t)', return_wfn=True, molecule=nh2)

    check_uhf(tnm=sys._getframe().f_code.co_name, fc=False)
    assert 'CC_PROGRAM           ICCPRO         VCC' in jrec['stdout']
    assert 'ABCDTYPE             IABCDT         STANDARD' in jrec['stdout']


@using_cfour
def test_sp_rohf_ccsd_t_(nh2):

    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options({
        'cfour_BASIS': 'qz2p',
        'cfour_reference': 'rohf',
        'cfour_occupation': [[3,1,1,0], [3,0,1,0]],
        'cfour_SCF_CONV': 12,
        'cfour_CC_CONV': 12,
        'cfour_print': 2,
        'cfour_cc_program': 'vcc'
    })

    e, jrec = qcdb.energy('c4-ccsd(t)', return_wfn=True, molecule=nh2)

    check_rohf(tnm=sys._getframe().f_code.co_name, fc=False, prog='vcc')
    assert 'CC_PROGRAM           ICCPRO         VCC' in jrec['stdout']
    assert 'ABCDTYPE             IABCDT         STANDARD' in jrec['stdout']


@using_cfour
def test_sp_rohf_ccsd_t_ao(nh2):

    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options({
        'cfour_BASIS': 'qz2p',
        'cfour_reference': 'rohf',
        'cfour_occupation': [[3,1,1,0], [3,0,1,0]],
        'cfour_SCF_CONV': 12,
        'cfour_CC_CONV': 12,
        'cfour_cc_program': 'vcc',
        'cfour_print': 2,
        'cfour_abcdtype': 'aobasis'
    })

    e, jrec = qcdb.energy('c4-ccsd(t)', return_wfn=True, molecule=nh2, prog='vcc')

    check_rohf(tnm=sys._getframe().f_code.co_name, fc=False, prog='vcc')
    assert 'CC_PROGRAM           ICCPRO         VCC' in jrec['stdout']
    assert 'ABCDTYPE             IABCDT         AOBASIS' in jrec['stdout']


@using_cfour
def test_sp_rohf_ccsd_t_ao_ecc(nh2):

    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options({
        'cfour_BASIS': 'qz2p',
        'cfour_reference': 'rohf',
        'cfour_occupation': [[3,1,1,0], [3,0,1,0]],
        'cfour_SCF_CONV': 12,
        'cfour_CC_CONV': 12,
        'cfour_cc_program': 'ecc',
        'cfour_abcdtype': 'aobasis'
    })

    e, jrec = qcdb.energy('c4-ccsd(t)', return_wfn=True, molecule=nh2)

    check_rohf(tnm=sys._getframe().f_code.co_name, fc=False, prog='ecc')
    assert 'CC_PROGRAM           ICCPRO         ECC' in jrec['stdout']
    assert 'ABCDTYPE             IABCDT         AOBASIS' in jrec['stdout']


@using_cfour
def test_sp_rohf_ccsd_t_fc(nh2):

    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options({
        'cfour_BASIS': 'qz2p',
        'cfour_reference': 'rohf',
        'cfour_dropmo': [1],
        'cfour_occupation': [[3,1,1,0], [3,0,1,0]],
        'cfour_SCF_CONV': 12,
        'cfour_CC_CONV': 12,
        'cfour_print': 2,
        'cfour_cc_program': 'ecc'
    })

    e, jrec = qcdb.energy('c4-ccsd(t)', return_wfn=True, molecule=nh2)

    check_rohf(tnm=sys._getframe().f_code.co_name, fc=True, prog='ecc')
    assert 'CC_PROGRAM           ICCPRO         ECC' in jrec['stdout']
    assert 'ABCDTYPE             IABCDT         STANDARD' in jrec['stdout']

