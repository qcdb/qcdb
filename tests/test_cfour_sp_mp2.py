import os
import sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))
from utils import *
from addons import *


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

@using_cfour
def test_sp_rhf_mp2(h2o):
    """cfour/sp-rhf-mp2/input.dat 
    #! single-point MP2/qz2p on water

    """
    h2o = qcdb.set_molecule(h2o)
    qcdb.set_options({
        'cfour_BASIS': 'qz2p',
        'd_convergence': 12
    })

    e, jrec = qcdb.energy('c4-mp2', return_wfn=True, molecule=h2o)

    scftot = -76.0627484601
    mp2corl = -0.270191667216
    mp2tot = -76.332940127333
    tnm = sys._getframe().f_code.co_name
    assert compare_values(scftot, qcdb.get_variable('scf total energy'), 6, tnm + ' SCF')
    assert compare_values(mp2corl, qcdb.get_variable('mp2 correlation energy'), 6, tnm + ' MP2 corl')
    assert compare_values(mp2tot, qcdb.get_variable('mp2 total energy'), 6, tnm + ' MP2')
    assert compare_values(scftot, jrec['qcvars']['HF TOTAL ENERGY'].data, 6, tnm)
    assert compare_values(mp2corl, jrec['qcvars']['MP2 CORRELATION ENERGY'].data, 6, tnm)
    assert compare_values(mp2tot, jrec['qcvars']['CURRENT ENERGY'].data, 6, tnm)


@using_cfour
def test_sp_uhf_mp2(nh2):

    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options({
        #'cfour_CALC_level': 'MP2',
        'cfour_BASIS': 'qz2p',
        'cfour_REFerence': 'UHF',
        'cfour_occupation': [[3,1,1,0],[3,0,1,0]],
        'cfour_SCF_CONV': 12,
    })

    qcdb.energy('c4-mp2', molecule=nh2)

    tnm = sys._getframe().f_code.co_name
    assert compare_values(-55.5893469688, qcdb.get_variable('scf total energy'), 6, tnm + ' SCF')
    assert compare_values(-55.784877360093, qcdb.get_variable('mp2 total energy'), 6, tnm + ' MP2')
    assert compare_values(-0.195530391306, qcdb.get_variable('mp2 correlation energy'), 6, tnm + ' MP2 corl')
    assert compare_values(-0.195530391306, qcdb.get_variable('mp2 doubles energy'), 6, tnm + ' MP2 corl')
    assert compare_values(0., qcdb.get_variable('mp2 singles energy'), 6, tnm + ' MP2 corl')
    assert compare_values(-0.0416164, qcdb.get_variable('mp2 same-spin correlation energy'), 6, tnm + ' MP2 SS corl')
    assert compare_values(-0.1539141, qcdb.get_variable('mp2 opposite-spin correlation energy'), 6, tnm + ' MP2 OS corl')


@using_cfour
def test_sp_rohf_mp2_sc(nh2):

    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options({
        #cfour_CALC_level=MP2
        'cfour_BASIS': 'qz2p',
        'cfour_REFerence': 'ROHF',
        'cfour_OCCUPATION': [[3,1,1,0],[3,0,1,0]],
        'cfour_SCF_CONV': 12,
        'cfour_CC_CONV': 12
    })

    qcdb.energy('c4-mp2', molecule=nh2)

    scftot = -55.5847372601
    scorl = -0.002983751786
    oscorl = -0.155770420921
    sscorl = -0.041785354569
    mp2corl = -0.200539527276
    mp2tot = -55.785276787341
    tnm = sys._getframe().f_code.co_name
    assert compare_values(scftot, qcdb.get_variable('scf total energy'), 6, tnm + ' SCF')
    assert compare_values(scorl, qcdb.get_variable('mp2 singles energy'), 6, tnm + ' MP2 singles')
    # non printed assert compare_values(oscorl, qcdb.get_variable('mp2 opposite-spin correlation energy'), 6, tnm + ' MP2 OS corl')
    # not printed assert compare_values(sscorl, qcdb.get_variable('mp2 same-spin correlation energy'), 6, tnm + ' MP2 SS corl')
    assert compare_values(mp2corl, qcdb.get_variable('mp2 correlation energy'), 6, tnm + ' MP2 corl')
    assert compare_values(mp2corl - scorl, qcdb.get_variable('mp2 doubles energy'), 6, tnm + ' MP2 corl')
    assert compare_values(mp2tot, qcdb.get_variable('mp2 total energy'), 6, tnm + ' MP2')



if __name__ == '__main__':
#    test_sp_rhf_mp2()
#    test_sp_uhf_mp2()
    test_sp_rohf_mp2_sc()

