import os
import sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))
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

@using_nwchem
def test_sp_rhf_ccsd(h2o):
    """cfour/sp-rhf-ccsd/input.dat 
    #! single point CCSD/qz2p on water

    """
    h2o = qcdb.set_molecule(h2o)
    qcdb.set_options({
        #'cfour_CALC_level': 'CCSD',
        'BASIS': 'cc-pvtz',
        'nwchem_scf__thresh': 12,
        'nwchem_ccsd__thresh': 12
    })

    e, jrec = qcdb.energy('nwc-ccsd', return_wfn=True, molecule=h2o)

    scftot = -76.057112730196  #qz2p -76.062748460117
    mp2tot =  -76.332242848821  #qz2p -76.332940127333
    ccsdcorl = -0.280877944529  #qz2p -0.275705491773
    ccsdtot = -76.337990674725  #qz2p-76.338453951890

    tnm = sys._getframe().f_code.co_name
    assert compare_values(scftot, qcdb.get_variable('hf total energy'), 6, tnm + ' SCF')
    assert compare_values(mp2tot, qcdb.get_variable('mp2 total energy'), 6, tnm + ' MP2')
    assert compare_values(ccsdcorl, qcdb.get_variable('ccsd correlation energy'), 6, tnm + ' CCSD corl')
    assert compare_values(ccsdtot, qcdb.get_variable('ccsd total energy'), 6, tnm + ' CCSD')

    # gradient gives every appearence of working but array should be tested and reorientation to input tested
    e, jrec = qcdb.gradient('nwc-ccsd', return_wfn=True, molecule=h2o)
    assert compare_values(scftot, qcdb.get_variable('hf total energy'), 6, tnm + ' SCF')
    assert compare_values(mp2tot, qcdb.get_variable('mp2 total energy'), 6, tnm + ' MP2')
    assert compare_values(ccsdcorl, qcdb.get_variable('ccsd correlation energy'), 6, tnm + ' CCSD corl')
    assert compare_values(ccsdtot, qcdb.get_variable('ccsd total energy'), 6, tnm + ' CCSD')


#@using_nwchem
#def test_sp_uhf_ccsd():
#    """cfour/sp-uhf-ccsd/input.dat
#    #! single-point CCSD/qz2p on NH2
#
#    """
#    qcdb.set_options({
#        #'cfour_CALC_level': 'CCSD',
#        'BASIS': 'cc-pvtz',
#        'REFerence': 'UHF',
#        'cfour_occupation': [[3,1,1,0],[3,0,1,0]],
#        'cfour_SCF_CONV': 12,
#        'cfour_CC_CONV': 12,
#    })
#
#    qcdb.energy('nwc-ccsd', molecule=nh2)
#
#    scftot = -55.586343126653
#    mp2tot = -55.783567346863
#    ccsdcorl = -0.214548265916
#    ccsdtot = -55.800891392570
#    tnm = sys._getframe().f_code.co_name
#    assert compare_values(scftot, qcdb.get_variable('scf total energy'), 6, tnm + 'SCF')
#    assert compare_values(mp2tot, qcdb.get_variable('mp2 total energy'), 6, tnm + 'MP2')
#    assert compare_values(ccsdcorl, qcdb.get_variable('ccsd correlation energy'), 6, tnm + 'CCSD')
#    assert compare_values(ccsdtot, qcdb.get_variable('ccsd total energy'), 6, tnm + 'CCSD')
#    assert compare_values(ccsdtot, qcdb.get_variable('current energy'), 6, tnm + 'CCSD')
#
#
#@using_nwchem
#def test_sp_rohf_ccsd():
#
#    qcdb.set_options({
#        #cfour_CALC_level': 'CCSD',
#        'cfour_BASIS': 'qz2p',
#        'cfour_REFerence': 'ROHF',
#        'cfour_OCCUPATION': [[3,1,1,0],[3,0,1,0]],
#        'cfour_SCF_CONV': 12,
#        'cfour_PRINT': 2
#    })
#
#    qcdb.energy('nwc-ccsd', molecule=nh2)
#
#    scftot = -55.5847372601
#    smp2corl = -0.002983751786
#    osmp2corl = -0.155770420921
#    ssmp2corl = -0.041785354569
#    mp2corl = -0.200539527276
#    mp2tot = -55.785276787341
#    osccsdcorl = -0.1745752
#    ssccsdcorl = -0.0432743
#    ccsdcorl = -0.217849506326
#    ccsdtot = -55.802586766392
#    tnm = sys._getframe().f_code.co_name
#    compare_values(scftot, qcdb.get_variable('scf total energy'), 6, tnm + 'SCF')
#    # not printed compare_values(smp2corl, qcdb.get_variable('mp2 singles energy'), 6, tnm + 'MP2 singles')
#    # not printed compare_values(osmp2corl, qcdb.get_variable('mp2 opposite-spin correlation energy'), 6, tnm + 'MP2 OS corl')
#    # not printed compare_values(ssmp2corl, qcdb.get_variable('mp2 same-spin correlation energy'), 6, tnm + 'MP2 SS corl')
#    # not printed compare_values(mp2corl, qcdb.get_variable('mp2 correlation energy'), 6, tnm + 'MP2 corl')
#    # not printed compare_values(mp2tot, qcdb.get_variable('mp2 total energy'), 6, tnm + 'MP2')
#    compare_values(osccsdcorl, qcdb.get_variable('ccsd opposite-spin correlation energy'), 6, tnm + 'CCSD OS corl')
#    compare_values(ssccsdcorl, qcdb.get_variable('ccsd same-spin correlation energy'), 6, tnm + 'CCSD SS corl')
#    compare_values(ccsdcorl, qcdb.get_variable('ccsd correlation energy'), 6, tnm + 'CCSD corl')
#    compare_values(ccsdtot, qcdb.get_variable('ccsd total energy'), 6, tnm + 'CCSD')
#    compare_values(ccsdcorl, qcdb.get_variable('current correlation energy'), 6, tnm + 'Current corl')
#    compare_values(ccsdtot, qcdb.get_variable('current energy'), 6, tnm + 'Current')

