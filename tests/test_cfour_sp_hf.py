import os
import sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))
from .utils import *
from .addons import *


@using_cfour
def test_sp_uhf_scf():
    """cfour/sp-uhf-scf/input.dat
    UHF-SCF energy calculation
    #global testing

    """
    nh2 = qcdb.set_molecule("""
        0 2
        N
        H 1 R
        H 1 R 2 A
        
        R=1.008
        A=105.0
    """)

    qcdb.set_options({
        #'cfour_CALC_level': 'HF',
        'cfour_BASIS': 'qz2p',
        'cfour_REFerence': 'UHF',
        'cfour_occupation': [[3,1,1,0],[3,0,1,0]],
        'cfour_SCF_CONV': 12
    })

    qcdb.energy('c4-hf')

    ans = -55.5893469688
    tnm = sys._getframe().f_code.co_name
    assert compare_values(ans, qcdb.get_variable('scf total energy'), 6, tnm + 'SCF')  #TEST
    assert compare_values(ans, qcdb.get_variable('current energy'), 6, tnm + 'Current')  #TEST
    assert compare_values(ans, qcdb.get_variable('current reference energy'), 6, tnm + 'Current ref')  #TEST


@using_cfour
def test_sp_rhf_scf_a():
    """cfour/sp-rhf-scf/input.dat 
    #! single-point HF/qz2p on water

    """
    h2o = qcdb.set_molecule("""
        O
        H 1 R
        H 1 R 2 A
        
        R=0.958
        A=104.5
    """)

    qcdb.set_options({
        'cfour_BASIS': 'qz2p',
        'd_convergence': 12
    })

    e, jrec = qcdb.energy('c4-scf', return_wfn=True)

    ans = -76.0627484601
    tnm = sys._getframe().f_code.co_name
    assert compare_values(ans, qcdb.get_variable('scf total energy'), 6, tnm + ' SCF')
    assert compare_values(ans, qcdb.get_variable('current energy'), 6, tnm + ' SCF')
    assert compare_values(ans, jrec['qcvars']['SCF TOTAL ENERGY'].data, 6, tnm)
    assert compare_values(ans, jrec['qcvars']['CURRENT ENERGY'].data, 6, tnm)


@using_cfour
def test_sp_rhf_scf_b():
    """cfour/sp-rhf-scf/input.dat 
    #! single-point HF/qz2p on water

    """
    h2o = qcdb.set_molecule("""
        O
        H 1 R
        H 1 R 2 A
        
        R=0.958
        A=104.5
    """)

    qcdb.set_options({
        'cfour_calc_level' : 'hf',
        'cfour_basis': 'qz2p',
        'cfour_scf_conv': 12
    })

    e, jrec = qcdb.energy('c4-cfour', return_wfn=True)

    ans = -76.0627484601
    tnm = sys._getframe().f_code.co_name
    assert compare_values(ans, qcdb.get_variable('scf total energy'), 6, tnm + ' SCF')
    assert compare_values(ans, qcdb.get_variable('current energy'), 6, tnm + ' SCF')
    assert compare_values(ans, jrec['qcvars']['SCF TOTAL ENERGY'].data, 6, tnm)
    assert compare_values(ans, jrec['qcvars']['CURRENT ENERGY'].data, 6, tnm)


if __name__ == '__main__':
    test_sp_uhf_scf()
