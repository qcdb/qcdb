import os
import sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))
from utils import *
from addons import *


h2o = qcdb.set_molecule("""
        O
        H 1 R
        H 1 R 2 A
        
        R=0.958
        A=104.5
""")

h2occsdans = -76.119378156918

hf_ccpvdz = -76.026760723833

mp2_sph_ccpvdz = -76.2307777348   
mp2_cart_ccpvdz = -76.2346763972  
mp2_cart_631gs = -76.1990355202   
mp2_sph_631gs = -76.1953470892    


@using_cfour
def test_bpo_1():
    """cfour/kw-1/input.dat"""

    qcdb.set_options({
        'basis': '6-31g',
        'cfour_CALC_level': 'CCSD',
    })

    e, jrec = qcdb.energy('c4-cfour', return_wfn=True, molecule=h2o)

    assert compare_values(h2occsdans, qcdb.get_variable('current energy'), 5, 'Total Energy')
    assert 'CC_PROGRAM           ICCPRO         VCC' in jrec['raw_output'], 'Cfour default VCC'

    qcdb.set_options({
        'cfour_cc_program': 'ecc',
    })

    e, jrec = qcdb.energy('c4-cfour', return_wfn=True, molecule=h2o)

    assert compare_values(h2occsdans, qcdb.get_variable('current energy'), 5, 'Total Energy')
    assert 'CC_PROGRAM           ICCPRO         ECC' in jrec['raw_output'], 'Cfour overwritten by ECC'


@using_cfour
def test_bpo_2():
    """cfour/kw-2/input.dat"""

    qcdb.set_options({
        'basis': '6-31g',
    })

    e, jrec = qcdb.energy('c4-ccsd', return_wfn=True, molecule=h2o)

    assert compare_values(h2occsdans, qcdb.get_variable('current energy'), 5, 'Total Energy')
    assert 'CC_PROGRAM           ICCPRO         ECC' in jrec['raw_output'], 'P4C4 default ECC'

    qcdb.set_options({
        'cfour_cc_program': 'vcc',
    })

    e, jrec = qcdb.energy('c4-ccsd', return_wfn=True, molecule=h2o)

    assert compare_values(h2occsdans, qcdb.get_variable('current energy'), 5, 'Total Energy')
    assert 'CC_PROGRAM           ICCPRO         VCC' in jrec['raw_output'], 'P4C4 overwritten by VCC'


@using_cfour
def test_sph_1():
    """cfour/kw-3/input.dat
    #! Basis set spherical/Cartesian with cfour_basis and cfour_spherical

    """
    qcdb.set_options({
        'cfour_basis': 'pvdz',
    })
   
    qcdb.energy('c4-mp2')
    assert compare_values(mp2_sph_ccpvdz, qcdb.get_variable('CURRENT ENERGY'), 6, 'cc-pVDZ Default Sph')
    
    qcdb.set_options({
        'cfour_basis': '6-31g*',
    })
    
    qcdb.energy('c4-mp2')
    assert compare_values(mp2_sph_631gs, qcdb.get_variable('CURRENT ENERGY'), 6, '6-31G* Default Sph')  
    
    qcdb.set_options({
        'cfour_basis': 'pvdz',
        'cfour_spherical': False,
    })
    
    qcdb.energy('c4-mp2')
    compare_values(mp2_cart_ccpvdz, qcdb.get_variable('CURRENT ENERGY'), 6, 'cc-pVDZ forced Cart')  
    
    qcdb.set_options({
        'cfour_basis': '6-31g*',
    })
    
    qcdb.energy('c4-mp2')
    assert compare_values(mp2_cart_631gs, qcdb.get_variable('CURRENT ENERGY'), 6, '6-31G* forced Cart')  #TEST


@using_cfour
def test_sph_2():

    qcdb.set_options({
        'cfour_basis': 'pvdz',
        'puream': False,
    })
    
    qcdb.energy('c4-mp2')
    compare_values(mp2_cart_ccpvdz, qcdb.get_variable('CURRENT ENERGY'), 6, 'cc-pVDZ forced Cart')  
    
    qcdb.set_options({
        'cfour_basis': '6-31g*',
    })
    
    qcdb.energy('c4-mp2')
    assert compare_values(mp2_cart_631gs, qcdb.get_variable('CURRENT ENERGY'), 6, '6-31G* forced Cart')  #TEST


@using_cfour
def test_sph_3():

    #qcdb.set_options({
    #    'basis': 'cc-pvdz',
    #})
   
    #qcdb.energy('c4-mp2')
    #assert compare_values(mp2_sph_ccpvdz, qcdb.get_variable('CURRENT ENERGY'), 6, 'cc-pVDZ Default Sph')
    
    qcdb.set_options({
        'basis': '6-31g*',
    })
    
    # I'm suspicious of this one
    qcdb.energy('c4-mp2')
    assert compare_values(mp2_sph_631gs, qcdb.get_variable('CURRENT ENERGY'), 6, '6-31G* Default Sph')  
    
    #qcdb.set_options({
    #    'basis': 'cc-pvdz',
    #    'cfour_spherical': False,
    #})
    #
    #qcdb.energy('c4-mp2')
    #compare_values(mp2_cart_ccpvdz, qcdb.get_variable('CURRENT ENERGY'), 6, 'cc-pVDZ forced Cart')  
    #
    #qcdb.set_options({
    #    'basis': '6-31g*',
    #})
    #
    #qcdb.energy('c4-mp2')
    #assert compare_values(mp2_cart_631gs, qcdb.get_variable('CURRENT ENERGY'), 6, '6-31G* forced Cart')  #TEST


@using_cfour
def test_sph_4lab():

    qcdb.set_options({
        'basis': 'cc-pvdz',
        'puream': False,
    })
    
    qcdb.energy('c4-mp2')
    compare_values(mp2_cart_ccpvdz, qcdb.get_variable('CURRENT ENERGY'), 6, 'cc-pVDZ forced Cart')  
    
    qcdb.set_options({
        'basis': '6-31g*',
    })
    
    qcdb.energy('c4-mp2')
    assert compare_values(mp2_cart_631gs, qcdb.get_variable('CURRENT ENERGY'), 6, '6-31G* forced Cart')  #TEST


@using_cfour
def test_sph_4():
    """cfour/kw-4/input.dat
    #! Basis set spherical/Cartesian with basis and puream

    """
    qcdb.set_options({
        'basis': 'cc-pvdz',
    })
    
    qcdb.energy('c4-mp2')
    assert compare_values(mp2_sph_ccpvdz, qcdb.get_variable('CURRENT ENERGY'), 6, 'cc-pVDZ Default Sph')
    
    qcdb.set_options({
        'basis': '6-31g*',
        'puream': False,
    })
    
    qcdb.energy('c4-mp2')
    assert compare_values(mp2_cart_631gs, qcdb.get_variable('CURRENT ENERGY'), 6, '6-31G* Default Cart')  
    
    qcdb.set_options({
        'basis': 'cc-pvdz',
        'puream': False,
    })
    
    qcdb.energy('c4-mp2')
    assert compare_values(mp2_cart_ccpvdz, qcdb.get_variable('CURRENT ENERGY'), 6, 'cc-pVDZ forced Cart') 
    
    qcdb.set_options({
        'basis': '6-31g*',
        'puream': True,
    })
    
    qcdb.energy('c4-mp2')
    assert compare_values(mp2_sph_631gs, qcdb.get_variable('CURRENT ENERGY'), 6, '6-31G* forced Sph')  


@using_cfour
def test_sph_5():
    """cfour/kw-5/input.dat
    #! Basis set spherical/Cartesian with basis and cfour_spherical

    """
    qcdb.set_options({
        'basis': 'cc-pvdz',
        'cfour_spherical': True,
    })

    qcdb.energy('c4-mp2')
    assert compare_values(mp2_sph_ccpvdz, qcdb.get_variable('CURRENT ENERGY'), 6, 'cc-pVDZ Default Sph')

    qcdb.set_options({
        'basis': '6-31g*',
        'cfour_spherical': False,
    })

    qcdb.energy('c4-mp2')
    assert compare_values(mp2_cart_631gs, qcdb.get_variable('CURRENT ENERGY'), 6, '6-31G* Default Cart')


@using_cfour
def test_conflict_1():
    """Conflicting option requirements btwn user (ANGSTROM) and driver (BOHR)"""

    qcdb.set_options({
        'cfour_units': 'angstrom',
        'cfour_basis': 'pvdz',
        'puream': True,
    })
    
    with pytest.raises(qcdb.OptionReconciliationError):
        qcdb.energy('c4-mp2')

@using_cfour
def test_sph_6():
    """cfour/kw-6/input.dat

    """
    qcdb.set_options({
        'cfour_units': 'boHR',
        'cfour_basis': 'pvdz',
        'puream': True,
    })
    
    qcdb.energy('c4-mp2')
    assert compare_values(mp2_sph_ccpvdz, qcdb.get_variable('CURRENT ENERGY'), 6, 'cc-pVDZ Default Sph')  #TEST
    
    qcdb.set_options({
        'cfour_basis': '6-31g*',
        'puream': True,
    })
    
    qcdb.energy('c4-mp2')
    assert compare_values(mp2_sph_631gs, qcdb.get_variable('CURRENT ENERGY'), 6, '6-31G* Default Sph')  #TEST
    
    qcdb.set_options({
        'cfour_basis': 'pvdz',
        'puream': False,
    })
    
    qcdb.energy('c4-mp2')
    assert compare_values(mp2_cart_ccpvdz, qcdb.get_variable('CURRENT ENERGY'), 6, 'cc-pVDZ forced Cart')  #TEST
    
    qcdb.set_options({
        'cfour_basis': '6-31g*',
        'puream': False,
    })
    
    qcdb.energy('c4-mp2')
    assert compare_values(mp2_cart_631gs, qcdb.get_variable('CURRENT ENERGY'), 6, '6-31G* forced Cart')  #TEST


@using_cfour
def test_conv_1():
    """cfour/kw-7/input.dat
    #! Translating psi4 options to cfour, part i

    """
    qcdb.set_options({
        'basis': 'cc-pvdz',
    })

    e, jrec = qcdb.energy('c4-scf', return_wfn=True)
    assert compare_values(hf_ccpvdz, qcdb.get_variable('current energy'), 5, '[a] Total Energy')
    assert 'SCF_CONV             ISCFCV          10D-  7' in jrec['raw_output'], 'SCF_CONV 7 default'

    qcdb.set_options({
        'cfour_scf_conv': 6,
    })

    e, jrec = qcdb.energy('c4-scf', return_wfn=True)
    assert compare_values(hf_ccpvdz, qcdb.get_variable('current energy'), 5, '[b] Total Energy')
    assert 'SCF_CONV             ISCFCV          10D-  6' in jrec['raw_output'], 'SCF_CONV 6 overwritten'

    qcdb.set_options({
        'cfour_scf_conv': 9,
        'd_convergence': 8,
    })

    e, jrec = qcdb.energy('c4-scf', return_wfn=True)
    assert compare_values(hf_ccpvdz, qcdb.get_variable('current energy'), 5, '[c] Total Energy')
    assert 'SCF_CONV             ISCFCV          10D-  9' in jrec['raw_output'], 'SCF_CONV 9 cfour kw trumps scf kw'


@using_cfour
def test_conv_2():
    """cfour/kw-8/input.dat
    #! Translating psi4 options to cfour, part ii

    """
    qcdb.set_options({
        'basis': 'cc-pvdz',
    })

    e, jrec = qcdb.energy('c4-scf', return_wfn=True)
    assert compare_values(hf_ccpvdz, qcdb.get_variable('current energy'), 5, '[a] Total Energy')
    assert 'SCF_CONV             ISCFCV          10D-  7' in jrec['raw_output'], 'SCF_CONV 7 default'

    qcdb.set_options({
        'd_convergence': 5,
    })

    e, jrec = qcdb.energy('c4-scf', return_wfn=True)
    assert compare_values(hf_ccpvdz, qcdb.get_variable('current energy'), 5, '[b] Total Energy')
    assert 'SCF_CONV             ISCFCV          10D-  5' in jrec['raw_output'], 'SCF_CONV 5 scf kw overwritten'

    qcdb.set_options({
        'scf__d_convergence': 5e-7,
    })

    e, jrec = qcdb.energy('c4-scf', return_wfn=True)
    assert compare_values(hf_ccpvdz, qcdb.get_variable('current energy'), 5, '[c] Total Energy')
    assert 'SCF_CONV             ISCFCV          10D-  7' in jrec['raw_output'], 'P4C4 default SCF_CONV 7'

    qcdb.set_options({
        'cfour_scf_conv': 6,
        'd_convergence': 8,
    })

    e, jrec = qcdb.energy('c4-scf', return_wfn=True)
    assert compare_values(hf_ccpvdz, qcdb.get_variable('current energy'), 5, '[d] Total Energy')
    assert 'SCF_CONV             ISCFCV          10D-  6' in jrec['raw_output'], 'SCF_CONV 6 cfour kw trumps scf kw'


def test_clash_mult():
    """Conflicting option requirements btwn user (3) and driver (1)"""

    qcdb.set_molecule("""
H
H 1 0.7
    """)

    qcdb.set_options({
        'basis': '6-31g',
        'cfour_multiplicity': 3,         # clash with implicit singlet in molecule {} above
    })

    with pytest.raises(qcdb.OptionReconciliationError):
        qcdb.energy('c4-scf')


def test_clash_puream():
    qcdb.set_options({
        'basis': 'cc-pvdz',
        'cfour_spherical': False,
    })

    qcdb.energy('c4-mp2')
    assert compare_values(cart_ccpvdz, qcdb.get_variable('CURRENT ENERGY'), 6, 'cc-pVDZ forced Cart')  #TEST


def test_clash_method():
    """Conflicting option requirements btwn user (CCSD) and driver (SCF"""

    qcdb.set_molecule("""
H
H 1 0.7
    """)

    qcdb.set_options({
        'basis': '6-31g',
        'cfour_calc_level': 'ccsd',        # clash with 'c4-scf' below
    })

    with pytest.raises(qcdb.OptionReconciliationError):
        qcdb.energy('c4-scf')

# memory ok indep but some problems in seq?

#def test_memory_1():
#    """driver mode default"""
#    qcdb.set_options({
#        'basis': '6-31g',
#    })
#
#    e, jrec = qcdb.energy('c4-scf', molecue=h2o, return_wfn=True)
#    assert 'MEMORY_SIZE          IMEMSZ          87500000           words' in jrec['raw_output']
#    assert 'MEM_UNIT             IMEMU          MB' in jrec['raw_output']
#
#
#def test_memory_2():
#    qcdb.set_options({
#        'memory': '300 mb',
#        'basis': '6-31g',
#    })
#
#    e, jrec = qcdb.energy('c4-scf', molecue=h2o, return_wfn=True)
#    assert 'MEMORY_SIZE          IMEMSZ          37500000           words' in jrec['raw_output']
#    assert 'MEM_UNIT             IMEMU          MB' in jrec['raw_output']
#
#
#def test_memory_3():
#    qcdb.set_options({
#        'cfour_memory_size': '300',
#        'cfour_mem_unit': 'mb',
#        'basis': '6-31g',
#    })
#
#    e, jrec = qcdb.energy('c4-scf', molecue=h2o, return_wfn=True)
#    assert 'MEMORY_SIZE          IMEMSZ          37500000           words' in jrec['raw_output']
#    assert 'MEM_UNIT             IMEMU          MB' in jrec['raw_output']
#
#
#def test_memory_4():
#    qcdb.set_options({
#        'memory': '600 mb',
#        'cfour_memory_size': '300',
#        'cfour_mem_unit': 'mb',
#        'basis': '6-31g',
#    })
#
#    e, jrec = qcdb.energy('c4-scf', molecue=h2o, return_wfn=True)
#    assert 'MEMORY_SIZE          IMEMSZ          37500000           words' in jrec['raw_output']
#    assert 'MEM_UNIT             IMEMU          MB' in jrec['raw_output']
#    # SHOULD FAIL
#
#
#def test_memory_5():
#    """runner mode default"""
#    qcdb.set_options({
#        'translate_qcdb': 'no',
#        'basis': '6-31g',
#    })
#
#    e, jrec = qcdb.energy('c4-scf', molecue=h2o, return_wfn=True)
#    assert 'MEMORY_SIZE          IMEMSZ         100000000' in jrec['raw_output']
#    assert 'MEM_UNIT             IMEMU          INTEGERWORDS' in jrec['raw_output']




def hide_test_memory_1():
    qcdb.set_molecule("""
H
H 1 0.7
    """)

    qcdb.set_options({
        'basis': '6-31g',
        #'memory': '300 mb',
        #'cfour_memory_size': 100000000  # clash with 300 mb above
    })


# cRAZY MEM AMT
#def test_memory_1():
#    qcdb.set_molecule("""
#H
#H 1 0.7
#    """)
#
#    qcdb.set_options({
#        'basis': '6-31g',
#        'memory': '300 mb',
#        'cfour_memory_size': 100000000  # clash with 300 mb above
#    })
#
#    qcdb.energy('c4-scf')

