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

mp2_sph_ccpvdz = -76.2307777348   
mp2_cart_ccpvdz = -76.2346763972  
mp2_cart_631gs = -76.1990355202   
mp2_sph_631gs = -76.1953470892    


@using_cfour
def test_sph_1():
    tnm = sys._getframe().f_code.co_name

    qcdb.set_options({
        'basis': 'cc-pvdz',
    })
   
    qcdb.energy('c4-mp2')
    assert compare_values(mp2_sph_ccpvdz, qcdb.get_variable('CURRENT ENERGY'), 6, tnm)


@using_cfour
def test_sph_2():
    tnm = sys._getframe().f_code.co_name

    qcdb.set_options({
        'basis': 'cc-pvdz',
        'puream': True,
    })
   
    qcdb.energy('c4-mp2')
    assert compare_values(mp2_sph_ccpvdz, qcdb.get_variable('CURRENT ENERGY'), 6, tnm)


@using_cfour
def test_sph_3():
    tnm = sys._getframe().f_code.co_name

    qcdb.set_options({
        'basis': 'cc-pvdz',
        'puream': False,
    })
   
    qcdb.energy('c4-mp2')
    assert compare_values(mp2_cart_ccpvdz, qcdb.get_variable('CURRENT ENERGY'), 6, tnm)


@using_cfour
def test_sph_4():
    tnm = sys._getframe().f_code.co_name

    qcdb.set_options({
        'basis': '6-31g*',
    })
   
    qcdb.energy('c4-mp2')
    assert compare_values(mp2_cart_631gs, qcdb.get_variable('CURRENT ENERGY'), 6, tnm)


@using_cfour
def test_sph_5():
    tnm = sys._getframe().f_code.co_name

    qcdb.set_options({
        'basis': '6-31g*',
        'puream': True,
    })
   
    qcdb.energy('c4-mp2')
    assert compare_values(mp2_sph_631gs, qcdb.get_variable('CURRENT ENERGY'), 6, tnm)


@using_cfour
def test_sph_6():
    tnm = sys._getframe().f_code.co_name

    qcdb.set_options({
        'basis': '6-31g*',
        'puream': False,
    })
   
    qcdb.energy('c4-mp2')
    assert compare_values(mp2_cart_631gs, qcdb.get_variable('CURRENT ENERGY'), 6, tnm)


@using_cfour
def test_sph_7():
    tnm = sys._getframe().f_code.co_name

    qcdb.set_options({
        'cfour_basis': '6-31g*',
    })
   
    qcdb.energy('c4-mp2')
    assert compare_values(mp2_sph_631gs, qcdb.get_variable('CURRENT ENERGY'), 6, tnm)


@using_cfour
def test_sph_8():
    tnm = sys._getframe().f_code.co_name

    qcdb.set_options({
        'cfour_basis': '6-31g*',
        'cfour_spherical': True,
    })
   
    qcdb.energy('c4-mp2')
    assert compare_values(mp2_sph_631gs, qcdb.get_variable('CURRENT ENERGY'), 6, tnm)


@using_cfour
def test_sph_9():
    tnm = sys._getframe().f_code.co_name

    qcdb.set_options({
        'cfour_basis': '6-31g*',
        'cfour_spherical': False,
    })
   
    qcdb.energy('c4-mp2')
    assert compare_values(mp2_cart_631gs, qcdb.get_variable('CURRENT ENERGY'), 6, tnm)


@using_cfour
def test_sph_10():
    tnm = sys._getframe().f_code.co_name

    qcdb.set_options({
        'cfour_basis': '6-31g*',
        'puream': True,
    })
   
    qcdb.energy('c4-mp2')
    assert compare_values(mp2_sph_631gs, qcdb.get_variable('CURRENT ENERGY'), 6, tnm)


@using_cfour
def test_sph_11():
    tnm = sys._getframe().f_code.co_name

    qcdb.set_options({
        'cfour_basis': '6-31g*',
        'puream': False,
    })
   
    qcdb.energy('c4-mp2')
    assert compare_values(mp2_cart_631gs, qcdb.get_variable('CURRENT ENERGY'), 6, tnm)


@using_cfour
def test_sph_12():
    tnm = sys._getframe().f_code.co_name

    qcdb.set_options({
        'cfour_basis': '6-31g*',
        'puream': True,
        'cfour_spherical': False,
    })
   
    qcdb.energy('c4-mp2')
    assert compare_values(mp2_cart_631gs, qcdb.get_variable('CURRENT ENERGY'), 6, tnm)


@using_cfour
def test_sph_13():
    tnm = sys._getframe().f_code.co_name

    qcdb.set_options({
        'cfour_basis': '6-31g*',
        'puream': False,
        'cfour_spherical': True,
    })
   
    qcdb.energy('c4-mp2')
    assert compare_values(mp2_sph_631gs, qcdb.get_variable('CURRENT ENERGY'), 6, tnm)
    

@using_cfour
def hide_test_sph_X():
    tnm = sys._getframe().f_code.co_name

    qcdb.set_options({
        'basis': 'cc-pvdz',
        'basis': '6-31g*',
        'cfour_basis': 'pvdz',
        'cfour_basis': '6-31g*',
        'puream': True,
        'puream': False,
        'cfour_spherical': True,
        'cfour_spherical': False,
    })
   
    qcdb.energy('c4-mp2')
    assert compare_values(mp2_sph_ccpvdz, qcdb.get_variable('CURRENT ENERGY'), 6, tnm)
    assert compare_values(mp2_cart_ccpvdz, qcdb.get_variable('CURRENT ENERGY'), 6, tnm)
    assert compare_values(mp2_sph_631gs, qcdb.get_variable('CURRENT ENERGY'), 6, tnm)
    assert compare_values(mp2_cart_631gs, qcdb.get_variable('CURRENT ENERGY'), 6, tnm)
    

