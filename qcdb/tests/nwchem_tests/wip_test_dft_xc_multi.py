#! pbe/sto-3g H2O DFT energy
import os
import sys

import pytest

import qcdb

from ..utils import *


@pytest.fixture
def h2o():
    h2o = qcdb.set_molecule('''
        O     0.00000000    0.000000000   0.11726921
        H     0.75698224    0.000000000   -0.46907685
        H     -0.75698224   0.000000000   -0.46907685
        ''')
    return h2o

def check_pw91(return_value):

    pw91    =   -76.391706600968

    assert compare_values(pw91, qcdb.variable('DFT TOTAL ENERGY'), 5, 'xc pw91')

@using("nwchem")
def test_1_pw91(h2o):
    qcdb.set_options({
        'basis' : 'cc-pvdz',
        'nwchem_dft__xc': 'xperdew91 perdew91'
        })
    val = qcdb.energy('nwc-pw91', molecule=h2o)
    check_pw91(val)

def check_pbe96(return_value):

    pbe96    =  -76.358574253086

    assert compare_values(pbe91, qcdb.variable('DFT TOTAL ENERGY'), 5, 'xc pbe96')

@using("nwchem")
def test_2_pbe96(h2o):
    qcdb.set_options({
        'basis' : 'cc-pvdz',
        'nwchem_dft__xc': 'xpbe96 perdew91'
        })
    val = qcdb.energy('nwc-pbe96', molecule=h2o)
    check_pbe96(val)

def check_bp86(return_value):

    bp86    =  -76.421591007635

    assert compare_values(bp86, qcdb.variable('DFT TOTAL ENERGY'), 5, 'xc bp86')

@using("nwchem")
def test_3_bp86(h2o):
    qcdb.set_options({
        'basis' : 'cc-pvdz',
        'nwchem_dft__xc': 'becke88 perdew86'
        })
    val = qcdb.energy('nwc-bp86', molecule=h2o)
    check_bp86(val)

def check_bp91(return_value):

    bp91    =  -76.413386121350

    assert compare_values(bp91, qcdb.variable('DFT TOTAL ENERGY'), 5, 'xc bp91')

@using("nwchem")
def test_4_bp91(h2o):
    qcdb.set_options({
        'basis' : 'cc-pvdz',
        'nwchem_dft__xc': 'becke88 perdew91'
        })
    val = qcdb.energy('nwc-bp91', molecule=h2o)
    check_bp91(val)

def check_blyp(return_value):

    blyp    =   -76.391706600968

    assert compare_values(blyp, qcdb.variable('DFT TOTAL ENERGY'), 5, 'xc blyp')

@using("nwchem")
def test_5_blyp(h2o):
    qcdb.set_options({
        'basis' : 'cc-pvdz',
        'nwchem_dft__xc': 'becke88 lyp'
        })
    val = qcdb.energy('nwc-blyp', molecule=h2o)
    check_blyp(val)

def check_b97(return_value):

    b97    =  -76.395774842062

    assert compare_values(bp97, qcdb.variable('DFT TOTAL ENERGY'), 5, 'xc b97')

@using("nwchem")
def test_6_b97(h2o):
    qcdb.set_options({
        'basis' : 'cc-pvdz',
        'nwchem_dft__xc': 'becke97 hfexch 0.1943'
        })
    val = qcdb.energy('nwc-b97', molecule=h2o)
    check_b97(val)

def check_mpw1pw(return_value):

    mpw1pw    =  -76.402259124777

    assert compare_values(mpw1pw, qcdb.variable('DFT TOTAL ENERGY'), 5, 'mpw1pw')

@using("nwchem")
def test_7_mpw1pw(h2o):
    qcdb.set_options({
        'basis' : 'cc-pvdz',
        'nwchem_dft__xc': 'mpw91 0.75 hfexch 0.25 perdew91'
        })
    val = qcdb.energy('nwc-mpw1pw', molecule=h2o)
    check_mpw1pw(val)

def check_mpwlyp1m(return_value):

    mpwlyp1m    = -76.845270838137

    assert compare_values(mpwlyp1m, qcdb.variable('DFT TOTAL ENERGY'), 5, 'mpwly1m')

@using("nwchem")
def test_8_mpwlyp1m(h2o):
    qcdb.set_options({
        'basis' : 'cc-pvdz',
        'nwchem_dft__xc': 'mpw91 0.95 hfexch 0.05 lyp'
        })
    val = qcdb.energy('nwc-mpwlyp1m', molecule=h2o)
    check_mpwlyp1m(val)

def check_mpwlyp1w(return_value):

    mpwlyp1w    =  -76.435413416144

    assert compare_values(mpwlyp1w, qcdb.variable('DFT TOTAL ENERGY'), 5, 'mpwlyp1w')

@using("nwchem")
def test_9_mpwlyp1w(h2o):
    qcdb.set_options({
        'basis' : 'cc-pvdz',
        'nwchem_dft__xc': 'mpw91 vwn_5 0.12 lyp 0.88'
        })
    val = qcdb.energy('nwc-mpwlyp1w', molecule=h2o)
    check_mpwlyp1w(val)

def check_b1lyp(return_value):

    b1lyp    =  -76.390361299310

    assert compare_values(b1lyp, qcdb.variable('DFT TOTAL ENERGY'), 5, 'b1lyp')

@using("nwchem")
def test_10_b1lyp(h2o):
    qcdb.set_options({
        'basis' : 'cc-pvdz',
        'nwchem_dft__xc': 'hfexch 0.25 becke88 0.75 lyp',
        })
    val = qcdb.energy('nwc-b1lyp', molecule=h2o)
    check_b1lyp(val)

def check_b1pw91(return_value):

    b1pw91    =  -76.404657807949

    assert compare_values(b1pw91, qcdb.variable('DFT TOTAL ENERGY'), 5, 'b1pw91')

@using("nwchem")
def test_11_b1pw91(h2o):
    qcdb.set_options({
        'basis' : 'cc-pvdz',
        'nwchem_dft__xc': 'hfexch 0.25 becke88 0.75 perdew91'
        })
    val = qcdb.energy('nwc-b1pw91', molecule=h2o)
    check_b1pw91(val)

def check_b2plyp(return_value):

    b2plyp    =  -76.289142410551

    assert compare_values(b2plyp, qcdb.variable('DFT TOTAL ENERGY'), 5, 'b2plyp')

@using("nwchem")
def test_12_b2plyp(h2o):
    qcdb.set_options({
        'basis' : 'cc-pvdz',
        'nwchem_dft__xc': 'hfexch 0.53 becke88 0.47 lyp 0.73'
        })
    val = qcdb.energy('nwc-b2plyp', molecule=h2o)
    check_b2plyp(val)

def check_b3lyp5(return_value):

    b3lyp5    =  -76.384402340966
    assert compare_values(b3lyp5, qcdb.variable('DFT TOTAL ENERGY'), 5, 'xc b97')

@using("nwchem")
def test_13_b3lyp5(h2o):
    qcdb.set_options({
        'basis' : 'cc-pvdz',
        'nwchem_dft__xc': 'hfexch 0.2 slater 0.8 becke88 0.72 vwn_5 0.190 lyp 0.81'
        })
    val = qcdb.energy('nwc-b3lyp5', molecule=h2o)
    check_b3lyp5(val)

def check_b86bpbe(return_value):

    b86bpbe    = -76.399128661846 

    assert compare_values(b86bpbe, qcdb.variable('DFT TOTAL ENERGY'), 5, 'xc b86cpbe')

@using("nwchem")
def test_14_b86bpbe(h2o):
    qcdb.set_options({
        'basis' : 'cc-pvdz',
        'nwchem_dft__xc': 'becke86 cpbe96'
        })
    val = qcdb.energy('nwc-b86bpbe', molecule=h2o)
    check_b86bpbe(val)

def check_b97_1(return_value):

    b97_1    = -85.701159230853

    assert compare_values(b97_1, qcdb.variable('DFT TOTAL ENERGY'), 5, 'xc b97-1')

@using("nwchem")
def test_15_b97_1(h2o):
    qcdb.set_options({
        'basis' : 'cc-pvdz',
        'nwchem_dft__xc': 'becke97gga1 hfexch'
        })
    val = qcdb.energy('nwc-b97-1', molecule=h2o)
    check_b97_1(val)

