#LCCD

import os
import sys

import qcdb

from ..utils import *


def check_lccd(return_value):
    hf      =   -74.962663062066
    lccd_tot=   -75.013238394202133
    lccd_corl=   -0.050575332136568

    assert compare_values(hf, qcdb.variable('HF TOTAL ENERGY'), 5, 'hf ref')
    assert compare_values(lccd_tot, qcdb.variable('LCCD TOTAL ENERGY'), 5, 'lccd tot')
    assert compare_values(lccd_corl, qcdb.variable('LCCD CORRELATION ENERGY'), 5, 'lccd corl')

@using("nwchem")
def test_1_lccd():
    h2o = qcdb.set_molecule('''
    H    0.000000000000000   1.079252144093028   1.474611055780858
    O    0.000000000000000   0.000000000000000   0.000000000000000
    H    0.000000000000000   1.079252144093028  -1.47461105578085
       units au ''')

    qcdb.set_options({
        'basis' : 'sto-3g',
        #'scf__e_convergence':   1.0e-10,
        'nwchem_scf__thresh':   1.0e-10,
        'nwchem_scf__tol2e' :   1.0e-10,
        'nwchem_scf__singlet':  True,
        'nwchem_scf__rhf'   :   True,
        'qc_module'         :   'TCE',
        'nwchem_tce__lccd'  :   True
        })
    val = qcdb.energy('nwc-lccd')
    check_lccd(val)

def check_lccsd(return_value):
    hf      =   -74.962663062066
    lccsd_tot   = -75.01355462478978
    lccsd_corl  =   -0.050891562718417

    assert compare_values(hf, qcdb.variable('HF TOTAL ENERGY'), 5, 'hf ref')
    assert compare_values(lccsd_tot, qcdb.variable('LCCSD TOTAL ENERGY'), 5, 'lccd tot')
    assert compare_values(lccsd_corl, qcdb.variable('LCCSD CORRELATION ENERGY'), 5, 'lccd corl')

@using("nwchem")
def test_2_lccsd():
    h2o = qcdb.set_molecule('''
    H    0.000000000000000   1.079252144093028   1.474611055780858
    O    0.000000000000000   0.000000000000000   0.000000000000000
    H    0.000000000000000   1.079252144093028  -1.47461105578085
       units au ''')

    qcdb.set_options({
        'basis' : 'sto-3g',
        #'scf__e_convergence':   1.0e-10,
        'nwchem_scf__thresh':   1.0e-10,
        'nwchem_scf__tol2e' :   1.0e-10,
        'nwchem_scf__singlet':  True,
        'nwchem_scf__rhf'   :   True,
        'qc_module'         :   'TCE',
        'nwchem_tce__lccsd'  :   True
        })
    val = qcdb.energy('nwc-lccsd')
    check_lccsd(val)
