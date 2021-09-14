#QCISD on water with sto-3g

import os
import sys

import qcdb

from ..utils import *


def check_qcisd(return_value):
    hf      =   -74.962663062066
    qcisd_tot=   -75.012808319270690
    qcisd_corl=   -0.050145257149164

    assert compare_values(hf, qcdb.variable('HF TOTAL ENERGY'), 5, 'hf ref')
    assert compare_values(qcisd_tot, qcdb.variable('QCISD TOTAL ENERGY'), 5, 'qcisd tot')
    assert compare_values(qcisd_corl, qcdb.variable('QCISD CORRELATION ENERGY'), 5, 'qcisd corl')

@using("nwchem")
def test_1_qcisd():
    h2o = qcdb.set_molecule('''
        H    0.000000000000000   1.079252144093028   1.474611055780858
        O    0.000000000000000   0.000000000000000   0.000000000000000
        H    0.000000000000000   1.079252144093028  -1.47461105578085
        units au ''')
    
    qcdb.set_options({
        'basis' : 'sto-3g',
        'memory': '1500 mb',
        #'scf__e_convergence':   1.0e-10,
        'nwchem_scf__thresh':   1.0e-10,
        'nwchem_scf__tol2e' :   1.0e-10,
        'nwchem_scf__singlet':  True,
        'nwchem_scf__rhf'   :   True,
        'qc_module'         :   'TCE',
        'nwchem_tce__qcisd'  :   True
        })
    val = qcdb.energy('nwc-qcisd')
    check_qcisd(val)
