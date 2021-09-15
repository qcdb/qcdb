#! single-point CCSD/cc-pvdz on water
import os
import sys

import qcdb
import qcelemental as qcel

from ..utils import *

import numpy as np


@using("nwchem")
def test_1_ccsd():
    h2o = qcdb.set_molecule('''
        O     0.000000000000    0.000000000000   0.22138519
        H     0.000000000000   -1.43013023      -0.88554075
        H     0.000000000000    1.43013023      -0.88554075
        ''')
    qcdb.set_options({
        'basis': 'cc-pvdz',
        'nwchem_scf__rhf': True,
        'nwchem_scf__thresh': 1.0e-12,
        'nwchem_tce__ccsd': True,
        'qc_module'  : 'tce',
        'nwchem_tce__dipole' : True,
        'nwchem_tce__thresh': 1.0e-12,
    })
    print('     Testing ccsd...')
    val = qcdb.energy('nwc-ccsd')

    ref = -75.645078266314
    ccsd_corl = -0.313405482575930
    ccsd_tot = -75.958483748890202
#   ccsd_dip_z   =   2.1073264  # ATL
    ccsd_dip = np.array([0.0000000, 0.000000, 1.1104046])

    assert compare_values(ref, qcdb.variable('HF TOTAL ENERGY'), 5, 'hf ref')
    assert compare_values(ccsd_corl, qcdb.variable('CCSD CORRELATION ENERGY'), 5, 'ccsd corl')
    assert compare_values(ccsd_tot, qcdb.variable('CCSD TOTAL ENERGY'), 5, 'ccsd total')
    assert compare_values(ccsd_dip, qcdb.variable('CCSD DIPOLE'), 5, 'ccsd dipole')
