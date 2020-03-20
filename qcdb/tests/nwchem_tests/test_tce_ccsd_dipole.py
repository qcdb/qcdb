#! single-point CCSD/cc-pvdz on water
import os
import sys

import qcdb

from ..addons import *
from ..utils import *

import numpy as np

h2o = qcdb.set_molecule('''
        O     0.000000000000    0.000000000000   0.22138519
        H     0.000000000000   -1.43013023      -0.88554075
        H     0.000000000000    1.43013023      -0.88554075
        ''')

print(h2o)


def check_ccsd(return_value):
    ref = -75.645078266314
    ccsd_corl = -0.313405482575930
    ccsd_tot = -75.958483748890202
    x   =   0.00000000
    y   =   0.00000000
    z   =   2.1073264
    #ccsd_dip = np.array[0.0000000, 0.000000, 2.1073264]

    assert compare_values(ref, qcdb.variable('HF TOTAL ENERGY'), 5, 'hf ref')
    assert compare_values(ccsd_corl, qcdb.variable('CCSD CORRELATION ENERGY'), 5, 'ccsd corl')
    assert compare_values(ccsd_tot, qcdb.variable('CCSD TOTAL ENERGY'), 5, 'ccsd total')
    assert compare_values(x, qcdb.variable('CCSD DIPOLE X'), 5, 'dipole x')
    assert compare_values(y, qcdb.variable('CCSD DIPOLE Y'), 5, 'dipole y')
    assert compare_values(z, qcdb.variable('CCSD DIPOLE Z'), 5, 'dipole z')

@using_nwchem
def test_1_ccsd():
    qcdb.set_options({
        'basis': 'cc-pvdz',
        'memory': '600 mb',
        'nwchem_scf__rhf': True,
        'nwchem_scf__thresh': 1.0e-12,
        'nwchem_tce__ccsd': True,
        'qc_module'  : 'tce',
        'nwchem_tce__dipole' : True,
        'nwchem_tce__thresh': 1.0e-12,
    })
    print('     Testing ccsd...')
    val = qcdb.energy('nwc-ccsd')
    check_ccsd(val)
