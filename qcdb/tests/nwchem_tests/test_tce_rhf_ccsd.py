#! single-point CCSD/cc-pvdz on water
import os
import sys
from ..utils import *
from ..addons import *
import qcdb


def check_ccsd(return_value):
    ref = -76.026760733967
    ccsd_corl = -0.213341251859034
    ccsd_tot = -76.240101985825859

    assert compare_values(ref, qcdb.variable('HF TOTAL ENERGY'), 5, 'hf ref')
    assert compare_values(ccsd_corl, qcdb.variable('CCSD CORRELATION ENERGY'), 5, 'ccsd corl')
    assert compare_values(ccsd_tot, qcdb.variable('CCSD TOTAL ENERGY'), 5, 'ccsd total')


@using_nwchem
def test_1_ccsd():
    h2o = qcdb.set_molecule('''
        O     0.000000000000    0.000000000000   -0.065638538099
        H     0.000000000000   -0.757480611647    0.520865616174
        H     0.000000000000    0.757480611647    0.520865616174
        ''')

    qcdb.set_options({
        'basis': 'cc-pvdz',
        'memory': '600 mb',
        'nwchem_scf__rhf': True,
        'nwchem_scf__thresh': 1.0e-12,
        'nwchem_tce__ccsd': True,
        'qc_module'  : 'tce',
        #'nwchem_tce__dipole' : True,
        'nwchem_tce__thresh': 1.0e-12,
    })
    print('     Testing ccsd...')
    val = qcdb.energy('nwc-ccsd')
    check_ccsd(val)
