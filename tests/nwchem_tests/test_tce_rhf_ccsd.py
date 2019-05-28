#! single-point CCSD/cc-pvdz on water
import os
import sys
from ..utils import *
from ..addons import *
import qcdb

h2o = qcdb.set_molecule('''
        O     0.000000000000    0.000000000000   -0.065638538099
        H     0.000000000000   -0.757480611647    0.520865616174
        H     0.000000000000    0.757480611647    0.520865616174
        ''')

print(h2o)


def check_ccsd_t_(return_value):
    ref = -76.026760733967
    nre = 9.187333574703
    ccsdcorl = -0.213341251859034
    ccsdtot = -76.240101985825859

    assert compare_values(ref, qcdb.get_variable('HF TOTAL ENERGY'), 5, 'hf ref')
    assert compare_values(ccsdcorl, qcdb.get_variable('CCSD CORRELATION ENERGY'), 5, 'ccsd corl')
    assert compare_values(ccsdtot, qcdb.get_variable('CCSD TOTAL ENERGY'), 5, 'ccsd total')
    assert compare_values(nre, qcdb.get_variable('NUCLEAR REPULSION ENERGY'), 5, 'nre')


@using_nwchem
def test_1_ccsd_t_():
    qcdb.set_options({
        'basis': 'cc-pvdz',
        'memory': '600 mb',
        'nwchem_scf__rhf': True,
        'nwchem_scf__thresh': 1.0e-12,
        'nwchem_tce__dft': False,
        'nwchem_tce__ccsd': True,
        'qc_module'  : 'tce',
        'nwchem_tce__thresh': 1.0e-12,
    })
    print('     Testing ccsd...')
    val = qcdb.energy('nwc-ccsd')
    check_ccsd_t_(val)
