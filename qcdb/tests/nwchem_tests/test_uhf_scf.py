#! single-point UHF/cc-pVDZ  on NH2
import os

from ..utils import *
from ..addons import *

import qcdb


def check_uhf_hf(return_value):
    ref = -55.566057523877
    nre = 7.580905897627

    assert compare_values(ref, qcdb.variable('HF TOTAL ENERGY'), 5, 'scf')
    assert compare_values(nre, qcdb.variable('NUCLEAR REPULSION ENERGY'), 5, 'nre')


@using_nwchem
def test_1_hf():
    nh2 = qcdb.set_molecule('''
           N        0.08546       -0.00020       -0.05091
           H       -0.25454       -0.62639        0.67895
           H       -0.25454       -0.31918       -0.95813
           ''')

    qcdb.set_options({
        'basis': 'cc-pvdz',
        'memory': '400 mb',
        'nwchem_scf__uhf': True,
        'nwchem_scf__nopen': 1,
        'scf__e_convergence': 1.0e-8
    })
    print('Testing hf...')
    val = qcdb.energy('nwc-hf')
    check_uhf_hf(val)
