#! single-point HF/cc-pVDZ (Cartesian) on NH2
import os
import sys
from ..utils import *
from ..addons import *
import qcdb

nh2 = qcdb.set_molecule('''
        N    0.00000000     0.00000000    -0.13636256
        H   -0.79970017     0.00000000     0.47726896
        H    0.79970017     0.00000000     0.47726896
        ''')
print(nh2)


def check_hf(return_value):
    ref = -55.562949656313
    nre = 7.680543797856

    assert compare_values(ref, qcdb.get_variable('HF TOTAL ENERGY'), 5, 'hf total')
    assert compare_values(nre, qcdb.get_variable('NUCLEAR REPULSION ENERGY'), 5, 'nre')


@using_nwchem
def test_1_rohf():
    qcdb.set_options({
        'basis': 'cc-pVDZ',
        'memory': '300 mb',
        'nwchem_scf': 'ROHF',
        'nwchem_scf_nopen': 1,
        'nwchem_scf_thresh': 1.0e-8
    })
    print('Testing HF energy ...')
    val = qcdb.energy('nwc-hf')
    check_hf(val)
