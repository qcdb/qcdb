import os
import sys
import pprint as pp

import numpy as py  # will keep for when QCDB transfer dipole from single var --> array
import pytest

import qcdb

from ..addons import *
from ..utils import *


def check_tce_ccsd(return_value):
    hf      =   -74.962663062078
    ccsd_tot=   -75.012790383079462
    ccsd_corl=   -0.050127321001234
    ccsd_dipole_X   =   0.0000000
    ccsd_dipole_Y   =   0.0000000
    ccsd_dipole_Z   =  -1.5819795
    ccsd_dipole = [ccsd_dipole_X, ccsd_dipole_Y, ccsd_dipole_Z]

    assert compare_values(hf, qcdb.variable('HF TOTAL ENERGY'), 5, 'hf ref')
    assert compare_values(ccsd_tot, qcdb.variable('CCSD TOTAL ENERGY'), 5, 'ccsd tot')
    assert compare_values(ccsd_corl, qcdb.variable('CCSD CORRELATION ENERGY'), 5, 'ccsd corl')
    assert compare_values(ccsd_dipole, qcdb.variable("CCSD DIPOLE"), 5, "ccsd dipole")

@using_nwchem
def test_1_ccsd_dipole():
    h2o = qcdb.set_molecule(
            '''H 0.0   0.5711156805885   0.7803306218431
               O 0.0   0.0   0.0
               H 0.0   0.5711156805885  -0.7803306218431'''
               )

    qcdb.set_options({
        'basis' :   'sto-3g',
        'qc_module' :   'TCE',
        'nwchem_scf__tol2e' :   1e-10,
        'nwchem_scf__thresh':   1e-14,
        'nwchem_scf__singlet':  True,
        'nwchem_scf__rhf'    :   True,
        'nwchem_tce__dipole':   True
        })

    val = qcdb.energy('nwc-ccsd')
    check_tce_ccsd(val)
