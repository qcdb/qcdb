#TCE CCD

import os
import sys

import qcdb

from ..addons import *
from ..utils import *


def tce_ccd(return_value):

    hf  =           -74.962663062074
    ccd_tot    =    -75.012515193677402
    ccd_corl   =     -0.049852131603433
    #x       = 0.0000000
    #y       = 0.0000000
    #z       = -1.6930528

    assert compare_values(hf, qcdb.variable('HF TOTAL ENERGY'), 5, 'hf ref')
    assert compare_values(ccd_tot, qcdb.variable('CCD TOTAL ENERGY'), 5, 'ccsd tot')
    assert compare_values(ccd_corl, qcdb.variable('CCD CORRELATION ENERGY'), 5, 'ccsd corl')
    #assert compare_values(x, qcdb.variable('CURRENT DIPOLE X'), 5, 'x dipole')
    #assert compare_values(y, qcdb.variable('CURRENT DIPOLE Y'), 5, 'y dipole')
    #assert compare_values(z, qcdb.variable('CURRENT DIPOLE Z'), 5, 'z dipole')

@using_nwchem
def test_1_ccd():
    h2o = qcdb.set_molecule('''
        H    0.000000000000000   1.079252144093028   1.474611055780858
        O    0.000000000000000   0.000000000000000   0.000000000000000
        H    0.000000000000000   1.079252144093028  -1.474611055780858
        units au''')

    qcdb.set_options({
        'basis'     :   'sto-3g',
        'scf__e_convergence'   :   1e-10,
        'nwchem_scf__rhf'  :   True,
        'nwchem_scf__maxiter': 50,
        'nwchem_scf__singlet': True,
        'nwchem_scf__tol2e':   1e-10,
        'nwchem_scf__thresh':  1e-10,
        'qc_module' :   'TCE',
        'nwchem_tce__scf' :   True,
        'nwchem_tce__ccd' :   True,
        'nwchem_tce__dipole':  True,
        })

    val = qcdb.energy('nwc-ccd')
    tce_ccd(val)
