#TCE CCD

import os
import sys
import qcdb
from ..addons import *
from ..utils import *


def tce_ccd(return_value):

    hf  =          -74.506112017705 
    ccd_tot    =   -74.789850325141188
    ccd_corl   =    -0.283738307435901

    assert compare_values(hf, qcdb.get_variable('HF TOTAL ENERGY'), 5, 'hf ref')
    assert compare_values(ccd_tot, qcdb.get_variable('CCD TOTAL ENERGY'), 5, 'ccsd tot')
    assert compare_values(ccd_corl, qcdb.get_variable('CCD CORRELATION ENERGY'), 5, 'ccsd corl')

@using_nwchem
def test_1_ccd():
    h2o = qcdb.set_molecule('''
        O                     0.000000000000     0.000000000000    -0.234154782060
        H                    -0.000000000000     2.702188571625     1.858103156322
        H                     0.000000000000    -2.702188571625     1.858103156322
        ''')

    qcdb.set_options({
        'basis'     :   'sto-3g',
        'scf__e_convergence'   :   1e-10,
        'nwchem_scf__rhf'  :   True,
        'nwchem_scf__singlet': True,
        'nwchem_scf__tol2e':   1e-10,
        'nwchem_scf__thresh':  1e-10,
        'qc_module' :   'TCE',
        'nwchem_tce__ccd' :   True,
        'nwchem_tce__dipole':  True,
        })

    val = qcdb.energy('nwc-ccd')
    tce_ccd(val)
