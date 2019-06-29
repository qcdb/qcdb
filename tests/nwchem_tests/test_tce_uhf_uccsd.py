#TCE UHF UCCSD

import os
import sys
import qcdb
from ..addons import *
from ..utils import *

def tce_uccsd(return_value):

    hf  =   -74.656470840577
    ccsd_tot    =   -74.695659393231864
    ccsd_corl   =    -0.039188552654924

    assert compare_values(hf, qcdb.get_variable('HF TOTAL ENERGY'), 5, 'hf ref')
    assert compare_values(ccsd_tot, qcdb.get_variable('CCSD TOTAL ENERGY'), 5, 'ccsd tot')
    assert compare_values(ccsd_corl, qcdb.get_variable('CCSD CORRELATION'), 5, 'ccsd corl')

@using_nwchem
def test_1_uccsd():
    h2o = qcdb.set_molecule('''
    1 2
        H
        O H 0.96
        H O 0.96 H 104.0
        ''')

    qcdb.set_options({
        'basis'     :   'sto-3g',
        'scf__e_convergence'   :   1e-10,
        'nwchem_charge'    : 1,
        'nwchem_scf__uhf'  :   True,
        'nwchem_scf__doublet': True,
        'nwchem_scf__tol2e':   1e-10,
        'qc_module' :   'TCE',
        'nwchem_tce__ccsd' :   True,
        })

    val = qcdb.energy('nwc-ccsd')
    tce_uccsd(val)


@using_nwchem
def test_mol_chg_error():
    h2o = qcdb.set_molecule('''
        H
        O H 0.96
        H O 0.96 H 104.0
        ''')

    qcdb.set_options({
        'basis'     :   'sto-3g',
        'nwchem_charge'    : 1,
        })

    with pytest.raises(qcdb.OptionReconciliationError):
        qcdb.energy('nwc-ccsd')
