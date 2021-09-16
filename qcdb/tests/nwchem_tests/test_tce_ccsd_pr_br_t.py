# TCE CCSD(T) and CCSD[T] calculations
import os
import sys

import qcdb

from ..utils import *


def check_ccsd_t_pr_br(return_value):
    ccsd_tot = -76.240077811301250
    ccsd_corl = -0.213269954065481
    t_br_corr = -0.003139909173705
    t_br_corl = -0.216409863239186
    ccsd_t_br = -76.243217720474960
    t_pr_corr = -0.003054718622142
    t_pr_corl = -0.216324672687623
    ccsd_t_pr = -76.243132529923390

    assert compare_values(ccsd_tot, qcdb.variable("CCSD TOTAL ENERGY"), 5, "ccsd total")
    assert compare_values(ccsd_corl, qcdb.variable("CCSD CORRELATION ENERGY"), 5, "ccsd corl")
    assert compare_values(t_br_corr, qcdb.variable("T(CCSD) CORRECTION ENERGY"), 5, "[t] corr")
    assert compare_values(t_br_corl, qcdb.variable("CCSD+T(CCSD) CORRELATION ENERGY"), 5, "ccsd[t] corl")
    assert compare_values(ccsd_t_br, qcdb.variable("CCSD+T(CCSD) TOTAL ENERGY"), 5, "ccsd[t] total")
    assert compare_values(t_pr_corr, qcdb.variable("(T) CORRECTION ENERGY"), 5, "(t) corr")
    assert compare_values(t_pr_corl, qcdb.variable("CCSD(T) CORRELATION ENERGY"), 5, "ccsd(t) corl")
    assert compare_values(ccsd_t_pr, qcdb.variable("CCSD(T) TOTAL ENERGY"), 5, "ccsd(t) tot")


@using("nwchem")
def test_1_ccsd_t():
    h2o = qcdb.set_molecule(
        """
        O     0.00000000     0.00000000     0.22138519
        H     0.00000000    -1.43013023    -0.88554075
        H     0.00000000     1.43013023    -0.88554075
        units au"""
    )

    qcdb.set_options(
        {
            "basis": "cc-pvdz",
            "nwchem_scf__rhf": True,
            "nwchem_scf__thresh": 1.0e-10,
            "nwchem_scf__tol2e": 1.0e-10,
            "nwchem_scf__singlet": True,
            "nwchem_tce__ccsd(t)": True,
            "qc_module": "TCE",
            "nwchem_tce__io": "ga",
        }
    )

    val = qcdb.energy("nwc-ccsd(t)")
    check_ccsd_t_pr_br(val)
