import os
import sys

import numpy as np

import qcdb

from ..utils import *


def check_eom_ccsd_grad(return_value):
    ref = -74.962905406245  # scf #sym a
    ccsd_corl = -0.049337013315753  # sym a
    ccsd_tot = -75.012242419560370  # sym a
    # may need to figure out how to pull out harvested items above
    eom_ccsd_grad = np.array(
        [[0.000000, -0.000091, 0.375029], [0.000000, -0.309350, -0.187479], [0.000000, 0.309441, -0.187550]]
    )

    assert compare_values(eom_ccsd_grad, qcdb.variable("CURRENT GRADIENT"), 5, "eom ccsd grad")


@using("nwchem")
def test_1_eom_ccsd_grad():
    h2o = qcdb.set_molecule(
        """
            O      0.000000000000     0.000000000000    -0.123909374404
            H      0.000000000000     1.429936611037     0.983265845431
            H      0.000000000000    -1.429936611037     0.983265845431
            symmetry c1
            units au
            """
    )

    qcdb.set_options(
        {
            "basis": "sto-3g",
            "memory": "1500 mb",
            "nwchem_scf__rhf": True,
            "scf__e_convergence": 1.0e-10,
            "nwchem_scf__tol2e": 1.0e-10,
            "nwchem_scf__thresh": 1.0e-10,
            "nwchem_scf__singlet": True,
            "qc_module": "TCE",
            "nwchem_tce__ccsd": True,
            "nwchem_tce__nroots": 1,
        }
    )

    print("Testing EOM-CCSD grad...")
    val = qcdb.gradient("nwc-eom-ccsd")
    check_eom_ccsd_grad(val)
