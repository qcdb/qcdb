import os
import sys

import numpy as np

import qcdb

from ..utils import *


@using("nwchem")
def test_grad():
    h2o = qcdb.set_molecule(
        """
        O      0.00000000    0.00000000    0.00000000
        H      0.00000000    1.93042809   -1.10715266
        H      0.00000000   -1.93042809   -1.10715266
        units au"""
    )

    qcdb.set_options(
        {
            "basis": "sto-3g",
            "scf__e_convergence": 1e-6,
            #'nwchem_driver__tight': True
        }
    )
    val = qcdb.gradient("nwc-scf")

    scf = -74.888142460799
    grads = np.array(
        [[0.000000, 0.000000, 0.058550], [0.000000, 0.140065, -0.029275], [0.000000, -0.140065, -0.029275]]
    )

    assert compare_values(scf, qcdb.variable("HF TOTAL ENERGY"), 5, "scf")
    assert compare_arrays(grads, qcdb.variable("CURRENT GRADIENT"), 5, "scf grad")
