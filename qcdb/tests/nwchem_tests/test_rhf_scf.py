# Geometry optimization HF/6-31g* on water
import os
import sys

import qcdb

from ..utils import *


def check_rhf(return_value):
    ref = -76.010496306999
    nre = 9.187334240165
    assert compare_values(ref, qcdb.variable("HF TOTAL ENERGY"), 5, "hf ref")
    assert compare_values(nre, qcdb.variable("NUCLEAR REPULSION ENERGY"), 5, "nre")


@using("nwchem")
def test_1_hf():
    h2o = qcdb.set_molecule(
        """
        O     0.000000000000    0.000000000000   -0.065638538099
        H     0.000000000000   -0.757480611647    0.520865616174
        H     0.000000000000    0.757480611647    0.520865616174
        """
    )

    qcdb.set_options(
        {
            "basis": "6-31g*",
            "scf__e_convergence": 1.0e-8,
            #'nwchem_geometry_center'    : False,
            #'nwchem_geometry_autosym'   : False,
            "nwchem_scf__rhf": True,
            #'nwchem_scf__thresh': 1.0e-8,
            "nwchem_scf__direct": True,
            #'nwchem_task_scf': 'optimize'
        }
    )
    print("Testing HF...")
    val = qcdb.energy("nwc-hf")
    check_rhf(val)
