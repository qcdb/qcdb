# single-point DFT-CCSDT/sto-3g on water

import os
import sys

import qcdb

from ..utils import *


def check_dft(return_value, is_dft):
    if is_dft:
        ref = -75.312572965120
        ccsdt_tot = -75.362578004794031
        ccsdt_corl = -0.050005039673792
        assert compare_values(ref, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft ref")
    else:  # change values for non-dft
        ref = -74.963048525888
        ccsdt_tot = -75.012605624645104
        ccsdt_corl = -0.049557089026546
        nre = 9.187333574703
        assert compare_values(ref, qcdb.variable("HF TOTAL ENERGY"), 5, "dft ref")
    assert compare_values(ccsdt_tot, qcdb.variable("CCSDT TOTAL ENERGY"), 5, "ccsdt total")
    assert compare_values(ccsdt_corl, qcdb.variable("CCSDT CORRELATION ENERGY"), 5, "ccsdt corl")


@using("nwchem")
def test_1_dft():
    h2o = qcdb.set_molecule(
        """
        O     0.000000000000    0.000000000000   -0.065638538099
        H     0.000000000000   -0.757480611647    0.520865616174
        H     0.000000000000    0.757480611647    0.520865616174
        """
    )

    qcdb.set_options(
        {
            "basis": "sto-3g",
            #'nwchem_dft__convergence__density': 1.0e-12,
            "nwchem_dft__xc": "b3lyp",
            # add grid options
            "nwchem_tce__dft": True,
            "qc_module": "tce",
            "nwchem_tce__ccsdt": True,
            "nwchem_tce__thresh": 1.0e-12,
        }
    )
    print("Testing CCSDT-DFT energy...")
    val = qcdb.energy("nwc-ccsdt", local_options={"memory": 6})
    check_dft(val, is_dft=True)


@using("nwchem")
def test_2_scf():
    h2o = qcdb.set_molecule(
        """
        O     0.000000000000    0.000000000000   -0.065638538099
        H     0.000000000000   -0.757480611647    0.520865616174
        H     0.000000000000    0.757480611647    0.520865616174
        """
    )

    qcdb.set_options(
        {
            "basis": "sto-3g",
            "nwchem_tce__scf": True,
            "qc_module": "tce",
            "nwchem_tce__ccsdt": True,
            "nwchem_tce__thresh": 1.0e-12,
        }
    )
    print("Test CCSDT-SCF energy ...")
    val = qcdb.energy("nwc-ccsdt")
    val = qcdb.energy("nwc-ccsdt", local_options={"memory": 6})
    check_dft(val, is_dft=False)
