# Tensor Contraction Engine Configuration Interaction (CI) energies: SD, SDT, SDTQ
import os
import sys

import pytest

import qcdb

from ..utils import *


@using("nwchem")
def test_1_cisd():
    h2o = qcdb.set_molecule(
        """
        O      0.000000000000     0.000000000000    -0.123909374404
        H      0.000000000000     1.429936611037     0.983265845431
        H      0.000000000000    -1.429936611037     0.983265845431
        """
    )

    qcdb.set_options(
        {
            "basis": "sto-3g",
            "qc_module": "TCE",
            "nwchem_tce__cisd": True,
        }
    )
    val = qcdb.energy("nwc-cisd")

    hf = -74.506112017320
    cisd_tot = -74.746025986067849
    cisd_corl = -0.239913968748276

    assert compare_values(hf, qcdb.variable("HF TOTAL ENERGY"), 5, "hf ref")
    assert compare_values(cisd_tot, qcdb.variable("CISD TOTAL ENERGY"), 5, "cisd tot")
    assert compare_values(cisd_corl, qcdb.variable("CISD CORRELATION ENERGY"), 5, "cisd corl")


@using("nwchem")
def test_2_cisdt():
    h2o = qcdb.set_molecule(
        """
        O      0.000000000000     0.000000000000    -0.123909374404
        H      0.000000000000     1.429936611037     0.983265845431
        H      0.000000000000    -1.429936611037     0.983265845431
        """
    )

    qcdb.set_options(
        {
            "basis": "sto-3g",
            "qc_module": "TCE",
            "nwchem_tce__cisdt": True,
        }
    )
    val = qcdb.energy("nwc-cisdt")

    hf = -74.506112017320
    cisdt_tot = -74.746791001337797
    cisdt_corl = -0.240678984018215

    assert compare_values(hf, qcdb.variable("HF TOTAL ENERGY"), 5, "hf ref")
    assert compare_values(cisdt_tot, qcdb.variable("CISDT TOTAL ENERGY"), 5, "cisdt tot")
    assert compare_values(cisdt_corl, qcdb.variable("CISDT CORRELATION ENERGY"), 5, "cisdt corl")


@pytest.mark.xfail(reason="cisdtq module not compiled")
@using("nwchem")
def test_3_cisdtq():
    h2o = qcdb.set_molecule(
        """
        O      0.000000000000     0.000000000000    -0.123909374404
        H      0.000000000000     1.429936611037     0.983265845431
        H      0.000000000000    -1.429936611037     0.983265845431
        """
    )

    qcdb.set_options(
        {
            "basis": "sto-3g",
            "qc_module": "TCE",
            "nwchem_tce__cisdtq": True,
        }
    )
    val = qcdb.energy("nwc-cisdtq", local_options={"memory": 1})

    hf = -74.506112017320
    cisdtq_tot = -74.788955327897597
    cisdtq_corl = -0.282843310578009

    assert compare_values(hf, qcdb.variable("HF TOTAL ENERGY"), 5, "hf ref")
    assert compare_values(cisdtq_tot, qcdb.variable("CISDTQ TOTAL ENERGY"), 5, "cisdtq tot")
    assert compare_values(cisdtq_corl, qcdb.variable("CISDTQ CORRELATION ENERGY"), 5, "cisdtq corl")
