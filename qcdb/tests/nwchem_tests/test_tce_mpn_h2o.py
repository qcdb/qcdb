# Tensor Contraction Engine Moller-Plesset perturbation to second order (MP2)
import os
import sys

import qcdb

from ..utils import *


def check_tce_mp2(return_value):
    hf = -75.645085110552
    mp2_tot = -75.934995820774219
    mp2_corl = -0.289910710222209

    assert compare_values(hf, qcdb.variable("HF TOTAL ENERGY"), 5, "hf ref")
    assert compare_values(mp2_tot, qcdb.variable("MP2 TOTAL ENERGY"), 5, "mbpt(2) tot")
    assert compare_values(mp2_corl, qcdb.variable("MP2 CORRELATION ENERGY"), 5, "mbpt(2) corl")


@using("nwchem")
def test_1_mp2():
    h2o = qcdb.set_molecule(
        """
        O      0.000000000000     0.000000000000    -0.123909374404
        H      0.000000000000     1.429936611037     0.983265845431
        H      0.000000000000    -1.429936611037     0.983265845431
        """
    )

    qcdb.set_options(
        {
            "basis": "cc-pvdz",
            "qc_module": "TCE",
            "nwchem_tce__mp2": True,
        }
    )
    val = qcdb.energy("nwc-mp2")
    check_tce_mp2(val)


def check_tce_mp3(return_value):
    hf = -75.645085110552
    mp2_tot = -75.934995820774219
    mp2_corl = -0.289910710222209
    mp3_tot = -75.924895885924158
    mp3_corr = 0.010099934850061

    assert compare_values(hf, qcdb.variable("HF TOTAL ENERGY"), 5, "hf ref")
    assert compare_values(mp2_tot, qcdb.variable("MP2 TOTAL ENERGY"), 5, "mbpt(2) tot")
    assert compare_values(mp2_corl, qcdb.variable("MP2 CORRELATION ENERGY"), 5, "mbpt(2) corl")
    assert compare_values(mp3_tot, qcdb.variable("MP3 TOTAL ENERGY"), 5, "mp3 tot")
    assert compare_values(mp3_corr, qcdb.variable("MP3 CORRECTION ENERGY"), 5, "mp3 corr")


# Check all corls to corr similar to qcng
@using("nwchem")
def test_2_mp3():
    h2o = qcdb.set_molecule(
        """
        O      0.000000000000     0.000000000000    -0.123909374404
        H      0.000000000000     1.429936611037     0.983265845431
        H      0.000000000000    -1.429936611037     0.983265845431
        """
    )

    qcdb.set_options(
        {
            "basis": "cc-pvdz",
            "qc_module": "TCE",
            "nwchem_tce__mp3": True,
        }
    )
    val = qcdb.energy("nwc-mp3", local_options={"memory": 1})
    check_tce_mp3(val)


def check_tce_mp4(return_value):
    hf = -75.645085110552
    mp2_tot = -75.934995820774219
    mp2_corl = -0.289910710222209
    mp3_tot = -75.924895885924158
    mp3_corr = 0.010099934850061
    mp4_tot = -75.964976635484248
    mp4_corr = -0.040080749560086

    assert compare_values(hf, qcdb.variable("HF TOTAL ENERGY"), 5, "hf ref")
    assert compare_values(mp2_tot, qcdb.variable("MP2 TOTAL ENERGY"), 5, "mbpt(2) tot")
    assert compare_values(mp2_corl, qcdb.variable("MP2 CORRELATION ENERGY"), 5, "mbpt(2) corl")
    assert compare_values(mp3_tot, qcdb.variable("MP3 TOTAL ENERGY"), 5, "mp3 tot")
    assert compare_values(mp3_corr, qcdb.variable("MP3 CORRECTION ENERGY"), 5, "mp3 corr")
    assert compare_values(mp4_tot, qcdb.variable("MP4 TOTAL ENERGY"), 5, "mp4 tot")
    assert compare_values(mp4_corr, qcdb.variable("MP4 CORRECTION ENERGY"), 5, "mp4 corr")


@using("nwchem")
def test_3_mp4():
    h2o = qcdb.set_molecule(
        """
        O      0.000000000000     0.000000000000    -0.123909374404
        H      0.000000000000     1.429936611037     0.983265845431
        H      0.000000000000    -1.429936611037     0.983265845431
        """
    )

    qcdb.set_options(
        {
            "basis": "cc-pvdz",
            "qc_module": "TCE",
            "nwchem_tce__mp4": True,
        }
    )
    val = qcdb.energy("nwc-mp4", local_options={"memory": 1})
    check_tce_mp4(val)
