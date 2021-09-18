#! pbe/sto-3g H2O DFT energy
import os
import sys

import pytest

import qcdb

from ..utils import *


@pytest.fixture
def h2o():
    h2o = qcdb.set_molecule(
        """
        O     0.00000000    0.000000000   0.11726921
        H     0.75698224    0.000000000   -0.46907685
        H     -0.75698224   0.000000000   -0.46907685
        """
    )
    return h2o


def check_pbe0(return_value):
    pbe0 = -76.338827415472
    assert compare_values(pbe0, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft pbe0")  # TEST


@using("nwchem")
def test_01_pbe0(h2o):
    qcdb.set_options(
        {
            "basis": "cc-pvdz",
            "nwchem_dft__xc": "pbe0",
        }
    )
    val = qcdb.energy("nwc-pbe0", molecule=h2o)
    check_pbe0(val)


@using("nwchem")
def test_02_b3lyp(h2o):
    qcdb.set_options({"basis": "cc-pvdz", "nwchem_dft__xc": "b3lyp"})
    print("Testing dft (b3lyp) energy...")
    val = qcdb.energy("nwc-b3lyp", molecule=h2o)
    b3lyp = -76.420359078705
    assert compare_values(b3lyp, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft b3lyp")


def check_b1b95(return_value):
    b1b95 = -76.391648839018
    assert compare_values(b1b95, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft b1b95")


@using("nwchem")
def test_03_b1b95(h2o):
    qcdb.set_options(
        {
            "basis": "cc-pvdz",
            "nwchem_dft__xc": "b1b95",
        }
    )
    print("Testing dft (b1b95) energy ...")
    val = qcdb.energy("nwc-b1b95", molecule=h2o)
    check_b1b95(val)


def check_b971(return_value):
    b971 = -76.396668942964
    assert compare_values(b971, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft b971")


@using("nwchem")
def test_04_b971(h2o):
    qcdb.set_options(
        {
            "basis": "cc-pvdz",
            "nwchem_dft__xc": "becke97-1",
        }
    )
    print("Testing dft (b971) energy ...")
    val = qcdb.energy("nwc-b97-1", molecule=h2o)
    check_b971(val)


def check_b972(return_value):
    b972 = -76.396836694369
    assert compare_values(b972, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft b972")


@using("nwchem")
def test_05_b972(h2o):
    qcdb.set_options(
        {
            "basis": "cc-pvdz",
            "nwchem_dft__xc": "becke97-2",
        }
    )
    print("Testing dft (b97-2) energy ...")
    val = qcdb.energy("nwc-b97-2", molecule=h2o)
    check_b972(val)


def check_b97gga1(return_value):
    b97gga1 = -76.406139772900
    assert compare_values(b97gga1, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft b97gga1")


@using("nwchem")
def test_06_b97gga1(h2o):
    qcdb.set_options(
        {
            "basis": "cc-pvdz",
            "nwchem_dft__xc": "becke97gga1",
        }
    )
    print("Testing dft (becke97gga1) energy ...")
    val = qcdb.energy("nwc-b97-gga1", molecule=h2o)
    check_b97gga1(val)


def check_bhandh(return_value):
    bhandh = -75.937855033020
    assert compare_values(bhandh, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft bhandh")


@using("nwchem")
def test_07_bhandh(h2o):
    qcdb.set_options({"basis": "cc-pvdz", "nwchem_dft__xc": "beckehandh"})
    print("Testing dft (bhandh) energy ...")
    val = qcdb.energy("nwc-bhandh", molecule=h2o)
    check_bhandh(val)


def check_bop(return_value):
    bop = -76.398931209278
    assert compare_values(bop, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft bop")


@using("nwchem")
def test_08_bop(h2o):
    qcdb.set_options({"basis": "cc-pvdz", "nwchem_dft__xc": "bop"})
    print("Testing dft (bop) energy...")
    val = qcdb.energy("nwc-bop", molecule=h2o)
    check_bop(val)


def check_dldf(return_value):
    dldf = -76.879179155507
    assert compare_values(dldf, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft dldf")


@using("nwchem")
def test_09_dldf(h2o):
    qcdb.set_options(
        {
            "basis": "cc-pvdz",
            "nwchem_dft__xc": "dldf",
        }
    )
    print("Testing dft (dldf) energy...")
    val = qcdb.energy("nwc-dldf", molecule=h2o)
    check_dldf(val)


def check_ft97(return_value):
    ft97 = -76.379498908163
    assert compare_values(ft97, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft ft97")


@using("nwchem")
def test_10_ft97(h2o):
    qcdb.set_options(
        {
            "basis": "cc-pvdz",
            "nwchem_dft__xc": "ft97",
        }
    )
    val = qcdb.energy("nwc-ft97", molecule=h2o)
    check_ft97(val)


def check_hcth(return_value):
    hcth = -76.405604490754
    assert compare_values(hcth, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft hcth")


@using("nwchem")
def test_11_hcth(h2o):
    qcdb.set_options({"basis": "cc-pvdz", "nwchem_dft__xc": "hcth"})
    val = qcdb.energy("nwc-hcth", molecule=h2o)
    check_hcth(val)


def check_hcth120(return_vaue):
    hcth120 = -76.411731222330
    assert compare_values(hcth120, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft hcth120")


@using("nwchem")
def test_12_hcth120(h2o):
    qcdb.set_options({"basis": "cc-pvdz", "nwchem_dft__xc": "hcth120"})


def check_hcth407p(return_value):
    hcth407p = -76.412591625926
    assert compare_values(hcth407p, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft hcth407p")


@using("nwchem")
def test_13_hcth407p(h2o):
    qcdb.set_options({"basis": "cc-pvdz", "nwchem_dft__xc": "hcth407p"})
    val = qcdb.energy("nwc-hcth407p", molecule=h2o)
    check_hcth407p(val)


def check_hcthp14(return_value):
    hcthp14 = -76.496959340155
    assert compare_values(hcthp14, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft hcthp14")


@using("nwchem")
def test_14_hcthp14(h2o):
    qcdb.set_options({"basis": "cc-pvdz", "nwchem_dft__xc": "hcthp14"})
    val = qcdb.energy("nwc-hcthp14", molecule=h2o)
    check_hcthp14(val)


def check_m05(return_value):
    m05 = -76.383962258113
    assert compare_values(m05, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft m05")


@using("nwchem")
def test_15_m05(h2o):
    qcdb.set_options(
        {
            "basis": "cc-pvdz",
            "nwchem_dft__xc": "m05",
        }
    )
    val = qcdb.energy("nwc-m05", molecule=h2o)
    check_m05(val)


def check_m05_2x(return_value):
    m05_2x = -76.408028315789
    assert compare_values(m05_2x, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft m05-2x")


@using("nwchem")
def test_16_m05_2x(h2o):
    qcdb.set_options({"basis": "cc-pvdz", "nwchem_dft__xc": "m05-2x"})
    val = qcdb.energy("nwc-m05-2x", molecule=h2o)
    check_m05_2x(val)


def check_m06(return_value):
    m06 = -76.386102783999
    assert compare_values(m06, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft m06")


@using("nwchem")
def test_17_m06(h2o):
    qcdb.set_options({"basis": "cc-pvdz", "nwchem_dft__xc": "m06"})
    val = qcdb.energy("nwc-m06", molecule=h2o)
    check_m06(val)


def check_m06_2x(return_value):
    m06_2x = -76.388644534921
    assert compare_values(m06_2x, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft m06-2x")


@using("nwchem")
def test_18_m06_2x(h2o):
    qcdb.set_options(
        {
            "basis": "cc-pvdz",
            "nwchem_dft__xc": "m06-2x",
        }
    )
    val = qcdb.energy("nwc-m06-2x", molecule=h2o)
    check_m06_2x(val)


def check_m06_hf(return_value):
    m06_hf = -76.386091998903
    assert compare_values(m06_hf, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft m06-hf")


@using("nwchem")
def test_19_m06_hf(h2o):
    qcdb.set_options({"basis": "cc-pvdz", "nwchem_dft__xc": "m06-hf"})
    val = qcdb.energy("nwc-m06-hf", molecule=h2o)
    check_m06_hf(val)


def check_m08_hx(return_value):
    m08_hx = -76.390114240676
    assert compare_values(m08_hx, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft m08-hx")


@using("nwchem")
def test_20_m08_hx(h2o):
    qcdb.set_options(
        {
            "basis": "cc-pvdz",
            "nwchem_dft__xc": "m08-hx",
        }
    )
    val = qcdb.energy("nwc-m08-hx", molecule=h2o)
    check_m08_hx(val)


def check_m08_so(return_value):
    m08_so = -76.370076256062
    assert compare_values(m08_so, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft m08-so")


@using("nwchem")
def test_21_m08_so(h2o):
    qcdb.set_options({"basis": "cc-pvdz", "nwchem_dft__xc": "m08-so"})
    val = qcdb.energy("nwc-m08-so", molecule=h2o)
    check_m08_so(val)


def check_m11(return_value):
    m11 = -76.392301723962
    assert compare_values(m11, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft m11")


@using("nwchem")
def test_22_m11(h2o):
    qcdb.set_options({"basis": "cc-pvdz", "nwchem_dft__xc": "m11"})
    val = qcdb.energy("nwc-m11", molecule=h2o)
    check_m11(val)


def check_m11_l(return_value):
    m11_l = -76.405634428063
    assert compare_values(m11_l, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft m11-l")


@using("nwchem")
def test_23_m11_l(h2o):
    qcdb.set_options(
        {
            "basis": "cc-pvdz",
            "nwchem_dft__xc": "m11-l",
        }
    )
    val = qcdb.energy("nwc-m11-l", molecule=h2o)
    check_m11_l(val)


def check_mpw1b95(return_value):
    mpw1b95 = -76.388487462370
    assert compare_values(mpw1b95, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft mpw1b95")


@using("nwchem")
def test_24_mpw1b95(h2o):
    qcdb.set_options({"basis": "cc-pvdz", "nwchem_dft__xc": "mpw1b95"})
    val = qcdb.energy("nwc-mpw1b95", molecule=h2o)
    check_mpw1b95(val)


def check_mpw1k(return_value):
    mpw1k = -76.396497325878
    assert compare_values(mpw1k, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft mpw1k")


@using("nwchem")
def test_25_mpw1k(h2o):
    qcdb.set_options({"basis": "cc-pvdz", "nwchem_dft__xc": "mpw1k"})
    val = qcdb.energy("nwc-mpw1k", molecule=h2o)
    check_mpw1k(val)


def check_pw6b95(return_value):
    pw6b95 = -76.500335097526
    assert compare_values(pw6b95, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft pw6b95")


@using("nwchem")
def test_26_pw6b95(h2o):
    qcdb.set_options(
        {
            "basis": "cc-pvdz",
            "nwchem_dft__xc": "pw6b95",
        }
    )
    val = qcdb.energy("nwc-pw6b95", molecule=h2o)
    check_pw6b95(val)


def check_pwb6k(return_value):
    pwb6k = -76.448685285204
    assert compare_values(pwb6k, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft pwb6k")


@using("nwchem")
def test_27_pwb6k(h2o):
    qcdb.set_options({"basis": "cc-pvdz", "nwchem_dft__xc": "pwb6k"})
    val = qcdb.energy("nwc-pwb6k", molecule=h2o)
    check_pwb6k(val)


def check_tpssh(return_value):
    tpssh = -76.416480345716
    assert compare_values(tpssh, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft tpss")


@using("nwchem")
def test_28_tpssh(h2o):
    qcdb.set_options({"basis": "cc-pvdz", "nwchem_dft__xc": "xctpssh"})
    val = qcdb.energy("nwc-tpssh", molecule=h2o)
    check_tpssh(val)


def check_bhlyp(return_value):
    bhlyp = -76.381422122691
    assert compare_values(bhlyp, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft bhlyp")


@using("nwchem")
def test_29_bhlyp(h2o):
    qcdb.set_options({"basis": "cc-pvdz", "nwchem_dft__xc": "bhlyp"})
    val = qcdb.energy("nwc-bhlyp", molecule=h2o)
    check_bhlyp(val)


def check_hcth407(return_value):
    hcth407 = -76.411625217878
    assert compare_values(hcth407, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft hcth407")


@using("nwchem")
def test_30_hcth407(h2o):
    qcdb.set_options({"basis": "cc-pvdz", "nwchem_dft__xc": "hcth407"})
    val = qcdb.energy("nwc-hcth407", molecule=h2o)
    check_hcth407(val)


def check_pbeop(return_value):
    pbeop = -76.344549764234
    assert compare_values(pbeop, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft pbeop")


@using("nwchem")
def test_31_pbeop(h2o):
    qcdb.set_options({"basis": "cc-pvdz", "nwchem_dft__xc": "pbeop"})
    val = qcdb.energy("nwc-pbeop", molecule=h2o)
    check_pbeop(val)


def check_b97d(return_value):
    b97d = -76.380347936365
    assert compare_values(b97d, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft b97d")


@using("nwchem")
def test_32_b97d(h2o):
    qcdb.set_options({"basis": "cc-pvdz", "nwchem_dft__xc": "becke97-d"})
    val = qcdb.energy("nwc-b97-d", molecule=h2o)
    check_b97d(val)


def check_cft97(return_value):
    cft97 = -67.648400200510
    assert compare_values(cft97, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft cft97")


@using("nwchem")
def test_33_cft97(h2o):
    qcdb.set_options({"basis": "cc-pvdz", "nwchem_dft__xc": "cft97"})
    val = qcdb.energy("nwc-cft97", molecule=h2o)
    check_cft97(val)


def check_acm(return_value):
    acm = -76.397101204907
    assert compare_values(acm, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft acm")


@using("nwchem")
def test_34_acm(h2o):
    qcdb.set_options({"basis": "cc-pvdz", "nwchem_dft__xc": "acm"})
    val = qcdb.energy("nwc-acm", molecule=h2o)
    check_acm(val)


def check_optx(return_value):
    optx = -76.055691019559
    assert compare_values(optx, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft optx")


@using("nwchem")
def test_35_optx(h2o):
    qcdb.set_options({"basis": "cc-pvdz", "nwchem_dft__xc": "optx"})
    val = qcdb.energy("nwc-optx", molecule=h2o)
    check_optx(val)


def check_b98(return_value):
    b98 = -76.394526199939
    assert compare_values(b98, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft b98")


@using("nwchem")
def test_36_b98(h2o):
    qcdb.set_options({"basis": "cc-pvdz", "nwchem_dft__xc": "becke98"})
    val = qcdb.energy("nwc-b98", molecule=h2o)
    check_b98(val)


def check_xtpss03(return_value):
    xtpss03 = -76.090088650605
    assert compare_values(xtpss03, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft xtpss03")


@using("nwchem")
def test_37_xtpss03(h2o):
    qcdb.set_options({"basis": "cc-pvdz", "nwchem_dft__xc": "xtpss03"})
    val = qcdb.energy("nwc-xtpss03", molecule=h2o)
    check_xtpss03(val)


def check_bb1k(return_value):
    bb1k = -76.387274033902
    assert compare_values(bb1k, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft bb1k")


@using("nwchem")
def test_38_bb1k(h2o):
    qcdb.set_options({"basis": "cc-pvdz", "nwchem_dft__xc": "bb1k"})
    val = qcdb.energy("nwc-bb1k", molecule=h2o)
    check_bb1k(val)


def check_vs98(return_value):
    vs98 = -76.441813521942
    assert compare_values(vs98, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft vs98")


@using("nwchem")
def test_39_vs98(h2o):
    qcdb.set_options({"basis": "cc-pvdz", "nwchem_dft__xc": "vs98"})
    val = qcdb.energy("nwc-vs98", molecule=h2o)
    check_vs98(val)


def check_m06_l(return_value):
    m06_l = -76.413915672071
    assert compare_values(m06_l, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft m06-l")


@using("nwchem")
def test_40_m06_l(h2o):
    qcdb.set_options({"basis": "cc-pvdz", "nwchem_dft__xc": "m06-l"})
    val = qcdb.energy("nwc-m06-l", molecule=h2o)
    check_m06_l(val)


def check_hcth147(return_value):
    hcth147 = -76.406989374352
    assert compare_values(hcth147, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft hcth147")


@using("nwchem")
def test_41_hcth147(h2o):
    qcdb.set_options({"basis": "cc-pvdz", "nwchem_dft__xc": "hcth147"})
    val = qcdb.energy("nwc-hcth147", molecule=h2o)
    check_hcth147(val)


def check_pw91(return_value):
    pw91 = -76.406989374352
    assert compare_values(pw91, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft hcth147")
