import numpy as np
import pytest
import qcengine as qcng
from qcengine.programs.tests.standard_suite_ref import std_molecules, std_refs
from qcengine.testing import using

import qcdb

from .standard_suite_runner import runner_asserter

# pytestmark = [pytest.mark.quick, pytest.mark.mp2]


_basis_keywords = ["cfour_basis"]
_trans_qcprog = {"c4": "cfour", "gms": "gamess", "nwc": "nwchem", "p4": "psi4"}


@pytest.fixture
def clsd_open_pmols():
    return {name: qcdb.Molecule.from_string(smol, name=name) for name, smol in std_molecules.items()}


# yapf: disable
_q1 = (qcng.exceptions.InputError, "unknown SCFTYPE")  # no rohf reference for nwc mp2
_q2 = (qcng.exceptions.InputError, "CCTYP IS PROGRAMMED ONLY FOR SCFTYP=RHF OR ROHF")
_q3 = (qcng.exceptions.InputError, "ccsd: nopen is not zero")
_q4 = (qcng.exceptions.InputError, "ROHF'S CCTYP MUST BE CCSD OR CR-CCL")  # No (T) w/ rohf in gms
_q5 = (qcdb.exceptions.ValidationError, "Derivative method 'name' (gms-ccsd) and derivative level 'dertype' (1) are not available.")  # no gms cc grad
_q6 = (qcng.exceptions.InputError, "Only RHF/UHF Hessians are currently implemented.")  # no rohf hess for psi4 hf
_q7 = (qcng.exceptions.UnknownError, "CALCLEVEL LCCD not implemented for ROHF references.")  # gms
_q8 = (qcng.exceptions.UnknownError, "select_lccd: Method 'lccd' with CC_TYPE 'CONV' and REFERENCE 'ROHF' not available")
_q9 = (KeyError, "gms-lccsd")  # no lccsd in gms
_q10 = (qcng.exceptions.InputError, "lccsd requires 'reference rhf'.")  # psi
_q11 = (qcng.exceptions.InputError, "Invalid type 'CONV' for DFOCC")  # psi only df/cd ccd
_q12 = (qcng.exceptions.UnknownError, "CCTYP HAS NO ANALYTIC NUCLEAR GRADIENT PROGRAMMED.")  # gms
_q13 = (qcng.exceptions.UnknownError, "UHF DF-CCSD has NOT been implemented yet!")  # psi only df/cd ccd
_q14 = (KeyError, "gms-mp3")  # no mp3 in gms
_q15 = (qcng.exceptions.UnknownError, "select_mp3: Method 'mp3' with MP_TYPE 'CONV' and REFERENCE 'ROHF' not available")  # only detci for conv rohf mp3 in psi4, and it's already peculiar for mp2
_q16 = (qcdb.exceptions.ValidationError, "Derivative method 'name' (gms-mp3) and derivative level 'dertype' (1) are not available.")  # no gms mp3 grad
_q17 = (qcdb.exceptions.ValidationError, "Derivative method 'name' (p4-ccd) and derivative level")  # no psi4 ccd
_q18 = (qcng.exceptions.InputError, "Frozen core is not available for the CC gradients.")  # psi ccenergy
_q19 = (qcng.exceptions.UnknownError, "UHF/ROHF gradients not implemented for this level of theory")  # cfour ccsdt
_q20 = (KeyError, "p4-ccsd+t(ccsd)")
_q21 = (KeyError, "gms-ccsdt")
_q22 = (KeyError, "p4-ccsdt")
_q23 = (KeyError, "ccsdt-3")  # gamess, nwchem, psi4
_q24 = (KeyError, "ccsdt-1a")  # gamess, nwchem, psi4
_q25 = (KeyError, "ccsdt-1b")  # gamess, nwchem, psi4
_q26 = (KeyError, "ccsdt-2")  # gamess, nwchem, psi4


_w1 = ("MP2 CORRELATION ENERGY", "nonstandard answer: NWChem TCE MP2 doesn't report singles (affects ROHF)")
_w2 = ("CCSD CORRELATION ENERGY", "nonstandard answer: GAMESS CCSD ROHF FC energy")
_w3 = ("(T) CORRECTION ENERGY", "nonstandard answer: NWChem CCSD(T) ROHF AE energy")
_w4 = ("MP2 CORRELATION ENERGY", "nonstandard answer: NWChem TCE MP3 doesn't report singles (affects ROHF), may be off by MP2 singles value")
_w5 = ("MP2 CORRELATION ENERGY", "nonstandard answer: GAMESS MP2 ROHF gradient ZAPT energies")
_w6 = ("CCSDTQ CORRELATION ENERGY", "misdirected calc: CFOUR NCC CCSDTQ gradient mixed fc/ae parts")
_w7 = ("CCSDT CORRELATION ENERGY", "misdirected calc: CFOUR NCC CCSDT & CCSDT(Q) gradient mixed fc/ae parts")
_w8 = ("CCSD CORRELATION ENERGY", "misdirected calc: CFOUR NCC CCSD gradient mixed fc/ae parts")
_w9 = ("CCD CORRELATION ENERGY", "misdirected calc: CFOUR NCC CCD gradient mixed fc/ae parts")
_w10 = ("MP2.5 CORRELATION ENERGY", "misdirected calc: CFOUR NCC MP3 gradient mixed fc/ae parts")
_w11 = ("MP3 TOTAL GRADIENT", "nonstandard answer: CFOUR MP3 RHF FC doesn't match findif")
_w12 = ("MP2 CORRELATION ENERGY", "nonstandard answer: CFOUR CCSD ROHF FC gradient right but energies wrong")
_w13 = ("CCSDT-3 CORRELATION ENERGY", "misdirected calc: CFOUR NCC CCSDT-3 gradient mixed fc/ae parts")
_w14 = ("CCSDT-3 TOTAL GRADIENT", "nonstandard answer: CFOUR CCSDT-3 FC doesn't match findif")
_w15 = ("CCSDT-3 TOTAL HESSIAN", "nonstandard answer: CFOUR CCSDT-3 FC doesn't match findif")
_w16 = ("CCSDT-1A TOTAL GRADIENT", "nonstandard answer: CFOUR CCSDT-1A FC doesn't match findif")
_w17 = ("CCSDT-1A TOTAL HESSIAN", "nonstandard answer: CFOUR CCSDT-1A FC doesn't match findif")
_w18 = ("CCSDT-1B TOTAL GRADIENT", "nonstandard answer: CFOUR CCSDT-1B FC doesn't match findif")
_w19 = ("CCSDT-1B TOTAL HESSIAN", "nonstandard answer: CFOUR CCSDT-1B FC doesn't match findif")
_w20 = ("CCSDT-2 TOTAL GRADIENT", "nonstandard answer: CFOUR CCSDT-2 FC doesn't match findif")
_w21 = ("CCSDT-2 TOTAL HESSIAN", "nonstandard answer: CFOUR CCSDT-2 FC doesn't match findif")
_w22 = ("CCSDT-1A CORRELATION ENERGY", "misdirected calc: CFOUR NCC CCSDT-1a gradient mixed fc/ae parts")
_w23 = ("CCSDT-1B CORRELATION ENERGY", "misdirected calc: CFOUR NCC CCSDT-1b gradient mixed fc/ae parts")
_w24 = ("CCSDT-2 CORRELATION ENERGY", "misdirected calc: CFOUR NCC CCSDT-2 gradient mixed fc/ae parts")
# yapf: enable

# <<<  Notes
#
# * `"xptd": {"fd": False}` is just a note that qcprog is doing internal findif, but the results are accurate
# * `"xptd": {"fd": True}` indicates internal findif, but probably 3-point, and results need a looser atol
#
# >>> Notes


def _trans_key(qc, bas, key):
    lkey = key.lower()

    # translate basis
    if lkey == "basis":
        return bas
    if bas.lower() == "cc-pvdz":
        if lkey == "cfour_basis":
            return "pvdz"
    elif bas.lower() == "aug-cc-pvdz":
        if lkey == "cfour_basis":
            return "aug-pvdz"
    elif bas.lower() == "cfour-qz2p":
        if lkey == "cfour_basis":
            return "qz2p"

    sys.exit(1)


# http://patorjk.com/software/taag/#p=display&c=bash&f=Soft&t=MP3


#  ,--.  ,--.,------.    ,------.
#  |  '--'  ||  .---'    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  .--.  ||  `--,     |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  |  |  ||  |`       |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `--'  `--'`--'        `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                      `---' `---'
#  <<<  HF Energy


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are all possible ways of computing <method> working?

        pytest.param({"call": "c4-hf",  "reference": "rhf",  "fcae": "ae", "keywords": {"basis": "<>", "cfour_scf_conv": 12},                                                                           }, id="hf  rhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-hf", "reference": "rhf",  "fcae": "ae", "keywords": {"basis": "<>"},                                                                                                 }, id="hf  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-hf", "reference": "rhf",  "fcae": "ae", "keywords": {"basis": "<>"},                                                                                                 }, id="hf  rhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-hf",  "reference": "rhf",  "fcae": "ae", "keywords": {"basis": "<>", "psi4_scf_type": "pk"},                                                                          }, id="hf  rhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-hf",  "reference": "uhf",  "fcae": "ae", "keywords": {"basis": "<>", "cfour_reference": "uhf", "cfour_scf_conv": 12},                                                 }, id="hf  uhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-hf", "reference": "uhf",  "fcae": "ae", "keywords": {"basis": "<>", "gamess_contrl__scftyp": "uhf"},                                                                 }, id="hf  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-hf", "reference": "uhf",  "fcae": "ae", "keywords": {"basis": "<>", "nwchem_scf__uhf": True},                                                                        }, id="hf  uhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-hf",  "reference": "uhf",  "fcae": "ae", "keywords": {"basis": "<>", "reference": "uhf", "psi4_scf_type": "pk"},                                                      }, id="hf  uhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-hf",  "reference": "rohf", "fcae": "ae", "keywords": {"basis": "<>", "cfour_reference": "rohf", "cfour_scf_conv": 12},                                                }, id="hf rohf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-hf", "reference": "rohf", "fcae": "ae", "keywords": {"basis": "<>", "gamess_contrl__scftyp": "rohf"},                                                                }, id="hf rohf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-hf", "reference": "rohf", "fcae": "ae", "keywords": {"basis": "<>", "nwchem_scf__rohf": True},                                                                       }, id="hf rohf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-hf",  "reference": "rohf", "fcae": "ae", "keywords": {"basis": "<>", "reference": "rohf", "psi4_scf_type": "pk"},                                                     }, id="hf rohf ae: psi4",       marks=using("psi4")),
        # yapf: enable
    ],
)
def test_hf_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "energy"))


#
#  ,--.  ,--.,------.     ,----.                     ,--.,--.                 ,--.
#  |  '--'  ||  .---'    '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  .--.  ||  `--,     |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  |  |  |  ||  |`       '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#  `--'  `--'`--'         `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  HF Gradient


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        # pytest.param(0, id="grd0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "c4-hf",  "reference": "rhf",  "fcae": "ae", "keywords": {"basis": "<>", "cfour_scf_conv": 12},                                                                           }, id="hf  rhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-hf", "reference": "rhf",  "fcae": "ae", "keywords": {"basis": "<>"},                                                                                                 }, id="hf  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-hf", "reference": "rhf",  "fcae": "ae", "keywords": {"basis": "<>"},                                                                                                 }, id="hf  rhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-hf",  "reference": "rhf",  "fcae": "ae", "keywords": {"basis": "<>", "psi4_scf_type": "pk"},                                                                          }, id="hf  rhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-hf",  "reference": "uhf",  "fcae": "ae", "keywords": {"basis": "<>", "cfour_reference": "uhf", "cfour_scf_conv": 12},                                                 }, id="hf  uhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-hf", "reference": "uhf",  "fcae": "ae", "keywords": {"basis": "<>", "gamess_contrl__scftyp": "uhf"},                                                                 }, id="hf  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-hf", "reference": "uhf",  "fcae": "ae", "keywords": {"basis": "<>", "nwchem_scf__uhf": True},                                                                        }, id="hf  uhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-hf",  "reference": "uhf",  "fcae": "ae", "keywords": {"basis": "<>", "reference": "uhf", "psi4_scf_type": "pk"},                                                      }, id="hf  uhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-hf",  "reference": "rohf", "fcae": "ae", "keywords": {"basis": "<>", "cfour_reference": "rohf", "cfour_scf_conv": 12},                                                }, id="hf rohf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-hf", "reference": "rohf", "fcae": "ae", "keywords": {"basis": "<>", "gamess_contrl__scftyp": "rohf"},                                                                }, id="hf rohf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-hf", "reference": "rohf", "fcae": "ae", "keywords": {"basis": "<>", "nwchem_scf__rohf": True},                                                                       }, id="hf rohf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-hf",  "reference": "rohf", "fcae": "ae", "keywords": {"basis": "<>", "reference": "rohf", "psi4_scf_type": "pk"},                                                     }, id="hf rohf ae: psi4",       marks=using("psi4")),
        # yapf: enable
    ],
)
def test_hf_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "gradient"))


#
#  ,--.  ,--.,------.    ,--.  ,--.                     ,--.
#  |  '--'  ||  .---'    |  '--'  | ,---.  ,---.  ,---. `--' ,--,--.,--,--,
#  |  .--.  ||  `--,     |  .--.  || .-. :(  .-' (  .-' ,--.' ,-.  ||      \
#  |  |  |  ||  |`       |  |  |  |\   --..-'  `).-'  `)|  |\ '-'  ||  ||  |
#  `--'  `--'`--'        `--'  `--' `----'`----' `----' `--' `--`--'`--''--'
#
#  <<<  HF Hessian


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(2, id="hes2"),
        # pytest.param(1, id="hes1", marks=pytest.mark.long),
        # pytest.param(0, id="hes0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "c4-hf",  "reference": "rhf",  "fcae": "ae", "keywords": {"basis": "<>", "cfour_scf_conv": 12},                                                                           }, id="hf  rhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-hf", "reference": "rhf",  "fcae": "ae", "keywords": {"basis": "<>"},                                                                                                 }, id="hf  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-hf", "reference": "rhf",  "fcae": "ae", "keywords": {"basis": "<>"},                                                                                                 }, id="hf  rhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-hf",  "reference": "rhf",  "fcae": "ae", "keywords": {"basis": "<>", "psi4_scf_type": "pk"},                                                                          }, id="hf  rhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-hf",  "reference": "uhf",  "fcae": "ae", "keywords": {"basis": "<>", "cfour_reference": "uhf", "cfour_scf_conv": 12},                                                 }, id="hf  uhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-hf", "reference": "uhf",  "fcae": "ae", "keywords": {"basis": "<>", "gamess_contrl__scftyp": "uhf"},                                                                 }, id="hf  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-hf", "reference": "uhf",  "fcae": "ae", "keywords": {"basis": "<>", "nwchem_scf__uhf": True},                                                                        }, id="hf  uhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-hf",  "reference": "uhf",  "fcae": "ae", "keywords": {"basis": "<>", "reference": "uhf", "psi4_scf_type": "pk"},                                                      }, id="hf  uhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-hf",  "reference": "rohf", "fcae": "ae", "keywords": {"basis": "<>", "cfour_reference": "rohf", "cfour_scf_conv": 12},                                                }, id="hf rohf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-hf", "reference": "rohf", "fcae": "ae", "keywords": {"basis": "<>", "gamess_contrl__scftyp": "rohf", "gamess_scf__dirscf": True},                                    }, id="hf rohf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-hf", "reference": "rohf", "fcae": "ae", "keywords": {"basis": "<>", "nwchem_scf__rohf": True},                                                                       }, id="hf rohf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-hf",  "reference": "rohf", "fcae": "ae", "keywords": {"basis": "<>", "reference": "rohf", "psi4_scf_type": "pk"},                                    "error": {2: _q6}}, id="hf rohf ae: psi4",       marks=using("psi4")),
        # yapf: enable
    ],
)
def test_hf_hessian_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "hessian"))


#
#
#  ,--.   ,--.,------.  ,---.     ,------.
#  |   `.'   ||  .--. ''.-.  \    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |'.'|  ||  '--' | .-' .'    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  |   |  ||  | --' /   '-.    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `--'   `--'`--'     '-----'    `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                               `---' `---'
#  <<<  MP2 Energy


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ]
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p"),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "c4-mp2",  "reference": "rhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", "cfour_scf_conv": 12},                                                                     }, id="mp2  rhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "c4-mp2",  "reference": "rhf",  "fcae": "ae", "keywords": {"basis": "<>", "cfour_scf_conv": 12},                                                                           }, id="mp2  rhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-mp2", "reference": "rhf",  "fcae": "ae", "keywords": {"basis": "<>", "gamess_mp2__nacore": 0},                                                                        }, id="mp2  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp2", "reference": "rhf",  "fcae": "ae", "keywords": {"basis": "<>", "qc_module": "tce"},                                                                             }, id="mp2  rhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-mp2", "reference": "rhf",  "fcae": "ae", "keywords": {"basis": "<>"},                                                                                                 }, id="mp2  rhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-mp2",  "reference": "rhf",  "fcae": "ae", "keywords": {"basis": "<>", "psi4_mp2_type": "conv"},                                                                        }, id="mp2  rhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-mp2",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", "cfour_dropmo": [1], "cfour_scf_conv": 12},                                                }, id="mp2  rhf fc: cfour",      marks=using("cfour")),
        pytest.param({"call": "c4-mp2",  "reference": "rhf",  "fcae": "fc", "keywords": {"basis": "<>", "cfour_dropmo": 1, "cfour_scf_conv": 12},                                                        }, id="mp2  rhf fc: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-mp2", "reference": "rhf",  "fcae": "fc", "keywords": {"basis": "<>"}                                                                                                  }, id="mp2  rhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp2", "reference": "rhf",  "fcae": "fc", "keywords": {"basis": "<>", "qc_module": "tce", "nwchem_tce__freeze": 1},                                                    }, id="mp2  rhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-mp2", "reference": "rhf",  "fcae": "fc", "keywords": {"basis": "<>", "nwchem_mp2__freeze": 1},                                                                        }, id="mp2  rhf fc: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "nwc-mp2", "reference": "rhf",  "fcae": "fc", "keywords": {"basis": "<>", "nwchem_mp2__freeze__core": 1},                                                                  }, id="mp2  rhf fc: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "nwc-mp2", "reference": "rhf",  "fcae": "fc", "keywords": {"basis": "<>", "nwchem_mp2__freeze__core__atomic": True},                                                       }, id="mp2  rhf fc: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "nwc-mp2", "reference": "rhf",  "fcae": "fc", "keywords": {"basis": "<>", "nwchem_mp2__freeze__atomic": {"O": 1}},                                                         }, id="mp2  rhf fc: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-mp2",  "reference": "rhf",  "fcae": "fc", "keywords": {"basis": "<>", "psi4_freeze_core": True, "psi4_mp2_type": "conv"},                                              }, id="mp2  rhf fc: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-mp2",  "reference": "uhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", "cfour_reference": "uhf", "cfour_scf_conv": 12},                                           }, id="mp2  uhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "c4-mp2",  "reference": "uhf",  "fcae": "ae", "keywords": {"basis": "<>", "cfour_reference": "uhf", "cfour_scf_conv": 12},                                                 }, id="mp2  uhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-mp2", "reference": "uhf",  "fcae": "ae", "keywords": {"basis": "<>", "gamess_contrl__scftyp": "uhf", "gamess_mp2__nacore": 0},                                        }, id="mp2  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp2", "reference": "uhf",  "fcae": "ae", "keywords": {"basis": "<>", "qc_module": "tce", "nwchem_scf__uhf": True},                                                    }, id="mp2  uhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-mp2", "reference": "uhf",  "fcae": "ae", "keywords": {"basis": "<>", "nwchem_scf__uhf": True},                                                                        }, id="mp2  uhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-mp2",  "reference": "uhf",  "fcae": "ae", "keywords": {"basis": "<>", "reference": "uhf", "psi4_mp2_type": "conv"},                                                    }, id="mp2  uhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-mp2",  "reference": "uhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", "cfour_reference": "uhf", "cfour_dropmo": [1], "cfour_scf_conv": 12},                      }, id="mp2  uhf fc: cfour",      marks=using("cfour")),
        pytest.param({"call": "c4-mp2",  "reference": "uhf",  "fcae": "fc", "keywords": {"basis": "<>", "cfour_reference": "uhf", "cfour_dropmo": 1, "cfour_scf_conv": 12},                              }, id="mp2  uhf fc: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-mp2", "reference": "uhf",  "fcae": "fc", "keywords": {"basis": "<>", "gamess_contrl__scftyp": "uhf"},                                                                 }, id="mp2  uhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp2", "reference": "uhf",  "fcae": "fc", "keywords": {"basis": "<>", "qc_module": "tce", "nwchem_tce__freeze": 1, "nwchem_scf__uhf": True},                           }, id="mp2  uhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-mp2", "reference": "uhf",  "fcae": "fc", "keywords": {"basis": "<>", "nwchem_scf__uhf": True, "nwchem_mp2__freeze": 1},                                               }, id="mp2  uhf fc: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-mp2",  "reference": "uhf",  "fcae": "fc", "keywords": {"basis": "<>", "reference": "uhf", "psi4_freeze_core": True, "psi4_mp2_type": "conv"},                          }, id="mp2  uhf fc: psi4",       marks=using("psi4")),

            # "cfour_occupation": [[3, 1, 1, 0], [3, 0, 1, 0]]
        pytest.param({"call": "c4-mp2",  "reference": "rohf", "fcae": "ae", "keywords": {"cfour_basis": "<>", "cfour_reference": "rohf", "cfour_scf_conv": 12},                                          }, id="mp2 rohf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "c4-mp2",  "reference": "rohf", "fcae": "ae", "keywords": {"basis": "<>", "cfour_reference": "rohf", "cfour_scf_conv": 12},                                                }, id="mp2 rohf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-mp2", "reference": "rohf", "fcae": "ae", "keywords": {"basis": "<>", "gamess_contrl__scftyp": "rohf", "gamess_mp2__nacore": 0, "gamess_mp2__ospt": "RMP"},            }, id="mp2 rohf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp2", "reference": "rohf", "fcae": "ae", "keywords": {"basis": "<>", "qc_module": "tce", "nwchem_scf__rohf": True, "nwchem_scf__thresh": 8,
                                                                                         "nwchem_tce__thresh": 8, "nwchem_tce__freeze": 0, "nwchem_scf__tol2e": 10},                         "wrong": _w1}, id="mp2 rohf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-mp2", "reference": "rohf", "fcae": "ae", "keywords": {"basis": "<>", "nwchem_scf__rohf": True},                                                           "error": _q1}, id="mp2 rohf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-mp2",  "reference": "rohf", "fcae": "ae", "keywords": {"basis": "<>", "reference": "rohf", "psi4_mp2_type": "conv"},                                                   }, id="mp2 rohf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-mp2",  "reference": "rohf", "fcae": "fc", "keywords": {"cfour_basis": "<>", "cfour_reference": "rohf", "cfour_dropmo": [1], "cfour_scf_conv": 12},                     }, id="mp2 rohf fc: cfour",      marks=using("cfour")),
        pytest.param({"call": "c4-mp2",  "reference": "rohf", "fcae": "fc", "keywords": {"basis": "<>", "cfour_reference": "rohf", "cfour_dropmo": 1, "cfour_scf_conv": 12},                             }, id="mp2 rohf fc: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-mp2", "reference": "rohf", "fcae": "fc", "keywords": {"basis": "<>", "gamess_contrl__scftyp": "rohf", "gamess_mp2__ospt": "RMP"},                                     }, id="mp2 rohf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp2", "reference": "rohf", "fcae": "fc", "keywords": {"basis": "<>", "qc_module": "tce", "nwchem_tce__freeze": 1, "nwchem_scf__rohf": True},              "wrong": _w1}, id="mp2 rohf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-mp2", "reference": "rohf", "fcae": "fc", "keywords": {"basis": "<>", "nwchem_scf__rohf": True, "nwchem_mp2__freeze": 1},                                  "error": _q1}, id="mp2 rohf fc: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-mp2",  "reference": "rohf", "fcae": "fc", "keywords": {"basis": "<>", "reference": "rohf", "psi4_freeze_core": True, "psi4_mp2_type": "conv"},                         }, id="mp2 rohf fc: psi4",       marks=using("psi4")),
        # yapf: enable
    ],
)
def test_mp2_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*energy_processor(inp, dertype, basis, subjects, clsd_open_pmols, request))


#
#  ,--.   ,--.,------.  ,---.      ,----.                     ,--.,--.                 ,--.
#  |   `.'   ||  .--. ''.-.  \    '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  |'.'|  ||  '--' | .-' .'    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  |  |   |  ||  | --' /   '-.    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#  `--'   `--'`--'     '-----'     `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  MP2 Gradient


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        # pytest.param(0, id="grd0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are all possible ways of computing <method> working?

        pytest.param({"call": "c4-mp2",  "reference": "rhf",  "fcae": "ae", "keywords": {"basis": "<>", "cfour_scf_conv": 12},                                                                                      }, id="mp2  rhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-mp2", "reference": "rhf",  "fcae": "ae", "keywords": {"basis": "<>", "gamess_mp2__nacore": 0},                                                                                   }, id="mp2  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp2", "reference": "rhf",  "fcae": "ae", "keywords": {"basis": "<>"},                                                                                                            }, id="mp2  rhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-mp2",  "reference": "rhf",  "fcae": "ae", "keywords": {"basis": "<>", "psi4_mp2_type": "conv"},                                                                                   }, id="mp2  rhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-mp2",  "reference": "rhf",  "fcae": "fc", "keywords": {"basis": "<>", "cfour_dropmo": 1, "cfour_scf_conv": 12},                                                                   }, id="mp2  rhf fc: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-mp2", "reference": "rhf",  "fcae": "fc", "keywords": {"basis": "<>"},                                                                                                            }, id="mp2  rhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp2", "reference": "rhf",  "fcae": "fc", "keywords": {"basis": "<>", "nwchem_mp2__freeze": 1},                                                                                   }, id="mp2  rhf fc: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-mp2",  "reference": "rhf",  "fcae": "fc", "keywords": {"basis": "<>", "psi4_freeze_core": True, "psi4_mp2_type": "conv", "psi4_points": 5},                                       }, id="mp2  rhf fc: psi4",       marks=using("psi4")),  # findif

        pytest.param({"call": "c4-mp2",  "reference": "uhf",  "fcae": "ae", "keywords": {"basis": "<>", "cfour_reference": "uhf", "cfour_scf_conv": 12},                                                            }, id="mp2  uhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-mp2", "reference": "uhf",  "fcae": "ae", "keywords": {"basis": "<>", "gamess_contrl__scftyp": "uhf", "gamess_mp2__nacore": 0},                                                   }, id="mp2  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp2", "reference": "uhf",  "fcae": "ae", "keywords": {"basis": "<>", "nwchem_scf__uhf": True},                                                                                   }, id="mp2  uhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-mp2",  "reference": "uhf",  "fcae": "ae", "keywords": {"basis": "<>", "reference": "uhf", "psi4_mp2_type": "conv"},                                                               }, id="mp2  uhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-mp2",  "reference": "uhf",  "fcae": "fc", "keywords": {"basis": "<>", "cfour_reference": "uhf", "cfour_dropmo": [1], "cfour_scf_conv": 12},                                       }, id="mp2  uhf fc: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-mp2", "reference": "uhf",  "fcae": "fc", "keywords": {"basis": "<>", "gamess_contrl__scftyp": "uhf"},                                                                            }, id="mp2  uhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp2", "reference": "uhf",  "fcae": "fc", "keywords": {"basis": "<>", "nwchem_scf__uhf": True, "nwchem_mp2__freeze": 1},                                                          }, id="mp2  uhf fc: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-mp2",  "reference": "uhf",  "fcae": "fc", "keywords": {"basis": "<>", "reference": "uhf", "psi4_mp2_type": "conv", "psi4_freeze_core": True, "psi4_points": 5},                   }, id="mp2  uhf fc: psi4",       marks=using("psi4")),  # findif

        pytest.param({"call": "c4-mp2",  "reference": "rohf", "fcae": "ae", "keywords": {"basis": "<>", "cfour_reference": "rohf", "cfour_scf_conv": 12},                                                           }, id="mp2 rohf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-mp2", "reference": "rohf", "fcae": "ae", "keywords": {"basis": "<>", "gamess_contrl__scftyp": "rohf", "gamess_mp2__nacore": 0, "gamess_system__memddi": 300},   "wrong": {1: _w5}}, id="mp2 rohf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp2", "reference": "rohf", "fcae": "ae", "keywords": {"basis": "<>", "nwchem_scf__rohf": True},                                                                 "error": {1: _q1}}, id="mp2 rohf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-mp2",  "reference": "rohf", "fcae": "ae", "keywords": {"basis": "<>", "reference": "rohf", "psi4_mp2_type": "conv", "psi4_points": 5},                                            }, id="mp2 rohf ae: psi4",       marks=using("psi4")),  # findif

        pytest.param({"call": "c4-mp2",  "reference": "rohf", "fcae": "fc", "keywords": {"basis": "<>", "cfour_reference": "rohf", "cfour_dropmo": [1], "cfour_scf_conv": 12},                                      }, id="mp2 rohf fc: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-mp2", "reference": "rohf", "fcae": "fc", "keywords": {"basis": "<>", "gamess_contrl__scftyp": "rohf", "gamess_system__memddi": 300},                            "wrong": {1: _w5}}, id="mp2 rohf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp2", "reference": "rohf", "fcae": "fc", "keywords": {"basis": "<>", "nwchem_scf__rohf": True, "nwchem_mp2__freeze": 1},                                        "error": {1: _q1}}, id="mp2 rohf fc: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-mp2",  "reference": "rohf", "fcae": "fc", "keywords": {"basis": "<>", "reference": "rohf", "psi4_mp2_type": "conv", "psi4_freeze_core": True, "psi4_points": 5},                  }, id="mp2 rohf fc: psi4",       marks=using("psi4")),  # findif
        # yapf: enable
    ],
)
def test_mp2_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*gradient_processor(inp, dertype, basis, subjects, clsd_open_pmols, request))


#
#  ,--.   ,--.,------. ,----.     ,------.
#  |   `.'   ||  .--. ''.-.  |    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |'.'|  ||  '--' |  .' <     |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  |   |  ||  | --' /'-'  |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `--'   `--'`--'     `----'     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                               `---' `---'
#  <<<  MP3 Energy


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ]
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p"),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # GAMESS doesn't do mp3. Psi4 only does rohf mp3 with detci, which is already non-std for mp2.
        # yapf: disable
        # ecc doesn't compute mp3
        pytest.param({"call": "c4-mp3",  "reference": "rhf",  "fcae": "ae", "keywords": {"basis": "<>", "cfour_scf_conv": 12, "cfour_cc_program": "vcc"},                                                }, id="mp3  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-mp3",  "reference": "rhf",  "fcae": "ae", "keywords": {"basis": "<>", "cfour_scf_conv": 12, "cfour_cc_program": "ncc"},                                                }, id="mp3  rhf ae: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-mp3", "reference": "rhf",  "fcae": "ae", "keywords": {"basis": "<>"},                                                                               "error": {0: _q14}}, id="mp3  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp3", "reference": "rhf",  "fcae": "ae", "keywords": {"basis": "<>", "qc_module": "tce"},                                                                             }, id="mp3  rhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-mp3",  "reference": "rhf",  "fcae": "ae", "keywords": {"basis": "<>", "psi4_mp_type": "conv"},                                                                         }, id="mp3  rhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-mp3",  "reference": "rhf",  "fcae": "fc", "keywords": {"basis": "<>", "cfour_dropmo": 1, "cfour_scf_conv": 12, "cfour_cc_program": "vcc"},                             }, id="mp3  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-mp3",  "reference": "rhf",  "fcae": "fc", "keywords": {"basis": "<>", "cfour_dropmo": 1, "cfour_scf_conv": 12, "cfour_cc_program": "ncc"},                             }, id="mp3  rhf fc: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-mp3", "reference": "rhf",  "fcae": "fc", "keywords": {"basis": "<>"},                                                                               "error": {0: _q14}}, id="mp3  rhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp3", "reference": "rhf",  "fcae": "fc", "keywords": {"basis": "<>", "qc_module": "tce", "nwchem_tce__freeze": 1},                                                    }, id="mp3  rhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-mp3",  "reference": "rhf",  "fcae": "fc", "keywords": {"basis": "<>", "psi4_freeze_core": True, "psi4_mp_type": "conv"},                                               }, id="mp3  rhf fc: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-mp3",  "reference": "uhf",  "fcae": "ae", "keywords": {"basis": "<>", "cfour_reference": "uhf", "cfour_scf_conv": 12},                                                 }, id="mp3  uhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-mp3", "reference": "uhf",  "fcae": "ae", "keywords": {"basis": "<>", "gamess_contrl__scftyp": "uhf"},                                               "error": {0: _q14}}, id="mp3  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp3", "reference": "uhf",  "fcae": "ae", "keywords": {"basis": "<>", "qc_module": "tce", "nwchem_scf__uhf": True},                                                    }, id="mp3  uhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-mp3",  "reference": "uhf",  "fcae": "ae", "keywords": {"basis": "<>", "reference": "uhf", "psi4_mp_type": "conv"},                                                     }, id="mp3  uhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-mp3",  "reference": "uhf",  "fcae": "fc", "keywords": {"basis": "<>", "cfour_reference": "uhf", "cfour_dropmo": 1, "cfour_scf_conv": 12},                              }, id="mp3  uhf fc: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-mp3", "reference": "uhf",  "fcae": "fc", "keywords": {"basis": "<>", "gamess_contrl__scftyp": "uhf"},                                               "error": {0: _q14}}, id="mp3  uhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp3", "reference": "uhf",  "fcae": "fc", "keywords": {"basis": "<>", "nwchem_tce__freeze": 1, "qc_module": "tce", "nwchem_scf__uhf": True},                           }, id="mp3  uhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-mp3",  "reference": "uhf",  "fcae": "fc", "keywords": {"basis": "<>", "reference": "uhf", "psi4_freeze_core": True, "psi4_mp_type": "conv"},                           }, id="mp3  uhf fc: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-mp3",  "reference": "rohf", "fcae": "ae", "keywords": {"basis": "<>", "cfour_reference": "rohf", "cfour_scf_conv": 12},                                                }, id="mp3 rohf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-mp3", "reference": "rohf", "fcae": "ae", "keywords": {"basis": "<>", "gamess_contrl__scftyp": "rohf"},                                              "error": {0: _q14}}, id="mp3 rohf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp3", "reference": "rohf", "fcae": "ae", "keywords": {"basis": "<>", "qc_module": "tce", "nwchem_scf__rohf": True, "nwchem_scf__thresh": 1.e-8},    "wrong": {0: _w4} }, id="mp3 rohf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-mp3",  "reference": "rohf", "fcae": "ae", "keywords": {"basis": "<>", "reference": "rohf", "psi4_mp_type": "conv"},                                  "error": {0: _q15}}, id="mp3 rohf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-mp3",  "reference": "rohf", "fcae": "fc", "keywords": {"basis": "<>", "cfour_reference": "rohf", "cfour_dropmo": 1, "cfour_scf_conv": 12},                             }, id="mp3 rohf fc: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-mp3", "reference": "rohf", "fcae": "fc", "keywords": {"basis": "<>", "gamess_contrl__scftyp": "rohf"},                                              "error": {0: _q14}}, id="mp3 rohf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp3", "reference": "rohf", "fcae": "fc", "keywords": {"basis": "<>", "nwchem_tce__freeze": 1, "qc_module": "tce", "nwchem_scf__rohf": True},        "wrong": {0: _w4} }, id="mp3 rohf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-mp3",  "reference": "rohf", "fcae": "fc", "keywords": {"basis": "<>", "reference": "rohf", "psi4_freeze_core": True, "psi4_mp_type": "conv"},        "error": {0: _q15}}, id="mp3 rohf fc: psi4",       marks=using("psi4")),
        # yapf: enable
    ],
)
def test_mp3_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "energy"))


#
#  ,--.   ,--.,------. ,----.      ,----.                     ,--.,--.                 ,--.
#  |   `.'   ||  .--. ''.-.  |    '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  |'.'|  ||  '--' |  .' <     |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  |  |   |  ||  | --' /'-'  |    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#  `--'   `--'`--'     `----'      `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  MP3 Gradient

@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        # pytest.param(0, id="grd0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "c4-mp3",  "reference": "rhf",  "fcae": "ae",                        "keywords": {"basis": "<>", "cfour_scf_conv": 12, "cfour_cc_program": "vcc"},                                                           }, id="mp3  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-mp3",  "reference": "rhf",  "fcae": "ae",                        "keywords": {"basis": "<>", "cfour_scf_conv": 12, "cfour_cc_program": "ncc"},                                                           }, id="mp3  rhf ae: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-mp3", "reference": "rhf",  "fcae": "ae",                        "keywords": {"basis": "<>"},                                                                                          "error": {1: _q16}}, id="mp3  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp3", "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {"basis": "<>", "qc_module": "tce"},                                                                                        }, id="mp3  rhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-mp3",  "reference": "rhf",  "fcae": "ae",                        "keywords": {"basis": "<>", "psi4_mp_type": "conv", },                                                                                  }, id="mp3  rhf ae: psi4",       marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-mp3", "reference": "rhf", "fcae": "ae", "keywords": {"psi4_dertype": "none", "psi4_cfour_cc_program": "vcc", "psi4_cfour_SCF_CONV": 12, "psi4_cfour_CC_CONV": 12, "psi4_cfour_LINEQ_CONV": 11, "psi4_points": 5, "psi4_fd_project": False}}, id="mp3  rhf ae: psi4-cfour-vcc"),
        # DEBUG pytest.param({"call": "p4-c4-mp3", "reference": "rhf", "fcae": "ae", "keywords": {"psi4_dertype": "none", "psi4_cfour_cc_program": "ncc", "psi4_cfour_SCF_CONV": 12, "psi4_cfour_CC_CONV": 12, "psi4_cfour_LINEQ_CONV": 11, "psi4_points": 5, "psi4_fd_project": False}}, id="mp3  rhf ae: psi4-cfour-ncc"),

        pytest.param({"call": "c4-mp3",  "reference": "rhf",  "fcae": "fc",                        "keywords": {"basis": "<>", "cfour_dropmo": 1, "cfour_scf_conv": 12, "cfour_cc_program": "vcc"},                      "wrong": {1: _w11}}, id="mp3  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-mp3",  "reference": "rhf",  "fcae": "fc",                        "keywords": {"basis": "<>", "cfour_dropmo": 1, "cfour_scf_conv": 12, "cfour_cc_program": "ncc"},                      "wrong": {1: _w10}}, id="mp3  rhf fc: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-mp3", "reference": "rhf",  "fcae": "fc",                        "keywords": {"basis": "<>"},                                                                                          "error": {1: _q16}}, id="mp3  rhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp3", "reference": "rhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {"basis": "<>", "qc_module": "tce", "nwchem_tce__freeze": 1},                                                               }, id="mp3  rhf fc: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-mp3",  "reference": "rhf",  "fcae": "fc", "xptd": {"fd": False}, "keywords": {"basis": "<>", "psi4_freeze_core": True, "psi4_mp_type": "conv", "psi4_points": 5, },                                      }, id="mp3  rhf fc: psi4",       marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-mp3", "reference": "rhf", "fcae": "fc", "keywords": {"psi4_cfour_dropmo": [1], "psi4_dertype": "none", "psi4_cfour_cc_program": "vcc", "psi4_cfour_SCF_CONV": 12, "psi4_cfour_CC_CONV": 12, "psi4_cfour_LINEQ_CONV": 11, "psi4_points": 5, "psi4_fd_project": False}}, id="mp3  rhf fc: psi4-cfour-vcc"),
        # DEBUG pytest.param({"call": "p4-c4-mp3", "reference": "rhf", "fcae": "fc", "keywords": {"psi4_cfour_dropmo": [1], "psi4_dertype": "none", "psi4_cfour_cc_program": "ncc", "psi4_cfour_SCF_CONV": 12, "psi4_cfour_CC_CONV": 12, "psi4_cfour_LINEQ_CONV": 11, "psi4_points": 5, "psi4_fd_project": False}}, id="mp3  rhf fc: psi4-cfour-ncc"),

        pytest.param({"call": "c4-mp3",  "reference": "uhf",  "fcae": "ae",                        "keywords": {"basis": "<>", "cfour_reference": "uhf", "cfour_scf_conv": 12, "cfour_cc_program": "vcc"},                                 }, id="mp3  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "gms-mp3", "reference": "uhf",  "fcae": "ae",                        "keywords": {"basis": "<>", "gamess_contrl__scftyp": "uhf"},                                                          "error": {1: _q16}}, id="mp3  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp3", "reference": "uhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {"basis": "<>", "qc_module": "tce", "nwchem_scf__uhf": True},                                                               }, id="mp3  uhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-mp3",  "reference": "uhf",  "fcae": "ae",                        "keywords": {"basis": "<>", "reference": "uhf", "psi4_mp_type": "conv"},                                                                }, id="mp3  uhf ae: psi4",       marks=using("psi4")),
        pytest.param({"call": "p4-mp3",  "reference": "uhf",  "fcae": "ae",                        "keywords": {"basis": "<>", "reference": "uhf", "psi4_mp_type": "conv", "psi4_points": 5, "psi4_dertype": "none"},                                                                }, id="mp3  uhf ae: psi4",       marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-mp3", "reference": "uhf", "fcae": "ae", "keywords": {"psi4_dertype": "none", "psi4_cfour_cc_program": "vcc", "psi4_cfour_SCF_CONV": 12, "psi4_cfour_CC_CONV": 12, "psi4_cfour_LINEQ_CONV": 11, "psi4_points": 5, "psi4_fd_project": False, "psi4_cfour_reference": "uhf"}}, id="mp3  uhf ae: psi4-cfour-vcc"),

        pytest.param({"call": "c4-mp3",  "reference": "uhf",  "fcae": "fc",                        "keywords": {"basis": "<>", "cfour_reference": "uhf", "cfour_dropmo": [1], "cfour_scf_conv": 12},                                       }, id="mp3  uhf fc: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-mp3", "reference": "uhf",  "fcae": "fc",                        "keywords": {"basis": "<>", "gamess_contrl__scftyp": "uhf"},                                                          "error": {1: _q16}}, id="mp3  uhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp3", "reference": "uhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {"basis": "<>", "qc_module": "tce", "nwchem_scf__uhf": True, "nwchem_tce__freeze": 1},                                      }, id="mp3  uhf fc: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-mp3",  "reference": "uhf",  "fcae": "fc", "xptd": {"fd": False}, "keywords": {"basis": "<>", "reference": "uhf", "psi4_mp_type": "conv", "psi4_freeze_core": True, "psi4_points": 5},                    }, id="mp3  uhf fc: psi4",       marks=using("psi4")),
        pytest.param({"call": "p4-mp3",  "reference": "uhf",  "fcae": "fc",                        "keywords": {"basis": "<>", "reference": "uhf", "psi4_mp_type": "conv", "psi4_points": 5, "psi4_freeze_core": True, "psi4_dertype": "none"},}, id="mp3  uhf fc: psi4",       marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-mp3", "reference": "uhf", "fcae": "fc", "keywords": {"psi4_cfour_dropmo": [1], "psi4_dertype": "none", "psi4_cfour_cc_program": "vcc", "psi4_cfour_SCF_CONV": 12, "psi4_cfour_CC_CONV": 12, "psi4_cfour_LINEQ_CONV": 11, "psi4_points": 5, "psi4_fd_project": False, "psi4_cfour_reference": "uhf"}}, id="mp3  uhf ae: psi4-cfour-vcc"),

        pytest.param({"call": "c4-mp3",  "reference": "rohf", "fcae": "ae", "keywords": {"basis": "<>", "cfour_reference": "rohf", "cfour_scf_conv": 12},                                                                                  }, id="mp3 rohf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-mp3", "reference": "rohf", "fcae": "ae", "keywords": {"basis": "<>", "gamess_contrl__scftyp": "rohf"},                                                                                "error": {1: _q16}}, id="mp3 rohf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp3", "reference": "rohf", "fcae": "ae", "keywords": {"basis": "<>", "qc_module": "tce", "nwchem_scf__rohf": True},                                                                   "wrong": {1: _w1} }, id="mp3 rohf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-mp3",  "reference": "rohf", "fcae": "ae", "keywords": {"basis": "<>", "reference": "rohf", "psi4_mp_type": "conv", "psi4_points": 5},                                                  "error": {1: _q15}}, id="mp3 rohf ae: psi4",       marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-mp3", "reference": "rohf", "fcae": "ae", "keywords": {"psi4_dertype": "none", "psi4_cfour_cc_program": "vcc", "psi4_cfour_SCF_CONV": 12, "psi4_cfour_CC_CONV": 12, "psi4_cfour_LINEQ_CONV": 11, "psi4_points": 5, "psi4_fd_project": False, "psi4_cfour_reference": "rohf"}}, id="mp3  rohf ae: psi4-cfour-vcc"),

        pytest.param({"call": "c4-mp3",  "reference": "rohf", "fcae": "fc", "keywords": {"basis": "<>", "cfour_reference": "rohf", "cfour_dropmo": [1], "cfour_scf_conv": 12},                                                             }, id="mp3 rohf fc: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-mp3", "reference": "rohf", "fcae": "fc", "keywords": {"basis": "<>", "gamess_contrl__scftyp": "rohf", "gamess_system__memddi": 300},                                                  "error": {1: _q16}}, id="mp3 rohf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp3", "reference": "rohf", "fcae": "fc", "keywords": {"basis": "<>", "qc_module": "tce", "nwchem_scf__rohf": True, "nwchem_tce__freeze": 1},                                          "wrong": {1: _w1} }, id="mp3 rohf fc: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-mp3",  "reference": "rohf", "fcae": "fc", "keywords": {"basis": "<>", "reference": "rohf", "psi4_mp_type": "conv", "psi4_freeze_core": True, "psi4_points": 5},                        "error": {1: _q15}}, id="mp3 rohf fc: psi4",       marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-mp3", "reference": "rohf", "fcae": "fc", "keywords": {"psi4_cfour_dropmo": [1], "psi4_dertype": "none", "psi4_cfour_cc_program": "vcc", "psi4_cfour_SCF_CONV": 12, "psi4_cfour_CC_CONV": 12, "psi4_cfour_LINEQ_CONV": 11, "psi4_points": 5, "psi4_fd_project": False, "psi4_cfour_reference": "rohf"}}, id="mp3  rohf fc: psi4-cfour-vcc"),
        # yapf: enable
    ],
)
def test_mp3_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "gradient"))


#
#   ,-----. ,-----.,------.      ,------.
#  '  .--./'  .--./|  .-.  \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |    |  |    |  |  \  :    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '--'\'  '--'\|  '--'  /    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----' `-----'`-------'     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                              `---' `---'
#  <<<  CCD Energy

@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p"),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # ecc doesn't run
        pytest.param({"call": "c4-ccd",  "reference": "rhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "vcc", "cfour_print": 2},                                                                  }, id="ccd  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccd",  "reference": "rhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "ecc",},                                                                                   }, id="ccd  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccd",  "reference": "rhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "ncc",},                                                                                   }, id="ccd  rhf ae: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccd", "reference": "rhf",  "fcae": "ae", "keywords": {"gamess_ccinp__ncore": 0},                                                                                                                                                     }, id="ccd  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccd", "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "tce"},                                                                                                                                                           }, id="ccd  rhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-ccd",  "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                                                           "error": {0: _q11}}, id="ccd  rhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-ccd",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", "cfour_dropmo": [1], "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "vcc", "cfour_print": 2},                                             }, id="ccd  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccd",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", "cfour_dropmo": [1], "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "ecc"},                                                               }, id="ccd  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccd",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", "cfour_dropmo": [1], "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "ncc"},                                                               }, id="ccd  rhf fc: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccd", "reference": "rhf",  "fcae": "fc", "keywords": {},                                                                                                                                                                             }, id="ccd  rhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccd", "reference": "rhf",  "fcae": "fc", "keywords": {"qc_module": "tce", "nwchem_tce__freeze": 1 },                                                                                                                                 }, id="ccd  rhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-ccd",  "reference": "rhf",  "fcae": "fc", "keywords": {"psi4_freeze_core": True},                                                                                                                                   "error": {0: _q11}}, id="ccd  rhf fc: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-ccd",  "reference": "uhf",  "fcae": "ae", "keywords": {"cfour_BASIS": "<>", "cfour_REFerence": "UHF", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "vcc", "cfour_print": 2},                                        }, id="ccd  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccd",  "reference": "uhf",  "fcae": "ae", "keywords": {"cfour_BASIS": "<>", "cfour_REFerence": "UHF", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "ecc"},                                                          }, id="ccd  uhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccd", "reference": "uhf",  "fcae": "ae", "keywords": {"gamess_contrl__scftyp": "uhf", "gamess_ccinp__ncore": 0},                                                                                                   "error": {0: _q2} }, id="ccd  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccd", "reference": "uhf",  "fcae": "ae", "keywords": {"qc_module": "tce", "nwchem_scf__uhf": True},                                                                                                                                  }, id="ccd  uhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-ccd",  "reference": "uhf",  "fcae": "ae", "keywords": {"psi4_qc_module": "occ"},                                                                                                                                    "error": {0: _q11}}, id="ccd  uhf ae: psi4-occ",   marks=using("psi4")),

        pytest.param({"call": "c4-ccd",  "reference": "uhf",  "fcae": "fc", "keywords": {"cfour_BASIS": "<>", "cfour_dropmo": [1], "cfour_REFerence": "UHF", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "vcc", "cfour_print": 2},                   }, id="ccd  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccd",  "reference": "uhf",  "fcae": "fc", "keywords": {"cfour_BASIS": "<>", "cfour_dropmo": [1], "cfour_REFerence": "UHF", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "ecc"},                                     }, id="ccd  uhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccd", "reference": "uhf",  "fcae": "fc", "keywords": {"gamess_contrl__scftyp": "uhf"},                                                                                                                             "error": {0: _q2} }, id="ccd  uhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccd", "reference": "uhf",  "fcae": "fc", "keywords": {"nwchem_tce__freeze": 1, "qc_module": "tce", "nwchem_scf__uhf": True},                                                                                                         }, id="ccd  uhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-ccd",  "reference": "uhf",  "fcae": "fc", "keywords": {"psi4_freeze_core": True, "reference": "uhf", "psi4_qc_module": "occ"},                                                                                      "error": {0: _q11}}, id="ccd  uhf fc: psi4-occ",   marks=using("psi4")),

        # rohf vcc = tce, but cfour paper disavows rohf, so I'm suspicious
        # pytest.param({"call": "c4-ccd",  "reference": "rohf", "fcae": "ae", "keywords": {"cfour_basis": "<>", "cfour_REFerence": "roHF", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "vcc", "cfour_print": 2},                                                        }, id="ccd rohf ae: cfour-vcc",  marks=using("cfour")),
        # pytest.param({"call": "gms-ccd", "reference": "rohf", "fcae": "ae", "keywords": {"gamess_contrl__scftyp": "rohf", "gamess_ccinp__ncore": 0, "gamess_ccinp__maxcc": 50},                                                                                        "error": {0: _q4} }, id="ccd rohf ae: gamess",     marks=using("gamess")),
        # pytest.param({"call": "nwc-ccd", "reference": "rohf", "fcae": "ae", "keywords": {"qc_module": "tce", "nwchem_scf__rohf": True},                                                                                                                                                  }, id="ccd rohf ae: nwchem-tce", marks=using("nwchem")),
        # pytest.param({"call": "p4-ccd",  "reference": "rohf", "fcae": "ae", "keywords": {"reference": "rohf"},                                                                                                                                                         "error": {0: _q11}}, id="ccd rohf ae: psi4",       marks=using("psi4")),

        # pytest.param({"call": "c4-ccd",  "reference": "rohf", "fcae": "fc", "keywords": {"cfour_BASIS": "<>", "cfour_dropmo": [1], "cfour_REFerence": "roHF", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_orbitals": 0, "cfour_cc_program": "vcc", "cfour_print": 2},              }, id="ccd rohf fc: cfour-vcc",  marks=using("cfour")),
        # pytest.param({"call": "gms-ccd", "reference": "rohf", "fcae": "fc", "keywords": {"gamess_contrl__scftyp": "rohf", "gamess_ccinp__iconv": 9, "gamess_scf__conv": 9},                                                                                                 "error": {0: _q4}}, id="ccd rohf fc: gamess",     marks=using("gamess")),
        # pytest.param({"call": "nwc-ccd", "reference": "rohf", "fcae": "fc", "keywords": {"nwchem_tce__freeze": 1, "qc_module": "tce", "nwchem_scf__rohf": True},                                                                                                                         }, id="ccd rohf fc: nwchem-tce", marks=using("nwchem")),
        # pytest.param({"call": "p4-ccd",  "reference": "rohf", "fcae": "fc", "keywords": {"psi4_e_convergence": 8, "psi4_r_convergence": 7, "psi4_freeze_core": True, "reference": "rohf"},                                                                                  "error": {0: _q11}}, id="ccd rohf fc: psi4",       marks=using("psi4")),
        # yapf: enable
    ],
)
def test_ccd_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "energy"))


#
#   ,-----. ,-----.,------.       ,----.                     ,--.,--.                 ,--.
#  '  .--./'  .--./|  .-.  \     '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  |    |  |    |  |  \  :    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  '  '--'\'  '--'\|  '--'  /    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#   `-----' `-----'`-------'      `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  CCD Gradient


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        # pytest.param(0, id="grd0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "c4-ccd",  "reference": "rhf",  "fcae": "ae",                        "keywords": {"basis": "<>", "cfour_scf_conv": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "vcc", "cfour_print": 2},                                                                    }, id="ccd  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccd",  "reference": "rhf",  "fcae": "ae",                        "keywords": {"basis": "<>", "cfour_scf_conv": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "ecc"},                                                                                      }, id="ccd  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccd",  "reference": "rhf",  "fcae": "ae",                        "keywords": {"basis": "<>", "cfour_scf_conv": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "ncc"},                                                                                      }, id="ccd  rhf ae: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccd", "reference": "rhf",  "fcae": "ae",                        "keywords": {"basis": "<>", "gamess_ccinp__ncore": 0},                                                                                                                "error": {1: _q12}}, id="ccd  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccd", "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {"basis": "<>", "qc_module": "tce"},                                                                                                                                        }, id="ccd  rhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-ccd",  "reference": "rhf",  "fcae": "ae",                        "keywords": {"basis": "<>", "psi4_cc_type": "conv"},                                                                                                                  "error": {1: _q11}}, id="ccd  rhf ae: psi4",       marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-ccd", "reference": "rhf", "fcae": "ae",                       "keywords": {"psi4_dertype": "none", "psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "vcc", **_p4c4_fd},                                                                        }, id="ccd  rhf fc: psi4-cfour-vcc"),
        # DEBUG pytest.param({"call": "p4-c4-ccd", "reference": "rhf", "fcae": "ae",                       "keywords": {"psi4_dertype": "none", "psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "ncc", **_p4c4_fd},                                                                        }, id="ccd  rhf fc: psi4-cfour-ncc"),

        pytest.param({"call": "c4-ccd",  "reference": "rhf",  "fcae": "fc",                        "keywords": {"basis": "<>", "cfour_dropmo": 1, "cfour_scf_conv": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "vcc", "cfour_print": 2},                                                 }, id="ccd  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccd",  "reference": "rhf",  "fcae": "fc",                        "keywords": {"basis": "<>", "cfour_dropmo": 1, "cfour_scf_conv": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "ecc"},                                                                   }, id="ccd  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccd",  "reference": "rhf",  "fcae": "fc",                        "keywords": {"basis": "<>", "cfour_dropmo": 1, "cfour_scf_conv": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "ncc"},                                                 "wrong": {1: _w9} }, id="ccd  rhf fc: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccd", "reference": "rhf",  "fcae": "fc",                        "keywords": {"basis": "<>"},                                                                                                                                          "error": {1: _q12}}, id="ccd  rhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccd", "reference": "rhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {"basis": "<>", "qc_module": "tce", "nwchem_tce__freeze": 1},                                                                                                               }, id="ccd  rhf fc: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-ccd",  "reference": "rhf",  "fcae": "fc",                        "keywords": {"basis": "<>", "psi4_freeze_core": True, "psi4_cc_type": "conv"},                                                                                        "error": {1: _q11}}, id="ccd  rhf fc: psi4",       marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-ccd", "reference": "rhf", "fcae": "fc",                       "keywords": {"psi4_cfour_dropmo": [1], "psi4_dertype": "none", "psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "vcc", **_p4c4_fd},                                              }, id="ccd  rhf fc: psi4-cfour-vcc"),

        pytest.param({"call": "c4-ccd",  "reference": "uhf",  "fcae": "ae",                        "keywords": {"basis": "<>", "cfour_reference": "uhf", "cfour_scf_conv": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "vcc", "cfour_print": 2},                                          }, id="ccd  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccd",  "reference": "uhf",  "fcae": "ae",                        "keywords": {"basis": "<>", "cfour_reference": "uhf", "cfour_scf_conv": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "ecc"},                                                            }, id="ccd  uhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccd", "reference": "uhf",  "fcae": "ae",                        "keywords": {"basis": "<>", "gamess_contrl__scftyp": "uhf", "gamess_ccinp__ncore": 0},                                                                                "error": {1: _q2} }, id="ccd  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccd", "reference": "uhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {"basis": "<>", "nwchem_scf__uhf": True, "qc_module": "tce"},                                                                                                               }, id="ccd  uhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-ccd",  "reference": "uhf",  "fcae": "ae",                        "keywords": {"basis": "<>", "reference": "uhf", "psi4_cc_type": "conv"},                                                                                              "error": {1: _q11}}, id="ccd  uhf ae: psi4",       marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-ccd", "reference": "uhf", "fcae": "ae",                       "keywords": {"psi4_dertype": "none", "psi4_cfour_reference": "uhf", "psi4_cfour_cc_program": "vcc", **_p4c4_fd},                                                                        }, id="ccd  uhf ae: psi4-cfour-vcc"),

        pytest.param({"call": "c4-ccd",  "reference": "uhf",  "fcae": "fc",                        "keywords": {"basis": "<>", "cfour_reference": "uhf", "cfour_dropmo": [1], "cfour_scf_conv": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "vcc", "cfour_print": 2},                     }, id="ccd  uhf fc: cfour-vcc", marks=using("cfour")),
        pytest.param({"call": "c4-ccd",  "reference": "uhf",  "fcae": "fc",                        "keywords": {"basis": "<>", "cfour_reference": "uhf", "cfour_dropmo": [1], "cfour_scf_conv": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "ecc"},                                       }, id="ccd  uhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccd", "reference": "uhf",  "fcae": "fc",                        "keywords": {"basis": "<>", "gamess_contrl__scftyp": "uhf"},                                                                                                          "error": {1: _q2} }, id="ccd  uhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccd", "reference": "uhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {"basis": "<>", "nwchem_scf__uhf": True, "nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                                      }, id="ccd  uhf fc: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-ccd",  "reference": "uhf",  "fcae": "fc",                        "keywords": {"basis": "<>", "reference": "uhf", "psi4_cc_type": "conv", "psi4_freeze_core": True},                                                                    "error": {1: _q11}}, id="ccd  uhf fc: psi4",       marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-ccd", "reference": "uhf", "fcae": "fc",                       "keywords": {"psi4_cfour_dropmo": [1], "psi4_dertype": "none", "psi4_cfour_reference": "uhf", "psi4_cfour_cc_program": "vcc", **_p4c4_fd},                                              }, id="ccd  uhf fc: psi4-cfour-vcc"),
        # yapf: enable
    ],
)
def test_ccd_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "gradient"))


#
#   ,-----. ,-----.,------.      ,--.  ,--.                     ,--.
#  '  .--./'  .--./|  .-.  \     |  '--'  | ,---.  ,---.  ,---. `--' ,--,--.,--,--,
#  |  |    |  |    |  |  \  :    |  .--.  || .-. :(  .-' (  .-' ,--.' ,-.  ||      \
#  '  '--'\'  '--'\|  '--'  /    |  |  |  |\   --..-'  `).-'  `)|  |\ '-'  ||  ||  |
#   `-----' `-----'`-------'     `--'  `--' `----'`----' `----' `--' `--`--'`--''--'
#
#  <<<  CCD Hessian

@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(2, id="hes2"),
        # pytest.param(1, id="hes1", marks=pytest.mark.long),
        # pytest.param(0, id="hes0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # ncc errors, gms complains anal/num not set up
        pytest.param({"call": "c4-ccd",  "reference": "rhf",  "fcae": "ae",                        "keywords": {"basis": "<>", "cfour_scf_conv": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "vcc", "cfour_print": 2},                                                                    }, id="ccd  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccd",  "reference": "rhf",  "fcae": "ae",                        "keywords": {"basis": "<>", "cfour_scf_conv": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "ecc"},                                                                                      }, id="ccd  rhf ae: cfour-ecc",  marks=using("cfour")),
       # pytest.param({"call": "c4-ccd",  "reference": "rhf",  "fcae": "ae",                        "keywords": {"basis": "<>", "cfour_scf_conv": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "ncc"},                                                                                      }, id="ccd  rhf ae: cfour-ncc",  marks=using("cfour")),
       # pytest.param({"call": "gms-ccd", "reference": "rhf",  "fcae": "ae",                        "keywords": {"basis": "<>", "gamess_ccinp__ncore": 0},                                                                                                                "error": {2: _q12}}, id="ccd  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccd", "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {"basis": "<>", "qc_module": "tce"},                                                                                                                                        }, id="ccd  rhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-ccd",  "reference": "rhf",  "fcae": "ae",                        "keywords": {"basis": "<>", "psi4_cc_type": "conv"},                                                                                                                  "error": {2: _q17}}, id="ccd  rhf ae: psi4",       marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-ccd", "reference": "rhf", "fcae": "ae",                       "keywords": {"psi4_dertype": "none", "psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "vcc", **_p4c4_fd},                                                                        }, id="ccd  rhf fc: psi4-cfour-vcc"),

        # FC: vcc errors for analytic hess
        # pytest.param({"call": "c4-ccd",  "reference": "rhf",  "fcae": "fc",                        "keywords": {"basis": "<>", "cfour_dropmo": 1, "cfour_scf_conv": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "vcc", "cfour_print": 2},                                                 }, id="ccd  rhf fc: cfour-vcc",  marks=using("cfour")),

        pytest.param({"call": "c4-ccd",  "reference": "uhf",  "fcae": "ae",                        "keywords": {"basis": "<>", "cfour_reference": "uhf", "cfour_scf_conv": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "vcc", "cfour_print": 2},                                          }, id="ccd  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccd",  "reference": "uhf",  "fcae": "ae",                        "keywords": {"basis": "<>", "cfour_reference": "uhf", "cfour_scf_conv": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "ecc"},                                                            }, id="ccd  uhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccd", "reference": "uhf",  "fcae": "ae",                        "keywords": {"basis": "<>", "gamess_contrl__scftyp": "uhf", "gamess_ccinp__ncore": 0},                                                                                "error": {2: _q2} }, id="ccd  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccd", "reference": "uhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {"basis": "<>", "nwchem_scf__uhf": True, "qc_module": "tce"},                                                                                                               }, id="ccd  uhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-ccd",  "reference": "uhf",  "fcae": "ae",                        "keywords": {"basis": "<>", "reference": "uhf", "psi4_cc_type": "conv"},                                                                                              "error": {2: _q17}}, id="ccd  uhf ae: psi4",       marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-ccd", "reference": "uhf", "fcae": "ae",                       "keywords": {"psi4_dertype": "none", "psi4_cfour_reference": "uhf", "psi4_cfour_cc_program": "vcc", **_p4c4_fd},                                                                        }, id="ccd  uhf ae: psi4-cfour-vcc"),

        # FC: vcc errors for analytic hess
        # pytest.param({"call": "c4-ccd",  "reference": "uhf",  "fcae": "fc",                        "keywords": {"basis": "<>", "cfour_reference": "uhf", "cfour_dropmo": [1], "cfour_scf_conv": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "vcc", "cfour_print": 2},                     }, id="ccd  uhf fc: cfour-vcc", marks=using("cfour")),
        # yapf: enable
    ],
)
def test_ccd_hessian_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "hessian"))


#
#   ,-----. ,-----. ,---.  ,------.      ,------.
#  '  .--./'  .--./'   .-' |  .-.  \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |    |  |    `.  `-. |  |  \  :    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '--'\'  '--'\.-'    ||  '--'  /    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----' `-----'`-----' `-------'     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                      `---' `---'
#  <<<  CCSD Energy


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ]
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p"),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "vcc", "cfour_print": 2},                                                                                       }, id="ccsd  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "ecc",},                                                                                                        }, id="ccsd  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "ncc",},                                                                                                        }, id="ccsd  rhf ae: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "rhf",  "fcae": "ae", "keywords": {"gamess_ccinp__ncore": 0},                                                                                                                                                                          }, id="ccsd  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "tce"},                                                                                                                                                                                }, id="ccsd  rhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "cc"},                                                                                                                                                                                 }, id="ccsd  rhf ae: nwchem-cc",  marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "ccenergy"},                                                                                                                                                                           }, id="ccsd  rhf ae: psi4-cc",    marks=using("psi4")),
        pytest.param({"call": "p4-ccsd",  "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "fnocc"},                                                                                                                                                                              }, id="ccsd  rhf ae: psi4-fnocc", marks=using("psi4")),

        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", "cfour_dropmo": [1], "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "vcc", "cfour_print": 2},                                                                  }, id="ccsd  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", "cfour_dropmo": [1], "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "ecc"},                                                                                    }, id="ccsd  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", "cfour_dropmo": [1], "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "ncc"},                                                                                    }, id="ccsd  rhf fc: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "rhf",  "fcae": "fc", "keywords": {},                                                                                                                                                                                                  }, id="ccsd  rhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "fc", "keywords": {"qc_module": "tce", "nwchem_tce__freeze": 1 },                                                                                                                                                      }, id="ccsd  rhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "fc", "keywords": {"nwchem_ccsd__freeze": 1},                                                                                                                                                                          }, id="ccsd  rhf fc: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "rhf",  "fcae": "fc", "keywords": {"psi4_freeze_core": True, "qc_module": "ccenergy"},                                                                                                                                                 }, id="ccsd  rhf fc: psi4-cc",    marks=using("psi4")),
        pytest.param({"call": "p4-ccsd",  "reference": "rhf",  "fcae": "fc", "keywords": {"psi4_freeze_core": True, "qc_module": "fnocc"},                                                                                                                                                    }, id="ccsd  rhf fc: psi4-fnocc", marks=using("psi4")),

        # "cfour_occupation": [[3, 1, 1, 0], [3, 0, 1, 0]]
        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "ae", "keywords": {"cfour_BASIS": "<>", "cfour_REFerence": "UHF", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "vcc", "cfour_print": 2},                                                             }, id="ccsd  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "ae", "keywords": {"cfour_BASIS": "<>", "cfour_REFerence": "UHF", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "ecc"},                                                                               }, id="ccsd  uhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "uhf",  "fcae": "ae", "keywords": {"gamess_contrl__scftyp": "uhf", "gamess_ccinp__ncore": 0},                                                                                                                         "error": {0: _q2}}, id="ccsd  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "uhf",  "fcae": "ae", "keywords": {"qc_module": "tce", "nwchem_scf__uhf": True},                                                                                                                                                       }, id="ccsd  uhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "uhf",  "fcae": "ae", "keywords": {"nwchem_scf__uhf": True},                                                                                                                                                          "error": {0: _q3}}, id="ccsd  uhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "uhf",  "fcae": "ae", "keywords": {"reference": "uhf", "qc_module": "ccenergy"},                                                                                                                                                       }, id="ccsd  uhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "fc", "keywords": {"cfour_BASIS": "<>", "cfour_dropmo": [1], "cfour_REFerence": "UHF", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "vcc", "cfour_print": 2},                                        }, id="ccsd  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "fc", "keywords": {"cfour_BASIS": "<>", "cfour_dropmo": [1], "cfour_REFerence": "UHF", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "ecc"},                                                          }, id="ccsd  uhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "uhf",  "fcae": "fc", "keywords": {"gamess_contrl__scftyp": "uhf"},                                                                                                                                                   "error": {0: _q2}}, id="ccsd  uhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "uhf",  "fcae": "fc", "keywords": {"nwchem_tce__freeze": 1, "qc_module": "tce", "nwchem_scf__uhf": True},                                                                                                                              }, id="ccsd  uhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "uhf",  "fcae": "fc", "keywords": {"nwchem_ccsd__freeze": 1, "nwchem_scf__uhf": True},                                                                                                                                "error": {0: _q3}}, id="ccsd  uhf fc: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "uhf",  "fcae": "fc", "keywords": {"psi4_freeze_core": True, "reference": "uhf", "qc_module": "ccenergy"},                                                                                                                             }, id="ccsd  uhf fc: psi4-cc",    marks=using("psi4")),

        pytest.param({"call": "c4-ccsd",  "reference": "rohf", "fcae": "ae", "keywords": {"cfour_basis": "<>", "cfour_REFerence": "roHF", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "vcc", "cfour_print": 2},                                                            }, id="ccsd rohf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rohf", "fcae": "ae", "keywords": {"cfour_basis": "<>", "cfour_REFerence": "roHF", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "ecc"},                                                                              }, id="ccsd rohf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "rohf", "fcae": "ae", "keywords": {"gamess_contrl__scftyp": "rohf", "gamess_ccinp__ncore": 0, "gamess_ccinp__maxcc": 50},                                                                                                              }, id="ccsd rohf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "rohf", "fcae": "ae", "keywords": {"qc_module": "tce", "nwchem_scf__rohf": True},                                                                                                                                                      }, id="ccsd rohf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "rohf", "fcae": "ae", "keywords": {"nwchem_scf__rohf": True},                                                                                                                                                         "error": {0: _q3}}, id="ccsd rohf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "rohf", "fcae": "ae", "keywords": {"reference": "rohf", "qc_module": "ccenergy"},                                                                                                                                                      }, id="ccsd rohf ae: psi4-cc",    marks=using("psi4")),

        pytest.param({"call": "c4-ccsd",  "reference": "rohf", "fcae": "fc", "keywords": {"cfour_BASIS": "<>", "cfour_dropmo": [1], "cfour_REFerence": "roHF", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_orbitals": 0, "cfour_cc_program": "vcc", "cfour_print": 2},                  }, id="ccsd rohf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rohf", "fcae": "fc", "keywords": {"cfour_BASIS": "<>", "cfour_dropmo": [1], "cfour_REFerence": "roHF", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_orbitals": 0, "cfour_cc_program": "ecc"},                                    }, id="ccsd rohf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "rohf", "fcae": "fc", "keywords": {"gamess_contrl__scftyp": "rohf", "gamess_ccinp__iconv": 9, "gamess_scf__conv": 9},                                                                                                 "wrong": {0: _w2}}, id="ccsd rohf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "rohf", "fcae": "fc", "keywords": {"nwchem_tce__freeze": 1, "qc_module": "tce", "nwchem_scf__rohf": True},                                                                                                                             }, id="ccsd rohf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "rohf", "fcae": "fc", "keywords": {"nwchem_ccsd__freeze": 1, "nwchem_scf__rohf": True},                                                                                                                               "error": {0: _q3}}, id="ccsd rohf fc: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "rohf", "fcae": "fc", "keywords": {"psi4_e_convergence": 8, "psi4_r_convergence": 7, "psi4_freeze_core": True, "reference": "rohf", "qc_module": "ccenergy"},                                                                          }, id="ccsd rohf fc: psi4-cc",    marks=using("psi4")),
        # yapf: enable
    ],
)
def test_ccsd_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "energy"))


# @pytest.mark.parametrize("mode", ["driver", "sandwich"])
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ]
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p"),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable

        ## local
        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "fc", "xptd": {"qc_module": "ecc"},      "keywords": {"cfour_dropmo": 1                                                                    },              }, id="ccsd  rhf fc: cfour dd1",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "fc", "xptd": {"qc_module": "ecc"},      "keywords": {"cfour_dropmo": 1,         "cfour_reference": "uhf"                                  },              }, id="ccsd  uhf fc: cfour dd1",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rohf", "fcae": "fc", "xptd": {"qc_module": "ecc"},      "keywords": {"cfour_dropmo": 1,         "cfour_reference": "rohf",       "cfour_orbitals": 0      },              }, id="ccsd rohf fc: cfour dd1",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "ae", "xptd": {"qc_module": "ecc"},      "keywords": {},                                                                                                   }, id="ccsd  rhf ae: cfour dd1",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "ae", "xptd": {"qc_module": "ecc"},      "keywords": {                           "cfour_reference": "uhf"                                  },              }, id="ccsd  uhf ae: cfour dd1",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rohf", "fcae": "ae", "xptd": {"qc_module": "ecc"},      "keywords": {                           "cfour_reference": "rohf"                                 },              }, id="ccsd rohf ae: cfour dd1",  marks=using("cfour")),

        pytest.param({"call": "gms-ccsd", "reference": "rhf",  "fcae": "fc", "xptd": {"qc_module": None},       "keywords": {},                                                                                                   }, id="ccsd  rhf fc: gamess dd1", marks=using("gamess")),
        pytest.param({"call": "gms-ccsd", "reference": "uhf",  "fcae": "fc", "xptd": {"qc_module": None},       "keywords": {                           "gamess_contrl__scftyp": "uhf"                            }, "error": _q2,}, id="ccsd  uhf fc: gamess dd1", marks=using("gamess")),
        pytest.param({"call": "gms-ccsd", "reference": "rohf", "fcae": "fc", "xptd": {"qc_module": None},       "keywords": {                           "gamess_contrl__scftyp": "rohf", "gamess_ccinp__maxcc": 50}, "wrong": _w2,}, id="ccsd rohf fc: gamess dd1", marks=using("gamess")),
        pytest.param({"call": "gms-ccsd", "reference": "rhf",  "fcae": "ae", "xptd": {"qc_module": None},       "keywords": {"gamess_ccinp__ncore": 0                                                             },              }, id="ccsd  rhf ae: gamess dd1", marks=using("gamess")),
        pytest.param({"call": "gms-ccsd", "reference": "uhf",  "fcae": "ae", "xptd": {"qc_module": None},       "keywords": {"gamess_ccinp__ncore": 0,  "gamess_contrl__scftyp": "uhf"                            }, "error": _q2,}, id="ccsd  uhf ae: gamess dd1", marks=using("gamess")),
        pytest.param({"call": "gms-ccsd", "reference": "rohf", "fcae": "ae", "xptd": {"qc_module": None},       "keywords": {"gamess_ccinp__ncore": 0,  "gamess_contrl__scftyp": "rohf", "gamess_ccinp__maxcc": 50},              }, id="ccsd rohf ae: gamess dd1", marks=using("gamess")),

        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "fc", "xptd": {"qc_module": "cc"},       "keywords": {"nwchem_ccsd__freeze": 1                                                             },              }, id="ccsd  rhf fc: nwchem dd1", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "uhf",  "fcae": "fc", "xptd": {"qc_module": "tce"},      "keywords": {"nwchem_tce__freeze": 1,   "nwchem_scf__uhf": True,         "qc_module": "tce"       },              }, id="ccsd  uhf fc: nwchem dd1", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "rohf", "fcae": "fc", "xptd": {"qc_module": "tce"},      "keywords": {"nwchem_tce__freeze": 1,   "nwchem_scf__rohf": True,        "qc_module": "tce"       },              }, id="ccsd rohf fc: nwchem dd1", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "ae", "xptd": {"qc_module": "cc"},       "keywords": {},                                                                                                   }, id="ccsd  rhf ae: nwchem dd1", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "uhf",  "fcae": "ae", "xptd": {"qc_module": "tce"},      "keywords": {                           "nwchem_scf__uhf": True,         "qc_module": "tce"       },              }, id="ccsd  uhf ae: nwchem dd1", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "rohf", "fcae": "ae", "xptd": {"qc_module": "tce"},      "keywords": {                           "nwchem_scf__rohf": True,        "qc_module": "tce"       },              }, id="ccsd rohf ae: nwchem dd1", marks=using("nwchem")),

        pytest.param({"call": "p4-ccsd",  "reference": "rhf",  "fcae": "fc", "xptd": {"qc_module": "ccenergy"}, "keywords": {"psi4_freeze_core": True                                                             },              }, id="ccsd  rhf fc: psi4 dd1",   marks=using("psi4")),
        pytest.param({"call": "p4-ccsd",  "reference": "uhf",  "fcae": "fc", "xptd": {"qc_module": "ccenergy"}, "keywords": {"psi4_freeze_core": True,  "psi4_reference": "uhf"                                   },              }, id="ccsd  uhf fc: psi4 dd1",   marks=using("psi4")),
        pytest.param({"call": "p4-ccsd",  "reference": "rohf", "fcae": "fc", "xptd": {"qc_module": "ccenergy"}, "keywords": {"psi4_freeze_core": True,  "psi4_reference": "rohf",                                 },              }, id="ccsd rohf fc: psi4 dd1",   marks=using("psi4")),
        pytest.param({"call": "p4-ccsd",  "reference": "rhf",  "fcae": "ae", "xptd": {"qc_module": "ccenergy"}, "keywords": {},                                                                                                   }, id="ccsd  rhf ae: psi4 dd1",   marks=using("psi4")),
        pytest.param({"call": "p4-ccsd",  "reference": "uhf",  "fcae": "ae", "xptd": {"qc_module": "ccenergy"}, "keywords": {                           "psi4_reference": "uhf"                                   },              }, id="ccsd  uhf ae: psi4 dd1",   marks=using("psi4")),
        pytest.param({"call": "p4-ccsd",  "reference": "rohf", "fcae": "ae", "xptd": {"qc_module": "ccenergy"}, "keywords": {                           "psi4_reference": "rohf"                                  },              }, id="ccsd rohf ae: psi4 dd1",   marks=using("psi4")),

        ## translated
        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "fc", "xptd": {"qc_module": "ecc"},      "keywords": {"freeze_core": True,       "reference": "rhf"                                        },              }, id="ccsd  rhf fc: cfour dd2",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "fc", "xptd": {"qc_module": "ecc"},      "keywords": {"freeze_core": True,       "reference": "uhf"                                        },              }, id="ccsd  uhf fc: cfour dd2",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rohf", "fcae": "fc", "xptd": {"qc_module": "ecc"},      "keywords": {"freeze_core": True,       "reference": "rohf",             "cfour_orbitals": 0      },              }, id="ccsd rohf fc: cfour dd2",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "ae", "xptd": {"qc_module": "ecc"},      "keywords": {"freeze_core": False,      "reference": "rhf"                                        },              }, id="ccsd  rhf ae: cfour dd2",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "ae", "xptd": {"qc_module": "ecc"},      "keywords": {"freeze_core": False,      "reference": "uhf"                                        },              }, id="ccsd  uhf ae: cfour dd2",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rohf", "fcae": "ae", "xptd": {"qc_module": "ecc"},      "keywords": {"freeze_core": False,      "reference": "rohf"                                       },              }, id="ccsd rohf ae: cfour dd2",  marks=using("cfour")),

        pytest.param({"call": "gms-ccsd", "reference": "rhf",  "fcae": "fc", "xptd": {"qc_module": None},       "keywords": {"freeze_core": True,       "reference": "rhf"                                        },              }, id="ccsd  rhf fc: gamess dd2", marks=using("gamess")),
        pytest.param({"call": "gms-ccsd", "reference": "uhf",  "fcae": "fc", "xptd": {"qc_module": None},       "keywords": {"freeze_core": True,       "reference": "uhf",                                       }, "error": _q2,}, id="ccsd  uhf fc: gamess dd2", marks=using("gamess")),
        pytest.param({"call": "gms-ccsd", "reference": "rohf", "fcae": "fc", "xptd": {"qc_module": None},       "keywords": {"freeze_core": True,       "reference": "rohf",             "gamess_ccinp__maxcc": 50}, "wrong": _w2,}, id="ccsd rohf fc: gamess dd2", marks=using("gamess")),
        pytest.param({"call": "gms-ccsd", "reference": "rhf",  "fcae": "ae", "xptd": {"qc_module": None},       "keywords": {"freeze_core": False,      "reference": "rhf"                                        },              }, id="ccsd  rhf ae: gamess dd2", marks=using("gamess")),
        pytest.param({"call": "gms-ccsd", "reference": "uhf",  "fcae": "ae", "xptd": {"qc_module": None},       "keywords": {"freeze_core": False,      "reference": "uhf"                                        }, "error": _q2,}, id="ccsd  uhf ae: gamess dd2", marks=using("gamess")),
        pytest.param({"call": "gms-ccsd", "reference": "rohf", "fcae": "ae", "xptd": {"qc_module": None},       "keywords": {"freeze_core": False,      "reference": "rohf",             "gamess_ccinp__maxcc": 50},              }, id="ccsd rohf ae: gamess dd2", marks=using("gamess")),

        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "fc", "xptd": {"qc_module": "cc"},       "keywords": {"freeze_core": True,       "reference": "rhf"                                        },              }, id="ccsd  rhf fc: nwchem dd2", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "uhf",  "fcae": "fc", "xptd": {"qc_module": "tce"},      "keywords": {"freeze_core": True,       "reference": "uhf",              "qc_module": "tce"       },              }, id="ccsd  uhf fc: nwchem dd2", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "rohf", "fcae": "fc", "xptd": {"qc_module": "tce"},      "keywords": {"freeze_core": True,       "reference": "rohf",             "qc_module": "tce"       },              }, id="ccsd rohf fc: nwchem dd2", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "ae", "xptd": {"qc_module": "cc"},       "keywords": {"freeze_core": False,      "reference": "rhf"                                        },              }, id="ccsd  rhf ae: nwchem dd2", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "uhf",  "fcae": "ae", "xptd": {"qc_module": "tce"},      "keywords": {"freeze_core": False,      "reference": "uhf",              "qc_module": "tce"       },              }, id="ccsd  uhf ae: nwchem dd2", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "rohf", "fcae": "ae", "xptd": {"qc_module": "tce"},      "keywords": {"freeze_core": False,      "reference": "rohf",             "qc_module": "tce"       },              }, id="ccsd rohf ae: nwchem dd2", marks=using("nwchem")),

        pytest.param({"call": "p4-ccsd",  "reference": "rhf",  "fcae": "fc", "xptd": {"qc_module": "ccenergy"}, "keywords": {"freeze_core": True,       "reference": "rhf"                                        },              }, id="ccsd  rhf fc: psi4 dd2",   marks=using("psi4")),
        pytest.param({"call": "p4-ccsd",  "reference": "uhf",  "fcae": "fc", "xptd": {"qc_module": "ccenergy"}, "keywords": {"freeze_core": True,       "reference": "uhf"                                        },              }, id="ccsd  uhf fc: psi4 dd2",   marks=using("psi4")),
        pytest.param({"call": "p4-ccsd",  "reference": "rohf", "fcae": "fc", "xptd": {"qc_module": "ccenergy"}, "keywords": {"freeze_core": True,       "reference": "rohf"                                       },              }, id="ccsd rohf fc: psi4 dd2",   marks=using("psi4")),
        pytest.param({"call": "p4-ccsd",  "reference": "rhf",  "fcae": "ae", "xptd": {"qc_module": "ccenergy"}, "keywords": {"freeze_core": False,      "reference": "rhf"                                        },              }, id="ccsd  rhf ae: psi4 dd2",   marks=using("psi4")),
        pytest.param({"call": "p4-ccsd",  "reference": "uhf",  "fcae": "ae", "xptd": {"qc_module": "ccenergy"}, "keywords": {"freeze_core": False,      "reference": "uhf"                                        },              }, id="ccsd  uhf ae: psi4 dd2",   marks=using("psi4")),
        pytest.param({"call": "p4-ccsd",  "reference": "rohf", "fcae": "ae", "xptd": {"qc_module": "ccenergy"}, "keywords": {"freeze_core": False,      "reference": "rohf"                                       },              }, id="ccsd rohf ae: psi4 dd2",   marks=using("psi4")),
        # yapf: enable
    ],
)
def test_ccsd_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request):
    qcprog, method = inp["call"].split("-", 1)
    qcprog = _trans_qcprog[qcprog.lower()]
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["reference"])]]

    inpcopy = {k: v for k, v in inp.items()}
    inpcopy["keywords"] = {
        k: (_trans_key(qcprog, basis, k) if v == "<>" else v) for k, v in inpcopy["keywords"].items()
    }

    inpcopy["driver"] = "energy"
    if not any([k.lower() in _basis_keywords for k in inpcopy["keywords"]]):
        inpcopy["keywords"]["basis"] = basis
    inpcopy["scf_type"] = "pk"
    inpcopy["corl_type"] = "conv"
    inpcopy["qc_module"] = "-".join(
        [qcprog, inp["keywords"].get("qc_module", inp["keywords"].get("cfour_cc_program", ""))]
    ).strip("-")
    print("INP", inpcopy)

    runner_asserter(inpcopy, subject, method, basis, tnm)


#
#   ,-----. ,-----. ,---.  ,------.       ,----.                     ,--.,--.                 ,--.
#  '  .--./'  .--./'   .-' |  .-.  \     '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  |    |  |    `.  `-. |  |  \  :    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  '  '--'\'  '--'\.-'    ||  '--'  /    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#   `-----' `-----'`-----' `-------'      `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  CCSD Gradient


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        # pytest.param(0, id="grd0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "ae",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "vcc", "cfour_print": 2},                                                            }, id="ccsd  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "ae",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ecc",},                                                                             }, id="ccsd  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "ae",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ncc",},                                                                             }, id="ccsd  rhf ae: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {"basis": "<>", "gamess_ccinp__ncore": 0, "gamess_contrl__numgrd": True,},                                                                  }, id="ccsd  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {"basis": "<>", "qc_module": "tce"},                                                                                                        }, id="ccsd  rhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {"basis": "<>"},                                                                                                                            }, id="ccsd  rhf ae: nwchem-cc",  marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "rhf",  "fcae": "ae",                        "keywords": {"basis": "<>", "psi4_cc_type": "conv", "qc_module": "ccenergy"},                                                                           }, id="ccsd  rhf ae: psi4-cc",    marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-ccsd", "reference": "rhf", "fcae": "ae", "keywords": {"psi4_dertype": "none", "psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "vcc", **_p4c4_fd},                                                      }, id="ccsd  rhf ae: psi4-cfour-vcc"),
        # DEBUG pytest.param({"call": "p4-c4-ccsd", "reference": "rhf", "fcae": "ae", "keywords": {"psi4_dertype": "none", "psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "ncc", **_p4c4_fd},                                                      }, id="ccsd  rhf ae: psi4-cfour-ncc"),

        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "fc",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "vcc", "cfour_dropmo": 1, "cfour_print": 2,},                                        }, id="ccsd  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "fc",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ecc", "cfour_dropmo": 1,},                                                          }, id="ccsd  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "fc",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ncc", "cfour_dropmo": 1,},                                        "wrong": {1: _w8} }, id="ccsd  rhf fc: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "rhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {"basis": "<>", "gamess_contrl__numgrd": True,},                                                                                            }, id="ccsd  rhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {"basis": "<>", "nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                               }, id="ccsd  rhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {"basis": "<>", "nwchem_ccsd__freeze": 1},                                                                                                  }, id="ccsd  rhf fc: nwchem-cc",  marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "rhf",  "fcae": "fc",                        "keywords": {"basis": "<>", "psi4_cc_type": "conv", "psi4_freeze_core": True, "qc_module": "ccenergy"},                               "error": {1: _q18}}, id="ccsd  rhf fc: psi4-cc",    marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-ccsd", "reference": "rhf", "fcae": "fc", "keywords": {"psi4_cfour_dropmo": [1], "psi4_dertype": "none", "psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "vcc", **_p4c4_fd},                            }, id="ccsd  rhf fc: psi4-cfour-vcc"),
        # DEBUG pytest.param({"call": "p4-c4-ccsd", "reference": "rhf", "fcae": "fc", "keywords": {"psi4_cfour_dropmo": [1], "psi4_dertype": "none", "psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "ncc", **_p4c4_fd},                            }, id="ccsd  rhf fc: psi4-cfour-ncc"),

        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "ae",                        "keywords": {"cfour_basis": "<>", "cfour_reference": "uhf", **_c4_tight, "cfour_cc_program": "vcc", "cfour_print": 2},                                  }, id="ccsd  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "ae",                        "keywords": {"cfour_basis": "<>", "cfour_reference": "uhf", **_c4_tight, "cfour_cc_program": "ecc",},                                                   }, id="ccsd  uhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "uhf",  "fcae": "ae",                        "keywords": {"basis": "<>", "gamess_ccinp__ncore": 0, "gamess_contrl__scftyp": "uhf", "gamess_contrl__numgrd": True,},                "error": {1: _q2} }, id="ccsd  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "uhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {"basis": "<>", "nwchem_scf__uhf": True, "qc_module": "tce"},                                                                               }, id="ccsd  uhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "uhf",  "fcae": "ae",                        "keywords": {"basis": "<>", "nwchem_scf__uhf": True},                                                                                 "error": {1: _q3} }, id="ccsd  uhf ae: nwchem-cc",  marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "uhf",  "fcae": "ae",                        "keywords": {"basis": "<>", "reference": "uhf", "psi4_cc_type": "conv"},                                                                                }, id="ccsd  uhf ae: psi4-cc",    marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-ccsd", "reference": "uhf", "fcae": "ae", "keywords": {"psi4_dertype": "none", "psi4_cfour_reference": "uhf", "psi4_cfour_cc_program": "vcc", **_p4c4_fd},                                                      }, id="ccsd  uhf ae: psi4-cfour-vcc"),

        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "fc",                        "keywords": {"cfour_basis": "<>", "cfour_reference": "uhf", **_c4_tight, "cfour_cc_program": "vcc", "cfour_dropmo": 1, "cfour_print": 2},               }, id="ccsd  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "fc",                        "keywords": {"cfour_basis": "<>", "cfour_reference": "uhf", **_c4_tight, "cfour_cc_program": "ecc", "cfour_dropmo": 1,},                                }, id="ccsd  uhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "uhf",  "fcae": "fc",                        "keywords": {"basis": "<>", "gamess_contrl__scftyp": "uhf", "gamess_contrl__numgrd": True,},                                          "error": {1: _q2} }, id="ccsd  uhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "uhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {"basis": "<>", "nwchem_scf__uhf": True, "nwchem_tce__freeze": 1, "qc_module": "tce"},                                                      }, id="ccsd  uhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "uhf",  "fcae": "fc",                        "keywords": {"basis": "<>", "nwchem_scf__uhf": True, "nwchem_ccsd__freeze": 1},                                                       "error": {1: _q3} }, id="ccsd  uhf fc: nwchem-cc",  marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "uhf",  "fcae": "fc",                        "keywords": {"basis": "<>", "reference": "uhf", "psi4_cc_type": "conv", "psi4_freeze_core": True},                                    "error": {1: _q18}}, id="ccsd  uhf fc: psi4",       marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-ccsd", "reference": "uhf", "fcae": "fc", "keywords": {"psi4_cfour_dropmo": [1], "psi4_dertype": "none", "psi4_cfour_reference": "uhf", "psi4_cfour_cc_program": "vcc", **_p4c4_fd},                            }, id="ccsd  uhf fc: psi4-cfour-vcc"),

        pytest.param({"call": "c4-ccsd",  "reference": "rohf", "fcae": "ae",                        "keywords": {"cfour_basis": "<>", "cfour_reference": "rohf", "cfour_cc_program": "vcc", **_c4_tight, "cfour_print": 2},                                 }, id="ccsd rohf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rohf", "fcae": "ae",                        "keywords": {"cfour_basis": "<>", "cfour_reference": "rohf", "cfour_cc_program": "ecc", **_c4_tight, "cfour_print": 2},                                 }, id="ccsd rohf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "rohf", "fcae": "ae", "xptd": {"fd": True},  "keywords": {"gamess_contrl__scftyp": "rohf", "gamess_ccinp__ncore": 0, "gamess_ccinp__maxcc": 50, "gamess_contrl__numgrd": True},                      }, id="ccsd rohf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "rohf", "fcae": "ae", "xptd": {"fd": True},  "keywords": {"qc_module": "tce", "nwchem_scf__rohf": True},                                                                                             }, id="ccsd rohf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "rohf", "fcae": "ae",                        "keywords": {"nwchem_scf__rohf": True},                                                                                               "error": {1: _q3} }, id="ccsd rohf ae: nwchem-cc",  marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "rohf", "fcae": "ae",                        "keywords": {"reference": "rohf", "qc_module": "ccenergy"},                                                                                             }, id="ccsd rohf ae: psi4-cc",    marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-ccsd", "reference": "rohf", "fcae": "ae", "keywords": {"psi4_dertype": "none", "psi4_cfour_reference": "rohf", "psi4_cfour_cc_program": "vcc", **_p4c4_fd},                                                    }, id="ccsd rohf ae: psi4-cfour-vcc"),

        # * vcc and ecc yield correct gradients w/pert_orb=0 but not at the same time with correct energies (orbitals=0)
        # * note that ecc not recc for rohf and no rohf gradients on list in cfour paper
        pytest.param({"call": "c4-ccsd",  "reference": "rohf", "fcae": "fc",                        "keywords": {"cfour_basis": "<>", "cfour_reference": "rohf", "cfour_cc_program": "vcc", "cfour_dropmo": 1, **_c4_tight, "cfour_print": 2, "cfour_pert_orb": 0}, "wrong": {1: _w12}}, id="ccsd rohf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rohf", "fcae": "fc",                        "keywords": {"cfour_basis": "<>", "cfour_reference": "rohf", "cfour_cc_program": "ecc", "cfour_dropmo": 1, **_c4_tight, "cfour_print": 2, "cfour_pert_orb": 0}, "wrong": {1: _w12}}, id="ccsd rohf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "rohf", "fcae": "fc",                        "keywords": {"gamess_contrl__scftyp": "rohf", "gamess_ccinp__iconv": 9, "gamess_scf__conv": 9, "gamess_contrl__numgrd": True},        "wrong": {1: _w2} }, id="ccsd rohf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "rohf", "fcae": "fc", "xptd": {"fd": True},  "keywords": {"nwchem_tce__freeze": 1, "qc_module": "tce", "nwchem_scf__rohf": True},                                                                    }, id="ccsd rohf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "rohf", "fcae": "fc",                        "keywords": {"nwchem_ccsd__freeze": 1, "nwchem_scf__rohf": True},                                                                     "error": {1: _q3} }, id="ccsd rohf fc: nwchem",     marks=using("nwchem")),
        # DEBUG pytest.param({"call": "p4-c4-ccsd", "reference": "rohf", "fcae": "fc", "keywords": {"psi4_cfour_dropmo": [1], "psi4_dertype": "none", "psi4_cfour_reference": "rohf", "psi4_cfour_cc_program": "vcc", "cfour_print": 2, **_p4c4_fd},        }, id="ccsd rohf fc: psi4-cfour-vcc"),
        # yapf: enable
    ],
)
def test_ccsd_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "gradient"))


# PERT_ORB
# specifies the type of perturbed orbitals used in energy derivative calculations. STANDARD means that the gradient formulation assumes that the perturbed orbitals are not those in which the (perturbed) Fock matrix is diagonal. CANONICAL means that the perturbed orbitals are assumed to be canonical. This keyword is set automatically to CANONICAL in derivative calculations with methods which include triple excitations (MBPT[4]/MP4, CCSD+T[CCSD], CCSD[T], QCISD[T] and all iterative schemes like CCSDT-n and CC3) apart from CCSDT. IJ_CANONICAL requests a canonical perturbed-orbital treatment only for the occupied-occupied block of the unperturbed density matrix in analytic derivative calculations.
# For testing purpose, it is possible to force the use standard perturbed orbitals even in case of iterative triple excitations via the option FORCE_STANDA.
# Note also that in case of unrelaxed derivatives standard orbitals must be used.
# Default: STANDARD for all methods without triples (except CCSDT), CANONICAL for all methods with triples in case of relaxed derivatives.


#
#   ,-----. ,-----. ,---.  ,------.      ,--.  ,--.                     ,--.
#  '  .--./'  .--./'   .-' |  .-.  \     |  '--'  | ,---.  ,---.  ,---. `--' ,--,--.,--,--,
#  |  |    |  |    `.  `-. |  |  \  :    |  .--.  || .-. :(  .-' (  .-' ,--.' ,-.  ||      \
#  '  '--'\'  '--'\.-'    ||  '--'  /    |  |  |  |\   --..-'  `).-'  `)|  |\ '-'  ||  ||  |
#   `-----' `-----'`-----' `-------'     `--'  `--' `----'`----' `----' `--' `--`--'`--''--'
#
#  <<<  CCSD Hessian

@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(2, id="hes2"),
        # pytest.param(1, id="hes1", marks=pytest.mark.long),
        # pytest.param(0, id="hes0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # * ncc errors 
        # * nwchem-cc Error in Communication
        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "ae",                        "keywords": {"basis": "<>", **_c4_tight, "cfour_cc_program": "vcc", "cfour_print": 2},                                                                    }, id="ccsd  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "ae",                        "keywords": {"basis": "<>", **_c4_tight, "cfour_cc_program": "ecc", "cfour_print": 2},                                                                    }, id="ccsd  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {"basis": "<>", "gamess_ccinp__ncore": 0, "gamess_force__method": "fullnum"},                                                                 }, id="ccsd  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {"basis": "<>", "qc_module": "tce"},                                                                                                          }, id="ccsd  rhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "rhf",  "fcae": "ae", "xptd": {"fd": False}, "keywords": {"basis": "<>", "psi4_cc_type": "conv", "psi4_points": 5, "psi4_fd_project": False},                                                          }, id="ccsd  rhf ae: psi4",       marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-ccsd", "reference": "rhf", "fcae": "ae", "keywords": {"psi4_dertype": "none", "psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "vcc", **_p4c4_fd},                                                        }, id="ccsd  rhf ae: psi4-cfour-vcc"),

        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "fc",                        "keywords": {"basis": "<>", **_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "vcc", "cfour_print": 2},                                                 }, id="ccsd  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "fc",                        "keywords": {"basis": "<>", **_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "ecc", "cfour_print": 2},                                                 }, id="ccsd  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "rhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {"basis": "<>", "gamess_force__method": "fullnum"},                                                                                           }, id="ccsd  rhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {"basis": "<>", "nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                                 }, id="ccsd  rhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "rhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {"basis": "<>", "psi4_cc_type": "conv", "psi4_freeze_core": True, "psi4_points": 5, "psi4_fd_project": False, "psi4_dertype": "none"},        }, id="ccsd  rhf fc: psi4",       marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-ccsd", "reference": "rhf", "fcae": "fc", "keywords": {"psi4_dertype": "none", "psi4_cfour_dropmo": [1], "psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "vcc", **_p4c4_fd},                              }, id="ccsd  rhf fc: psi4-cfour-vcc"),

        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "ae",                        "keywords": {"basis": "<>", **_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "vcc", "cfour_print": 2},                                          }, id="ccsd  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "ae",                        "keywords": {"basis": "<>", **_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "ecc"},                                                            }, id="ccsd  uhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "uhf",  "fcae": "ae",                        "keywords": {"basis": "<>", "gamess_contrl__scftyp": "uhf", "gamess_ccinp__ncore": 0},                                                  "error": {2: _q2} }, id="ccsd  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "uhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {"basis": "<>", "nwchem_scf__uhf": True, "qc_module": "tce"},                                                                                 }, id="ccsd  uhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "uhf",  "fcae": "ae",                        "keywords": {"basis": "<>", "reference": "uhf", "psi4_cc_type": "conv", "psi4_points": 5, "psi4_fd_project": False},                                      }, id="ccsd  uhf ae: psi4",       marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-ccsd", "reference": "uhf", "fcae": "ae", "keywords": {"psi4_dertype": "none", "psi4_cfour_reference": "uhf", "psi4_cfour_cc_program": "vcc", **_p4c4_fd},                                                        }, id="ccsd  uhf ae: psi4-cfour-vcc"),

        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "fc",                        "keywords": {"basis": "<>", **_c4_tight, "cfour_reference": "uhf", "cfour_dropmo": 1, "cfour_cc_program": "vcc", "cfour_print": 2},                       }, id="ccsd  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "fc",                        "keywords": {"basis": "<>", **_c4_tight, "cfour_reference": "uhf", "cfour_dropmo": 1, "cfour_cc_program": "ecc", "cfour_print": 2},                       }, id="ccsd  uhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "uhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {"basis": "<>", "gamess_contrl__scftyp": "uhf", "gamess_force__method": "fullnum"},                                         "error": {2: _q2} }, id="ccsd  uhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "uhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {"basis": "<>", "nwchem_scf__uhf": True, "nwchem_tce__freeze": 1, "qc_module": "tce"},                                                        }, id="ccsd  uhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "uhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {"basis": "<>", "psi4_reference": "uhf", "psi4_cc_type": "conv", "psi4_freeze_core": True, "psi4_points": 5, "psi4_fd_project": False, "psi4_dertype": "none"}}, id="ccsd  uhf fc: psi4", marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-ccsd", "reference": "uhf", "fcae": "fc", "keywords": {"psi4_dertype": "none", "psi4_cfour_dropmo": [1], "psi4_cfour_reference": "uhf", "psi4_cfour_cc_program": "vcc", **_p4c4_fd},                              }, id="ccsd  uhf fc: psi4-cfour-vcc"),
        # yapf: enable
    ],
)
def test_ccsd_hessian_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "hessian"))


#                                        ,--.                                                                                                          
#   ,-----. ,-----. ,---.  ,------.      |  |    ,--------. ,-. ,-----. ,-----. ,---.  ,------. ,-.      ,------.                                      
#  '  .--./'  .--./'   .-' |  .-.  \ ,---    ---.'--.  .--'/ .''  .--./'  .--./'   .-' |  .-.  \'. \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--. 
#  |  |    |  |    `.  `-. |  |  \  :'---    ---'   |  |  |  | |  |    |  |    `.  `-. |  |  \  :|  |    |  `--, |      \| .-. :|  .--'| .-. |\  '  /  
#  '  '--'\'  '--'\.-'    ||  '--'  /    |  |       |  |  |  | '  '--'\'  '--'\.-'    ||  '--'  /|  |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '   
#   `-----' `-----'`-----' `-------'     `--'       `--'   \ '. `-----' `-----'`-----' `-------'.' /     `------'`--''--' `----'`--'   .`-  /.-'  /    
#                                                           `-'                                 `-'                                    `---' `---'     
#  <<<  CCSD+T(CCSD) Energy aka CCSD[T] Energy

@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ]
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p"),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "c4-ccsd+t(ccsd)",  "reference": "rhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "vcc", "cfour_print": 2},                                                                 }, id="ccsd_t_ccsd_  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd+t(ccsd)",  "reference": "rhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ecc",},                                                                                  }, id="ccsd_t_ccsd_  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd+t(ccsd)",  "reference": "rhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ncc",},                                                                                  }, id="ccsd_t_ccsd_  rhf ae: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd+t(ccsd)", "reference": "rhf",  "fcae": "ae", "keywords": {"gamess_ccinp__ncore": 0},                                                                                                                      }, id="ccsd_t_ccsd_  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd+t(ccsd)", "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "tce"},                                                                                                                            }, id="ccsd_t_ccsd_  rhf ae: nwchem-tce",  marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd+t(ccsd)", "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "cc"},                                                                                                                             }, id="ccsd_t_ccsd_  rhf ae: nwchem-cc",  marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd+t(ccsd)",  "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                            "error": {0: _q20}}, id="ccsd_t_ccsd_  rhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-ccsd+t(ccsd)",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_print": 2},                                            }, id="ccsd_t_ccsd_  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd+t(ccsd)",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_dropmo": [1], "cfour_cc_program": "ecc"},                                                              }, id="ccsd_t_ccsd_  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd+t(ccsd)",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_dropmo": [1], "cfour_cc_program": "ncc"},                                                              }, id="ccsd_t_ccsd_  rhf fc: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd+t(ccsd)", "reference": "rhf",  "fcae": "fc", "keywords": {},                                                                                                                                              }, id="ccsd_t_ccsd_  rhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd+t(ccsd)", "reference": "rhf",  "fcae": "fc", "keywords": {"qc_module": "tce", "nwchem_tce__freeze": 1 },                                                                                                  }, id="ccsd_t_ccsd_  rhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd+t(ccsd)", "reference": "rhf",  "fcae": "fc", "keywords": {"nwchem_ccsd__freeze": 1, "qc_module": "cc"},                                                                                                   }, id="ccsd_t_ccsd_  rhf fc: nwchem-cc",  marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd+t(ccsd)",  "reference": "rhf",  "fcae": "fc", "keywords": {"psi4_freeze_core": True},                                                                                                    "error": {0: _q20}}, id="ccsd_t_ccsd_  rhf fc: psi4-cc",    marks=using("psi4")),

        pytest.param({"call": "c4-ccsd+t(ccsd)",  "reference": "uhf",  "fcae": "ae", "keywords": {"cfour_BASIS": "<>", **_c4_tight, "cfour_REFerence": "UHF", "cfour_cc_program": "vcc", "cfour_print": 2},                                       }, id="ccsd_t_ccsd_  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd+t(ccsd)",  "reference": "uhf",  "fcae": "ae", "keywords": {"cfour_BASIS": "<>", **_c4_tight, "cfour_REFerence": "UHF", "cfour_cc_program": "ecc"},                                                         }, id="ccsd_t_ccsd_  uhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd+t(ccsd)", "reference": "uhf",  "fcae": "ae", "keywords": {"gamess_contrl__scftyp": "uhf", "gamess_ccinp__ncore": 0},                                                                    "error": {0: _q2} }, id="ccsd_t_ccsd_  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd+t(ccsd)", "reference": "uhf",  "fcae": "ae", "keywords": {"qc_module": "tce", "nwchem_scf__uhf": True},                                                                                                   }, id="ccsd_t_ccsd_  uhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd+t(ccsd)", "reference": "uhf",  "fcae": "ae", "keywords": {"nwchem_scf__uhf": True},                                                                                                     "error": {0: _q3} }, id="ccsd_t_ccsd_  uhf ae: nwchem-cc",  marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd+t(ccsd)",  "reference": "uhf",  "fcae": "ae", "keywords": {"reference": "uhf", "qc_module": "ccenergy"},                                                                                 "error": {0: _q20}}, id="ccsd_t_ccsd_  uhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-ccsd+t(ccsd)",  "reference": "uhf",  "fcae": "fc", "keywords": {"cfour_BASIS": "<>", **_c4_tight, "cfour_dropmo": [1], "cfour_REFerence": "UHF", "cfour_cc_program": "vcc", "cfour_print": 2},                  }, id="ccsd_t_ccsd_  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd+t(ccsd)",  "reference": "uhf",  "fcae": "fc", "keywords": {"cfour_BASIS": "<>", **_c4_tight, "cfour_dropmo": [1], "cfour_REFerence": "UHF", "cfour_cc_program": "ecc"},                                    }, id="ccsd_t_ccsd_  uhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd+t(ccsd)", "reference": "uhf",  "fcae": "fc", "keywords": {"gamess_contrl__scftyp": "uhf"},                                                                                              "error": {0: _q2} }, id="ccsd_t_ccsd_  uhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd+t(ccsd)", "reference": "uhf",  "fcae": "fc", "keywords": {"nwchem_tce__freeze": 1, "qc_module": "tce", "nwchem_scf__uhf": True},                                                                          }, id="ccsd_t_ccsd_  uhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd+t(ccsd)", "reference": "uhf",  "fcae": "fc", "keywords": {"nwchem_ccsd__freeze": 1, "nwchem_scf__uhf": True},                                                                           "error": {0: _q3} }, id="ccsd_t_ccsd_  uhf fc: nwchem-cc",  marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd+t(ccsd)",  "reference": "uhf",  "fcae": "fc", "keywords": {"psi4_freeze_core": True, "reference": "uhf"},                                                                                "error": {0: _q20}}, id="ccsd_t_ccsd_  uhf fc: psi4-cc",    marks=using("psi4")),
        # yapf: enable
    ],
)
def test_ccsdpt_prccsd_pr_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "energy"))



#                                                                                                                                 
#   ,-----. ,-----. ,---.  ,------. ,--------.        ,--.            ,------.                                                    
#  '  .--./'  .--./'   .-' |  .-.  \'--.  .--',-----./   | ,--,--.    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.               
#  |  |    |  |    `.  `-. |  |  \  :  |  |   '-----'`|  |' ,-.  |    |  `--, |      \| .-. :|  .--'| .-. |\  '  /                
#  '  '--'\'  '--'\.-'    ||  '--'  /  |  |           |  |\ '-'  |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '                 
#   `-----' `-----'`-----' `-------'   `--'           `--' `--`--'    `------'`--''--' `----'`--'   .`-  /.-'  /                  
#                                                                                                   `---' `---'                   
#  <<<  CCSDT-1a Energy

@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ]
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "c4-ccsdt-1a",  "reference": "rhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "vcc",},                                                                                                               }, id="ccsdt-1a  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-1a",  "reference": "rhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ecc",},                                                                                                               }, id="ccsdt-1a  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-1a",  "reference": "rhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ncc",},                                                                                                               }, id="ccsdt-1a  rhf ae: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsdt-1a", "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                                                         "error": {0: _q24}}, id="ccsdt-1a  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsdt-1a", "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                                                         "error": {0: _q24}}, id="ccsdt-1a  rhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-ccsdt-1a",  "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                                                         "error": {0: _q24}}, id="ccsdt-1a  rhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-ccsdt-1a",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "vcc", "cfour_dropmo": 1,},                                                                                            }, id="ccsdt-1a  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-1a",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ecc", "cfour_dropmo": 1,},                                                                                            }, id="ccsdt-1a  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-1a",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ncc", "cfour_dropmo": 1,},                                                                                            }, id="ccsdt-1a  rhf fc: cfour-ncc",  marks=using("cfour")),

        pytest.param({"call": "c4-ccsdt-1a",  "reference": "uhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "vcc",},                                                                                     }, id="ccsdt-1a  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-1a",  "reference": "uhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "ecc",},                                                                                     }, id="ccsdt-1a  uhf ae: cfour-ecc",  marks=using("cfour")),

        pytest.param({"call": "c4-ccsdt-1a",  "reference": "uhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "vcc", "cfour_dropmo": 1},                                                                   }, id="ccsdt-1a  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-1a",  "reference": "uhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "ecc", "cfour_dropmo": 1},                                                                   }, id="ccsdt-1a  uhf fc: cfour-ecc",  marks=using("cfour")),
        # yapf: enable
    ],
)
def test_ccsdt1a_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "energy"))


#
#   ,-----. ,-----. ,---.  ,------. ,--------.        ,--.             ,----.                     ,--.,--.                 ,--.   
#  '  .--./'  .--./'   .-' |  .-.  \'--.  .--',-----./   | ,--,--.    '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-. 
#  |  |    |  |    `.  `-. |  |  \  :  |  |   '-----'`|  |' ,-.  |    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-' 
#  '  '--'\'  '--'\.-'    ||  '--'  /  |  |           |  |\ '-'  |    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |   
#   `-----' `-----'`-----' `-------'   `--'           `--' `--`--'     `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'   
#
#  <<<  CCSDT-1a Gradient
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        # pytest.param(0, id="grd0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # * vcc turns into ecc
        pytest.param({"call": "c4-ccsdt-1a",  "reference": "rhf",  "fcae": "ae",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ecc",},                                                     }, id="ccsdt-1a  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-1a",  "reference": "rhf",  "fcae": "ae",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ncc",},                                                     }, id="ccsdt-1a  rhf ae: cfour-ncc",  marks=using("cfour")),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt-1a", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "ecc"},                                                      }, id="ccsdt-1a  rhf ae: psi4-cfour-ecc"),
        # DEBUG pytest.param({"call": "p4-mrccsdt-1a", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_reference": "rhf"},                                                                                             }, id="ccsdt-1a  rhf ae: psi4-mrcc"),

        pytest.param({"call": "c4-ccsdt-1a",  "reference": "rhf",  "fcae": "fc",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ecc", "cfour_dropmo": 1,},                "wrong": {1: _w16}}, id="ccsdt-1a  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-1a",  "reference": "rhf",  "fcae": "fc",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ncc", "cfour_dropmo": 1,},                "wrong": {1: _w22}}, id="ccsdt-1a  rhf fc: cfour-ncc",  marks=using("cfour")),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt-1a",  "reference": "rhf",  "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_cfour_cc_program": "ecc", "psi4_cfour_dropmo": [1],},                                                        }, id="ccsdt-1a  rhf fc: psi4-cfour-ecc"),
        # DEBUG pytest.param({"call": "p4-mrccsdt-1a", "reference": "rhf", "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_reference": "rhf", "psi4_freeze_core": True},                                                                   }, id="ccsdt-1a  rhf fc: psi4-mrcc"),

        pytest.param({"call": "c4-ccsdt-1a",  "reference": "uhf",  "fcae": "ae",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_reference": "uhf"},                                     "error": {1: _q19}}, id="ccsdt-1a  uhf ae: cfour",      marks=using("cfour")),

        pytest.param({"call": "c4-ccsdt-1a",  "reference": "rohf", "fcae": "ae",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_reference": "rohf"},                                    "error": {1: _q19}}, id="ccsdt-1a rohf ae: cfour",      marks=using("cfour")),
        # yapf: enable
    ],
)
def test_ccsdt1a_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "gradient"))


#
#   ,-----. ,-----. ,---.  ,------. ,--------.        ,--.            ,--.  ,--.                     ,--.                         
#  '  .--./'  .--./'   .-' |  .-.  \'--.  .--',-----./   | ,--,--.    |  '--'  | ,---.  ,---.  ,---. `--' ,--,--.,--,--,          
#  |  |    |  |    `.  `-. |  |  \  :  |  |   '-----'`|  |' ,-.  |    |  .--.  || .-. :(  .-' (  .-' ,--.' ,-.  ||      \         
#  '  '--'\'  '--'\.-'    ||  '--'  /  |  |           |  |\ '-'  |    |  |  |  |\   --..-'  `).-'  `)|  |\ '-'  ||  ||  |         
#   `-----' `-----'`-----' `-------'   `--'           `--' `--`--'    `--'  `--' `----'`----' `----' `--' `--`--'`--''--'         
#                               
#  <<<  CCSDT-1a Hessian

@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(2, id="hes2"),
        # pytest.param(1, id="hes1", marks=pytest.mark.long),
        # pytest.param(0, id="hes0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "c4-ccsdt-1a",  "reference": "rhf",  "fcae": "ae",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ecc",},                                                     }, id="ccsdt-1a  rhf ae: cfour-ecc",  marks=using("cfour")),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt-1a", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "ecc"},                                                      }, id="ccsdt-1a  rhf ae: psi4-cfour-ecc"),
        # DEBUG pytest.param({"call": "p4-mrccsdt-1a", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_reference": "rhf"},                                                                                             }, id="ccsdt-1a  rhf ae: psi4-mrcc"),

        pytest.param({"call": "c4-ccsdt-1a",  "reference": "rhf",  "fcae": "fc",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ecc", "cfour_dropmo": 1,},                "wrong": {2: _w17}}, id="ccsdt-1a  rhf fc: cfour-ecc",  marks=using("cfour")),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt-1a", "reference": "rhf", "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_dropmo": [1], "psi4_cfour_cc_program": "ecc"},                            }, id="ccsdt-1a  rhf fc: psi4-cfour-ecc"),
        # DEBUG pytest.param({"call": "p4-mrccsdt-1a", "reference": "rhf", "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_reference": "rhf", "psi4_freeze_core": True},                                                                   }, id="ccsdt-1a  rhf fc: psi4-mrcc"),

        pytest.param({"call": "c4-ccsdt-1a",  "reference": "uhf",  "fcae": "ae",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_reference": "uhf"},                                     "error": {2: _q19}}, id="ccsdt-1a  uhf ae: cfour",      marks=using("cfour")),

        pytest.param({"call": "c4-ccsdt-1a",  "reference": "rohf", "fcae": "ae",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_reference": "rohf"},                                    "error": {2: _q19}}, id="ccsdt-1a rohf ae: cfour",      marks=using("cfour")),
        # yapf: enable
    ],
)
def test_ccsdt1a_hessian_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "hessian"))


#                                                                                                                                
#   ,-----. ,-----. ,---.  ,------. ,--------.        ,--.,--.       ,------.                                                    
#  '  .--./'  .--./'   .-' |  .-.  \'--.  .--',-----./   ||  |-.     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.               
#  |  |    |  |    `.  `-. |  |  \  :  |  |   '-----'`|  || .-. '    |  `--, |      \| .-. :|  .--'| .-. |\  '  /                
#  '  '--'\'  '--'\.-'    ||  '--'  /  |  |           |  || `-' |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '                 
#   `-----' `-----'`-----' `-------'   `--'           `--' `---'     `------'`--''--' `----'`--'   .`-  /.-'  /                  
#                                                                                                  `---' `---'                   
#  <<<  CCSDT-1b Energy

@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ]
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "c4-ccsdt-1b",  "reference": "rhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "vcc",},                                                                                                               }, id="ccsdt-1b  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-1b",  "reference": "rhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ecc",},                                                                                                               }, id="ccsdt-1b  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-1b",  "reference": "rhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ncc",},                                                                                                               }, id="ccsdt-1b  rhf ae: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsdt-1b", "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                                                         "error": {0: _q25}}, id="ccsdt-1b  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsdt-1b", "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                                                         "error": {0: _q25}}, id="ccsdt-1b  rhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-ccsdt-1b",  "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                                                         "error": {0: _q25}}, id="ccsdt-1b  rhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-ccsdt-1b",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "vcc", "cfour_dropmo": 1,},                                                                                            }, id="ccsdt-1b  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-1b",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ecc", "cfour_dropmo": 1,},                                                                                            }, id="ccsdt-1b  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-1b",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ncc", "cfour_dropmo": 1,},                                                                                            }, id="ccsdt-1b  rhf fc: cfour-ncc",  marks=using("cfour")),

        pytest.param({"call": "c4-ccsdt-1b",  "reference": "uhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "vcc",},                                                                                     }, id="ccsdt-1b  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-1b",  "reference": "uhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "ecc",},                                                                                     }, id="ccsdt-1b  uhf ae: cfour-ecc",  marks=using("cfour")),

        pytest.param({"call": "c4-ccsdt-1b",  "reference": "uhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "vcc", "cfour_dropmo": 1},                                                                   }, id="ccsdt-1b  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-1b",  "reference": "uhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "ecc", "cfour_dropmo": 1},                                                                   }, id="ccsdt-1b  uhf fc: cfour-ecc",  marks=using("cfour")),
        # yapf: enable
    ],
)
def test_ccsdt1b_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "energy"))


#
#   ,-----. ,-----. ,---.  ,------. ,--------.        ,--.,--.        ,----.                     ,--.,--.                 ,--.   
#  '  .--./'  .--./'   .-' |  .-.  \'--.  .--',-----./   ||  |-.     '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-. 
#  |  |    |  |    `.  `-. |  |  \  :  |  |   '-----'`|  || .-. '    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-' 
#  '  '--'\'  '--'\.-'    ||  '--'  /  |  |           |  || `-' |    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |   
#   `-----' `-----'`-----' `-------'   `--'           `--' `---'      `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'   
#
#  <<<  CCSDT-1b Gradient

@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        # pytest.param(0, id="grd0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # * vcc turns into ecc
        pytest.param({"call": "c4-ccsdt-1b",  "reference": "rhf",  "fcae": "ae",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ecc",},                                                     }, id="ccsdt-1b  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-1b",  "reference": "rhf",  "fcae": "ae",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ncc",},                                                     }, id="ccsdt-1b  rhf ae: cfour-ncc",  marks=using("cfour")),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt-1b", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "ecc"},                                                      }, id="ccsdt-1b  rhf ae: psi4-cfour-ecc"),
        # DEBUG pytest.param({"call": "p4-mrccsdt-1b", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_reference": "rhf"},                                                                                             }, id="ccsdt-1b  rhf ae: psi4-mrcc"),

        pytest.param({"call": "c4-ccsdt-1b",  "reference": "rhf",  "fcae": "fc",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ecc", "cfour_dropmo": 1,},                "wrong": {1: _w18}}, id="ccsdt-1b  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-1b",  "reference": "rhf",  "fcae": "fc",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ncc", "cfour_dropmo": 1,},                "wrong": {1: _w23}}, id="ccsdt-1b  rhf fc: cfour-ncc",  marks=using("cfour")),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt-1b",  "reference": "rhf",  "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_cfour_cc_program": "ecc", "psi4_cfour_dropmo": [1],},                                                        }, id="ccsdt-1b  rhf fc: psi4-cfour-ecc"),
        # DEBUG pytest.param({"call": "p4-mrccsdt-1b", "reference": "rhf", "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_reference": "rhf", "psi4_freeze_core": True},                                                                   }, id="ccsdt-1b  rhf fc: psi4-mrcc"),

        pytest.param({"call": "c4-ccsdt-1b",  "reference": "uhf",  "fcae": "ae",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_reference": "uhf"},                                     "error": {1: _q19}}, id="ccsdt-1b  uhf ae: cfour",      marks=using("cfour")),

        pytest.param({"call": "c4-ccsdt-1b",  "reference": "rohf", "fcae": "ae",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_reference": "rohf"},                                    "error": {1: _q19}}, id="ccsdt-1b rohf ae: cfour",      marks=using("cfour")),
        # yapf: enable
    ],
)
def test_ccsdt1b_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "gradient"))


#
#   ,-----. ,-----. ,---.  ,------. ,--------.        ,--.,--.       ,--.  ,--.                     ,--.                         
#  '  .--./'  .--./'   .-' |  .-.  \'--.  .--',-----./   ||  |-.     |  '--'  | ,---.  ,---.  ,---. `--' ,--,--.,--,--,          
#  |  |    |  |    `.  `-. |  |  \  :  |  |   '-----'`|  || .-. '    |  .--.  || .-. :(  .-' (  .-' ,--.' ,-.  ||      \         
#  '  '--'\'  '--'\.-'    ||  '--'  /  |  |           |  || `-' |    |  |  |  |\   --..-'  `).-'  `)|  |\ '-'  ||  ||  |         
#   `-----' `-----'`-----' `-------'   `--'           `--' `---'     `--'  `--' `----'`----' `----' `--' `--`--'`--''--'         
#                                                                             
#  <<<  CCSDT-1b Hessian

@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(2, id="hes2"),
        # pytest.param(1, id="hes1", marks=pytest.mark.long),
        # pytest.param(0, id="hes0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # * adz ae the two debugs don't match to 1e-5 - went with p4c4 to pass
        pytest.param({"call": "c4-ccsdt-1b",  "reference": "rhf",  "fcae": "ae",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ecc",},                                                     }, id="ccsdt-1b  rhf ae: cfour-ecc",  marks=using("cfour")),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt-1b", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "ecc"},                                                      }, id="ccsdt-1b  rhf ae: psi4-cfour-ecc"),
        # DEBUG pytest.param({"call": "p4-mrccsdt-1b", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_reference": "rhf"},                                                                                             }, id="ccsdt-1b  rhf ae: psi4-mrcc"),

        # * fc is wrong wrt the two debugs for adz - commented b/c not filling in ref for other bases
        # pytest.param({"call": "c4-ccsdt-1b",  "reference": "rhf",  "fcae": "fc",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ecc", "cfour_dropmo": 1,},                "wrong": {2: _w17}}, id="ccsdt-1b  rhf fc: cfour-ecc",  marks=using("cfour")),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt-1b", "reference": "rhf", "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_dropmo": [1], "psi4_cfour_cc_program": "ecc"},                            }, id="ccsdt-1b  rhf fc: psi4-cfour-ecc"),
        # DEBUG pytest.param({"call": "p4-mrccsdt-1b", "reference": "rhf", "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_reference": "rhf", "psi4_freeze_core": True},                                                                   }, id="ccsdt-1b  rhf fc: psi4-mrcc"),

        pytest.param({"call": "c4-ccsdt-1b",  "reference": "uhf",  "fcae": "ae",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_reference": "uhf"},                                     "error": {2: _q19}}, id="ccsdt-1b  uhf ae: cfour",      marks=using("cfour")),

        pytest.param({"call": "c4-ccsdt-1b",  "reference": "rohf", "fcae": "ae",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_reference": "rohf"},                                    "error": {2: _q19}}, id="ccsdt-1b rohf ae: cfour",      marks=using("cfour")),
        # yapf: enable
    ],
)
def test_ccsdt1b_hessian_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "hessian"))


#                                                                                                                           
#   ,-----. ,-----. ,---.  ,------. ,--------.        ,---.     ,------.                                                    
#  '  .--./'  .--./'   .-' |  .-.  \'--.  .--',-----.'.-.  \    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.               
#  |  |    |  |    `.  `-. |  |  \  :  |  |   '-----' .-' .'    |  `--, |      \| .-. :|  .--'| .-. |\  '  /                
#  '  '--'\'  '--'\.-'    ||  '--'  /  |  |          /   '-.    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '                 
#   `-----' `-----'`-----' `-------'   `--'          '-----'    `------'`--''--' `----'`--'   .`-  /.-'  /                  
#                                                                                             `---' `---'                   
#  <<<  CCSDT-2 Energy

@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ]
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "c4-ccsdt-2",  "reference": "rhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "vcc",},                                                                                                               }, id="ccsdt-2  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-2",  "reference": "rhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ecc",},                                                                                                               }, id="ccsdt-2  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-2",  "reference": "rhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ncc",},                                                                                                               }, id="ccsdt-2  rhf ae: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsdt-2", "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                                                         "error": {0: _q26}}, id="ccsdt-2  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsdt-2", "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                                                         "error": {0: _q26}}, id="ccsdt-2  rhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-ccsdt-2",  "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                                                         "error": {0: _q26}}, id="ccsdt-2  rhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-ccsdt-2",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "vcc", "cfour_dropmo": 1,},                                                                                            }, id="ccsdt-2  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-2",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ecc", "cfour_dropmo": 1,},                                                                                            }, id="ccsdt-2  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-2",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ncc", "cfour_dropmo": 1,},                                                                                            }, id="ccsdt-2  rhf fc: cfour-ncc",  marks=using("cfour")),

        pytest.param({"call": "c4-ccsdt-2",  "reference": "uhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "vcc",},                                                                                     }, id="ccsdt-2  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-2",  "reference": "uhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "ecc",},                                                                                     }, id="ccsdt-2  uhf ae: cfour-ecc",  marks=using("cfour")),

        pytest.param({"call": "c4-ccsdt-2",  "reference": "uhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "vcc", "cfour_dropmo": 1},                                                                   }, id="ccsdt-2  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-2",  "reference": "uhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "ecc", "cfour_dropmo": 1},                                                                   }, id="ccsdt-2  uhf fc: cfour-ecc",  marks=using("cfour")),
        # yapf: enable
    ],
)
def test_ccsdt2_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "energy"))


#
#   ,-----. ,-----. ,---.  ,------. ,--------.        ,---.      ,----.                     ,--.,--.                 ,--.   
#  '  .--./'  .--./'   .-' |  .-.  \'--.  .--',-----.'.-.  \    '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-. 
#  |  |    |  |    `.  `-. |  |  \  :  |  |   '-----' .-' .'    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-' 
#  '  '--'\'  '--'\.-'    ||  '--'  /  |  |          /   '-.    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |   
#   `-----' `-----'`-----' `-------'   `--'          '-----'     `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'   
#
#  <<<  CCSDT-2 Gradient

@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        # pytest.param(0, id="grd0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # * ccsdt-2 not available in mrcc
        # * vcc turns into ecc
        pytest.param({"call": "c4-ccsdt-2",  "reference": "rhf",  "fcae": "ae",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ecc",},                                                     }, id="ccsdt-2  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-2",  "reference": "rhf",  "fcae": "ae",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ncc",},                                                     }, id="ccsdt-2  rhf ae: cfour-ncc",  marks=using("cfour")),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt-2", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "ecc"},                                                      }, id="ccsdt-2  rhf ae: psi4-cfour-ecc"),

        pytest.param({"call": "c4-ccsdt-2",  "reference": "rhf",  "fcae": "fc",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ecc", "cfour_dropmo": 1,},                "wrong": {1: _w20}}, id="ccsdt-2  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-2",  "reference": "rhf",  "fcae": "fc",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ncc", "cfour_dropmo": 1,},                "wrong": {1: _w24}}, id="ccsdt-2  rhf fc: cfour-ncc",  marks=using("cfour")),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt-2",  "reference": "rhf",  "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_cfour_cc_program": "ecc", "psi4_cfour_dropmo": [1],},                                                        }, id="ccsdt-2  rhf fc: psi4-cfour-ecc"),

        pytest.param({"call": "c4-ccsdt-2",  "reference": "uhf",  "fcae": "ae",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_reference": "uhf"},                                     "error": {1: _q19}}, id="ccsdt-2  uhf ae: cfour",      marks=using("cfour")),

        pytest.param({"call": "c4-ccsdt-2",  "reference": "rohf", "fcae": "ae",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_reference": "rohf"},                                    "error": {1: _q19}}, id="ccsdt-2 rohf ae: cfour",      marks=using("cfour")),
        # yapf: enable
    ],
)
def test_ccsdt2_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "gradient"))



#
#   ,-----. ,-----. ,---.  ,------. ,--------.        ,---.     ,--.  ,--.                     ,--.                         
#  '  .--./'  .--./'   .-' |  .-.  \'--.  .--',-----.'.-.  \    |  '--'  | ,---.  ,---.  ,---. `--' ,--,--.,--,--,          
#  |  |    |  |    `.  `-. |  |  \  :  |  |   '-----' .-' .'    |  .--.  || .-. :(  .-' (  .-' ,--.' ,-.  ||      \         
#  '  '--'\'  '--'\.-'    ||  '--'  /  |  |          /   '-.    |  |  |  |\   --..-'  `).-'  `)|  |\ '-'  ||  ||  |         
#   `-----' `-----'`-----' `-------'   `--'          '-----'    `--'  `--' `----'`----' `----' `--' `--`--'`--''--'         
#                                                                                       
#  <<<  CCSDT-2 Hessian

@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(2, id="hes2"),
        # pytest.param(1, id="hes1", marks=pytest.mark.long),
        # pytest.param(0, id="hes0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # * ccsdt-2 not available in mrcc
        # * adz has 7e-6 sized discrepancies from debug - adding false fd marking to pass it
        pytest.param({"call": "c4-ccsdt-2",  "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ecc",},                                                     }, id="ccsdt-2  rhf ae: cfour-ecc",  marks=using("cfour")),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt-2", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "ecc"},                                                      }, id="ccsdt-2  rhf ae: psi4-cfour-ecc"),

        # * no use testing hes if grad wrong
        # pytest.param({"call": "c4-ccsdt-2",  "reference": "rhf",  "fcae": "fc",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ecc", "cfour_dropmo": 1,},              "wrong": {2: _w17}}, id="ccsdt-2  rhf fc: cfour-ecc",  marks=using("cfour")),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt-2", "reference": "rhf", "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_dropmo": [1], "psi4_cfour_cc_program": "ecc"},                            }, id="ccsdt-2  rhf fc: psi4-cfour-ecc"),

        pytest.param({"call": "c4-ccsdt-2",  "reference": "uhf",  "fcae": "ae",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_reference": "uhf"},                                     "error": {2: _q19}}, id="ccsdt-2  uhf ae: cfour",      marks=using("cfour")),

        pytest.param({"call": "c4-ccsdt-2",  "reference": "rohf", "fcae": "ae",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_reference": "rohf"},                                    "error": {2: _q19}}, id="ccsdt-2 rohf ae: cfour",      marks=using("cfour")),
        # yapf: enable
    ],
)
def test_ccsdt2_hessian_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "hessian"))


#                                                                                                                           
#   ,-----. ,-----. ,---.  ,------. ,--------.       ,----.     ,------.                                                    
#  '  .--./'  .--./'   .-' |  .-.  \'--.  .--',-----.'.-.  |    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.               
#  |  |    |  |    `.  `-. |  |  \  :  |  |   '-----'  .' <     |  `--, |      \| .-. :|  .--'| .-. |\  '  /                
#  '  '--'\'  '--'\.-'    ||  '--'  /  |  |          /'-'  |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '                 
#   `-----' `-----'`-----' `-------'   `--'          `----'     `------'`--''--' `----'`--'   .`-  /.-'  /                  
#                                                                                             `---' `---'                   
#  <<<  CCSDT-3 Energy

@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ]
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "c4-ccsdt-3",  "reference": "rhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "vcc",},                                                                                                               }, id="ccsdt-3  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-3",  "reference": "rhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ecc",},                                                                                                               }, id="ccsdt-3  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-3",  "reference": "rhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ncc",},                                                                                                               }, id="ccsdt-3  rhf ae: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsdt-3", "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                                                         "error": {0: _q23}}, id="ccsdt-3  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsdt-3", "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                                                         "error": {0: _q23}}, id="ccsdt-3  rhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-ccsdt-3",  "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                                                         "error": {0: _q23}}, id="ccsdt-3  rhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-ccsdt-3",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "vcc", "cfour_dropmo": 1,},                                                                                            }, id="ccsdt-3  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-3",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ecc", "cfour_dropmo": 1,},                                                                                            }, id="ccsdt-3  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-3",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ncc", "cfour_dropmo": 1,},                                                                                            }, id="ccsdt-3  rhf fc: cfour-ncc",  marks=using("cfour")),

        pytest.param({"call": "c4-ccsdt-3",  "reference": "uhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "vcc",},                                                                                     }, id="ccsdt-3  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-3",  "reference": "uhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "ecc",},                                                                                     }, id="ccsdt-3  uhf ae: cfour-ecc",  marks=using("cfour")),

        pytest.param({"call": "c4-ccsdt-3",  "reference": "uhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "vcc", "cfour_dropmo": 1},                                                                   }, id="ccsdt-3  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-3",  "reference": "uhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "ecc", "cfour_dropmo": 1},                                                                   }, id="ccsdt-3  uhf fc: cfour-ecc",  marks=using("cfour")),
        # yapf: enable
    ],
)
def test_ccsdt3_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "energy"))


#
#   ,-----. ,-----. ,---.  ,------. ,--------.       ,----.      ,----.                     ,--.,--.                 ,--.   
#  '  .--./'  .--./'   .-' |  .-.  \'--.  .--',-----.'.-.  |    '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-. 
#  |  |    |  |    `.  `-. |  |  \  :  |  |   '-----'  .' <     |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-' 
#  '  '--'\'  '--'\.-'    ||  '--'  /  |  |          /'-'  |    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |   
#   `-----' `-----'`-----' `-------'   `--'          `----'      `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'   
#
#  <<<  CCSDT-3 Gradient

@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        # pytest.param(0, id="grd0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # * vcc turns into ecc
        pytest.param({"call": "c4-ccsdt-3",  "reference": "rhf",  "fcae": "ae",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ecc",},                                                     }, id="ccsdt-3  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-3",  "reference": "rhf",  "fcae": "ae",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ncc",},                                                     }, id="ccsdt-3  rhf ae: cfour-ncc",  marks=using("cfour")),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt-3", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "ecc"},                                                      }, id="ccsdt-3  rhf ae: psi4-cfour-ecc"),
        # DEBUG pytest.param({"call": "p4-mrccsdt-3", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_reference": "rhf"},                                                                                             }, id="ccsdt-3  rhf ae: psi4-mrcc"),

        # p4c4 by psi findif or mrcc by psi findif differ from ecc analytic by 1e-6
        pytest.param({"call": "c4-ccsdt-3",  "reference": "rhf",  "fcae": "fc",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ecc", "cfour_dropmo": 1,},                "wrong": {1: _w14}}, id="ccsdt-3  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-3",  "reference": "rhf",  "fcae": "fc",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ncc", "cfour_dropmo": 1,},                "wrong": {1: _w13}}, id="ccsdt-3  rhf fc: cfour-ncc",  marks=using("cfour")),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt-3",  "reference": "rhf",  "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_cfour_cc_program": "ecc", "psi4_cfour_dropmo": [1],},                                                        }, id="ccsdt-3  rhf fc: psi4-cfour-ecc"),
        # DEBUG pytest.param({"call": "p4-mrccsdt-3", "reference": "rhf", "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_reference": "rhf", "psi4_freeze_core": True},                                                                   }, id="ccsdt-3  rhf fc: psi4-mrcc"),

        pytest.param({"call": "c4-ccsdt-3",  "reference": "uhf",  "fcae": "ae",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_reference": "uhf"},                                     "error": {1: _q19}}, id="ccsdt-3  uhf ae: cfour",      marks=using("cfour")),

        pytest.param({"call": "c4-ccsdt-3",  "reference": "rohf", "fcae": "ae",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_reference": "rohf"},                                    "error": {1: _q19}}, id="ccsdt-3 rohf ae: cfour",      marks=using("cfour")),
        # yapf: enable
    ],
)
def test_ccsdt3_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "gradient"))


#
#   ,-----. ,-----. ,---.  ,------. ,--------.       ,----.     ,--.  ,--.                     ,--.                         
#  '  .--./'  .--./'   .-' |  .-.  \'--.  .--',-----.'.-.  |    |  '--'  | ,---.  ,---.  ,---. `--' ,--,--.,--,--,          
#  |  |    |  |    `.  `-. |  |  \  :  |  |   '-----'  .' <     |  .--.  || .-. :(  .-' (  .-' ,--.' ,-.  ||      \         
#  '  '--'\'  '--'\.-'    ||  '--'  /  |  |          /'-'  |    |  |  |  |\   --..-'  `).-'  `)|  |\ '-'  ||  ||  |         
#   `-----' `-----'`-----' `-------'   `--'          `----'     `--'  `--' `----'`----' `----' `--' `--`--'`--''--'         
#                                                                                                                           
#  <<<  CCSDT-3 Hessian

@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(2, id="hes2"),
        # pytest.param(1, id="hes1", marks=pytest.mark.long),
        # pytest.param(0, id="hes0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "c4-ccsdt-3",  "reference": "rhf",  "fcae": "ae",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ecc",},                                                     }, id="ccsdt-3  rhf ae: cfour-ecc",  marks=using("cfour")),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt-3", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "ecc"},                                                      }, id="ccsdt-3  rhf ae: psi4-cfour-ecc"),
        # DEBUG pytest.param({"call": "p4-mrccsdt-3", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_reference": "rhf"},                                                                                             }, id="ccsdt-3  rhf ae: psi4-mrcc"),

        # p4c4 by psi findif or mrcc by psi findif differ from ecc analytic gradient by 1e-4
        pytest.param({"call": "c4-ccsdt-3",  "reference": "rhf",  "fcae": "fc",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ecc", "cfour_dropmo": 1,},                "wrong": {2: _w15}}, id="ccsdt-3  rhf fc: cfour-ecc",  marks=using("cfour")),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt-3", "reference": "rhf", "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_dropmo": [1], "psi4_cfour_cc_program": "ecc"},                            }, id="ccsdt-3  rhf fc: psi4-cfour-ecc"),
        # DEBUG pytest.param({"call": "p4-mrccsdt-3", "reference": "rhf", "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_reference": "rhf", "psi4_freeze_core": True},                                                                   }, id="ccsdt-3  rhf fc: psi4-mrcc"),

        pytest.param({"call": "c4-ccsdt-3",  "reference": "uhf",  "fcae": "ae",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_reference": "uhf"},                                     "error": {2: _q19}}, id="ccsdt-3  uhf ae: cfour",      marks=using("cfour")),

        pytest.param({"call": "c4-ccsdt-3",  "reference": "rohf", "fcae": "ae",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_reference": "rohf"},                                    "error": {2: _q19}}, id="ccsdt-3 rohf ae: cfour",      marks=using("cfour")),
        # yapf: enable
    ],
)
def test_ccsdt3_hessian_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "hessian"))


#
#   ,-----. ,-----. ,---.  ,------. ,--------.    ,------.
#  '  .--./'  .--./'   .-' |  .-.  \'--.  .--'    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |    |  |    `.  `-. |  |  \  :  |  |       |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '--'\'  '--'\.-'    ||  '--'  /  |  |       |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----' `-----'`-----' `-------'   `--'       `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                               `---' `---'
#  <<<  CCSDT Energy


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ]
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "c4-ccsdt",  "reference": "rhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", "cfour_cc_program": "vcc",},                                                                                                                            }, id="ccsdt  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt",  "reference": "rhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", "cfour_cc_program": "ecc",},                                                                                                                            }, id="ccsdt  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt",  "reference": "rhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", "cfour_cc_program": "ncc",},                                                                                                                            }, id="ccsdt  rhf ae: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsdt", "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                                                         "error": {0: _q21}}, id="ccsdt  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsdt", "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "tce"},                                                                                                                                                         }, id="ccsdt  rhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-ccsdt",  "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                                                         "error": {0: _q22}}, id="ccsdt  rhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-ccsdt",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", "cfour_cc_program": "vcc", "cfour_dropmo": 1,},                                                                                                         }, id="ccsdt  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", "cfour_cc_program": "ecc", "cfour_dropmo": 1,},                                                                                                         }, id="ccsdt  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", "cfour_cc_program": "ncc", "cfour_dropmo": 1,},                                                                                                         }, id="ccsdt  rhf fc: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "nwc-ccsdt", "reference": "rhf",  "fcae": "fc", "keywords": {"qc_module": "tce", "nwchem_tce__freeze": 1 },                                                                                                                               }, id="ccsdt  rhf fc: nwchem-tce", marks=using("nwchem")),

        pytest.param({"call": "c4-ccsdt",  "reference": "uhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", "cfour_reference": "uhf", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "vcc",},                                                       }, id="ccsdt  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "nwc-ccsdt", "reference": "uhf",  "fcae": "ae", "keywords": {"qc_module": "tce", "nwchem_scf__uhf": True},                                                                                                                                }, id="ccsdt  uhf ae: nwchem-tce", marks=using("nwchem")),

        # cfour uhf/rohf ecc does not converge, ncc does not run
        pytest.param({"call": "c4-ccsdt",  "reference": "uhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", "cfour_reference": "uhf", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "vcc", "cfour_dropmo": 1,},                                    }, id="ccsdt  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "nwc-ccsdt", "reference": "uhf",  "fcae": "fc", "keywords": {"nwchem_tce__freeze": 1, "qc_module": "tce", "nwchem_scf__uhf": True},                                                                                                       }, id="ccsdt  uhf fc: nwchem-tce", marks=using("nwchem")),

        pytest.param({"call": "c4-ccsdt",  "reference": "rohf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", "cfour_reference": "rohf", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "vcc",},                                                     }, id="ccsdt rohf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "nwc-ccsdt", "reference": "rohf", "fcae": "ae", "keywords": {"qc_module": "tce", "nwchem_scf__rohf": True},                                                                                                                               }, id="ccsdt rohf ae: nwchem-tce", marks=using("nwchem")),

        pytest.param({"call": "c4-ccsdt",  "reference": "rohf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", "cfour_reference": "rohf", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "vcc", "cfour_dropmo": 1, "cfour_orbitals": 0},              }, id="ccsdt rohf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "nwc-ccsdt", "reference": "rohf", "fcae": "fc", "keywords": {"nwchem_tce__freeze": 1, "qc_module": "tce", "nwchem_scf__rohf": True},                                                                                                      }, id="ccsdt rohf fc: nwchem-tce", marks=using("nwchem")),
        # yapf: enable
    ],
)
def test_ccsdt_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "energy"))


#   ,-----. ,-----. ,---.  ,------. ,--------.     ,----.                     ,--.,--.                 ,--.
#  '  .--./'  .--./'   .-' |  .-.  \'--.  .--'    '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  |    |  |    `.  `-. |  |  \  :  |  |       |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  '  '--'\'  '--'\.-'    ||  '--'  /  |  |       '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#   `-----' `-----'`-----' `-------'   `--'        `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  CCSDT Gradient


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        # pytest.param(0, id="grd0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # * vcc turns into ecc
        pytest.param({"call": "c4-ccsdt",  "reference": "rhf",  "fcae": "ae",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ecc",},                                                     }, id="ccsdt  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt",  "reference": "rhf",  "fcae": "ae",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ncc",},                                                     }, id="ccsdt  rhf ae: cfour-ncc",  marks=using("cfour")),
        # VLONG pytest.param({"call": "nwc-ccsdt", "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {"basis": "<>", "qc_module": "tce"},                                                                        }, id="ccsdt  rhf ae: nwchem-tce", marks=using("nwchem")),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "ecc"},                                                      }, id="ccsdt  rhf ae: psi4-cfour-ecc"),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "ncc"},                                                      }, id="ccsdt  rhf ae: psi4-cfour-ncc"),

        pytest.param({"call": "c4-ccsdt",  "reference": "rhf",  "fcae": "fc",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ecc", "cfour_dropmo": 1,},                                  }, id="ccsdt  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt",  "reference": "rhf",  "fcae": "fc",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ncc", "cfour_dropmo": 1,},                "wrong": {1: _w7} }, id="ccsdt  rhf fc: cfour-ncc",  marks=using("cfour")),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt", "reference": "rhf", "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_dropmo": [1], "psi4_cfour_cc_program": "ecc"},                            }, id="ccsdt  rhf fc: psi4-cfour-ecc"),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt", "reference": "rhf", "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_dropmo": [1], "psi4_cfour_cc_program": "ncc"},                            }, id="ccsdt  rhf fc: psi4-cfour-ncc"),

        pytest.param({"call": "c4-ccsdt",  "reference": "uhf",  "fcae": "ae",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_reference": "uhf"},                                     "error": {1: _q19}}, id="ccsdt  uhf ae: cfour",      marks=using("cfour")),

        pytest.param({"call": "c4-ccsdt",  "reference": "rohf", "fcae": "ae",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_reference": "rohf"},                                    "error": {1: _q19}}, id="ccsdt rohf ae: cfour",      marks=using("cfour")),
        # yapf: enable
    ],
)
def test_ccsdt_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "gradient"))


#                                                                                                     
#   ,-----. ,-----. ,---.  ,------. ,--------.    ,--.  ,--.                     ,--.                 
#  '  .--./'  .--./'   .-' |  .-.  \'--.  .--'    |  '--'  | ,---.  ,---.  ,---. `--' ,--,--.,--,--,  
#  |  |    |  |    `.  `-. |  |  \  :  |  |       |  .--.  || .-. :(  .-' (  .-' ,--.' ,-.  ||      \ 
#  '  '--'\'  '--'\.-'    ||  '--'  /  |  |       |  |  |  |\   --..-'  `).-'  `)|  |\ '-'  ||  ||  | 
#   `-----' `-----'`-----' `-------'   `--'       `--'  `--' `----'`----' `----' `--' `--`--'`--''--' 
# 
#  <<<  CCSDT Hessian

@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(2, id="hes2"),
        # pytest.param(1, id="hes1", marks=pytest.mark.long),
        # pytest.param(0, id="hes0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # * ncc errors
        pytest.param({"call": "c4-ccsdt",  "reference": "rhf",  "fcae": "ae",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ecc",},                                                     }, id="ccsdt  rhf ae: cfour-ecc",  marks=using("cfour")),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "ecc"},                                                      }, id="ccsdt  rhf ae: psi4-cfour-ecc"),

        pytest.param({"call": "c4-ccsdt",  "reference": "rhf",  "fcae": "fc",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ecc", "cfour_dropmo": 1,},                                  }, id="ccsdt  rhf fc: cfour-ecc",  marks=using("cfour")),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt", "reference": "rhf", "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_dropmo": [1], "psi4_cfour_cc_program": "ecc"},                            }, id="ccsdt  rhf fc: psi4-cfour-ecc"),

        pytest.param({"call": "c4-ccsdt",  "reference": "uhf",  "fcae": "ae",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_reference": "uhf"},                                     "error": {2: _q19}}, id="ccsdt  uhf ae: cfour",      marks=using("cfour")),

        pytest.param({"call": "c4-ccsdt",  "reference": "rohf", "fcae": "ae",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_reference": "rohf"},                                    "error": {2: _q19}}, id="ccsdt rohf ae: cfour",      marks=using("cfour")),
        # yapf: enable
    ],
)
def test_ccsdt_hessian_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "hessian"))


#
#   ,-----. ,-----. ,---.  ,------. ,--------. ,-. ,-----.   ,-.      ,------.
#  '  .--./'  .--./'   .-' |  .-.  \'--.  .--'/ .''  .-.  '  '. \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |    |  |    `.  `-. |  |  \  :  |  |  |  | |  | |  |   |  |    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '--'\'  '--'\.-'    ||  '--'  /  |  |  |  | '  '-'  '-. |  |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----' `-----'`-----' `-------'   `--'   \ '. `-----'--'.' /     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                              `-'           `-'                                    `---' `---'
#  <<<  CCSDT(Q) Energy


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ]
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "c4-ccsdt(q)",  "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                                                                                            }, id="ccsdt_q_  rhf ae: cfour",      marks=using("cfour")),

        pytest.param({"call": "c4-ccsdt(q)",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_dropmo": 1},                                                                                                                                                                           }, id="ccsdt_q_  rhf fc: cfour",      marks=using("cfour")),
        # yapf: enable
    ],
)
def test_ccsdt_prq_pr_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "energy"))


#   ,-----. ,-----. ,---.  ,------. ,--------. ,-. ,-----.   ,-.       ,----.                     ,--.,--.                 ,--.
#  '  .--./'  .--./'   .-' |  .-.  \'--.  .--'/ .''  .-.  '  '. \     '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  |    |  |    `.  `-. |  |  \  :  |  |  |  | |  | |  |   |  |    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  '  '--'\'  '--'\.-'    ||  '--'  /  |  |  |  | '  '-'  '-. |  |    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#   `-----' `-----'`-----' `-------'   `--'   \ '. `-----'--'.' /      `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#                                              `-'           `-'
#  <<<  CCSDT(Q) Gradient


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        # pytest.param(0, id="grd0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are all possible ways of computing <method> working?
        pytest.param({"call": "c4-ccsdt(q)",  "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                                                                                               }, id="ccsdt_q_  rhf ae: cfour",      marks=using("cfour")),

        pytest.param({"call": "c4-ccsdt(q)",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_dropmo": 1},                                                                                                                                                             "wrong": {1: _w7}}, id="ccsdt_q_  rhf fc: cfour",      marks=using("cfour")),
        # yapf: enable
    ],
)
def test_ccsdt_prq_pr_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "gradient"))


#
#   ,-----. ,-----. ,---.  ,------. ,--------. ,-----.       ,------.
#  '  .--./'  .--./'   .-' |  .-.  \'--.  .--''  .-.  '      |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |    |  |    `.  `-. |  |  \  :  |  |   |  | |  |      |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '--'\'  '--'\.-'    ||  '--'  /  |  |   '  '-'  '-.    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----' `-----'`-----' `-------'   `--'    `-----'--'    `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                                          `---' `---'
#  <<<  CCSDTQ Energy


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ]
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        # pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # ecc/vcc sent to cfour, but it switches to ncc
        # adz only converges to cc_conv=9
        pytest.param({"call": "c4-ccsdtq",  "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                                                                                            }, id="ccsdtq  rhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "c4-ccsdtq",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_dropmo": 1},                                                                                                                                                                           }, id="ccsdtq  rhf fc: cfour",      marks=using("cfour")),
        # yapf: enable
    ],
)
def test_ccsdtq_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "energy"))


#
#   ,-----. ,-----. ,---.  ,------. ,--------. ,-----.        ,----.                     ,--.,--.                 ,--.
#  '  .--./'  .--./'   .-' |  .-.  \'--.  .--''  .-.  '      '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  |    |  |    `.  `-. |  |  \  :  |  |   |  | |  |      |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  '  '--'\'  '--'\.-'    ||  '--'  /  |  |   '  '-'  '-.    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#   `-----' `-----'`-----' `-------'   `--'    `-----'--'     `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  CCSDTQ Gradient


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        # pytest.param(0, id="grd0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are all possible ways of computing <method> working?

        # vcc/ecc error out
        pytest.param({"call": "c4-ccsdtq",  "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                                                                                               }, id="ccsdtq  rhf ae: cfour",      marks=using("cfour")),

        pytest.param({"call": "c4-ccsdtq",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_dropmo": 1},                                                                                                                                                             "wrong": {1: _w6}}, id="ccsdtq  rhf fc: cfour",      marks=using("cfour")),
        # yapf: enable
    ],
)
def test_ccsdtq_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "gradient"))


def energy_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, driver):
    qcprog, method = inp["call"].split("-", 1)
    qcprog = _trans_qcprog[qcprog.lower()]
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["reference"])]]

    inpcopy = {k: v for k, v in inp.items()}
    inpcopy["keywords"] = {
        k: (_trans_key(qcprog, basis, k) if v == "<>" else v) for k, v in inpcopy["keywords"].items()
    }

    inpcopy["driver"] = driver
    if not any([k.lower() in _basis_keywords for k in inpcopy["keywords"]]):
        inpcopy["keywords"]["basis"] = basis
    inpcopy["scf_type"] = "pk"
    inpcopy["corl_type"] = "conv"
    inpcopy["qc_module"] = "-".join(
        [qcprog, inp["keywords"].get("qc_module", inp["keywords"].get("cfour_cc_program", ""))]
    ).strip("-")
    print("INP", inpcopy)

    return inpcopy, subject, method, basis, tnm


def _processor(inp, dertype, basis, subjects, clsd_open_pmols, request, driver):
    qcprog, method = inp["call"].split("-", 1)
    qcprog = _trans_qcprog[qcprog.lower()]
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["reference"])]]

    inpcopy = {k: v for k, v in inp.items() if k != "error"}
    if inp.get("error", False) and inp["error"].get(dertype, False):
        inpcopy["error"] = inp["error"][dertype]
    if inp.get("wrong", False) and inp["wrong"].get(dertype, False):
        inpcopy["wrong"] = inp["wrong"][dertype]
    if inp.get("marks", False) and inp["marks"].get(dertype, False):
        pytest.xfail(inp["marks"][dertype])
        # request.node.add_marker(inp["marks"][dertype])

    inpcopy["keywords"] = {
        k: (_trans_key(qcprog, basis, k) if v == "<>" else v) for k, v in inpcopy["keywords"].items()
    }
    inpcopy["driver"] = driver
    if not any([k.lower() in _basis_keywords for k in inpcopy["keywords"]]):
        inpcopy["keywords"]["basis"] = basis
    inpcopy["scf_type"] = "pk"
    inpcopy["corl_type"] = "conv"
    inpcopy["qc_module"] = "-".join(
        [qcprog, inp["keywords"].get("qc_module", inp["keywords"].get("cfour_cc_program", ""))]
    ).strip("-")
    inpcopy["keywords"]["function_kwargs"] = {"dertype": dertype}
    print("INP", inpcopy)

    return inpcopy, subject, method, basis, tnm
