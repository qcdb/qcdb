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


_w1 = ("MP2 CORRELATION ENERGY", "nonstandard answer: NWChem TCE MP2 doesn't report singles (affects ROHF)")
_w2 = ("CCSD CORRELATION ENERGY", "nonstandard answer: GAMESS CCSD ROHF FC energy")
_w3 = ("(T) CORRECTION ENERGY", "nonstandard answer: NWChem CCSD(T) ROHF AE energy")
_w4 = ("MP2 CORRELATION ENERGY", "nonstandard answer: NWChem TCE MP3 doesn't report singles (affects ROHF), may be off by MP2 singles value")
_w5 = ("MP2 CORRELATION ENERGY", "nonstandard answer: GAMESS MP2 ROHF gradient ZAPT energies")
_w6 = ("CCSDTQ CORRELATION ENERGY", "misdirected calc: CFOUR NCC CCSDTQ gradient mixed fc/ae parts")
_w7 = ("CCSDT CORRELATION ENERGY", "misdirected calc: CFOUR NCC CCSDT & CCSDT(Q) gradient mixed fc/ae parts")
# yapf: enable


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
    runner_asserter(*energy_processor(inp, dertype, basis, subjects, clsd_open_pmols, request))


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
        #    pytest.param(0, id="grd0", marks=pytest.mark.long),
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
def test_hf_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*gradient_processor(inp, dertype, basis, subjects, clsd_open_pmols, request))



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
        0,
    ],
    ids=["ene0"],
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
        #    pytest.param(0, id="grd0", marks=pytest.mark.long),
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
        0,
    ],
    ids=["ene0"],
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
        pytest.param({"call": "c4-mp3",  "reference": "rhf",  "fcae": "ae", "keywords": {"basis": "<>", "cfour_scf_conv": 12},                                                                           }, id="mp3  rhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "nwc-mp3", "reference": "rhf",  "fcae": "ae", "keywords": {"basis": "<>", "qc_module": "tce"},                                                                             }, id="mp3  rhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-mp3",  "reference": "rhf",  "fcae": "ae", "keywords": {"basis": "<>", "psi4_mp_type": "conv"},                                                                         }, id="mp3  rhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-mp3",  "reference": "rhf",  "fcae": "fc", "keywords": {"basis": "<>", "cfour_dropmo": 1, "cfour_scf_conv": 12},                                                        }, id="mp3  rhf fc: cfour",      marks=using("cfour")),
        pytest.param({"call": "nwc-mp3", "reference": "rhf",  "fcae": "fc", "keywords": {"basis": "<>", "qc_module": "tce", "nwchem_tce__freeze": 1},                                                    }, id="mp3  rhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-mp3",  "reference": "rhf",  "fcae": "fc", "keywords": {"basis": "<>", "psi4_freeze_core": True, "psi4_mp_type": "conv"},                                               }, id="mp3  rhf fc: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-mp3",  "reference": "uhf",  "fcae": "ae", "keywords": {"basis": "<>", "cfour_reference": "uhf", "cfour_scf_conv": 12},                                                 }, id="mp3  uhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "nwc-mp3", "reference": "uhf",  "fcae": "ae", "keywords": {"basis": "<>", "qc_module": "tce", "nwchem_scf__uhf": True},                                                    }, id="mp3  uhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-mp3",  "reference": "uhf",  "fcae": "ae", "keywords": {"basis": "<>", "reference": "uhf", "psi4_mp_type": "conv"},                                                     }, id="mp3  uhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-mp3",  "reference": "uhf",  "fcae": "fc", "keywords": {"basis": "<>", "cfour_reference": "uhf", "cfour_dropmo": 1, "cfour_scf_conv": 12},                              }, id="mp3  uhf fc: cfour",      marks=using("cfour")),
        pytest.param({"call": "nwc-mp3", "reference": "uhf",  "fcae": "fc", "keywords": {"basis": "<>", "nwchem_tce__freeze": 1, "qc_module": "tce", "nwchem_scf__uhf": True},                           }, id="mp3  uhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-mp3",  "reference": "uhf",  "fcae": "fc", "keywords": {"basis": "<>", "reference": "uhf", "psi4_freeze_core": True, "psi4_mp_type": "conv"},                           }, id="mp3  uhf fc: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-mp3",  "reference": "rohf", "fcae": "ae", "keywords": {"basis": "<>", "cfour_reference": "rohf", "cfour_scf_conv": 12},                                                }, id="mp3 rohf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "nwc-mp3", "reference": "rohf", "fcae": "ae", "keywords": {"basis": "<>", "qc_module": "tce", "nwchem_scf__rohf": True, "nwchem_scf__thresh": 1.e-8},          "wrong": _w4}, id="mp3 rohf ae: nwchem-tce", marks=using("nwchem")),

        pytest.param({"call": "c4-mp3",  "reference": "rohf", "fcae": "fc", "keywords": {"basis": "<>", "cfour_reference": "rohf", "cfour_dropmo": 1, "cfour_scf_conv": 12},                             }, id="mp3 rohf fc: cfour",      marks=using("cfour")),
        pytest.param({"call": "nwc-mp3", "reference": "rohf", "fcae": "fc", "keywords": {"basis": "<>", "nwchem_tce__freeze": 1, "qc_module": "tce", "nwchem_scf__rohf": True},              "wrong": _w4}, id="mp3 rohf fc: nwchem-tce", marks=using("nwchem")),
        # yapf: enable
    ],
)
def test_mp3_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*energy_processor(inp, dertype, basis, subjects, clsd_open_pmols, request))


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
        0,
    ],
    ids=["ene0"],
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
        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "vcc", "cfour_print": 2},                                                                                   }, id="ccsd  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "ecc",},                                                                                                    }, id="ccsd  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "ncc",},                                                                                                    }, id="ccsd  rhf ae: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                                                                                              }, id="ccsd  rhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "rhf",  "fcae": "ae", "keywords": {"gamess_ccinp__ncore": 0},                                                                                                                                                                      }, id="ccsd  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "tce"},                                                                                                                                                                            }, id="ccsd  rhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "cc"},                                                                                                                                                                             }, id="ccsd  rhf ae: nwchem-cc",  marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                                                                                              }, id="ccsd  rhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                                                                                              }, id="ccsd  rhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", "cfour_dropmo": [1], "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "vcc", "cfour_print": 2},                                                              }, id="ccsd  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", "cfour_dropmo": [1], "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "ecc"},                                                                                }, id="ccsd  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", "cfour_dropmo": [1], "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "ncc"},                                                                                }, id="ccsd  rhf fc: cfour-ncc",  marks=using("cfour")),
        # pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_dropmo": 1},                                                                                                                                                                           }, id="ccsd  rhf fc: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "rhf",  "fcae": "fc", "keywords": {},                                                                                                                                                                                              }, id="ccsd  rhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "fc", "keywords": {"qc_module": "tce", "nwchem_tce__freeze": 1 },                                                                                                                                                  }, id="ccsd  rhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "fc", "keywords": {"nwchem_ccsd__freeze": 1},                                                                                                                                                                      }, id="ccsd  rhf fc: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "rhf",  "fcae": "fc", "keywords": {"psi4_freeze_core": True},                                                                                                                                                                      }, id="ccsd  rhf fc: psi4",       marks=using("psi4")),

        # "cfour_occupation": [[3, 1, 1, 0], [3, 0, 1, 0]]
        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "ae", "keywords": {"cfour_BASIS": "<>", "cfour_REFerence": "UHF", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "vcc", "cfour_print": 2},                                                         }, id="ccsd  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "ae", "keywords": {"cfour_BASIS": "<>", "cfour_REFerence": "UHF", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "ecc"},                                                                           }, id="ccsd  uhf ae: cfour-ecc",  marks=using("cfour")),
        # pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "ae", "keywords": {"cfour_reference": "uhf"},                                                                                                                                                                    }, id="ccsd  uhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "uhf",  "fcae": "ae", "keywords": {"gamess_contrl__scftyp": "uhf", "gamess_ccinp__ncore": 0},                                                                                                                         "error": _q2,}, id="ccsd  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "uhf",  "fcae": "ae", "keywords": {"qc_module": "tce", "nwchem_scf__uhf": True},                                                                                                                                                   }, id="ccsd  uhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "uhf",  "fcae": "ae", "keywords": {"nwchem_scf__uhf": True},                                                                                                                                                          "error": _q3,}, id="ccsd  uhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "uhf",  "fcae": "ae", "keywords": {"reference": "uhf"},                                                                                                                                                                            }, id="ccsd  uhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "fc", "keywords": {"cfour_BASIS": "<>", "cfour_dropmo": [1], "cfour_REFerence": "UHF", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "vcc", "cfour_print": 2},                                    }, id="ccsd  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "fc", "keywords": {"cfour_BASIS": "<>", "cfour_dropmo": [1], "cfour_REFerence": "UHF", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "ecc"},                                                      }, id="ccsd  uhf fc: cfour-ecc",  marks=using("cfour")),
        # pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "fc", "keywords": {"cfour_dropmo": 1, "cfour_reference": "uhf"},                                                                                                                                                 }, id="ccsd  uhf fc: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "uhf",  "fcae": "fc", "keywords": {"gamess_contrl__scftyp": "uhf"},                                                                                                                                                   "error": _q2,}, id="ccsd  uhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "uhf",  "fcae": "fc", "keywords": {"nwchem_tce__freeze": 1, "qc_module": "tce", "nwchem_scf__uhf": True},                                                                                                                          }, id="ccsd  uhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "uhf",  "fcae": "fc", "keywords": {"nwchem_ccsd__freeze": 1, "nwchem_scf__uhf": True},                                                                                                                                "error": _q3,}, id="ccsd  uhf fc: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "uhf",  "fcae": "fc", "keywords": {"psi4_freeze_core": True, "reference": "uhf"},                                                                                                                                                  }, id="ccsd  uhf fc: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-ccsd",  "reference": "rohf", "fcae": "ae", "keywords": {"cfour_basis": "<>", "cfour_REFerence": "roHF", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "vcc", "cfour_print": 2},                                                        }, id="ccsd rohf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rohf", "fcae": "ae", "keywords": {"cfour_basis": "<>", "cfour_REFerence": "roHF", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "ecc"},                                                                          }, id="ccsd rohf ae: cfour-ecc",  marks=using("cfour")),
        # pytest.param({"call": "c4-ccsd",  "reference": "rohf", "fcae": "ae", "keywords": {"cfour_reference": "rohf"},                                                                                                                                                                   }, id="ccsd rohf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "rohf", "fcae": "ae", "keywords": {"gamess_contrl__scftyp": "rohf", "gamess_ccinp__ncore": 0, "gamess_ccinp__maxcc": 50},                                                                                                          }, id="ccsd rohf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "rohf", "fcae": "ae", "keywords": {"qc_module": "tce", "nwchem_scf__rohf": True},                                                                                                                                                  }, id="ccsd rohf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "rohf", "fcae": "ae", "keywords": {"basis": "cfour-qz2p", "nwchem_scf__rohf": True},                                                                                                                                  "error": _q3,}, id="ccsd rohf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "rohf", "fcae": "ae", "keywords": {"reference": "rohf", "qc_module": "ccenergy"},                                                                                                                                                  }, id="ccsd rohf ae: psi4",       marks=using("psi4")),  # TODO another way for ccenergy? (fc, too)

        pytest.param({"call": "c4-ccsd",  "reference": "rohf", "fcae": "fc", "keywords": {"cfour_BASIS": "<>", "cfour_dropmo": [1], "cfour_REFerence": "roHF", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_orbitals": 0, "cfour_cc_program": "vcc", "cfour_print": 2},              }, id="ccsd rohf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rohf", "fcae": "fc", "keywords": {"cfour_BASIS": "<>", "cfour_dropmo": [1], "cfour_REFerence": "roHF", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_orbitals": 0, "cfour_cc_program": "ecc"},                                }, id="ccsd rohf fc: cfour-ecc",  marks=using("cfour")),
        # pytest.param({"call": "c4-ccsd",  "reference": "rohf", "fcae": "fc", "keywords": {"cfour_dropmo": 1, "cfour_reference": "rohf", "cfour_orbitals": 0},                                                                                                                           }, id="ccsd rohf fc: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "rohf", "fcae": "fc", "keywords": {"gamess_contrl__scftyp": "rohf", "gamess_ccinp__iconv": 9, "gamess_scf__conv": 9},                                                                                                 "wrong": _w2,}, id="ccsd rohf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "rohf", "fcae": "fc", "keywords": {"nwchem_tce__freeze": 1, "qc_module": "tce", "nwchem_scf__rohf": True},                                                                                                                         }, id="ccsd rohf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "rohf", "fcae": "fc", "keywords": {"nwchem_ccsd__freeze": 1, "nwchem_scf__rohf": True},                                                                                                                               "error": _q3,}, id="ccsd rohf fc: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "rohf", "fcae": "fc", "keywords": {"psi4_e_convergence": 8, "psi4_r_convergence": 7, "psi4_freeze_core": True, "reference": "rohf", "qc_module": "ccenergy"},                                                                      }, id="ccsd rohf fc: psi4",       marks=using("psi4")),
        # yapf: enable
    ],
)
def test_ccsd_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*energy_processor(inp, dertype, basis, subjects, clsd_open_pmols, request))


#@pytest.mark.parametrize("mode", ["driver", "sandwich"])
@pytest.mark.parametrize(
    "dertype",
    [
        0,
    ],
    ids=["ene0"],
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
        0,
    ],
    ids=["ene0"],
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
        pytest.param({"call": "nwc-ccsdt", "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "tce"},                                                                                                                                                         }, id="ccsdt  rhf ae: nwchem-tce", marks=using("nwchem")),

        pytest.param({"call": "c4-ccsdt",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", "cfour_cc_program": "vcc", "cfour_dropmo": 1,},                                                                                                         }, id="ccsdt  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", "cfour_cc_program": "ecc", "cfour_dropmo": 1,},                                                                                                         }, id="ccsdt  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", "cfour_cc_program": "ncc", "cfour_dropmo": 1,},                                                                                                         }, id="ccsdt  rhf fc: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "nwc-ccsdt", "reference": "rhf",  "fcae": "fc", "keywords": {"qc_module": "tce", "nwchem_tce__freeze": 1 },                                                                                                                               }, id="ccsdt  rhf fc: nwchem-tce", marks=using("nwchem")),

        pytest.param({"call": "c4-ccsdt",  "reference": "uhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", "cfour_reference": "uhf", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "vcc",},                                                       }, id="ccsdt  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "nwc-ccsdt", "reference": "uhf",  "fcae": "ae", "keywords": {"qc_module": "tce", "nwchem_scf__uhf": True},                                                                                                                                }, id="ccsdt  uhf ae: nwchem-tce", marks=using("nwchem")),

        pytest.param({"call": "c4-ccsdt",  "reference": "uhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", "cfour_reference": "uhf", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "vcc", "cfour_dropmo": 1,},                                    }, id="ccsdt  uhf fc: cfour-vcc",  marks=using("cfour")),
        # cfour uhf/rohf ecc does not converge, ncc does not run
        pytest.param({"call": "nwc-ccsdt", "reference": "uhf",  "fcae": "fc", "keywords": {"nwchem_tce__freeze": 1, "qc_module": "tce", "nwchem_scf__uhf": True},                                                                                                       }, id="ccsdt  uhf fc: nwchem-tce", marks=using("nwchem")),

        pytest.param({"call": "c4-ccsdt",  "reference": "rohf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", "cfour_reference": "rohf", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "vcc",},                                                     }, id="ccsdt rohf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "nwc-ccsdt", "reference": "rohf", "fcae": "ae", "keywords": {"qc_module": "tce", "nwchem_scf__rohf": True},                                                                                                                               }, id="ccsdt rohf ae: nwchem-tce", marks=using("nwchem")),

        pytest.param({"call": "c4-ccsdt",  "reference": "rohf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", "cfour_reference": "rohf", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "vcc", "cfour_dropmo": 1, "cfour_orbitals": 0},              }, id="ccsdt rohf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "nwc-ccsdt", "reference": "rohf", "fcae": "fc", "keywords": {"nwchem_tce__freeze": 1, "qc_module": "tce", "nwchem_scf__rohf": True},                                                                                                      }, id="ccsdt rohf fc: nwchem-tce", marks=using("nwchem")),
        # yapf: enable
    ],
)
def test_ccsdt_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*energy_processor(inp, dertype, basis, subjects, clsd_open_pmols, request))


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
        #    pytest.param(0, id="grd0", marks=pytest.mark.long),
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
        pytest.param({"call": "c4-ccsdt",  "reference": "rhf",  "fcae": "ae", "keywords": {"cfour_cc_program": "ecc",},                                                                                                                                                                 }, id="ccsdt  rhf ae: cfour-ecc",      marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt",  "reference": "rhf",  "fcae": "ae", "keywords": {"cfour_cc_program": "ncc",},                                                                                                                                                                 }, id="ccsdt  rhf ae: cfour-ncc",      marks=using("cfour")),

        pytest.param({"call": "c4-ccsdt",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_cc_program": "ecc","cfour_dropmo": 1},                                                                                                                                                }, id="ccsdt  rhf fc: cfour-ecc",      marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_cc_program": "ncc","cfour_dropmo": 1},                                                                                                                               "wrong": {1: _w7}}, id="ccsdt  rhf fc: cfour-ncc",      marks=using("cfour")),

        # vcc turns into ecc
#        pytest.param({"call": "c4-ccsdt",  "reference": "rhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "ecc",},                                                                      }, id="ccsdt  rhf ae: cfour-ecc",  marks=using("cfour")),
#        pytest.param({"call": "c4-ccsdt",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "ecc", "cfour_dropmo": 1,},                                                   }, id="ccsdt  rhf fc: cfour-ecc",  marks=using("cfour")),
#        pytest.param({"call": "c4-ccsdt",  "reference": "rhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "ncc",},                                                                      }, id="ccsdt  rhf ae: cfour-ncc",  marks=using("cfour")),
#        pytest.param({"call": "c4-ccsdt",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_basis": "<>", "cfour_SCF_CONV": 12, "cfour_CC_CONV": 12, "cfour_cc_program": "ncc", "cfour_dropmo": 1,},                                                   }, id="ccsdt  rhf fc: cfour-ncc",  marks=using("cfour")),
        # yapf: enable
    ],
)
def test_ccsdt_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*gradient_processor(inp, dertype, basis, subjects, clsd_open_pmols, request))


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
        0,
    ],
    ids=["ene0"],
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
    runner_asserter(*energy_processor(inp, dertype, basis, subjects, clsd_open_pmols, request))


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
        #    pytest.param(0, id="grd0", marks=pytest.mark.long),
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
    runner_asserter(*gradient_processor(inp, dertype, basis, subjects, clsd_open_pmols, request))


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
        0,
    ],
    ids=["ene0"],
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
    runner_asserter(*energy_processor(inp, dertype, basis, subjects, clsd_open_pmols, request))


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
        #    pytest.param(0, id="grd0", marks=pytest.mark.long),
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
    runner_asserter(*gradient_processor(inp, dertype, basis, subjects, clsd_open_pmols, request))


def energy_processor(inp, dertype, basis, subjects, clsd_open_pmols, request):
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

    return inpcopy, subject, method, basis, tnm


def gradient_processor(inp, dertype, basis, subjects, clsd_open_pmols, request):
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
    inpcopy["driver"] = "gradient"
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
