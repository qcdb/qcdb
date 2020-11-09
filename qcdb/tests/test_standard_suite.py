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


_w1 = ("MP2 CORRELATION ENERGY", "nonstandard answer: NWChem TCE MP2 doesn't report singles (affects ROHF)")
_w2 = ("CCSD CORRELATION ENERGY", "nonstandard answer: GAMESS CCSD ROHF FC energy")
_w3 = ("(T) CORRECTION ENERGY", "nonstandard answer: NWChem CCSD(T) ROHF AE energy")
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


#
#  ,--.   ,--.,------.  ,---.     ,------.
#  |   `.'   ||  .--. ''.-.  \    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |'.'|  ||  '--' | .-' .'    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  |   |  ||  | --' /   '-.    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `--'   `--'`--'     '-----'    `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                               `---' `---'


@pytest.mark.parametrize("dertype", [0,], ids=["ene0"])
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
        pytest.param({"call": "gms-mp2", "reference": "rhf",  "fcae": "fc", "keywords": {"basis": "<>", "nwchem_mp2__freeze": 1},                                                                        }, id="mp2  rhf fc: gamess",     marks=using("gamess")),
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
    qcprog, method = inp["call"].split("-", 1)
    qcprog = _trans_qcprog[qcprog.lower()]
    tnm = request.node.name
    subject = clsd_open_pmols[subjects[std_refs.index(inp["reference"])]]

    inpcopy = {k: v for k, v in inp.items()}
    inpcopy["keywords"] = {
        k: (_trans_key(qcprog, basis, k) if v == "<>" else v) for k, v in inpcopy["keywords"].items()
    }

    inpcopy["driver"] = "energy"
    inpcopy["scf_type"] = "pk"
    inpcopy["corl_type"] = "conv"
    inpcopy["qc_module"] = "-".join([qcprog, inp["keywords"].get("qc_module", "")]).strip("-")
    print("INP", inpcopy)

    runner_asserter(inpcopy, subject, method, basis, tnm)


@pytest.mark.parametrize("dertype", [1, pytest.param(0, marks=pytest.mark.long),], ids=["grd1", "grd0"])
@pytest.mark.parametrize(
    "basis, subjects", [pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # pytest.param({"driver": "gradient", "subject": "bh3p", "options": {"reference": "uhf",  "mp2_type": "cd",   "qc_module": "occ", "freeze_core": "true",                     }, "error": {1: _e8f},}, id="mp2  uhf    cd   fc: * dfocc",),
        # yapf: enable
    ],
)
def hide_test_mp2_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    inpcopy = {k: v for k, v in inp.items() if k != "error"}
    if inp.get("error", False) and inp["error"].get(dertype, False):
        inpcopy["error"] = inp["error"][dertype]
    if inp.get("marks", False) and inp["marks"].get(dertype, False):
        pytest.xfail(inp["marks"][dertype])


@pytest.mark.parametrize("dertype", [0,], ids=["ene0"])
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
        # pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                                                                                            }, id="ccsd  rhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "rhf",  "fcae": "ae", "keywords": {"gamess_ccinp__ncore": 0},                                                                                                                                                                      }, id="ccsd  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "tce"},                                                                                                                                                                            }, id="ccsd  rhf ae: nwchem-tce", marks=using("nwchem")),
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
