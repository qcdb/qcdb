import pytest
from qcengine.programs.tests.standard_suite_ref import std_molecules, std_refs

import qcdb

from .standard_suite_runner import runner_asserter
from .test_standard_suite import _processor, clsd_open_pmols
from .utils import using

from .test_standard_suite import _q2, _q3

#@pytest.fixture
#def clsd_open_pmols():
#    frame_not_important = {
#        name[:-4]: qcdb.Molecule.from_string(smol, name=name[:-4])
#        for name, smol in std_molecules.items()
#        if name.endswith("-xyz")
#    }
#    frame_part_of_spec = {
#        name[:-4] + "-fixed": qcdb.Molecule.from_string(smol + "\nno_com\nno_reorient\n", name=name[:-4])
#        for name, smol in std_molecules.items()
#        if name.endswith("-xyz")
#    }
#
#    return {**frame_not_important, **frame_part_of_spec}


#
#   ,-----. ,-----. ,---.  ,------.      ,--.   ,--.                                                                 ,--.
#  '  .--./'  .--./'   .-' |  .-.  \     |   `.'   | ,--,--.,--,--,  ,--,--. ,---.  ,---. ,--,--,--. ,---. ,--,--, ,-'  '-.
#  |  |    |  |    `.  `-. |  |  \  :    |  |'.'|  |' ,-.  ||      \' ,-.  || .-. || .-. :|        || .-. :|      \'-.  .-'
#  '  '--'\'  '--'\.-'    ||  '--'  /    |  |   |  |\ '-'  ||  ||  |\ '-'  |' '-' '\   --.|  |  |  |\   --.|  ||  |  |  |
#   `-----' `-----'`-----' `-------'     `--'   `--' `--`--'`--''--' `--`--'.`-  /  `----'`--`--`--' `----'`--''--'  `--'
#                                                                           `---'
#  <<<  CCSD Management


#@pytest.mark.parametrize(
#    "scramble",
#    # * this parameter alters the input molecule by Cartesian frame and/or atom ordering to imitate real-world inputs.
#    # * scramble dictionary are arguments to qcdb.Molecule.scramble() or qcel.models.Molecule.scramble() that computes a shifted/rotated/atom-mapped input molecule from `subjects` below and the transformations to be applied to the reference data.
#    # * arguments do_shift=T/F and other boolean values generate random perturbations, so multiple lines or runs makes a fuller test.
#    # * specific shifts etc or problematic runs can be reproduced by specifying arrays like the commented example below. note that b/c do_resort is natom-dependent, may need to exclude subjects while debugging.
#    # * `id`s below are numbered since `pytest -k` doesn't recognize case and so duplicate entries can be added. the 0, 1, 2, 3, 9 progression is handy for debugging since perturbations are added systematically
#    [
#        # pytest.param({"do_shift": [ 2.660760432055,  1.477336796939, -2.098045335573], "do_rotate": [[ 0.321861140022,  0.445246880671,  0.835560064745], [-0.447874157747,  0.849136107328, -0.279958229125], [-0.834154749052, -0.28411808546 ,  0.472718487209]], "do_resort": [2, 0, 1]}, id="mTTT3"),
#        # pytest.param({"do_shift": False, "do_rotate": False, "do_resort": False}, id="srm0"),
#        # pytest.param({"do_shift": True, "do_rotate": False, "do_resort": False}, id="Srm1"),
#        # pytest.param({"do_shift": True, "do_rotate": True, "do_resort": False}, id="SRm2"),
#        # pytest.param({"do_shift": False, "do_rotate": False, "do_resort": True}, id="srM3"),
#        # pytest.param({"do_shift": True, "do_rotate": True, "do_resort": True}, id="SRM8"),
#        pytest.param({"do_shift": True, "do_rotate": True, "do_resort": True}, id="SRM9"),
#    ],
#)
#@pytest.mark.parametrize(
#    "frame",
#    # * this parameter alters the input molecule by fix_com and fix_orientation to imitate user signalling frame matters or not.
#    [
#        pytest.param("fixed"),  # fix_=True (no_com/no_reorient); atres.mol.geom = atin.mol.geom aka scrambled
#        pytest.param("free"),  # fix_=False (def)              ; atres.mol.geom = qcdb (aka psi4) interna orientation
#    ],
#)
@pytest.mark.parametrize(
    "driver,dertype",
    [
        pytest.param("energy", 0, id="ene0"),
#        pytest.param("gradient", 1, id="grd1"),
#        pytest.param("hessian", 2, id="hes2"),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
#    # this parameter, along with rhf/uhf below covers four different molecules, with nat=2-4.
    [
#        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
    ],
)
#@pytest.mark.parametrize(
#    "inp",
#    [
#        # yapf: disable
#        pytest.param({"call": "c4-hf",  "reference": "rhf",  "fcae": "ae", "keywords": {"basis": "<>", "cfour_scf_conv": 12},                                                                           }, id="hf  rhf ae: cfour",      marks=using("cfour")),
#        pytest.param({"call": "gms-hf", "reference": "rhf",  "fcae": "ae", "keywords": {"basis": "<>"},                                                                                                 }, id="hf  rhf ae: gamess",     marks=using("gamess")),
#        pytest.param({"call": "nwc-hf", "reference": "rhf",  "fcae": "ae", "keywords": {"basis": "<>"},                                                                                                 }, id="hf  rhf ae: nwchem",     marks=using("nwchem")),
#        pytest.param({"call": "p4-hf",  "reference": "rhf",  "fcae": "ae", "keywords": {"basis": "<>", "psi4_scf_type": "pk"},                                                                          }, id="hf  rhf ae: psi4",       marks=using("psi4")),
#
#        pytest.param({"call": "c4-hf",  "reference": "uhf",  "fcae": "ae", "keywords": {"basis": "<>", "cfour_reference": "uhf", "cfour_scf_conv": 12},                                                 }, id="hf  uhf ae: cfour",      marks=using("cfour")),
#        pytest.param({"call": "gms-hf", "reference": "uhf",  "fcae": "ae", "keywords": {"basis": "<>", "gamess_contrl__scftyp": "uhf"},                                                                 }, id="hf  uhf ae: gamess",     marks=using("gamess")),
#        pytest.param({"call": "nwc-hf", "reference": "uhf",  "fcae": "ae", "keywords": {"basis": "<>", "nwchem_scf__uhf": True},                                                                        }, id="hf  uhf ae: nwchem",     marks=using("nwchem")),
#        pytest.param({"call": "p4-hf",  "reference": "uhf",  "fcae": "ae", "keywords": {"basis": "<>", "reference": "uhf", "psi4_scf_type": "pk"},                                                      }, id="hf  uhf ae: psi4",       marks=using("psi4")),
#        # yapf: enable
#    ],
#)
#@pytest.mark.parametrize(
#    "qcprog",
#    [
#        pytest.param("cfour", marks=using("cfour")),
#        pytest.param("gamess", marks=using("gamess")),
#        pytest.param("nwchem", marks=using("nwchem")),
#        pytest.param("psi4", marks=using("psi4")),
#    ],
#)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
#        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "vcc", "cfour_print": 2},                                                                                        }, id="ccsd  rhf ae: cfour-vcc",  marks=using("cfour")),
#        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "ecc",},                                                                                                         }, id="ccsd  rhf ae: cfour-ecc",  marks=using("cfour")),
#        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "ncc",},                                                                                                         }, id="ccsd  rhf ae: cfour-ncc",  marks=using("cfour")),
#        pytest.param({"call": "gms-ccsd", "reference": "rhf",  "fcae": "ae", "keywords": {"gamess_ccinp__ncore": 0},                                                                                                                        }, id="ccsd  rhf ae: gamess",     marks=using("gamess")),
#        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "tce"},                                                                                                                              }, id="ccsd  rhf ae: nwchem-tce", marks=using("nwchem")),
#        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "cc"},                                                                                                                               }, id="ccsd  rhf ae: nwchem-cc",  marks=using("nwchem")),
#        pytest.param({"call": "p4-ccsd",  "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "ccenergy"},                                                                                                                         }, id="ccsd  rhf ae: psi4-cc",    marks=using("psi4")),
#        pytest.param({"call": "p4-ccsd",  "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "fnocc"},                                                                                                                            }, id="ccsd  rhf ae: psi4-fnocc", marks=using("psi4")),

#        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_print": 2},                                                                   }, id="ccsd  rhf fc: cfour-vcc",  marks=using("cfour")),
#        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": [1], "cfour_cc_program": "ecc"},                                                                                     }, id="ccsd  rhf fc: cfour-ecc",  marks=using("cfour")),
#        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": [1], "cfour_cc_program": "ncc"},                                                                                     }, id="ccsd  rhf fc: cfour-ncc",  marks=using("cfour")),
#        pytest.param({"call": "gms-ccsd", "reference": "rhf",  "fcae": "fc", "keywords": {},                                                                                                                                                }, id="ccsd  rhf fc: gamess",     marks=using("gamess")),
#        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "fc", "keywords": {"nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                                                     }, id="ccsd  rhf fc: nwchem-tce", marks=using("nwchem")),
#        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "fc", "keywords": {"nwchem_ccsd__freeze": 1, "qc_module": "cc"},                                                                                                     }, id="ccsd  rhf fc: nwchem-cc",  marks=using("nwchem")),
#        pytest.param({"call": "p4-ccsd",  "reference": "rhf",  "fcae": "fc", "keywords": {"psi4_freeze_core": True, "qc_module": "ccenergy"},                                                                                               }, id="ccsd  rhf fc: psi4-cc",    marks=using("psi4")),
#        pytest.param({"call": "p4-ccsd",  "reference": "rhf",  "fcae": "fc", "keywords": {"psi4_freeze_core": True, "qc_module": "fnocc"},                                                                                                  }, id="ccsd  rhf fc: psi4-fnocc", marks=using("psi4")),

        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "ae", "cfg": {"module_fallback": False}, "keywords": {"reference": "rhf", "freeze_core": False},                  }, id="mfF ccsd  rhf ae: cfour", marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "rhf",  "fcae": "ae", "cfg": {"module_fallback": False}, "keywords": {"reference": "rhf", "freeze_core": False},                  }, id="mfF ccsd  rhf ae: gamess", marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "ae", "cfg": {"module_fallback": False}, "keywords": {"reference": "rhf", "freeze_core": False},                  }, id="mfF ccsd  rhf ae: nwchem", marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "rhf",  "fcae": "ae", "cfg": {"module_fallback": False}, "keywords": {"reference": "rhf", "freeze_core": False},                  }, id="mfF ccsd  rhf ae: psi4", marks=using("psi4")),

        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "fc", "cfg": {"module_fallback": False}, "keywords": {"reference": "rhf", "freeze_core": True},                   }, id="mfF ccsd  rhf fc: cfour", marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "rhf",  "fcae": "fc", "cfg": {"module_fallback": False}, "keywords": {"reference": "rhf", "freeze_core": True},                   }, id="mfF ccsd  rhf fc: gamess", marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "fc", "cfg": {"module_fallback": False}, "keywords": {"reference": "rhf", "freeze_core": True},                   }, id="mfF ccsd  rhf fc: nwchem", marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "rhf",  "fcae": "fc", "cfg": {"module_fallback": False}, "keywords": {"reference": "rhf", "freeze_core": True},                   }, id="mfF ccsd  rhf fc: psi4", marks=using("psi4")),

        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "ae", "cfg": {"module_fallback": False}, "keywords": {"reference": "uhf", "freeze_core": False},                  }, id="mfF ccsd  uhf ae: cfour", marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "uhf",  "fcae": "ae", "cfg": {"module_fallback": False}, "keywords": {"reference": "uhf", "freeze_core": False}, "error": {0: _q2}}, id="mfF ccsd  uhf ae: gamess", marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "uhf",  "fcae": "ae", "cfg": {"module_fallback": False}, "keywords": {"reference": "uhf", "freeze_core": False}, "error": {0: _q3}}, id="mfF ccsd  uhf ae: nwchem", marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "uhf",  "fcae": "ae", "cfg": {"module_fallback": False}, "keywords": {"reference": "uhf", "freeze_core": False},                  }, id="mfF ccsd  uhf ae: psi4", marks=using("psi4")),


        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "ae", "cfg": {"module_fallback": True}, "keywords": {"reference": "rhf", "freeze_core": False},                   }, id="mfT ccsd  rhf ae: cfour", marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "rhf",  "fcae": "ae", "cfg": {"module_fallback": True}, "keywords": {"reference": "rhf", "freeze_core": False},                   }, id="mfT ccsd  rhf ae: gamess", marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "ae", "cfg": {"module_fallback": True}, "keywords": {"reference": "rhf", "freeze_core": False},                   }, id="mfT ccsd  rhf ae: nwchem", marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "rhf",  "fcae": "ae", "cfg": {"module_fallback": True}, "keywords": {"reference": "rhf", "freeze_core": False},                   }, id="mfT ccsd  rhf ae: psi4", marks=using("psi4")),

        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "fc", "cfg": {"module_fallback": True}, "keywords": {"reference": "rhf", "freeze_core": True},                    }, id="mfT ccsd  rhf fc: cfour", marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "rhf",  "fcae": "fc", "cfg": {"module_fallback": True}, "keywords": {"reference": "rhf", "freeze_core": True},                    }, id="mfT ccsd  rhf fc: gamess", marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "fc", "cfg": {"module_fallback": True}, "keywords": {"reference": "rhf", "freeze_core": True},                    }, id="mfT ccsd  rhf fc: nwchem", marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "rhf",  "fcae": "fc", "cfg": {"module_fallback": True}, "keywords": {"reference": "rhf", "freeze_core": True},                    }, id="mfT ccsd  rhf fc: psi4", marks=using("psi4")),

        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "ae", "cfg": {"module_fallback": True}, "keywords": {"reference": "uhf", "freeze_core": False},                   }, id="mfT ccsd  uhf ae: cfour", marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "uhf",  "fcae": "ae", "cfg": {"module_fallback": True}, "keywords": {"reference": "uhf", "freeze_core": False}, "error": {0: _q2} }, id="mfT ccsd  uhf ae: gamess", marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "uhf",  "fcae": "ae", "cfg": {"module_fallback": True}, "keywords": {"reference": "uhf", "freeze_core": False},                   }, id="mfT ccsd  uhf ae: nwchem", marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "uhf",  "fcae": "ae", "cfg": {"module_fallback": True}, "keywords": {"reference": "uhf", "freeze_core": False},                   }, id="mfT ccsd  uhf ae: psi4", marks=using("psi4")),

#        # "cfour_occupation": [[3, 1, 1, 0], [3, 0, 1, 0]]
#        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_REFerence": "UHF", "cfour_cc_program": "vcc", "cfour_print": 2},                                                              }, id="ccsd  uhf ae: cfour-vcc",  marks=using("cfour")),
#        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_REFerence": "UHF", "cfour_cc_program": "ecc"},                                                                                }, id="ccsd  uhf ae: cfour-ecc",  marks=using("cfour")),
#        pytest.param({"call": "gms-ccsd", "reference": "uhf",  "fcae": "ae", "keywords": {"gamess_contrl__scftyp": "uhf", "gamess_ccinp__ncore": 0},                                                                      "error": {0: _q2 }}, id="ccsd  uhf ae: gamess",     marks=using("gamess")),
#        pytest.param({"call": "nwc-ccsd", "reference": "uhf",  "fcae": "ae", "keywords": {"nwchem_scf__uhf": True, "qc_module": "tce"},                                                                                                     }, id="ccsd  uhf ae: nwchem-tce", marks=using("nwchem")),
#        pytest.param({"call": "nwc-ccsd", "reference": "uhf",  "fcae": "ae", "keywords": {"nwchem_scf__uhf": True, "qc_module": "cc"},                                                                                    "error": {0: _q3 }}, id="ccsd  uhf ae: nwchem-cc",  marks=using("nwchem")),
#        pytest.param({"call": "p4-ccsd",  "reference": "uhf",  "fcae": "ae", "keywords": {"reference": "uhf", "qc_module": "ccenergy"},                                                                                                     }, id="ccsd  uhf ae: psi4",       marks=using("psi4")),
#
#        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight,  "cfour_REFerence": "UHF", "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_print": 2},                                        }, id="ccsd  uhf fc: cfour-vcc",  marks=using("cfour")),
#        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight,  "cfour_REFerence": "UHF", "cfour_dropmo": [1], "cfour_cc_program": "ecc"},                                                          }, id="ccsd  uhf fc: cfour-ecc",  marks=using("cfour")),
#        pytest.param({"call": "gms-ccsd", "reference": "uhf",  "fcae": "fc", "keywords": {"gamess_contrl__scftyp": "uhf"},                                                                                                "error": {0: _q2 }}, id="ccsd  uhf fc: gamess",     marks=using("gamess")),
#        pytest.param({"call": "nwc-ccsd", "reference": "uhf",  "fcae": "fc", "keywords": {"nwchem_scf__uhf": True, "nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                            }, id="ccsd  uhf fc: nwchem-tce", marks=using("nwchem")),
#        pytest.param({"call": "nwc-ccsd", "reference": "uhf",  "fcae": "fc", "keywords": {"nwchem_scf__uhf": True, "nwchem_ccsd__freeze": 1, "qc_module": "cc"},                                                          "error": {0: _q3 }}, id="ccsd  uhf fc: nwchem-cc",  marks=using("nwchem")),
#        pytest.param({"call": "p4-ccsd",  "reference": "uhf",  "fcae": "fc", "keywords": {"psi4_freeze_core": True, "reference": "uhf", "qc_module": "ccenergy"},                                                                           }, id="ccsd  uhf fc: psi4-cc",    marks=using("psi4")),
#
#        pytest.param({"call": "c4-ccsd",  "reference": "rohf", "fcae": "ae", "keywords": {**_c4_tight, "cfour_REFerence": "roHF", "cfour_cc_program": "vcc", "cfour_print": 2},                                                             }, id="ccsd rohf ae: cfour-vcc",  marks=using("cfour")),
#        pytest.param({"call": "c4-ccsd",  "reference": "rohf", "fcae": "ae", "keywords": {**_c4_tight, "cfour_REFerence": "roHF", "cfour_cc_program": "ecc"},                                                                               }, id="ccsd rohf ae: cfour-ecc",  marks=using("cfour")),
#        pytest.param({"call": "gms-ccsd", "reference": "rohf", "fcae": "ae", "keywords": {"gamess_contrl__scftyp": "rohf", "gamess_ccinp__ncore": 0, "gamess_ccinp__maxcc": 50},                                                            }, id="ccsd rohf ae: gamess",     marks=using("gamess")),
#        pytest.param({"call": "nwc-ccsd", "reference": "rohf", "fcae": "ae", "keywords": {"nwchem_scf__rohf": True, "qc_module": "tce"},                                                                                                    }, id="ccsd rohf ae: nwchem-tce", marks=using("nwchem")),
#        pytest.param({"call": "nwc-ccsd", "reference": "rohf", "fcae": "ae", "keywords": {"nwchem_scf__rohf": True, "qc_module": "cc"},                                                                                   "error": {0: _q3 }}, id="ccsd rohf ae: nwchem",     marks=using("nwchem")),
#        pytest.param({"call": "p4-ccsd",  "reference": "rohf", "fcae": "ae", "keywords": {"reference": "rohf", "qc_module": "ccenergy"},                                                                                                    }, id="ccsd rohf ae: psi4-cc",    marks=using("psi4")),
#
#        pytest.param({"call": "c4-ccsd",  "reference": "rohf", "fcae": "fc", "keywords": {**_c4_tight, "cfour_REFerence": "roHF", "cfour_dropmo": [1], "cfour_orbitals": 0, "cfour_cc_program": "vcc", "cfour_print": 2},                   }, id="ccsd rohf fc: 0cfour-vcc",  marks=using("cfour")),
#        pytest.param({"call": "c4-ccsd",  "reference": "rohf", "fcae": "fc", "keywords": {**_c4_tight, "cfour_REFerence": "roHF", "cfour_dropmo": [1], "cfour_orbitals": 0, "cfour_cc_program": "ecc"},                                     }, id="ccsd rohf fc: 0cfour-ecc",  marks=using("cfour")),
#        pytest.param({"call": "gms-ccsd", "reference": "rohf", "fcae": "fc", "keywords": {"gamess_contrl__scftyp": "rohf", "gamess_ccinp__iconv": 9, "gamess_scf__conv": 9},                                              "wrong": {0: _w2 }}, id="ccsd rohf fc: gamess",     marks=using("gamess")),
#        pytest.param({"call": "nwc-ccsd", "reference": "rohf", "fcae": "fc", "keywords": {"nwchem_scf__rohf": True, "nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                           }, id="ccsd rohf fc: nwchem-tce", marks=using("nwchem")),
#        pytest.param({"call": "nwc-ccsd", "reference": "rohf", "fcae": "fc", "keywords": {"nwchem_scf__rohf": True, "nwchem_ccsd__freeze": 1, "qc_module": "cc"},                                                         "error": {0: _q3 }}, id="ccsd rohf fc: nwchem-cc",  marks=using("nwchem")),
#        pytest.param({"call": "p4-ccsd",  "reference": "rohf", "fcae": "fc", "keywords": {"reference": "rohf", "psi4_freeze_core": True,  "qc_module": "ccenergy", "psi4_e_convergence": 8, "psi4_r_convergence": 7},                       }, id="ccsd rohf fc: psi4-cc",    marks=using("psi4")),
        # yapf: enable
    ],
)
def test_ccsd_energy_management(inp, driver, dertype, basis, subjects, clsd_open_pmols, request):
#def test_ccsd_energy_management(inp, scramble, frame, driver, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, driver) #, scramble=scramble, frame=frame)
    )
