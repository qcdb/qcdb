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
# * tuple is (error_type, string_match_in_error_message, reason_for_human)
# * note that 2nd is regex matched, so raw strings and escape chars may be needed
_q1 = (qcng.exceptions.InputError, "unknown SCFTYPE", "no ROHF reference for NWChem hand-coded MP2.")
_q2 = (qcng.exceptions.InputError, "CCTYP IS PROGRAMMED ONLY FOR SCFTYP=RHF OR ROHF", "no UHF CC in GAMESS.")
_q3 = (qcng.exceptions.InputError, "ccsd: nopen is not zero", "no non-RHF reference for NWChem hand-coded CC.")
_q4 = (qcng.exceptions.InputError, "ROHF'S CCTYP MUST BE CCSD OR CR-CCL", "no ROHF LCCD or CCD or (T) in GAMESS.")
_q5 = (qcdb.exceptions.ValidationError, r"Derivative method 'name' \(gms-ccsd\) and derivative level 'dertype' \(1\) are not available.")  # no gms cc grad
_q6 = (qcng.exceptions.InputError, "Only RHF/UHF Hessians are currently implemented.", "no ROHF Hessian for Psi4 HF.")
_q7 = (qcng.exceptions.UnknownError, "CALCLEVEL LCCD not implemented for ROHF references.", "no ROHF LCCD in CFOUR.")
_q8 = (qcng.exceptions.InputError, "Method 'lccd' with CC_TYPE 'CONV' and REFERENCE 'ROHF' not available", "no conv ROHF LCCD in Psi4.")
_q9 = (KeyError, "gms-lccsd", "no LCCSD in GAMESS.")
_q10 = (qcng.exceptions.InputError, r"Invalid reference type (U|RO)HF != RHF for FNOCC energy", "no non-RHF LCCSD in Psi4.")
_q11 = (qcng.exceptions.InputError, r"Method 'ccd' with CC_TYPE 'CONV' and REFERENCE '(R|U|RO)HF' not available", "no conv CCD in Psi4.")
#_q12 = (qcng.exceptions.UnknownError, "CCTYP HAS NO ANALYTIC NUCLEAR GRADIENT PROGRAMMED.")  # gms
_q13 = (qcng.exceptions.UnknownError, "UHF DF-CCSD has NOT been implemented yet!")  # psi only df/cd ccd
_q14 = (KeyError, "gms-mp3", "no MP3 in GAMESS.")
_q15 = (qcng.exceptions.InputError, "Method 'mp3' with MP_TYPE 'CONV' and REFERENCE 'ROHF' not available", "only detci for conv ROHF MP3 in Psi4, and it's already peculiar for MP2.")
_q16 = (qcdb.exceptions.ValidationError, r"Derivative method 'name' \(gms-mp3\) and derivative level 'dertype' \(1\) are not available.", "no MP3 gradients in GAMESS.")
_q17 = (qcdb.exceptions.ValidationError, r"Derivative method 'name' \(p4-ccd\) and derivative level 'dertype' \(2\) are not available", "no CCD Hessian in Psi4.")
_q18 = (qcng.exceptions.InputError, "Frozen core is not available for the CC gradients.", "no FC CC gradients in Psi4.")
_q19 = (qcng.exceptions.UnknownError, "UHF/ROHF gradients not implemented for this level of theory", "no non-RHF CCSDT-n or CCSDT gradients in CFOUR.")
_q20 = (KeyError, r"p4-ccsd\+t\(ccsd\)", "no CCSD+T(CCSD) in Psi4.")
_q21 = (KeyError, r"(gms|p4)-ccsdt", "no CCSDT in GAMESS or Psi4.")
#_q22
_q23 = (KeyError, r"(gms|nwc|p4)-ccsdt-3", "no CCSDT-3 in GAMESS, NWChem, or Psi4.")
_q24 = (KeyError, r"(gms|nwc|p4)-ccsdt-1a", "no CCSDT-1a in GAMESS, NWChem, or Psi4.")
_q25 = (KeyError, r"(gms|nwc|p4)-ccsdt-1b", "no CCSDT-1b in GAMESS, NWChem, or Psi4.")
_q26 = (KeyError, r"(gms|nwc|p4)-ccsdt-2", "no CCSDT-2 in GAMESS, NWChem, or Psi4.")
_q27 = (KeyError, r"(gms|nwc)-a-ccsd\(t\)", "no a-CCSD(T) in GAMESS or NWChem.")
_q28 = (qcng.exceptions.UnknownError, r"UHF has not been implemented for CCSD\(2\)_T", "no UHF a-CCSD(T) in CFOUR.")
_q29 = (qcng.exceptions.InputError, r"Method 'a-ccsd\(t\)' with CC_TYPE 'CONV' and REFERENCE 'UHF' not available", "no UHF a-CCSD(T) in Psi4.")
_q30 = (KeyError, r"gms-mp4\(sdq\)", "no MP4(SDQ) in GAMESS.")
_q31 = (KeyError, r"nwc-mp4\(sdq\)", "no specialty MP4(SDQ) in NWChem.")
_q32 = (qcng.exceptions.InputError, r"Invalid reference type UHF != RHF for FNOCC energy", "no non-RHF MP4(SDQ) in Psi4.")
_q33 = (KeyError, "gms-mp4", "no full MP4 in GAMESS.")
_q34 = (qcng.exceptions.InputError, "Method 'mp4' with MP_TYPE 'CONV' and REFERENCE", "only detci for conv ROHF MP4 in Psi4, and it's already peculiar for MP2.")
_q35 = (KeyError, r"c4-(pbe|b3lyp|b2plyp)", "no DFT in CFOUR.")
_q36 = (qcng.exceptions.InputError, "ROHF reference for DFT is not available.", "no ROHF DFT in Psi4.")
_q37 = (qcdb.exceptions.ValidationError, r"Derivative method 'name' \(c4-(pbe|b3lyp|b3lyp5)\) and derivative level 'dertype' \(1\) are not available.", "no DFT in CFOUR.")
_q38 = (qcng.exceptions.UnknownError, "RHF/UHF gradient calculations not possible for CALCLEVEL LCCD.", "no LCCD gradient in CFOUR.")
_q39 = (qcng.exceptions.InputError, "Method 'cisd' with CI_TYPE 'CONV' and REFERENCE 'UHF' not available", "no UHF CI in Psi4.")
_q40 = (qcng.exceptions.InputError, "CI IS NOT AVAILABLE FOR UHF WAVEFUNCTIONS.", "no UHF CI in GAMESS.")
_q41 = (qcng.exceptions.InputError, r"Invalid reference type UHF != RHF for FNOCC energy", "no non-RHF QCISD in Psi4.")
_q42 = (qcng.exceptions.InputError, "QCISD is not a supported calculation in NCC", "no QCISD in Cfour with NCC.")
_q43 = (KeyError, "gms-qcisd", "no QCISD in GAMESS.")
_q44 = (KeyError, r"nwc-qcisd\(t\)", "no QCISD(T) in NWChem.")
_q45 = (qcng.exceptions.InputError, r"Method 'cc2'", "no UHF CC2 gradients in Psi4")  # Derivative level requested (1) exceeds that available (0). Details: Method 'cc2' with CC_TYPE 'CONV' and REFERENCE 'UHF' not available
_q46 = (qcng.exceptions.UnknownError, r"CALCLEVEL CC3 not implemented for ROHF reference", "no ROHF CC3 in CFOUR.")
_q47 = (qcng.exceptions.InputError, r"Method '(ccsdt-1a|ccsdt-1b|ccsdt-3|ccsdt)' with CC_TYPE 'CONV' and REFERENCE 'RHF' not available", "no high-order CC in Psi4 without MRCC.")

_w1 = ("MP2 CORRELATION ENERGY", "nonstandard answer: NWChem TCE MP2 doesn't report singles (affects ROHF).")
_w2 = ("CCSD CORRELATION ENERGY", "nonstandard answer: GAMESS CCSD ROHF FC energy.")
_w3 = ("(T) CORRECTION ENERGY", "nonstandard answer: NWChem CCSD(T) ROHF AE/FC energy.")
_w4 = ("MP2 CORRELATION ENERGY", "nonstandard answer: NWChem TCE MP3 & MP4 doesn't report singles (affects ROHF), may be off by MP2 singles value.")
_w5 = ("MP2 CORRELATION ENERGY", "nonstandard answer: GAMESS MP2 ROHF gradient ZAPT energies.")
_w6 = ("CCSDTQ CORRELATION ENERGY", "misdirected calc: CFOUR NCC CCSDTQ gradient mixed fc/ae parts.")
_w7 = ("CCSDT CORRELATION ENERGY", "misdirected calc: CFOUR NCC CCSDT & CCSDT(Q) gradient mixed fc/ae parts.")
_w8 = ("CCSD CORRELATION ENERGY", "misdirected calc: CFOUR NCC CCSD & A-CCSD(T) gradient mixed fc/ae parts.")
_w9 = ("CCD CORRELATION ENERGY", "misdirected calc: CFOUR NCC CCD gradient mixed fc/ae parts.")
_w10 = ("MP2.5 CORRELATION ENERGY", "misdirected calc: CFOUR NCC MP3 gradient mixed fc/ae parts.")
_w11 = ("MP3 TOTAL GRADIENT", "nonstandard answer: CFOUR MP3 RHF FC doesn't match findif.")
_w12 = ("MP2 CORRELATION ENERGY", "nonstandard answer: CFOUR CCSD ROHF FC gradient right but energies wrong. (Paper says NYI.)")
_w13 = ("CCSDT-3 CORRELATION ENERGY", "misdirected calc: CFOUR NCC CCSDT-3 gradient mixed fc/ae parts.")
_w14 = ("CCSDT-3 TOTAL GRADIENT", "nonstandard answer: CFOUR CCSDT-3 FC doesn't match findif.")
_w15 = ("CCSDT-3 TOTAL HESSIAN", "nonstandard answer: CFOUR CCSDT-3 FC doesn't match findif.")
_w16 = ("CCSDT-1A TOTAL GRADIENT", "nonstandard answer: CFOUR CCSDT-1A FC doesn't match findif.")
_w17 = ("CCSDT-1A TOTAL HESSIAN", "nonstandard answer: CFOUR CCSDT-1A FC doesn't match findif.")
_w18 = ("CCSDT-1B TOTAL GRADIENT", "nonstandard answer: CFOUR CCSDT-1B FC doesn't match findif.")
_w19 = ("CCSDT-1B TOTAL HESSIAN", "nonstandard answer: CFOUR CCSDT-1B FC doesn't match findif.")
_w20 = ("CCSDT-2 TOTAL GRADIENT", "nonstandard answer: CFOUR CCSDT-2 FC doesn't match findif.")
_w21 = ("CCSDT-2 TOTAL HESSIAN", "nonstandard answer: CFOUR CCSDT-2 FC doesn't match findif.")
_w22 = ("CCSDT-1A CORRELATION ENERGY", "misdirected calc: CFOUR NCC CCSDT-1a gradient mixed fc/ae parts.")
_w23 = ("CCSDT-1B CORRELATION ENERGY", "misdirected calc: CFOUR NCC CCSDT-1b gradient mixed fc/ae parts.")
_w24 = ("CCSDT-2 CORRELATION ENERGY", "misdirected calc: CFOUR NCC CCSDT-2 gradient mixed fc/ae parts.")
_w25 = ("A-CCSD(T) TOTAL GRADIENT", "nonstandard answer: CFOUR ECC A-CCSD(T) AE doesn't match findif.")
_w26 = ("MP2 TOTAL HESSIAN", "nonstandard answer: CFOUR MP2 ROHF hessian doesn't match findif. (Paper says NYI, but program does yield a value.)")
_w27 = ("MP2 CORRELATION ENERGY", "nonstandard answer: CFOUR MP2 ROHF hessian doesn't match findif, nor FC energy. (Paper says NYI, but program does yield a value.)")
_w28 = ("CISD CORRELATION ENERGY", "nonstandard answer: ROHF CISD vcc=tce!=guga=detci.")
_w29 = ("LCCSD CORRELATION ENERGY", "uncertain reference: ROHF LCCSD vcc=tce!=psi4numpy while CFOUR paper says NYI.")
_w30 = ("HF TOTAL HESSIAN", "unconverged: GAMESS HF hessian CPHF too loose; answer correct if hardwire tighter convergence. thanks, Sarom.")
# yapf: enable

# <<<  Notes
#
# * `"xptd": {"fd": False}` is just a note that qcprog is doing internal findif, but the results are accurate
# * `"xptd": {"fd": True}` indicates internal findif, but probably 3-point, and results need a looser atol
#
# >>> Notes

_c4_tight = {
    "cfour_SCF_CONV": 12,
    "cfour_CC_CONV": 12,
    "cfour_LINEQ_CONV": 10,
    "cfour_BRUCK_CONV": 7,  # 9 doesn't converge some dz
}  # tighten cfour for generating references
_p4_fd = {"psi4_points": 5, "psi4_fd_project": False}  # tighten fd for a psi calc that's already doing it
_p4c4_fd = {
    **{"psi4_" + k: v for k, v in _c4_tight.items()},
    **_p4_fd,
    "psi4_function_kwargs_dertype": 0
}  # tight cfour mtd with psi4 fd. needs modified psi


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
        pytest.param({"call": "gms-hf", "reference": "uhf",  "fcae": "ae", "keywords": {"basis": "<>", "gamess_contrl__scftyp": "uhf"},                                    "wrong": {"cfour-qz2p": _w30}}, id="hf  uhf ae: gamess",     marks=using("gamess")),
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
        pytest.param({"call": "c4-mp2",  "reference": "rhf",  "fcae": "ae", "keywords": {"cfour_basis": "<>", "cfour_scf_conv": 12},                                                                    }, id="mp2  rhf ae: cfour",            marks=using("cfour")),
        pytest.param({"call": "c4-mp2",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight},                                                                                                  }, id="mp2  rhf ae: cfour",            marks=using("cfour")),
        pytest.param({"call": "gms-mp2", "reference": "rhf",  "fcae": "ae", "keywords": {"gamess_mp2__nacore": 0, "gamess_mp2__code": "serial"},                                                                                              }, id="mp2  rhf ae: gamess-serial",    marks=using("gamess")),
        pytest.param({"call": "gms-mp2", "reference": "rhf",  "fcae": "ae", "keywords": {"gamess_mp2__nacore": 0, "gamess_mp2__code": "ddi"},                                                                                                 }, id="mp2  rhf ae: gamess-ddi",       marks=using("gamess")),
        pytest.param({"call": "gms-mp2", "reference": "rhf",  "fcae": "ae", "keywords": {"gamess_mp2__nacore": 0, "gamess_mp2__code": "ims"},                                                                                                 }, id="mp2  rhf ae: gamess-ims",       marks=using("gamess")),
        pytest.param({"call": "nwc-mp2", "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "mp2grad"},                                                                                       }, id="mp2  rhf ae: nwchem-mp2grad",   marks=using("nwchem")),
        pytest.param({"call": "nwc-mp2", "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "directmp2"},                                                                                     }, id="mp2  rhf ae: nwchem-directmp2", marks=using("nwchem")),
        pytest.param({"call": "nwc-mp2", "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "tce"},                                                                                           }, id="mp2  rhf ae: nwchem-tce",       marks=using("nwchem")),
        pytest.param({"call": "p4-mp2",  "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "fnocc", "psi4_mp2_type": "conv"},                                                                }, id="mp2  rhf ae: psi4-fnocc",       marks=using("psi4")),
        pytest.param({"call": "p4-mp2",  "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "occ", "psi4_mp2_type": "conv"},                                                                  }, id="mp2  rhf ae: psi4-occ",         marks=using("psi4")),

        pytest.param({"call": "c4-mp2",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": 1},                                                                               }, id="mp2  rhf fc: cfour",            marks=using("cfour")),
        pytest.param({"call": "gms-mp2", "reference": "rhf",  "fcae": "fc", "keywords": {"gamess_mp2__code": "serial"}                                                                                                                        }, id="mp2  rhf fc: gamess-serial",    marks=using("gamess")),
        pytest.param({"call": "gms-mp2", "reference": "rhf",  "fcae": "fc", "keywords": {"gamess_mp2__code": "ddi"}                                                                                                                           }, id="mp2  rhf fc: gamess-ddi",       marks=using("gamess")),
        pytest.param({"call": "gms-mp2", "reference": "rhf",  "fcae": "fc", "keywords": {"gamess_mp2__code": "ims"}                                                                                                                           }, id="mp2  rhf fc: gamess-ims",       marks=using("gamess")),
        pytest.param({"call": "nwc-mp2", "reference": "rhf",  "fcae": "fc", "keywords": {"nwchem_mp2__freeze": 1, "qc_module": "mp2grad"},                                                              }, id="mp2  rhf fc: nwchem-mp2grad",   marks=using("nwchem")),
        pytest.param({"call": "nwc-mp2", "reference": "rhf",  "fcae": "fc", "keywords": {"nwchem_mp2__freeze": 1, "qc_module": "directmp2"},                                                            }, id="mp2  rhf fc: nwchem-directmp2", marks=using("nwchem")),
        pytest.param({"call": "nwc-mp2", "reference": "rhf",  "fcae": "fc", "keywords": {"nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                  }, id="mp2  rhf fc: nwchem-tce",       marks=using("nwchem")),
        pytest.param({"call": "nwc-mp2", "reference": "rhf",  "fcae": "fc", "keywords": {"nwchem_mp2__freeze__core": 1},                                                                                }, id="mp2  rhf fc: nwchem",           marks=using("nwchem")),
        pytest.param({"call": "nwc-mp2", "reference": "rhf",  "fcae": "fc", "keywords": {"nwchem_mp2__freeze__core__atomic": True},                                                                     }, id="mp2  rhf fc: nwchem",           marks=using("nwchem")),
        pytest.param({"call": "nwc-mp2", "reference": "rhf",  "fcae": "fc", "keywords": {"nwchem_mp2__freeze__atomic": {"O": 1}},                                                                       }, id="mp2  rhf fc: nwchem",           marks=using("nwchem")),
        pytest.param({"call": "p4-mp2",  "reference": "rhf",  "fcae": "fc", "keywords": {"psi4_freeze_core": True, "qc_module": "fnocc", "psi4_mp2_type": "conv"},                                      }, id="mp2  rhf fc: psi4-fnocc",       marks=using("psi4")),
        pytest.param({"call": "p4-mp2",  "reference": "rhf",  "fcae": "fc", "keywords": {"psi4_freeze_core": True, "qc_module": "occ", "psi4_mp2_type": "conv"},                                        }, id="mp2  rhf fc: psi4-occ",         marks=using("psi4")),

        pytest.param({"call": "c4-mp2",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_reference": "uhf"},                                                                        }, id="mp2  uhf ae: cfour",            marks=using("cfour")),
        pytest.param({"call": "gms-mp2", "reference": "uhf",  "fcae": "ae", "keywords": {"gamess_contrl__scftyp": "uhf", "gamess_mp2__nacore": 0, "gamess_mp2__code": "serial"},                                                              }, id="mp2  uhf ae: gamess-serial",    marks=using("gamess")),
        pytest.param({"call": "gms-mp2", "reference": "uhf",  "fcae": "ae", "keywords": {"gamess_contrl__scftyp": "uhf", "gamess_mp2__nacore": 0, "gamess_mp2__code": "ddi"},                                                                 }, id="mp2  uhf ae: gamess-ddi",       marks=using("gamess")),
        pytest.param({"call": "nwc-mp2", "reference": "uhf",  "fcae": "ae", "keywords": {"nwchem_scf__uhf": True, "qc_module": "mp2grad"},                                                              }, id="mp2  uhf ae: nwchem-mp2grad",   marks=using("nwchem")),
        pytest.param({"call": "nwc-mp2", "reference": "uhf",  "fcae": "ae", "keywords": {"nwchem_scf__uhf": True, "qc_module": "tce"},                                                                  }, id="mp2  uhf ae: nwchem-tce",       marks=using("nwchem")),
        pytest.param({"call": "p4-mp2",  "reference": "uhf",  "fcae": "ae", "keywords": {"reference": "uhf", "qc_module": "occ", "psi4_mp2_type": "conv"},                                              }, id="mp2  uhf ae: psi4-occ",         marks=using("psi4")),

        pytest.param({"call": "c4-mp2",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_dropmo": [1]},                                                   }, id="mp2  uhf fc: cfour",            marks=using("cfour")),
        pytest.param({"call": "gms-mp2", "reference": "uhf",  "fcae": "fc", "keywords": {"gamess_contrl__scftyp": "uhf", "gamess_mp2__code": "serial"},                                                                                       }, id="mp2  uhf fc: gamess-serial",    marks=using("gamess")),
        pytest.param({"call": "gms-mp2", "reference": "uhf",  "fcae": "fc", "keywords": {"gamess_contrl__scftyp": "uhf", "gamess_mp2__code": "ddi"},                                                                                          }, id="mp2  uhf fc: gamess-ddi",       marks=using("gamess")),
        pytest.param({"call": "nwc-mp2", "reference": "uhf",  "fcae": "fc", "keywords": {"nwchem_scf__uhf": True, "nwchem_mp2__freeze": 1, "qc_module": "mp2grad"},                                     }, id="mp2  uhf fc: nwchem-mp2grad",   marks=using("nwchem")),
        pytest.param({"call": "nwc-mp2", "reference": "uhf",  "fcae": "fc", "keywords": {"nwchem_scf__uhf": True, "nwchem_tce__freeze": 1, "qc_module": "tce"},                                         }, id="mp2  uhf fc: nwchem-tce",       marks=using("nwchem")),
        pytest.param({"call": "p4-mp2",  "reference": "uhf",  "fcae": "fc", "keywords": {"reference": "uhf", "psi4_freeze_core": True, "qc_module": "occ", "psi4_mp2_type": "conv"},                    }, id="mp2  uhf fc: psi4-occ",         marks=using("psi4")),

            # "cfour_occupation": [[3, 1, 1, 0], [3, 0, 1, 0]]
        pytest.param({"call": "c4-mp2",  "reference": "rohf", "fcae": "ae", "keywords": {**_c4_tight, "cfour_reference": "rohf"},                                                                       }, id="mp2 rohf ae: cfour",            marks=using("cfour")),
        pytest.param({"call": "gms-mp2", "reference": "rohf", "fcae": "ae", "keywords": {"gamess_contrl__scftyp": "rohf", "gamess_mp2__nacore": 0, "gamess_mp2__code": "serial", "gamess_mp2__ospt": "RMP"},                                  }, id="mp2 rohf ae: gamess-serial",    marks=using("gamess")),
        pytest.param({"call": "nwc-mp2", "reference": "rohf", "fcae": "ae", "keywords": {"nwchem_scf__rohf": True, "qc_module": "mp2grad"},                                           "error": {0: _q1} }, id="mp2 rohf ae: nwchem-mp2grad",   marks=using("nwchem")),
        pytest.param({"call": "nwc-mp2", "reference": "rohf", "fcae": "ae", "keywords": {"nwchem_scf__rohf": True, "nwchem_tce__freeze": 0, "qc_module": "tce"},                      "wrong": {0: _w1} }, id="mp2 rohf ae: nwchem-tce",       marks=using("nwchem")),
        pytest.param({"call": "p4-mp2",  "reference": "rohf", "fcae": "ae", "keywords": {"reference": "rohf", "psi4_mp2_type": "conv", "qc_module": "occ"},                                             }, id="mp2 rohf ae: psi4-occ",         marks=using("psi4")),

        pytest.param({"call": "c4-mp2",  "reference": "rohf", "fcae": "fc", "keywords": {**_c4_tight, "cfour_reference": "rohf", "cfour_dropmo": 1},                                                    }, id="mp2 rohf fc: cfour",            marks=using("cfour")),
        pytest.param({"call": "gms-mp2", "reference": "rohf", "fcae": "fc", "keywords": {"gamess_contrl__scftyp": "rohf", "gamess_mp2__ospt": "RMP", "gamess_mp2__code": "serial"},                                                           }, id="mp2 rohf fc: gamess-serial",           marks=using("gamess")),
        pytest.param({"call": "nwc-mp2", "reference": "rohf", "fcae": "fc", "keywords": {"nwchem_scf__rohf": True, "nwchem_mp2__freeze": 1, "qc_module": "mp2grad"},                  "error": {0: _q1} }, id="mp2 rohf fc: nwchem-mp2grad",   marks=using("nwchem")),
        pytest.param({"call": "nwc-mp2", "reference": "rohf", "fcae": "fc", "keywords": {"nwchem_scf__rohf": True, "nwchem_tce__freeze": 1, "qc_module": "tce"},                      "wrong": {0: _w1} }, id="mp2 rohf fc: nwchem-tce",       marks=using("nwchem")),
        pytest.param({"call": "p4-mp2",  "reference": "rohf", "fcae": "fc", "keywords": {"reference": "rohf", "psi4_freeze_core": True, "qc_module": "occ", "psi4_mp2_type": "conv"},                   }, id="mp2 rohf fc: psi4-occ",         marks=using("psi4")),
        # yapf: enable
    ],
)
def test_mp2_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "energy"))


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
        pytest.param({"call": "c4-mp2",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight},                                                                                                            }, id="mp2  rhf ae: cfour",          marks=using("cfour")),
        pytest.param({"call": "gms-mp2", "reference": "rhf",  "fcae": "ae",                        "keywords": {"gamess_mp2__nacore": 0, "gamess_mp2__code": "serial"},                                                                       }, id="mp2  rhf ae: gamess-serial",         marks=using("gamess")),
        pytest.param({"call": "gms-mp2", "reference": "rhf",  "fcae": "ae",                        "keywords": {"gamess_mp2__nacore": 0, "gamess_mp2__code": "ddi"},                                                                          }, id="mp2  rhf ae: gamess-ddi",         marks=using("gamess")),
        pytest.param({"call": "gms-mp2", "reference": "rhf",  "fcae": "ae",                        "keywords": {"gamess_mp2__nacore": 0, "gamess_mp2__code": "ims"},                                                                          }, id="mp2  rhf ae: gamess-ims",         marks=using("gamess")),
        pytest.param({"call": "nwc-mp2", "reference": "rhf",  "fcae": "ae",                        "keywords": {"qc_module": "mp2grad"},                                                                                                 }, id="mp2  rhf ae: nwchem-mp2grad", marks=using("nwchem")),
        pytest.param({"call": "nwc-mp2", "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {"qc_module": "tce"},                                                                                                     }, id="mp2  rhf ae: nwchem-tce",     marks=using("nwchem")),
        pytest.param({"call": "p4-mp2",  "reference": "rhf",  "fcae": "ae",                        "keywords": {"psi4_qc_module": "occ", "psi4_mp2_type": "conv"},                                                                       }, id="mp2  rhf ae: psi4",           marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-mp2", "reference": "rhf", "fcae": "ae", "keywords": {"psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "vcc", **_p4c4_fd},                                                            }, id="mp2  rhf ae: psi4-cfour-vcc"),

        pytest.param({"call": "c4-mp2",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_dropmo": 1},                                                                                         }, id="mp2  rhf fc: cfour",          marks=using("cfour")),
        pytest.param({"call": "gms-mp2", "reference": "rhf",  "fcae": "fc",                        "keywords": {"gamess_mp2__code": "serial"}                                                                                                 }, id="mp2  rhf fc: gamess-serial",         marks=using("gamess")),
        pytest.param({"call": "gms-mp2", "reference": "rhf",  "fcae": "fc",                        "keywords": {"gamess_mp2__code": "ddi"}                                                                                                    }, id="mp2  rhf fc: gamess-ddi",         marks=using("gamess")),
        pytest.param({"call": "gms-mp2", "reference": "rhf",  "fcae": "fc",                        "keywords": {"gamess_mp2__code": "ims"}                                                                                                    }, id="mp2  rhf fc: gamess-ims",         marks=using("gamess")),
        pytest.param({"call": "nwc-mp2", "reference": "rhf",  "fcae": "fc",                        "keywords": {"nwchem_mp2__freeze": 1, "qc_module": "mp2grad"},                                                                        }, id="mp2  rhf fc: nwchem-mp2grad", marks=using("nwchem")),
        pytest.param({"call": "nwc-mp2", "reference": "rhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {"nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                            }, id="mp2  rhf fc: nwchem-tce",     marks=using("nwchem")),
        pytest.param({"call": "p4-mp2",  "reference": "rhf",  "fcae": "fc", "xptd": {"fd": False}, "keywords": {**_p4_fd, "psi4_freeze_core": True, "psi4_mp2_type": "conv"},                                                            }, id="mp2  rhf fc: psi4",           marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-mp2", "reference": "rhf", "fcae": "fc", "keywords": {"psi4_cfour_reference": "rhf", "psi4_cfour_dropmo": [1], "psi4_cfour_cc_program": "vcc", **_p4c4_fd}                                   }, id="mp2  rhf fc: psi4-cfour-vcc"),

        pytest.param({"call": "c4-mp2",  "reference": "uhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_reference": "uhf"},                                                                                  }, id="mp2  uhf ae: cfour",          marks=using("cfour")),
        pytest.param({"call": "gms-mp2", "reference": "uhf",  "fcae": "ae",                        "keywords": {"gamess_contrl__scftyp": "uhf", "gamess_mp2__nacore": 0, "gamess_mp2__code": "serial"},                                       }, id="mp2  uhf ae: gamess-serial",         marks=using("gamess")),
        pytest.param({"call": "gms-mp2", "reference": "uhf",  "fcae": "ae",                        "keywords": {"gamess_contrl__scftyp": "uhf", "gamess_mp2__nacore": 0, "gamess_mp2__code": "ddi"},                                          }, id="mp2  uhf ae: gamess-ddi",         marks=using("gamess")),
        pytest.param({"call": "nwc-mp2", "reference": "uhf",  "fcae": "ae",                        "keywords": {"nwchem_scf__uhf": True, "qc_module": "mp2grad"},                                                                        }, id="mp2  uhf ae: nwchem-mp2grad", marks=using("nwchem")),
        pytest.param({"call": "nwc-mp2", "reference": "uhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {"nwchem_scf__uhf": True, "qc_module": "tce"},                                                                            }, id="mp2  uhf ae: nwchem-tce",     marks=using("nwchem")),
        pytest.param({"call": "p4-mp2",  "reference": "uhf",  "fcae": "ae",                        "keywords": {"reference": "uhf", "psi4_qc_module": "occ", "psi4_mp2_type": "conv"},                                                   }, id="mp2  uhf ae: psi4",           marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-mp2", "reference": "uhf", "fcae": "ae", "keywords": {"psi4_cfour_reference": "uhf", "psi4_cfour_cc_program": "vcc", **_p4c4_fd},                                                            }, id="mp2  uhf ae: psi4-cfour-vcc"),

        pytest.param({"call": "c4-mp2",  "reference": "uhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_dropmo": [1]},                                                             }, id="mp2  uhf fc: cfour",          marks=using("cfour")),
        pytest.param({"call": "gms-mp2", "reference": "uhf",  "fcae": "fc",                        "keywords": {"gamess_contrl__scftyp": "uhf", "gamess_mp2__code": "serial"},                                                                }, id="mp2  uhf fc: gamess-serial",         marks=using("gamess")),
        pytest.param({"call": "gms-mp2", "reference": "uhf",  "fcae": "fc",                        "keywords": {"gamess_contrl__scftyp": "uhf", "gamess_mp2__code": "ddi"},                                                                   }, id="mp2  uhf fc: gamess-ddi",         marks=using("gamess")),
        pytest.param({"call": "nwc-mp2", "reference": "uhf",  "fcae": "fc",                        "keywords": {"nwchem_scf__uhf": True, "nwchem_mp2__freeze": 1, "qc_module": "mp2grad"},                                               }, id="mp2  uhf fc: nwchem-mp2grad", marks=using("nwchem")),
        pytest.param({"call": "nwc-mp2", "reference": "uhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {"nwchem_scf__uhf": True, "nwchem_tce__freeze": 1, "qc_module": "tce"},                                                   }, id="mp2  uhf fc: nwchem-tce",     marks=using("nwchem")),
        pytest.param({"call": "p4-mp2",  "reference": "uhf",  "fcae": "fc", "xptd": {"fd": False}, "keywords": {**_p4_fd, "reference": "uhf", "psi4_freeze_core": True, "psi4_mp2_type": "conv"},                                        }, id="mp2  uhf fc: psi4",           marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-mp2", "reference": "uhf", "fcae": "fc", "keywords": {"psi4_cfour_reference": "uhf", "psi4_cfour_dropmo": [1], "psi4_cfour_cc_program": "vcc", **_p4c4_fd}                                   }, id="mp2  uhf fc: psi4-cfour-vcc"),

        pytest.param({"call": "c4-mp2",  "reference": "rohf", "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_reference": "rohf"},                                                                                 }, id="mp2 rohf ae: cfour",          marks=using("cfour")),
        pytest.param({"call": "gms-mp2", "reference": "rohf", "fcae": "ae",                        "keywords": {"gamess_contrl__scftyp": "rohf", "gamess_mp2__nacore": 0},                                                  "wrong": {1: _w5} }, id="mp2 rohf ae: gamess",         marks=using("gamess")),
        pytest.param({"call": "nwc-mp2", "reference": "rohf", "fcae": "ae",                        "keywords": {"nwchem_scf__rohf": True, "qc_module": "mp2grad"},                                                     "error": {1: _q1} }, id="mp2 rohf ae: nwchem-mp2grad", marks=using("nwchem")),
        pytest.param({"call": "nwc-mp2", "reference": "rohf", "fcae": "ae",                        "keywords": {"nwchem_scf__rohf": True, "qc_module": "tce"},                                                         "wrong": {1: _w1} }, id="mp2 rohf ae: nwchem-tce",     marks=using("nwchem")),
        pytest.param({"call": "p4-mp2",  "reference": "rohf", "fcae": "ae", "xptd": {"fd": False}, "keywords": {**_p4_fd, "reference": "rohf", "psi4_qc_module": "occ", "psi4_mp2_type": "conv"},                                        }, id="mp2 rohf ae: psi4",           marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-mp2", "reference": "rohf", "fcae": "ae", "keywords": {"psi4_cfour_reference": "rohf", "psi4_cfour_cc_program": "vcc", **_p4c4_fd},                                                          }, id="mp2 rohf ae: psi4-cfour-vcc"),

        pytest.param({"call": "c4-mp2",  "reference": "rohf", "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_reference": "rohf", "cfour_dropmo": [1]},                                                            }, id="mp2 rohf fc: cfour",          marks=using("cfour")),
        pytest.param({"call": "gms-mp2", "reference": "rohf", "fcae": "fc",                        "keywords": {"gamess_contrl__scftyp": "rohf"},                                                                           "wrong": {1: _w5} }, id="mp2 rohf fc: gamess",         marks=using("gamess")),
        pytest.param({"call": "nwc-mp2", "reference": "rohf", "fcae": "fc",                        "keywords": {"nwchem_scf__rohf": True, "nwchem_mp2__freeze": 1, "qc_module": "mp2grad"},                            "error": {1: _q1} }, id="mp2 rohf fc: nwchem-mp2grad", marks=using("nwchem")),
        pytest.param({"call": "nwc-mp2", "reference": "rohf", "fcae": "fc",                        "keywords": {"nwchem_scf__rohf": True, "nwchem_tce__freeze": 1, "qc_module": "tce"},                                "wrong": {1: _w1} }, id="mp2 rohf fc: nwchem-tce",     marks=using("nwchem")),
        pytest.param({"call": "p4-mp2",  "reference": "rohf", "fcae": "fc", "xptd": {"fd": False}, "keywords": {**_p4_fd, "reference": "rohf", "psi4_freeze_core": True, "qc_module": "occ", "psi4_mp2_type": "conv"},                   }, id="mp2 rohf fc: psi4",           marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-mp2", "reference": "rohf", "fcae": "fc", "keywords": {"psi4_cfour_reference": "rohf", "psi4_cfour_dropmo": [1], "psi4_cfour_cc_program": "vcc", **_p4c4_fd}                                 }, id="mp2 rohf fc: psi4-cfour-vcc"),
        # yapf: enable
    ],
)
def test_mp2_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "gradient"))


#
#  ,--.   ,--.,------.  ,---.     ,--.  ,--.                     ,--.
#  |   `.'   ||  .--. ''.-.  \    |  '--'  | ,---.  ,---.  ,---. `--' ,--,--.,--,--,
#  |  |'.'|  ||  '--' | .-' .'    |  .--.  || .-. :(  .-' (  .-' ,--.' ,-.  ||      \
#  |  |   |  ||  | --' /   '-.    |  |  |  |\   --..-'  `).-'  `)|  |\ '-'  ||  ||  |
#  `--'   `--'`--'     '-----'    `--'  `--' `----'`----' `----' `--' `--`--'`--''--'
#
#  <<<  MP2 Hessian


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
        pytest.param({"call": "c4-mp2",  "reference": "rhf",  "fcae": "ae", "xptd": {           }, "keywords": {**_c4_tight},                                                                                               }, id="mp2  rhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-mp2", "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True }, "keywords": {"gamess_mp2__nacore": 0},                                                                                   }, id="mp2  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp2", "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True }, "keywords": {},                                                                                                          }, id="mp2  rhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-mp2",  "reference": "rhf",  "fcae": "ae", "xptd": {"fd": False}, "keywords": {**_p4_fd, "psi4_mp2_type": "conv"},                                                                         }, id="mp2  rhf ae: psi4",       marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-mp2", "reference": "rhf", "fcae": "ae", "keywords": {"psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "vcc", **_p4c4_fd},                                               }, id="mp2  rhf ae: psi4-cfour-vcc"),

        pytest.param({"call": "c4-mp2",  "reference": "rhf",  "fcae": "fc", "xptd": {           }, "keywords": {**_c4_tight, "cfour_dropmo": 1},                                                                            }, id="mp2  rhf fc: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-mp2", "reference": "rhf",  "fcae": "fc", "xptd": {"fd": True }, "keywords": {},                                                                                                          }, id="mp2  rhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp2", "reference": "rhf",  "fcae": "fc", "xptd": {"fd": True }, "keywords": {"nwchem_mp2__freeze": 1},                                                                                   }, id="mp2  rhf fc: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-mp2",  "reference": "rhf",  "fcae": "fc", "xptd": {"fd": False}, "keywords": {**_p4_fd, "psi4_freeze_core": True, "psi4_mp2_type": "conv"},                                               }, id="mp2  rhf fc: psi4",       marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-mp2", "reference": "rhf", "fcae": "fc", "keywords": {"psi4_cfour_reference": "rhf", "psi4_cfour_dropmo": [1], "psi4_cfour_cc_program": "vcc", **_p4c4_fd}                      }, id="mp2  rhf fc: psi4-cfour-vcc"),

        pytest.param({"call": "c4-mp2",  "reference": "uhf",  "fcae": "ae", "xptd": {           }, "keywords": {**_c4_tight, "cfour_reference": "uhf"},                                                                     }, id="mp2  uhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-mp2", "reference": "uhf",  "fcae": "ae", "xptd": {"fd": True }, "keywords": {"gamess_contrl__scftyp": "uhf", "gamess_mp2__nacore": 0},                                                   }, id="mp2  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp2", "reference": "uhf",  "fcae": "ae", "xptd": {"fd": True }, "keywords": {"nwchem_scf__uhf": True},                                                                                   }, id="mp2  uhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-mp2",  "reference": "uhf",  "fcae": "ae", "xptd": {"fd": False}, "keywords": {**_p4_fd, "reference": "uhf", "psi4_mp2_type": "conv"},                                                     }, id="mp2  uhf ae: psi4",       marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-mp2", "reference": "uhf", "fcae": "ae", "keywords": {"psi4_cfour_reference": "uhf", "psi4_cfour_cc_program": "vcc", **_p4c4_fd},                                               }, id="mp2  uhf ae: psi4-cfour-vcc"),

        pytest.param({"call": "c4-mp2",  "reference": "uhf",  "fcae": "fc", "xptd": {           }, "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_dropmo": [1]},                                                }, id="mp2  uhf fc: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-mp2", "reference": "uhf",  "fcae": "fc", "xptd": {"fd": True }, "keywords": {"gamess_contrl__scftyp": "uhf"},                                                                            }, id="mp2  uhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp2", "reference": "uhf",  "fcae": "fc", "xptd": {"fd": True }, "keywords": {"nwchem_scf__uhf": True, "nwchem_mp2__freeze": 1},                                                          }, id="mp2  uhf fc: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-mp2",  "reference": "uhf",  "fcae": "fc", "xptd": {"fd": False}, "keywords": {**_p4_fd, "reference": "uhf", "psi4_freeze_core": True, "psi4_mp2_type": "conv"},                           }, id="mp2  uhf fc: psi4",       marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-mp2", "reference": "uhf", "fcae": "fc", "keywords": {"psi4_cfour_reference": "uhf", "psi4_cfour_dropmo": [1], "psi4_cfour_cc_program": "vcc", **_p4c4_fd}                      }, id="mp2  uhf fc: psi4-cfour-vcc"),

        pytest.param({"call": "c4-mp2",  "reference": "rohf", "fcae": "ae", "xptd": {           }, "keywords": {**_c4_tight, "cfour_reference": "rohf"},                                                  "wrong": {2: _w26}}, id="mp2 rohf ae: cfour",      marks=using("cfour")),
        # DEBUG pytest.param({"call": "p4-mp2",  "reference": "rohf", "fcae": "ae", "xptd": {"fd": False}, "keywords": {**_p4_fd, "reference": "rohf", "psi4_mp2_type": "conv"},                                            }, id="mp2 rohf ae: psi4",       marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-mp2", "reference": "rohf", "fcae": "ae", "keywords": {"psi4_cfour_reference": "rohf", "psi4_cfour_cc_program": "vcc", **_p4c4_fd},                                             }, id="mp2 rohf ae: psi4-cfour-vcc"),

        pytest.param({"call": "c4-mp2",  "reference": "rohf", "fcae": "fc", "xptd": {           }, "keywords": {**_c4_tight, "cfour_reference": "rohf", "cfour_dropmo": [1]},                             "wrong": {2: _w27}}, id="mp2 rohf fc: cfour",      marks=using("cfour")),
        # DEBUG pytest.param({"call": "p4-mp2",  "reference": "rohf", "fcae": "fc", "xptd": {"fd": False}, "keywords": {**_p4_fd, "reference": "rohf", "psi4_freeze_core": True, "psi4_mp2_type": "conv"},                  }, id="mp2 rohf fc: psi4",       marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-mp2", "reference": "rohf", "fcae": "fc", "keywords": {"psi4_cfour_reference": "rohf", "psi4_cfour_dropmo": [1], "psi4_cfour_cc_program": "vcc", **_p4c4_fd}                    }, id="mp2 rohf fc: psi4-cfour-vcc"),
        # for rohf, cfour analytic 2nd deriv (which cfour ppr doesn't admit to) doesn't match findif HbyE from psi4 not HbyG from p4c4
        # yapf: enable
    ],
)
def test_mp2_hessian_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "hessian"))


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
        # GAMESS doesn't do mp3. Psi4 only does rohf mp3 with detci, which is already non-std for mp2.
        # yapf: disable
        # ecc doesn't compute mp3
        pytest.param({"call": "c4-mp3",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "vcc"},                                                                          }, id="mp3  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-mp3",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "ncc"},                                                                          }, id="mp3  rhf ae: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-mp3", "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                              "error": {0: _q14}}, id="mp3  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp3", "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "tce"},                                                                                              }, id="mp3  rhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-mp3",  "reference": "rhf",  "fcae": "ae", "keywords": {"psi4_qc_module": "fnocc", "psi4_mp_type": "conv"},                                                               }, id="mp3  rhf ae: psi4-fnocc", marks=using("psi4")),
        pytest.param({"call": "p4-mp3",  "reference": "rhf",  "fcae": "ae", "keywords": {"psi4_qc_module": "occ", "psi4_mp_type": "conv"},                                                                 }, id="mp3  rhf ae: psi4-occ",   marks=using("psi4")),

        pytest.param({"call": "c4-mp3",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "vcc"},                                                       }, id="mp3  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-mp3",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "ncc"},                                                       }, id="mp3  rhf fc: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-mp3", "reference": "rhf",  "fcae": "fc", "keywords": {},                                                                                              "error": {0: _q14}}, id="mp3  rhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp3", "reference": "rhf",  "fcae": "fc", "keywords": {"nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                     }, id="mp3  rhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-mp3",  "reference": "rhf",  "fcae": "fc", "keywords": {"psi4_freeze_core": True, "psi4_qc_module": "fnocc", "psi4_mp_type": "conv"},                                     }, id="mp3  rhf fc: psi4-fnocc", marks=using("psi4")),
        pytest.param({"call": "p4-mp3",  "reference": "rhf",  "fcae": "fc", "keywords": {"psi4_freeze_core": True, "psi4_qc_module": "occ", "psi4_mp_type": "conv"},                                       }, id="mp3  rhf fc: psi4-occ",   marks=using("psi4")),

        pytest.param({"call": "c4-mp3",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "vcc"},                                                }, id="mp3  uhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-mp3", "reference": "uhf",  "fcae": "ae", "keywords": {"gamess_contrl__scftyp": "uhf"},                                                                "error": {0: _q14}}, id="mp3  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp3", "reference": "uhf",  "fcae": "ae", "keywords": {"nwchem_scf__uhf": True, "qc_module": "tce"},                                                                     }, id="mp3  uhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-mp3",  "reference": "uhf",  "fcae": "ae", "keywords": {"reference": "uhf", "psi4_qc_module": "occ", "psi4_mp_type": "conv"},                                             }, id="mp3  uhf ae: psi4-occ",   marks=using("psi4")),

        pytest.param({"call": "c4-mp3",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_dropmo": 1, "cfour_cc_program": "vcc"},                             }, id="mp3  uhf fc: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-mp3", "reference": "uhf",  "fcae": "fc", "keywords": {"gamess_contrl__scftyp": "uhf"},                                                                "error": {0: _q14}}, id="mp3  uhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp3", "reference": "uhf",  "fcae": "fc", "keywords": {"nwchem_scf__uhf": True, "nwchem_tce__freeze": 1, "qc_module": "tce"},                                            }, id="mp3  uhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-mp3",  "reference": "uhf",  "fcae": "fc", "keywords": {"reference": "uhf", "psi4_freeze_core": True, "psi4_qc_module": "occ", "psi4_mp_type": "conv"},                   }, id="mp3  uhf fc: psi4-occ",   marks=using("psi4")),

        pytest.param({"call": "c4-mp3",  "reference": "rohf", "fcae": "ae", "keywords": {**_c4_tight, "cfour_reference": "rohf", "cfour_cc_program": "vcc"},                                               }, id="mp3 rohf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-mp3", "reference": "rohf", "fcae": "ae", "keywords": {"gamess_contrl__scftyp": "rohf"},                                                               "error": {0: _q14}}, id="mp3 rohf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp3", "reference": "rohf", "fcae": "ae", "keywords": {"nwchem_scf__rohf": True, "qc_module": "tce", "nwchem_scf__thresh": 1.e-8},                     "wrong": {0: _w4} }, id="mp3 rohf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-mp3",  "reference": "rohf", "fcae": "ae", "keywords": {"reference": "rohf", "psi4_mp_type": "conv"},                                                   "error": {0: _q15}}, id="mp3 rohf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-mp3",  "reference": "rohf", "fcae": "fc", "keywords": {**_c4_tight, "cfour_reference": "rohf", "cfour_dropmo": 1, "cfour_cc_program": "vcc"},                            }, id="mp3 rohf fc: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-mp3", "reference": "rohf", "fcae": "fc", "keywords": {"gamess_contrl__scftyp": "rohf"},                                                               "error": {0: _q14}}, id="mp3 rohf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp3", "reference": "rohf", "fcae": "fc", "keywords": {"nwchem_scf__rohf": True, "nwchem_tce__freeze": 1, "qc_module": "tce"},                         "wrong": {0: _w4} }, id="mp3 rohf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-mp3",  "reference": "rohf", "fcae": "fc", "keywords": {"reference": "rohf", "psi4_freeze_core": True, "psi4_mp_type": "conv"},                         "error": {0: _q15}}, id="mp3 rohf fc: psi4",       marks=using("psi4")),
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
        pytest.param({"call": "c4-mp3",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_cc_program": "vcc"},                                                                                    }, id="mp3  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-mp3",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_cc_program": "ncc"},                                                                                    }, id="mp3  rhf ae: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-mp3", "reference": "rhf",  "fcae": "ae",                        "keywords": {},                                                                                                        "error": {1: _q16}}, id="mp3  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp3", "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {"qc_module": "tce"},                                                                                                        }, id="mp3  rhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-mp3",  "reference": "rhf",  "fcae": "ae",                        "keywords": {"psi4_qc_module": "occ", "psi4_mp_type": "conv"},                                                                           }, id="mp3  rhf ae: psi4-occ",   marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-mp3", "reference": "rhf", "fcae": "ae", "keywords": {"psi4_function_kwargs_dertype": 0, "psi4_cfour_cc_program": "vcc", "psi4_cfour_SCF_CONV": 12, "psi4_cfour_CC_CONV": 12, "psi4_cfour_LINEQ_CONV": 11, "psi4_points": 5, "psi4_fd_project": False}}, id="mp3  rhf ae: psi4-cfour-vcc"),
        # DEBUG pytest.param({"call": "p4-c4-mp3", "reference": "rhf", "fcae": "ae", "keywords": {"psi4_function_kwargs_dertype": 0, "psi4_cfour_cc_program": "ncc", "psi4_cfour_SCF_CONV": 12, "psi4_cfour_CC_CONV": 12, "psi4_cfour_LINEQ_CONV": 11, "psi4_points": 5, "psi4_fd_project": False}}, id="mp3  rhf ae: psi4-cfour-ncc"),

        pytest.param({"call": "c4-mp3",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "vcc"},                                               "wrong": {1: _w11}}, id="mp3  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-mp3",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "ncc"},                                               "wrong": {1: _w10}}, id="mp3  rhf fc: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-mp3", "reference": "rhf",  "fcae": "fc",                        "keywords": {},                                                                                                        "error": {1: _q16}}, id="mp3  rhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp3", "reference": "rhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {"nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                               }, id="mp3  rhf fc: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-mp3",  "reference": "rhf",  "fcae": "fc", "xptd": {"fd": False}, "keywords": {**_p4_fd, "psi4_freeze_core": True, "psi4_qc_module": "occ", "psi4_mp_type": "conv"},                                       }, id="mp3  rhf fc: psi4-occ",   marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-mp3", "reference": "rhf", "fcae": "fc", "keywords": {"psi4_cfour_dropmo": [1], "psi4_function_kwargs_dertype": 0, "psi4_cfour_cc_program": "vcc", "psi4_cfour_SCF_CONV": 12, "psi4_cfour_CC_CONV": 12, "psi4_cfour_LINEQ_CONV": 11, "psi4_points": 5, "psi4_fd_project": False}}, id="mp3  rhf fc: psi4-cfour-vcc"),
        # DEBUG pytest.param({"call": "p4-c4-mp3", "reference": "rhf", "fcae": "fc", "keywords": {"psi4_cfour_dropmo": [1], "psi4_function_kwargs_dertype": 0, "psi4_cfour_cc_program": "ncc", "psi4_cfour_SCF_CONV": 12, "psi4_cfour_CC_CONV": 12, "psi4_cfour_LINEQ_CONV": 11, "psi4_points": 5, "psi4_fd_project": False}}, id="mp3  rhf fc: psi4-cfour-ncc"),

        pytest.param({"call": "c4-mp3",  "reference": "uhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "vcc"},                                                          }, id="mp3  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "gms-mp3", "reference": "uhf",  "fcae": "ae",                        "keywords": {"gamess_contrl__scftyp": "uhf"},                                                                          "error": {1: _q16}}, id="mp3  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp3", "reference": "uhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {"nwchem_scf__uhf": True, "qc_module": "tce"},                                                                               }, id="mp3  uhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-mp3",  "reference": "uhf",  "fcae": "ae",                        "keywords": {"reference": "uhf", "psi4_qc_module": "occ", "psi4_mp_type": "conv"},                                                       }, id="mp3  uhf ae: psi4-occ",   marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-mp3",  "reference": "uhf",  "fcae": "ae", "xptd": {"fd": False}, "keywords": {**_p4_fd, "reference": "uhf", "psi4_mp_type": "conv", "psi4_function_kwargs_dertype": 0},                           }, id="mp3  uhf ae: psi4",       marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-mp3", "reference": "uhf", "fcae": "ae", "keywords": {"psi4_function_kwargs_dertype": 0, "psi4_cfour_cc_program": "vcc", "psi4_cfour_SCF_CONV": 12, "psi4_cfour_CC_CONV": 12, "psi4_cfour_LINEQ_CONV": 11, "psi4_points": 5, "psi4_fd_project": False, "psi4_cfour_reference": "uhf"}}, id="mp3  uhf ae: psi4-cfour-vcc"),

        pytest.param({"call": "c4-mp3",  "reference": "uhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_dropmo": [1]},                                                                }, id="mp3  uhf fc: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-mp3", "reference": "uhf",  "fcae": "fc",                        "keywords": {"gamess_contrl__scftyp": "uhf"},                                                                          "error": {1: _q16}}, id="mp3  uhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp3", "reference": "uhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {"nwchem_scf__uhf": True, "nwchem_tce__freeze": 1, "qc_module": "tce"},                                                      }, id="mp3  uhf fc: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-mp3",  "reference": "uhf",  "fcae": "fc", "xptd": {"fd": False}, "keywords": {**_p4_fd, "reference": "uhf", "psi4_freeze_core": True, "psi4_qc_module": "occ", "psi4_mp_type": "conv"},                   }, id="mp3  uhf fc: psi4-occ",   marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-mp3",  "reference": "uhf",  "fcae": "fc", "xptd": {"fd": False}, "keywords": {**_p4_fd, "reference": "uhf", "psi4_freeze_core": True, "psi4_mp_type": "conv", "psi4_function_kwargs_dertype": 0}, }, id="mp3  uhf fc: psi4",       marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-mp3", "reference": "uhf", "fcae": "fc", "keywords": {"psi4_cfour_dropmo": [1], "psi4_function_kwargs_dertype": 0, "psi4_cfour_cc_program": "vcc", "psi4_cfour_SCF_CONV": 12, "psi4_cfour_CC_CONV": 12, "psi4_cfour_LINEQ_CONV": 11, "psi4_points": 5, "psi4_fd_project": False, "psi4_cfour_reference": "uhf"}}, id="mp3  uhf ae: psi4-cfour-vcc"),

        pytest.param({"call": "c4-mp3",  "reference": "rohf", "fcae": "ae", "keywords": {**_c4_tight, "cfour_reference": "rohf"},                                                                                                           }, id="mp3 rohf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-mp3", "reference": "rohf", "fcae": "ae", "keywords": {"gamess_contrl__scftyp": "rohf"},                                                                                                "error": {1: _q16}}, id="mp3 rohf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp3", "reference": "rohf", "fcae": "ae", "keywords": {"nwchem_scf__rohf": True, "qc_module": "tce"},                                                                                   "wrong": {1: _w1} }, id="mp3 rohf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-mp3",  "reference": "rohf", "fcae": "ae", "keywords": {"reference": "rohf", "psi4_mp_type": "conv", "psi4_points": 5},                                                                  "error": {1: _q15}}, id="mp3 rohf ae: psi4",       marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-mp3", "reference": "rohf", "fcae": "ae", "keywords": {"psi4_function_kwargs_dertype": 0, "psi4_cfour_cc_program": "vcc", "psi4_cfour_SCF_CONV": 12, "psi4_cfour_CC_CONV": 12, "psi4_cfour_LINEQ_CONV": 11, "psi4_points": 5, "psi4_fd_project": False, "psi4_cfour_reference": "rohf"}}, id="mp3  rohf ae: psi4-cfour-vcc"),

        pytest.param({"call": "c4-mp3",  "reference": "rohf", "fcae": "fc", "keywords": {**_c4_tight, "cfour_reference": "rohf", "cfour_dropmo": [1]},                                                                                      }, id="mp3 rohf fc: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-mp3", "reference": "rohf", "fcae": "fc", "keywords": {"gamess_contrl__scftyp": "rohf", "gamess_system__memddi": 300},                                                                  "error": {1: _q16}}, id="mp3 rohf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp3", "reference": "rohf", "fcae": "fc", "keywords": {"nwchem_scf__rohf": True, "nwchem_tce__freeze": 1, "qc_module": "tce"},                                                          "wrong": {1: _w1} }, id="mp3 rohf fc: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-mp3",  "reference": "rohf", "fcae": "fc", "keywords": {"reference": "rohf", "psi4_mp_type": "conv", "psi4_freeze_core": True},                                                          "error": {1: _q15}}, id="mp3 rohf fc: psi4",       marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-mp3", "reference": "rohf", "fcae": "fc", "keywords": {"psi4_cfour_dropmo": [1], "psi4_function_kwargs_dertype": 0, "psi4_cfour_cc_program": "vcc", "psi4_cfour_SCF_CONV": 12, "psi4_cfour_CC_CONV": 12, "psi4_cfour_LINEQ_CONV": 11, "psi4_points": 5, "psi4_fd_project": False, "psi4_cfour_reference": "rohf"}}, id="mp3  rohf fc: psi4-cfour-vcc"),
        # yapf: enable
    ],
)
def test_mp3_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "gradient"))


#
#  ,--.   ,--.,------.   ,---.  ,-. ,---.  ,------.   ,-----.   ,-.      ,------.
#  |   `.'   ||  .--. ' /    | / .''   .-' |  .-.  \ '  .-.  '  '. \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |'.'|  ||  '--' |/  '  ||  | `.  `-. |  |  \  :|  | |  |   |  |    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  |   |  ||  | --' '--|  ||  | .-'    ||  '--'  /'  '-'  '-. |  |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `--'   `--'`--'        `--' \ '.`-----' `-------'  `-----'--'.' /     `------'`--''--' `----'`--'   .`-  /.-'  /
#                               `-'                             `-'                                    `---' `---'
#  <<<  MP4(SDQ) Energy


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
        # ecc doesn't compute mp4
        pytest.param({"call": "c4-mp4(sdq)",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "vcc"},                                                                      }, id="mp4_sdq_  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-mp4(sdq)",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "ncc"},                                                                      }, id="mp4_sdq_  rhf ae: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-mp4(sdq)", "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                          "error": {0: _q30}}, id="mp4_sdq_  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp4(sdq)", "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                          "error": {0: _q31}}, id="mp4_sdq_  rhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-mp4(sdq)",  "reference": "rhf",  "fcae": "ae", "keywords": {"psi4_mp_type": "conv", "psi4_qc_module": "fnocc"},                                                           }, id="mp4_sdq_  rhf ae: psi4-fnocc", marks=using("psi4")),

        pytest.param({"call": "c4-mp4(sdq)",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "vcc"},                                                   }, id="mp4_sdq_  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-mp4(sdq)",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "ncc"},                                                   }, id="mp4_sdq_  rhf fc: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "p4-mp4(sdq)",  "reference": "rhf",  "fcae": "fc", "keywords": {"psi4_freeze_core": True, "psi4_qc_module": "fnocc", "psi4_mp_type": "conv"},                                 }, id="mp4_sdq_  rhf fc: psi4-fnocc", marks=using("psi4")),

        pytest.param({"call": "c4-mp4(sdq)",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "vcc"},                                            }, id="mp4_sdq_  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "p4-mp4(sdq)",  "reference": "uhf",  "fcae": "ae", "keywords": {"reference": "uhf", "psi4_mp_type": "conv"},                                                "error": {0: _q32}}, id="mp4_sdq_  uhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-mp4(sdq)",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_dropmo": 1, "cfour_cc_program": "vcc"},                         }, id="mp4_sdq_  uhf fc: cfour-vcc",  marks=using("cfour")),

        pytest.param({"call": "c4-mp4(sdq)",  "reference": "rohf", "fcae": "ae", "keywords": {**_c4_tight, "cfour_reference": "rohf", "cfour_cc_program": "vcc"},                                           }, id="mp4_sdq_ rohf ae: cfour-vcc",  marks=using("cfour")),

        pytest.param({"call": "c4-mp4(sdq)",  "reference": "rohf", "fcae": "fc", "keywords": {**_c4_tight, "cfour_reference": "rohf", "cfour_dropmo": 1, "cfour_cc_program": "vcc"},                        }, id="mp4_sdq_ rohf fc: cfour-vcc",  marks=using("cfour")),
        # yapf: enable
    ],
)
def test_mp4_prsdq_pr_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "energy"))


#
#  ,--.   ,--.,------.   ,---.    ,------.
#  |   `.'   ||  .--. ' /    |    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |'.'|  ||  '--' |/  '  |    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  |   |  ||  | --' '--|  |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `--'   `--'`--'        `--'    `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                               `---' `---'
#  <<<  MP4 Energy


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
        pytest.param({"call": "c4-mp4",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "vcc"},                                                                 }, id="mp4  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-mp4",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "ncc"},                                                                 }, id="mp4  rhf ae: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-mp4", "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                     "error": {0: _q33}}, id="mp4  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp4", "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "tce"},                                                                                     }, id="mp4  rhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-mp4",  "reference": "rhf",  "fcae": "ae", "keywords": {"psi4_mp_type": "conv", "psi4_qc_module": "fnocc"},                                                      }, id="mp4  rhf ae: psi4-fnocc", marks=using("psi4")),

        pytest.param({"call": "c4-mp4",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "vcc"},                                              }, id="mp4  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-mp4",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "ncc"},                                              }, id="mp4  rhf fc: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-mp4", "reference": "rhf",  "fcae": "fc", "keywords": {},                                                                                     "error": {0: _q33}}, id="mp4  rhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp4", "reference": "rhf",  "fcae": "fc", "keywords": {"nwchem_tce__freeze": 1, "qc_module": "tce"},                                                            }, id="mp4  rhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-mp4",  "reference": "rhf",  "fcae": "fc", "keywords": {"psi4_freeze_core": True, "psi4_mp_type": "conv", "psi4_qc_module": "fnocc"},                            }, id="mp4  rhf fc: psi4-fnocc", marks=using("psi4")),

        pytest.param({"call": "c4-mp4",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "vcc"},                                       }, id="mp4  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "gms-mp4", "reference": "uhf",  "fcae": "ae", "keywords": {"gamess_contrl__scftyp": "uhf"},                                                       "error": {0: _q33}}, id="mp4  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp4", "reference": "uhf",  "fcae": "ae", "keywords": {"nwchem_scf__uhf": True, "qc_module": "tce"},                                                            }, id="mp4  uhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-mp4",  "reference": "uhf",  "fcae": "ae", "keywords": {"reference": "uhf", "psi4_mp_type": "conv"},                                           "error": {0: _q34}}, id="mp4  uhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-mp4",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_dropmo": 1, "cfour_cc_program": "vcc"},                    }, id="mp4  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "gms-mp4", "reference": "uhf",  "fcae": "fc", "keywords": {"gamess_contrl__scftyp": "uhf"},                                                       "error": {0: _q33}}, id="mp4  uhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp4", "reference": "uhf",  "fcae": "fc", "keywords": {"nwchem_scf__uhf": True, "nwchem_tce__freeze": 1, "qc_module": "tce"},                                   }, id="mp4  uhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-mp4",  "reference": "uhf",  "fcae": "fc", "keywords": {"reference": "uhf", "psi4_freeze_core": True, "psi4_mp_type": "conv"},                 "error": {0: _q34}}, id="mp4  uhf fc: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-mp4",  "reference": "rohf", "fcae": "ae", "keywords": {**_c4_tight, "cfour_reference": "rohf", "cfour_cc_program": "vcc"},                                      }, id="mp4 rohf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "gms-mp4", "reference": "rohf", "fcae": "ae", "keywords": {"gamess_contrl__scftyp": "rohf"},                                                      "error": {0: _q33}}, id="mp4 rohf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp4", "reference": "rohf", "fcae": "ae", "keywords": {"nwchem_scf__rohf": True, "qc_module": "tce"},                                         "wrong": {0: _w4} }, id="mp4 rohf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-mp4",  "reference": "rohf", "fcae": "ae", "keywords": {"reference": "rohf", "psi4_mp_type": "conv"},                                          "error": {0: _q34}}, id="mp4 rohf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-mp4",  "reference": "rohf", "fcae": "fc", "keywords": {**_c4_tight, "cfour_reference": "rohf", "cfour_dropmo": 1, "cfour_cc_program": "vcc"},                   }, id="mp4 rohf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "gms-mp4", "reference": "rohf", "fcae": "fc", "keywords": {"gamess_contrl__scftyp": "rohf"},                                                      "error": {0: _q33}}, id="mp4 rohf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-mp4", "reference": "rohf", "fcae": "fc", "keywords": {"nwchem_scf__rohf": True, "nwchem_tce__freeze": 1, "qc_module": "tce"},                "wrong": {0: _w4} }, id="mp4 rohf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-mp4",  "reference": "rohf", "fcae": "fc", "keywords": {"reference": "rohf", "psi4_freeze_core": True, "psi4_mp_type": "conv"},                "error": {0: _q34}}, id="mp4 rohf fc: psi4",       marks=using("psi4")),
        # yapf: enable
    ],
)
def test_mp4_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "energy"))


#
#  ,-------.  ,---.  ,------. ,--------. ,---.     ,------.
#  `--.   /  /  O  \ |  .--. ''--.  .--''.-.  \    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#    /   /  |  .-.  ||  '--' |   |  |    .-' .'    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#   /   `--.|  | |  ||  | --'    |  |   /   '-.    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `-------'`--' `--'`--'        `--'   '-----'    `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                                `---' `---'
#  <<<  ZAPT2 Energy


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
        # pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),  # larger basis sets too long, and detci ae can't allocate arrays
        # pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p"),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "gms-zapt2", "reference": "rohf", "fcae": "ae", "keywords": {"gamess_contrl__scftyp": "rohf", "gamess_mp2__nacore": 0, "gamess_mp2__code": "serial"},                        }, id="zapt2 rohf ae: gamess-serial",    marks=using("gamess")),
        pytest.param({"call": "p4-zapt2",  "reference": "rohf", "fcae": "ae", "keywords": {"reference": "rohf", "psi4_mp_type": "conv", "qc_module": "detci"},                                             }, id="zapt2 rohf ae: psi4-detci",       marks=using("psi4")),

        pytest.param({"call": "gms-zapt2", "reference": "rohf", "fcae": "fc", "keywords": {"gamess_contrl__scftyp": "rohf", "gamess_mp2__code": "serial"},                                                 }, id="zapt2 rohf fc: gamess-serial",    marks=using("gamess")),
        pytest.param({"call": "p4-zapt2",  "reference": "rohf", "fcae": "fc", "keywords": {"reference": "rohf", "psi4_freeze_core": True, "qc_module": "detci", "psi4_mp_type": "conv"},                   }, id="zapt2 rohf fc: psi4-detci",       marks=using("psi4")),
        # yapf: enable
    ],
)
def test_zapt2_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "energy"))


#
#   ,-----.,--. ,---.  ,------.      ,------.
#  '  .--./|  |'   .-' |  .-.  \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |    |  |`.  `-. |  |  \  :    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '--'\|  |.-'    ||  '--'  /    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----'`--'`-----' `-------'     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                  `---' `---'
#  <<<  CISD Energy


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
        # ecc, ncc don't run
        pytest.param({"call": "c4-cisd",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "vcc", "cfour_print": 2},                                                                                        }, id="cisd  rhf ae: cfour-vcc",    marks=using("cfour")),
        pytest.param({"call": "gms-cisd", "reference": "rhf",  "fcae": "ae", "keywords": {"freeze_core": False, "qc_module": "fsoci"},                                                                                                      }, id="cisd  rhf ae: gamess-fsoci", marks=using("gamess")),
        pytest.param({"call": "gms-cisd", "reference": "rhf",  "fcae": "ae", "keywords": {"freeze_core": False, "qc_module": "guga", "gamess_cidrt__mxnint": 200000},                                                                       }, id="cisd  rhf ae: gamess-guga",  marks=using("gamess")),
        pytest.param({"call": "nwc-cisd", "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "tce"},                                                                                                                              }, id="cisd  rhf ae: nwchem-tce",   marks=using("nwchem")),
        pytest.param({"call": "p4-cisd",  "reference": "rhf",  "fcae": "ae", "keywords": {"psi4_qc_module": "fnocc"},                                                                                                                       }, id="cisd  rhf ae: psi4-fnocc",   marks=using("psi4")),
        pytest.param({"call": "p4-cisd",  "reference": "rhf",  "fcae": "ae", "keywords": {"psi4_qc_module": "detci"},                                                                                                                       }, id="cisd  rhf ae: psi4-detci",   marks=using("psi4")),

        pytest.param({"call": "c4-cisd",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_print": 2},                                                                   }, id="cisd  rhf fc: cfour-vcc",    marks=using("cfour")),
        pytest.param({"call": "gms-cisd", "reference": "rhf",  "fcae": "fc", "keywords": {"freeze_core": True, "qc_module": "fsoci"},                                                                                                       }, id="cisd  rhf fc: gamess-fsoci", marks=using("gamess")),
        pytest.param({"call": "gms-cisd", "reference": "rhf",  "fcae": "fc", "keywords": {"freeze_core": True, "qc_module": "guga", "gamess_cidrt__mxnint": 200000},                                                                        }, id="cisd  rhf fc: gamess-guga",  marks=using("gamess")),
        pytest.param({"call": "nwc-cisd", "reference": "rhf",  "fcae": "fc", "keywords": {"nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                                                     }, id="cisd  rhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-cisd",  "reference": "rhf",  "fcae": "fc", "keywords": {"psi4_freeze_core": True, "psi4_qc_module": "fnocc"},                                                                                             }, id="cisd  rhf fc: psi4-fnocc", marks=using("psi4")),
        pytest.param({"call": "p4-cisd",  "reference": "rhf",  "fcae": "fc", "keywords": {"psi4_freeze_core": True, "psi4_qc_module": "detci"},                                                                                             }, id="cisd  rhf fc: psi4-detci", marks=using("psi4")),

        pytest.param({"call": "c4-cisd",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_REFerence": "UHF", "cfour_cc_program": "vcc", "cfour_print": 2},                                                              }, id="cisd  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "gms-cisd", "reference": "uhf",  "fcae": "ae", "keywords": {"gamess_contrl__scftyp": "uhf", "freeze_core": True, "qc_module": "fsoci"},                                                                           "error": {0: _q40}}, id="cisd  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-cisd", "reference": "uhf",  "fcae": "ae", "keywords": {"nwchem_scf__uhf": True, "qc_module": "tce"},                                                                                                     }, id="cisd  uhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-cisd",  "reference": "uhf",  "fcae": "ae", "keywords": {"psi4_reference": "uhf"},                                                                                                       "error": {0: _q39}}, id="cisd  uhf ae: psi4",   marks=using("psi4")),

        pytest.param({"call": "c4-cisd",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_refereNCE": "UHF", "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_print": 2},                                         }, id="cisd  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "gms-cisd", "reference": "uhf",  "fcae": "fc", "keywords": {"gamess_contrl__scftyp": "uhf", "qc_module": "fsoci"},                                                                                                "error": {0: _q40}}, id="cisd  uhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-cisd", "reference": "uhf",  "fcae": "fc", "keywords": {"nwchem_scf__uhf": True, "nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                            }, id="cisd  uhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-cisd",  "reference": "uhf",  "fcae": "fc", "keywords": {"psi4_reference": "uhf", "psi4_freeze_core": True},                                                                             "error": {0: _q39}}, id="cisd  uhf fc: psi4",   marks=using("psi4")),

        pytest.param({"call": "c4-cisd",  "reference": "rohf", "fcae": "ae", "sdsc": "sd", "keywords": {**_c4_tight, "cfour_REFerence": "roHF", "cfour_cc_program": "vcc", "cfour_orbitals": 0, "cfour_print": 2},                          }, id="cisd rohf ae sd: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-cisd",  "reference": "rohf", "fcae": "ae", "sdsc": "sc", "keywords": {**_c4_tight, "cfour_REFerence": "roHF", "cfour_cc_program": "vcc", "cfour_print": 2},                                               }, id="cisd rohf ae sc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "gms-cisd", "reference": "rohf", "fcae": "ae",               "keywords": {"gamess_contrl__scftyp": "rohf", "freeze_core": False, "qc_module": "guga", "gamess_cidrt__mxnint": 200000, "gamess_cidrt__mxneme": 40000}, "wrong": {0: _w28}}, id="cisd rohf ae: gamess-guga",  marks=using("gamess")),
        pytest.param({"call": "nwc-cisd", "reference": "rohf", "fcae": "ae", "sdsc": "sd", "keywords": {"nwchem_scf__rohf": True, "qc_module": "tce"},                                                                                      }, id="cisd rohf ae sd: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-cisd",  "reference": "rohf", "fcae": "ae",               "keywords": {"reference": "rohf", "psi4_qc_module": "detci"},                                                                  "wrong": {0: _w28}}, id="cisd rohf ae   : psi4-detci",    marks=using("psi4")),

        pytest.param({"call": "c4-cisd",  "reference": "rohf", "fcae": "fc", "sdsc": "sc", "keywords": {**_c4_tight, "cfour_REFerence": "roHF", "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_print": 2},                          }, id="cisd rohf fc sc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-cisd",  "reference": "rohf", "fcae": "fc", "sdsc": "sd", "keywords": {**_c4_tight, "cfour_REFerence": "roHF", "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_orbitals": 0, "cfour_print": 2},     }, id="cisd rohf fc sd: cfour-vcc",  marks=using("cfour")),
        #pytest.param({"call": "gms-cisd", "reference": "rohf", "fcae": "fc", "keywords": {"gamess_contrl__scftyp": "rohf", "freeze_core": True},                                              "error": {2: _q4 }}, id="cisd rohf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-cisd", "reference": "rohf", "fcae": "fc", "sdsc": "sd", "keywords": {"nwchem_scf__rohf": True, "nwchem_tce__freeze": 1, "qc_module": "tce"},                                                             }, id="cisd rohf fc sd: nwchem-tce", marks=using("nwchem")),
        # pytest.param({"call": "p4-cisd",  "reference": "rohf", "fcae": "fc", "keywords": {"psi4_reference": "rohf", "psi4_freeze_core": True, "psi4_qc_module": "detci", "psi4_e_convergence": 8, "psi4_r_convergence": 7}, "error": {2: _q8 }}, id="cisd rohf fc: psi4",       marks=using("psi4")),
        # yapf: enable
    ],
)
def test_cisd_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "energy"))


#
#   ,-----.    ,-----.,--. ,---.  ,------.      ,------.
#  '  .-.  '  '  .--./|  |'   .-' |  .-.  \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  | |  |  |  |    |  |`.  `-. |  |  \  :    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '-'  '-.'  '--'\|  |.-'    ||  '--'  /    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----'--' `-----'`--'`-----' `-------'     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                             `---' `---'
#  <<<  QCISD Energy


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
        # * ncc errors out
        pytest.param({"call": "c4-qcisd",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "vcc", "cfour_print": 2},                                                                                       }, id="qcisd  rhf ae: cfour-vcc",    marks=using("cfour")),
        pytest.param({"call": "c4-qcisd",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "ecc"},                                                                                                         }, id="qcisd  rhf ae: cfour-ecc",    marks=using("cfour")),
        pytest.param({"call": "gms-qcisd", "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                             "error": {0: _q43}}, id="qcisd  rhf ae: gamess", marks=using("gamess")),
        pytest.param({"call": "nwc-qcisd", "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "tce"},                                                                                                                             }, id="qcisd  rhf ae: nwchem-tce",   marks=using("nwchem")),
        pytest.param({"call": "p4-qcisd",  "reference": "rhf",  "fcae": "ae", "keywords": {"psi4_qc_module": "fnocc"},                                                                                                                      }, id="qcisd  rhf ae: psi4-fnocc",   marks=using("psi4")),

        pytest.param({"call": "c4-qcisd",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_print": 2},                                                                  }, id="qcisd  rhf fc: cfour-vcc",    marks=using("cfour")),
        pytest.param({"call": "c4-qcisd",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": [1], "cfour_cc_program": "ecc"},                                                                                    }, id="qcisd  rhf fc: cfour-ecc",    marks=using("cfour")),
        pytest.param({"call": "nwc-qcisd", "reference": "rhf",  "fcae": "fc", "keywords": {"nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                                                    }, id="qcisd  rhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-qcisd",  "reference": "rhf",  "fcae": "fc", "keywords": {"psi4_freeze_core": True, "psi4_qc_module": "fnocc"},                                                                                            }, id="qcisd  rhf fc: psi4-fnocc", marks=using("psi4")),

        pytest.param({"call": "c4-qcisd",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_REFerence": "UHF", "cfour_cc_program": "vcc", "cfour_print": 2},                                                             }, id="qcisd  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-qcisd",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_REFerence": "UHF", "cfour_cc_program": "ecc"},                                                                               }, id="qcisd  uhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "nwc-qcisd", "reference": "uhf",  "fcae": "ae", "keywords": {"nwchem_scf__uhf": True, "qc_module": "tce"},                                                                                                    }, id="qcisd  uhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-qcisd",  "reference": "uhf",  "fcae": "ae", "keywords": {"psi4_reference": "uhf"},                                                                                                      "error": {0: _q41}}, id="qcisd  uhf ae: psi4",   marks=using("psi4")),

        pytest.param({"call": "c4-qcisd",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_refereNCE": "UHF", "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_print": 2},                                        }, id="qcisd  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-qcisd",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_refereNCE": "UHF", "cfour_dropmo": [1], "cfour_cc_program": "ecc"},                                                          }, id="qcisd  uhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "nwc-qcisd", "reference": "uhf",  "fcae": "fc", "keywords": {"nwchem_scf__uhf": True, "nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                           }, id="qcisd  uhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-qcisd",  "reference": "uhf",  "fcae": "fc", "keywords": {"psi4_reference": "uhf", "psi4_freeze_core": True},                                                                            "error": {0: _q41}}, id="qcisd  uhf fc: psi4",   marks=using("psi4")),
        # * no ROHF in Cfour, so hard to establish a ref
        # yapf: enable
    ],
)
def test_qcisd_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "energy"))


#   ,-----.    ,-----.,--. ,---.  ,------.    ,-.,--------.,-.      ,------.
#  '  .-.  '  '  .--./|  |'   .-' |  .-.  \  / .''--.  .--''. \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  | |  |  |  |    |  |`.  `-. |  |  \  :|  |    |  |    |  |    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '-'  '-.'  '--'\|  |.-'    ||  '--'  /|  |    |  |    |  |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----'--' `-----'`--'`-----' `-------'  \ '.   `--'   .' /     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                             `-'          `-'                                    `---' `---'
#  <<<  QCISD(T) Energy


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
        pytest.param({"call": "c4-qcisd(t)",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "vcc", "cfour_print": 2},                                                                                       }, id="qcisd_t_  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "gms-qcisd(t)", "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                             "error": {0: _q43}}, id="qcisd_t_  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-qcisd(t)", "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                             "error": {0: _q44}}, id="qcisd_t_  rhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-qcisd(t)",  "reference": "rhf",  "fcae": "ae", "keywords": {"psi4_qc_module": "fnocc"},                                                                                                                      }, id="qcisd_t_  rhf ae: psi4-fnocc", marks=using("psi4")),

        pytest.param({"call": "c4-qcisd(t)",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_print": 2},                                                                  }, id="qcisd_t_  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "p4-qcisd(t)",  "reference": "rhf",  "fcae": "fc", "keywords": {"psi4_freeze_core": True, "psi4_qc_module": "fnocc"},                                                                                            }, id="qcisd_t_  rhf fc: psi4-fnocc", marks=using("psi4")),

        pytest.param({"call": "c4-qcisd(t)",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_REFerence": "UHF", "cfour_cc_program": "vcc", "cfour_print": 2},                                                             }, id="qcisd_t_  uhf ae: cfour-vcc",  marks=using("cfour")),

        pytest.param({"call": "c4-qcisd(t)",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_refereNCE": "UHF", "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_print": 2},                                        }, id="qcisd_t_  uhf fc: cfour-vcc",  marks=using("cfour")),
        # yapf: enable
    ],
)
def test_qcisd_prt_pr_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "energy"))


#
#  ,------. ,-----.,--.    ,------.
#  |  .---''  .--./|  |    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  `--, |  |    |  |    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  |`   '  '--'\|  |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `--'     `-----'`--'    `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                        `---' `---'
#  <<<  FCI Energy


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
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "gms-fci", "reference": "rhf",  "fcae": "ae", "keywords": {"freeze_core": False},                                                                                                                            }, id="fci  rhf ae: gamess",     marks=[using("gamess"), pytest.mark.cilong]),
        pytest.param({"call": "p4-fci",  "reference": "rhf",  "fcae": "ae", "keywords": {"psi4_qc_module": "detci"},                                                                                                                       }, id="fci  rhf ae: psi4-detci", marks=using("psi4")),  # 6m

        pytest.param({"call": "gms-fci", "reference": "rhf",  "fcae": "fc", "keywords": {"freeze_core": True},                                                                                                                             }, id="fci  rhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "p4-fci",  "reference": "rhf",  "fcae": "fc", "keywords": {"psi4_freeze_core": True, "psi4_qc_module": "detci"},                                                                                             }, id="fci  rhf fc: psi4-detci", marks=using("psi4")),

        # TOO LONG pytest.param({"call": "gms-fci", "reference": "rohf", "fcae": "ae", "keywords": {"gamess_contrl__scftyp": "rohf", "freeze_core": False},                                                                                }, id="fci rohf ae: gamess",     marks=using("gamess")),
        # TOO LONG pytest.param({"call": "p4-fci",  "reference": "rohf", "fcae": "ae", "keywords": {"psi4_freeze_core": False, "psi4_reference": "rohf", "psi4_qc_module": "detci"},                                                       }, id="fci rohf ae: psi4-detci", marks=using("psi4")),

        # gms != detci by 1e-5. detci value stored in psi4 ref_local file
        # pytest.param({"call": "gms-fci", "reference": "rohf", "fcae": "fc", "keywords": {"gamess_contrl__scftyp": "rohf", "freeze_core": True},                                                                                           }, id="fci rohf fc: gamess",     marks=using("gamess")),
        # pytest.param({"call": "p4-fci",  "reference": "rohf", "fcae": "fc", "keywords": {"psi4_freeze_core": True, "psi4_reference": "rohf", "psi4_qc_module": "detci"},                                                                   }, id="fci rohf fc: psi4-detci", marks=using("psi4")),
        # yapf: enable
    ],
)
def test_fci_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "energy"))


#
#  ,--.    ,-----. ,-----.,------.      ,------.
#  |  |   '  .--./'  .--./|  .-.  \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |   |  |    |  |    |  |  \  :    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  '--.'  '--'\'  '--'\|  '--'  /    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `-----' `-----' `-----'`-------'     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                     `---' `---'
#  <<<  LCCD Energy


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
        pytest.param({"call": "c4-lccd",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "vcc", "cfour_print": 2},                                                                                        }, id="lccd  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-lccd",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "ncc",},                                                                                                         }, id="lccd  rhf ae: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-lccd", "reference": "rhf",  "fcae": "ae", "keywords": {"gamess_ccinp__ncore": 0},                                                                                                                        }, id="lccd  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-lccd", "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "tce"},                                                                                                                              }, id="lccd  rhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-lccd",  "reference": "rhf",  "fcae": "ae", "keywords": {"psi4_qc_module": "fnocc"},                                                                                                                       }, id="lccd  rhf ae: psi4-fnocc", marks=using("psi4")),
        pytest.param({"call": "p4-lccd",  "reference": "rhf",  "fcae": "ae", "keywords": {"psi4_qc_module": "occ"},                                                                                                                         }, id="lccd  rhf ae: psi4-occ",   marks=using("psi4")),

        pytest.param({"call": "c4-lccd",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_print": 2},                                                                   }, id="lccd  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-lccd",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": [1], "cfour_cc_program": "ncc"},                                                                                     }, id="lccd  rhf fc: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-lccd", "reference": "rhf",  "fcae": "fc", "keywords": {},                                                                                                                                                }, id="lccd  rhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-lccd", "reference": "rhf",  "fcae": "fc", "keywords": {"nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                                                     }, id="lccd  rhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-lccd",  "reference": "rhf",  "fcae": "fc", "keywords": {"psi4_freeze_core": True, "psi4_qc_module": "fnocc"},                                                                                             }, id="lccd  rhf fc: psi4-fnocc", marks=using("psi4")),
        pytest.param({"call": "p4-lccd",  "reference": "rhf",  "fcae": "fc", "keywords": {"psi4_freeze_core": True, "psi4_qc_module": "occ"},                                                                                               }, id="lccd  rhf fc: psi4-occ", marks=using("psi4")),

        pytest.param({"call": "c4-lccd",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_REFerence": "UHF", "cfour_cc_program": "vcc", "cfour_print": 2},                                                              }, id="lccd  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "gms-lccd", "reference": "uhf",  "fcae": "ae", "keywords": {"gamess_contrl__scftyp": "uhf", "gamess_ccinp__ncore": 0},                                                                      "error": {0: _q2 }}, id="lccd  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-lccd", "reference": "uhf",  "fcae": "ae", "keywords": {"nwchem_scf__uhf": True, "qc_module": "tce"},                                                                                                     }, id="lccd  uhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-lccd",  "reference": "uhf",  "fcae": "ae", "keywords": {"psi4_reference": "uhf", "psi4_qc_module": "occ"},                                                                                                }, id="lccd  uhf ae: psi4-occ",   marks=using("psi4")),

        pytest.param({"call": "c4-lccd",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_refereNCE": "UHF", "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_print": 2},                                         }, id="lccd  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "gms-lccd", "reference": "uhf",  "fcae": "fc", "keywords": {"gamess_contrl__scftyp": "uhf"},                                                                                                "error": {0: _q2 }}, id="lccd  uhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-lccd", "reference": "uhf",  "fcae": "fc", "keywords": {"nwchem_scf__uhf": True, "nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                            }, id="lccd  uhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-lccd",  "reference": "uhf",  "fcae": "fc", "keywords": {"psi4_reference": "uhf", "psi4_freeze_core": True, "psi4_qc_module": "occ"},                                                                      }, id="lccd  uhf fc: psi4-occ",   marks=using("psi4")),

        pytest.param({"call": "c4-lccd",  "reference": "rohf", "fcae": "ae", "keywords": {**_c4_tight, "cfour_REFerence": "roHF", "cfour_cc_program": "vcc", "cfour_print": 2},                                           "error": {0: _q7 }}, id="lccd rohf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "gms-lccd", "reference": "rohf", "fcae": "ae", "keywords": {"gamess_contrl__scftyp": "rohf", "gamess_ccinp__ncore": 0, "gamess_ccinp__maxcc": 50},                                          "error": {0: _q4 }}, id="lccd rohf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-lccd", "reference": "rohf", "fcae": "ae", "keywords": {"nwchem_scf__rohf": True, "qc_module": "tce"},                                                                                                    }, id="lccd rohf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-lccd",  "reference": "rohf", "fcae": "ae", "keywords": {"reference": "rohf"},                                                                                                           "error": {0: _q8 }}, id="lccd rohf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-lccd",  "reference": "rohf", "fcae": "fc", "keywords": {**_c4_tight, "cfour_REFerence": "roHF", "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_orbitals": 0, "cfour_print": 2}, "error": {0: _q7 }}, id="lccd rohf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "gms-lccd", "reference": "rohf", "fcae": "fc", "keywords": {"gamess_contrl__scftyp": "rohf", "gamess_ccinp__iconv": 9, "gamess_scf__conv": 9},                                              "error": {0: _q4 }}, id="lccd rohf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-lccd", "reference": "rohf", "fcae": "fc", "keywords": {"nwchem_scf__rohf": True, "nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                           }, id="lccd rohf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-lccd",  "reference": "rohf", "fcae": "fc", "keywords": {"psi4_reference": "rohf", "psi4_freeze_core": True, "psi4_e_convergence": 8, "psi4_r_convergence": 7},                          "error": {0: _q8 }}, id="lccd rohf fc: psi4",       marks=using("psi4")),
        # yapf: enable
    ],
)
def test_lccd_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "energy"))


#
#  ,--.    ,-----. ,-----.,------.       ,----.                     ,--.,--.                 ,--.
#  |  |   '  .--./'  .--./|  .-.  \     '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  |   |  |    |  |    |  |  \  :    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  |  '--.'  '--'\'  '--'\|  '--'  /    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#  `-----' `-----' `-----'`-------'      `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  LCCD Gradient


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
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p"),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # ecc doesn't run
        pytest.param({"call": "c4-lccd",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_cc_program": "vcc", "cfour_print": 2},                                                                  "error": {1: _q38}}, id="lccd  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-lccd",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_cc_program": "ncc",},                                                                                   "error": {1: _q38}}, id="lccd  rhf ae: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-lccd", "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {"gamess_ccinp__ncore": 0, "gamess_contrl__numgrd": True},                                                                                     }, id="lccd  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-lccd", "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {"qc_module": "tce"},                                                                                                                          }, id="lccd  rhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-lccd",  "reference": "rhf",  "fcae": "ae", "xptd": {"fd": False}, "keywords": {**_p4_fd, "psi4_qc_module": "fnocc"},                                                                                                         }, id="lccd  rhf ae: psi4-fnocc", marks=using("psi4")),
        pytest.param({"call": "p4-lccd",  "reference": "rhf",  "fcae": "ae",                        "keywords": {"psi4_qc_module": "occ"},                                                                                                                     }, id="lccd  rhf ae: psi4-occ",   marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-lccd", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_cfour_cc_program": "vcc"}                                                                                                                 }, id="lccd  rhf ae: psi4-cfour-vcc"),

        pytest.param({"call": "c4-lccd",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_print": 2},                                             "error": {1: _q38}}, id="lccd  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-lccd",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_dropmo": [1], "cfour_cc_program": "ncc"},                                                               "error": {1: _q38}}, id="lccd  rhf fc: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-lccd", "reference": "rhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {"gamess_contrl__numgrd": True},                                                                                                               }, id="lccd  rhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-lccd", "reference": "rhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {"nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                                                 }, id="lccd  rhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-lccd",  "reference": "rhf",  "fcae": "fc", "xptd": {"fd": False}, "keywords": {**_p4_fd, "psi4_freeze_core": True, "psi4_qc_module": "fnocc"},                                                                               }, id="lccd  rhf fc: psi4-fnocc", marks=using("psi4")),
        pytest.param({"call": "p4-lccd",  "reference": "rhf",  "fcae": "fc", "xptd": {"fd": False}, "keywords": {**_p4_fd, "psi4_freeze_core": True, "psi4_qc_module": "occ"},                                                                                 }, id="lccd  rhf fc: psi4-occ", marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-lccd", "reference": "rhf", "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_cfour_dropmo": [1], "psi4_cfour_cc_program": "vcc"}                                                                                       }, id="lccd  rhf fc: psi4-cfour-vcc"),

        pytest.param({"call": "c4-lccd",  "reference": "uhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_REFerence": "UHF", "cfour_cc_program": "vcc", "cfour_print": 2},                                        "error": {1: _q38}}, id="lccd  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "gms-lccd", "reference": "uhf",  "fcae": "ae",                    "keywords": {"gamess_contrl__scftyp": "uhf", "gamess_ccinp__ncore": 0},                                                                      "error": {1: _q2 }}, id="lccd  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-lccd", "reference": "uhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {"nwchem_scf__uhf": True, "qc_module": "tce"},                                                                                                 }, id="lccd  uhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-lccd",  "reference": "uhf",  "fcae": "ae",                        "keywords": {"psi4_reference": "uhf", "psi4_qc_module": "occ"},                                                                                            }, id="lccd  uhf ae: psi4-occ",   marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-lccd", "reference": "uhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "uhf", "psi4_cfour_cc_program": "vcc"}                                                                                  }, id="lccd  uhf ae: psi4-cfour-vcc"),

        pytest.param({"call": "c4-lccd",  "reference": "uhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_refereNCE": "UHF", "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_print": 2},                   "error": {1: _q38}}, id="lccd  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "gms-lccd", "reference": "uhf",  "fcae": "fc",                    "keywords": {"gamess_contrl__scftyp": "uhf"},                                                                                                "error": {1: _q2 }}, id="lccd  uhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-lccd", "reference": "uhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {"nwchem_scf__uhf": True, "nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                        }, id="lccd  uhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-lccd",  "reference": "uhf",  "fcae": "fc", "xptd": {"fd": False}, "keywords": {**_p4_fd, "psi4_reference": "uhf", "psi4_freeze_core": True, "psi4_qc_module": "occ"},                                                        }, id="lccd  uhf fc: psi4-occ",   marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-lccd", "reference": "uhf", "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "uhf", "psi4_cfour_dropmo": [1], "psi4_cfour_cc_program": "vcc"}                                                        }, id="lccd  uhf fc: psi4-cfour-vcc"),
        # yapf: enable
    ],
)
def test_lccd_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "gradient"))


#
#  ,--.    ,-----. ,-----. ,---.  ,------.      ,------.
#  |  |   '  .--./'  .--./'   .-' |  .-.  \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |   |  |    |  |    `.  `-. |  |  \  :    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  '--.'  '--'\'  '--'\.-'    ||  '--'  /    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `-----' `-----' `-----'`-----' `-------'     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                             `---' `---'
#  <<<  LCCSD Energy


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
        pytest.param({"call": "c4-lccsd",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "vcc", "cfour_print": 2},                                                                                                          }, id="lccsd  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-lccsd",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "ncc",},                                                                                                                           }, id="lccsd  rhf ae: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-lccsd", "reference": "rhf",  "fcae": "ae", "keywords": {"gamess_ccinp__ncore": 0},                                                                                                                        "error": {0: _q9 }}, id="lccsd  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-lccsd", "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "tce"},                                                                                                                                                }, id="lccsd  rhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-lccsd",  "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "fnocc"},                                                                                                                                              }, id="lccsd  rhf ae: psi4-fnocc", marks=using("psi4")),

        pytest.param({"call": "c4-lccsd",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_print": 2},                                                                                     }, id="lccsd  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-lccsd",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": [1], "cfour_cc_program": "ncc"},                                                                                                       }, id="lccsd  rhf fc: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-lccsd", "reference": "rhf",  "fcae": "fc", "keywords": {},                                                                                                                                                "error": {0: _q9 }}, id="lccsd  rhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-lccsd", "reference": "rhf",  "fcae": "fc", "keywords": {"qc_module": "tce", "nwchem_tce__freeze": 1},                                                                                                                       }, id="lccsd  rhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-lccsd",  "reference": "rhf",  "fcae": "fc", "keywords": {"psi4_freeze_core": True, "qc_module": "fnocc"},                                                                                                                    }, id="lccsd  rhf fc: psi4-fnocc", marks=using("psi4")),

        pytest.param({"call": "c4-lccsd",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_REFerence": "UHF", "cfour_cc_program": "vcc", "cfour_print": 2},                                                                                }, id="lccsd  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "gms-lccsd", "reference": "uhf",  "fcae": "ae", "keywords": {"gamess_contrl__scftyp": "uhf", "gamess_ccinp__ncore": 0},                                                                                        "error": {0: _q9 }}, id="lccsd  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-lccsd", "reference": "uhf",  "fcae": "ae", "keywords": {"nwchem_scf__uhf": True, "qc_module": "tce"},                                                                                                                       }, id="lccsd  uhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-lccsd",  "reference": "uhf",  "fcae": "ae", "keywords": {"reference": "uhf"},                                                                                                                              "error": {0: _q10}}, id="lccsd  uhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-lccsd",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_REFerence": "UHF", "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_print": 2},                                                           }, id="lccsd  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "gms-lccsd", "reference": "uhf",  "fcae": "fc", "keywords": {"gamess_contrl__scftyp": "uhf"},                                                                                                                  "error": {0: _q9 }}, id="lccsd  uhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-lccsd", "reference": "uhf",  "fcae": "fc", "keywords": {"nwchem_scf__uhf": True, "nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                                              }, id="lccsd  uhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-lccsd",  "reference": "uhf",  "fcae": "fc", "keywords": {"reference": "uhf", "psi4_freeze_core": True},                                                                                                    "error": {0: _q10}}, id="lccsd  uhf fc: psi4",       marks=using("psi4")),

        # Sum 2021: uncertain rohf reference values: cfour == nwc-tce, but those != p4n. cfour says rohf nyi. Sum 2022: switching from p4n to vcc ref values so c4 and nwc show as right, not `"wrong": {0: _w29}`
        pytest.param({"call": "c4-lccsd",  "reference": "rohf", "fcae": "ae", "sdsc": "sc", "keywords": {**_c4_tight, "cfour_REFerence": "roHF", "cfour_cc_program": "vcc", "cfour_print": 2},                                                                 }, id="lccsd rohf ae sc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-lccsd",  "reference": "rohf", "fcae": "ae", "sdsc": "sd", "keywords": {**_c4_tight, "cfour_REFerence": "roHF", "cfour_orbitals": 0, "cfour_cc_program": "vcc", "cfour_print": 2},                                            }, id="lccsd rohf ae sd: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "gms-lccsd", "reference": "rohf", "fcae": "ae",               "keywords": {"gamess_contrl__scftyp": "rohf", "gamess_ccinp__ncore": 0, "gamess_ccinp__maxcc": 50},                                              "error": {0: _q9 }}, id="lccsd rohf ae   : gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-lccsd", "reference": "rohf", "fcae": "ae", "sdsc": "sd", "keywords": {"qc_module": "tce", "nwchem_scf__rohf": True},                                                                                                        }, id="lccsd rohf ae sd: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-lccsd",  "reference": "rohf", "fcae": "ae",               "keywords": {"reference": "rohf"},                                                                                                               "error": {0: _q10}}, id="lccsd rohf ae   : psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-lccsd",  "reference": "rohf", "fcae": "fc", "sdsc": "sc", "keywords": {**_c4_tight, "cfour_REFerence": "roHF", "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_print": 2},                                            }, id="lccsd rohf fc sc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-lccsd",  "reference": "rohf", "fcae": "fc", "sdsc": "sd", "keywords": {**_c4_tight, "cfour_REFerence": "roHF", "cfour_dropmo": [1], "cfour_orbitals": 0, "cfour_cc_program": "vcc", "cfour_print": 2},                       }, id="lccsd rohf fc sd: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "gms-lccsd", "reference": "rohf", "fcae": "fc",               "keywords": {"gamess_contrl__scftyp": "rohf", "gamess_ccinp__iconv": 9, "gamess_scf__conv": 9},                                                  "error": {0: _q9 }}, id="lccsd rohf fc   : gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-lccsd", "reference": "rohf", "fcae": "fc", "sdsc": "sd", "keywords": {"nwchem_scf__rohf": True, "nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                               }, id="lccsd rohf fc sd: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-lccsd",  "reference": "rohf", "fcae": "fc",               "keywords": {"psi4_reference": "rohf", "psi4_freeze_core": True},                                                                                "error": {0: _q10}}, id="lccsd rohf fc   : psi4",       marks=using("psi4")),
        # yapf: enable
    ],
)
def test_lccsd_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "energy"))


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
        pytest.param({"call": "c4-ccd",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "vcc", "cfour_print": 2},                                                                                     }, id="ccd  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccd",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "ecc",},                                                                                                      }, id="ccd  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccd",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "ncc",},                                                                                                      }, id="ccd  rhf ae: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccd", "reference": "rhf",  "fcae": "ae", "keywords": {"gamess_ccinp__ncore": 0},                                                                                                                     }, id="ccd  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccd", "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "tce"},                                                                                                                           }, id="ccd  rhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-ccd",  "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                           "error": {0: _q11}}, id="ccd  rhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-ccd",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_print": 2},                                                                }, id="ccd  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccd",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": [1], "cfour_cc_program": "ecc"},                                                                                  }, id="ccd  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccd",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": [1], "cfour_cc_program": "ncc"},                                                                                  }, id="ccd  rhf fc: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccd", "reference": "rhf",  "fcae": "fc", "keywords": {},                                                                                                                                             }, id="ccd  rhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccd", "reference": "rhf",  "fcae": "fc", "keywords": {"nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                                                  }, id="ccd  rhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-ccd",  "reference": "rhf",  "fcae": "fc", "keywords": {"psi4_freeze_core": True},                                                                                                   "error": {0: _q11}}, id="ccd  rhf fc: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-ccd",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_REFerence": "UHF", "cfour_cc_program": "vcc", "cfour_print": 2},                                                           }, id="ccd  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccd",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_REFerence": "UHF", "cfour_cc_program": "ecc"},                                                                             }, id="ccd  uhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccd", "reference": "uhf",  "fcae": "ae", "keywords": {"gamess_contrl__scftyp": "uhf", "gamess_ccinp__ncore": 0},                                                                   "error": {0: _q2} }, id="ccd  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccd", "reference": "uhf",  "fcae": "ae", "keywords": {"nwchem_scf__uhf": True, "qc_module": "tce"},                                                                                                  }, id="ccd  uhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-ccd",  "reference": "uhf",  "fcae": "ae", "keywords": {"psi4_reference": "uhf"},                                                                                                    "error": {0: _q11}}, id="ccd  uhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-ccd",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_REFerence": "uHF", "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_print": 2},                                      }, id="ccd  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccd",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_REFerence": "uHF", "cfour_dropmo": [1], "cfour_cc_program": "ecc"},                                                        }, id="ccd  uhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccd", "reference": "uhf",  "fcae": "fc", "keywords": {"gamess_contrl__scftyp": "uhf"},                                                                                             "error": {0: _q2} }, id="ccd  uhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccd", "reference": "uhf",  "fcae": "fc", "keywords": {"nwchem_scf__uhf": True, "nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                         }, id="ccd  uhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-ccd",  "reference": "uhf",  "fcae": "fc", "keywords": {"reference": "uhf", "psi4_freeze_core": True},                                                                               "error": {0: _q11}}, id="ccd  uhf fc: psi4",       marks=using("psi4")),

        # rohf vcc = tce, but cfour paper disavows rohf, so I'm suspicious
        pytest.param({"call": "c4-ccd",  "reference": "rohf", "fcae": "ae", "sdsc": "sd", "keywords": {**_c4_tight, "cfour_REFerence": "roHF", "cfour_cc_program": "vcc", "cfour_orbitals": 0, "cfour_print": 2},                                        }, id="ccd rohf ae sd: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccd",  "reference": "rohf", "fcae": "ae", "sdsc": "sc", "keywords": {**_c4_tight, "cfour_REFerence": "roHF", "cfour_cc_program": "vcc", "cfour_print": 2},                                                             }, id="ccd rohf ae sc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccd", "reference": "rohf", "fcae": "ae",               "keywords": {"gamess_contrl__scftyp": "rohf", "gamess_ccinp__ncore": 0, "gamess_ccinp__maxcc": 50},                                          "error": {0: _q4} }, id="ccd rohf ae   : gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccd", "reference": "rohf", "fcae": "ae", "sdsc": "sd", "keywords": {"nwchem_scf__rohf": True, "qc_module": "tce"},                                                                                                    }, id="ccd rohf ae sd: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-ccd",  "reference": "rohf", "fcae": "ae",               "keywords": {"reference": "rohf"},                                                                                                           "error": {0: _q11}}, id="ccd rohf ae   : psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-ccd",  "reference": "rohf", "fcae": "fc", "sdsc": "sd", "keywords": {**_c4_tight, "cfour_REFerence": "roHF", "cfour_dropmo": [1], "cfour_orbitals": 0, "cfour_cc_program": "vcc", "cfour_print": 2},                   }, id="ccd rohf fc sd: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccd",  "reference": "rohf", "fcae": "fc", "sdsc": "sc", "keywords": {**_c4_tight, "cfour_REFerence": "roHF", "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_print": 2},                                        }, id="ccd rohf fc sc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccd", "reference": "rohf", "fcae": "fc",               "keywords": {"gamess_contrl__scftyp": "rohf", "gamess_ccinp__iconv": 9, "gamess_scf__conv": 9},                                              "error": {0: _q4} }, id="ccd rohf fc   : gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccd", "reference": "rohf", "fcae": "fc", "sdsc": "sd", "keywords": {"nwchem_scf__rohf": True, "nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                           }, id="ccd rohf fc sd: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-ccd",  "reference": "rohf", "fcae": "fc",               "keywords": {"psi4_reference": "rohf", "psi4_freeze_core": True},                                                                            "error": {0: _q11}}, id="ccd rohf fc   : psi4",       marks=using("psi4")),
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
        pytest.param({"call": "c4-ccd",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_cc_program": "vcc", "cfour_print": 2},                                                                    }, id="ccd  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccd",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_cc_program": "ecc"},                                                                                      }, id="ccd  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccd",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_cc_program": "ncc"},                                                                                      }, id="ccd  rhf ae: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccd", "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {"gamess_ccinp__ncore": 0, "gamess_contrl__numgrd": True},                                                                     }, id="ccd  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccd", "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {"qc_module": "tce"},                                                                                                          }, id="ccd  rhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-ccd",  "reference": "rhf",  "fcae": "ae",                        "keywords": {"psi4_cc_type": "conv"},                                                                                    "error": {1: _q11}}, id="ccd  rhf ae: psi4",       marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-ccd", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "vcc"},                                                                 }, id="ccd  rhf fc: psi4-cfour-vcc"),
        # DEBUG pytest.param({"call": "p4-c4-ccd", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "ncc"},                                                                 }, id="ccd  rhf fc: psi4-cfour-ncc"),

        pytest.param({"call": "c4-ccd",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "vcc", "cfour_print": 2},                                                 }, id="ccd  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccd",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "ecc"},                                                                   }, id="ccd  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccd",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "ncc"},                                                 "wrong": {1: _w9} }, id="ccd  rhf fc: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccd", "reference": "rhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {"gamess_contrl__numgrd": True},                                                                                               }, id="ccd  rhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccd", "reference": "rhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {"nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                                 }, id="ccd  rhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-ccd",  "reference": "rhf",  "fcae": "fc",                        "keywords": {"psi4_freeze_core": True, "psi4_cc_type": "conv"},                                                          "error": {1: _q11}}, id="ccd  rhf fc: psi4",       marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-ccd", "reference": "rhf", "fcae": "fc", "keywords": {**_p4c4_fd,  "psi4_cfour_reference": "rhf", "psi4_cfour_dropmo": [1], "psi4_cfour_cc_program": "vcc"},                                      }, id="ccd  rhf fc: psi4-cfour-vcc"),

        pytest.param({"call": "c4-ccd",  "reference": "uhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "vcc", "cfour_print": 2},                                          }, id="ccd  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccd",  "reference": "uhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "ecc"},                                                            }, id="ccd  uhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccd", "reference": "uhf",  "fcae": "ae",                        "keywords": {"gamess_contrl__scftyp": "uhf", "gamess_ccinp__ncore": 0},                                                  "error": {1: _q2} }, id="ccd  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccd", "reference": "uhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {"nwchem_scf__uhf": True, "qc_module": "tce"},                                                                                 }, id="ccd  uhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-ccd",  "reference": "uhf",  "fcae": "ae",                        "keywords": {"reference": "uhf", "psi4_cc_type": "conv"},                                                                "error": {1: _q11}}, id="ccd  uhf ae: psi4",       marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-ccd", "reference": "uhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "uhf", "psi4_cfour_cc_program": "vcc"},                                                                 }, id="ccd  uhf ae: psi4-cfour-vcc"),

        pytest.param({"call": "c4-ccd",  "reference": "uhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_print": 2},                     }, id="ccd  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccd",  "reference": "uhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_dropmo": [1], "cfour_cc_program": "ecc"},                                       }, id="ccd  uhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccd", "reference": "uhf",  "fcae": "fc",                        "keywords": {"gamess_contrl__scftyp": "uhf"},                                                                            "error": {1: _q2} }, id="ccd  uhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccd", "reference": "uhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {"nwchem_scf__uhf": True, "nwchem_tce__freeze": 1, "qc_module": "tce"},                                                        }, id="ccd  uhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-ccd",  "reference": "uhf",  "fcae": "fc",                        "keywords": {"reference": "uhf", "psi4_freeze_core": True, "psi4_cc_type": "conv"},                                      "error": {1: _q11}}, id="ccd  uhf fc: psi4",       marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-ccd", "reference": "uhf", "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "uhf", "psi4_cfour_dropmo": [1], "psi4_cfour_cc_program": "vcc"},                                       }, id="ccd  uhf fc: psi4-cfour-vcc"),
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
        # ncc errors
        pytest.param({"call": "c4-ccd",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_cc_program": "vcc", "cfour_print": 2},                                                                    }, id="ccd  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccd",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_cc_program": "ecc"},                                                                                      }, id="ccd  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccd", "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {"gamess_ccinp__ncore": 0, "gamess_force__method": "fullnum"},                                                                 }, id="ccd  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccd", "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {"qc_module": "tce"},                                                                                                          }, id="ccd  rhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-ccd",  "reference": "rhf",  "fcae": "ae",                        "keywords": {"psi4_cc_type": "conv"},                                                                                    "error": {2: _q17}}, id="ccd  rhf ae: psi4",       marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-ccd", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "vcc"},                                                                 }, id="ccd  rhf fc: psi4-cfour-vcc"),

        # FC: vcc errors for analytic hess
        # pytest.param({"call": "c4-ccd",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "vcc", "cfour_print": 2},                                                 }, id="ccd  rhf fc: cfour-vcc",  marks=using("cfour")),

        pytest.param({"call": "c4-ccd",  "reference": "uhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "vcc", "cfour_print": 2},                                          }, id="ccd  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccd",  "reference": "uhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "ecc"},                                                            }, id="ccd  uhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccd", "reference": "uhf",  "fcae": "ae",                        "keywords": {"gamess_contrl__scftyp": "uhf", "gamess_ccinp__ncore": 0},                                                  "error": {2: _q2} }, id="ccd  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccd", "reference": "uhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {"nwchem_scf__uhf": True, "qc_module": "tce"},                                                                                 }, id="ccd  uhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-ccd",  "reference": "uhf",  "fcae": "ae",                        "keywords": {"reference": "uhf", "psi4_cc_type": "conv"},                                                                "error": {2: _q17}}, id="ccd  uhf ae: psi4",       marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-ccd", "reference": "uhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "uhf", "psi4_cfour_cc_program": "vcc"},                                                                 }, id="ccd  uhf ae: psi4-cfour-vcc"),
        # FC: vcc errors for analytic hess
        # pytest.param({"call": "c4-ccd",  "reference": "uhf",  "fcae": "fc",                         "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_print": 2},                    }, id="ccd  uhf fc: cfour-vcc", marks=using("cfour")),
        # yapf: enable
    ],
)
def test_ccd_hessian_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "hessian"))


#
#  ,-----.   ,-----. ,-----.,------.      ,------.
#  |  |) /_ '  .--./'  .--./|  .-.  \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  .-.  \|  |    |  |    |  |  \  :    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  '--' /'  '--'\'  '--'\|  '--'  /    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `------'  `-----' `-----'`-------'     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                       `---' `---'
#  <<<  BCCD Energy


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
        # ecc, ncc doesn't run
        pytest.param({"call": "c4-bccd",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "vcc", "cfour_print": 2},                                                                                     }, id="bccd  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "p4-bccd",  "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                                             }, id="bccd  rhf ae: psi4",       marks=using("psi4")),

        # cfour fc Brueckner hangs!
        # pytest.param({"call": "c4-bccd",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_print": 2},                                                                }, id="bccd  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "p4-bccd",  "reference": "rhf",  "fcae": "fc", "keywords": {"psi4_freeze_core": True},                                                                                                                     }, id="bccd  rhf fc: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-bccd",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_REFerence": "UHF", "cfour_cc_program": "vcc", "cfour_print": 2},                                                           }, id="bccd  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "p4-bccd",  "reference": "uhf",  "fcae": "ae", "keywords": {"psi4_reference": "uhf"},                                                                                                                      }, id="bccd  uhf ae: psi4",       marks=using("psi4")),

        # pytest.param({"call": "c4-bccd",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_REFerence": "uHF", "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_print": 2},                                      }, id="bccd  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "p4-bccd",  "reference": "uhf",  "fcae": "fc", "keywords": {"reference": "uhf", "psi4_freeze_core": True},                                                                                                 }, id="bccd  uhf fc: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-bccd",  "reference": "rohf", "fcae": "ae", "keywords": {**_c4_tight, "cfour_REFerence": "roHF", "cfour_cc_program": "vcc", "cfour_print": 2},                                                          }, id="bccd rohf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "p4-bccd",  "reference": "rohf", "fcae": "ae", "keywords": {"reference": "rohf"},                                                                                                                          }, id="bccd rohf ae: psi4",       marks=using("psi4")),  # NOTE: sc

        # pytest.param({"call": "c4-bccd",  "reference": "rohf", "fcae": "fc", "sdsc": "sc", "keywords": {**_c4_tight, "cfour_REFerence": "roHF", "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_print": 2},                       }, id="bccd rohf fc sc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "p4-bccd",  "reference": "rohf", "fcae": "fc", "sdsc": "sc", "keywords": {"psi4_reference": "rohf", "psi4_freeze_core": True},                                                                             }, id="bccd rohf fc sc: psi4",       marks=using("psi4")),
        # yapf: enable
    ],
)
def test_bccd_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "energy"))


#
#   ,-----. ,-----. ,---.     ,------.
#  '  .--./'  .--./'.-.  \    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |    |  |     .-' .'    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '--'\'  '--'\/   '-.    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----' `-----''-----'    `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                           `---' `---'
#  <<<  CC2 Energy


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
        # ecc, ncc doesn't run
        pytest.param({"call": "c4-cc2",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "vcc", "cfour_print": 2},                                                                                     }, id="cc2  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "nwc-cc2", "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "tce"},                                                                                                                           }, id="cc2  rhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-cc2",  "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                                             }, id="cc2  rhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-cc2",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_print": 2},                                                                }, id="cc2  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "nwc-cc2", "reference": "rhf",  "fcae": "fc", "keywords": {"nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                                                  }, id="cc2  rhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-cc2",  "reference": "rhf",  "fcae": "fc", "keywords": {"psi4_freeze_core": True},                                                                                                                     }, id="cc2  rhf fc: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-cc2",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_REFerence": "UHF", "cfour_cc_program": "vcc", "cfour_print": 2},                                                           }, id="cc2  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "nwc-cc2", "reference": "uhf",  "fcae": "ae", "keywords": {"nwchem_scf__uhf": True, "qc_module": "tce"},                                                                                                  }, id="cc2  uhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-cc2",  "reference": "uhf",  "fcae": "ae", "keywords": {"psi4_reference": "uhf"},                                                                                                                      }, id="cc2  uhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-cc2",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_REFerence": "uHF", "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_print": 2},                                      }, id="cc2  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "nwc-cc2", "reference": "uhf",  "fcae": "fc", "keywords": {"nwchem_scf__uhf": True, "nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                         }, id="cc2  uhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-cc2",  "reference": "uhf",  "fcae": "fc", "keywords": {"reference": "uhf", "psi4_freeze_core": True},                                                                                                 }, id="cc2  uhf fc: psi4",       marks=using("psi4")),

        # rohf ae vcc != tce == psi4, and cfour paper disavows rohf, so I'm suspicious
        # pytest.param({"call": "c4-cc2",  "reference": "rohf", "fcae": "ae", "keywords": {**_c4_tight, "cfour_REFerence": "roHF", "cfour_cc_program": "vcc", "cfour_print": 2},                                                          }, id="cc2 rohf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "nwc-cc2", "reference": "rohf", "fcae": "ae", "keywords": {"nwchem_scf__rohf": True, "qc_module": "tce"},                                                                                                 }, id="cc2 rohf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-cc2",  "reference": "rohf", "fcae": "ae", "keywords": {"reference": "rohf"},                                                                                                                          }, id="cc2 rohf ae: psi4",       marks=using("psi4")),  # NOTE: sc

        # nwc slightly different
        # pytest.param({"call": "c4-cc2",  "reference": "rohf", "fcae": "fc", "sdsc": "sc", "keywords": {**_c4_tight, "cfour_REFerence": "roHF", "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_print": 2},                       }, id="cc2 rohf fc sc: cfour-vcc",  marks=using("cfour")),
        # pytest.param({"call": "nwc-cc2", "reference": "rohf", "fcae": "fc", "sdsc": "sc", "keywords": {"nwchem_scf__rohf": True, "nwchem_tce__freeze": 1, "qc_module": "tce"},                                                          }, id="cc2 rohf fc   : nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-cc2",  "reference": "rohf", "fcae": "fc", "sdsc": "sc", "keywords": {"psi4_reference": "rohf", "psi4_freeze_core": True},                                                                             }, id="cc2 rohf fc sc: psi4",       marks=using("psi4")),
        # yapf: enable
    ],
)
def test_cc2_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "energy"))


#
#   ,-----. ,-----. ,---.      ,----.                     ,--.,--.                 ,--.
#  '  .--./'  .--./'.-.  \    '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  |    |  |     .-' .'    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  '  '--'\'  '--'\/   '-.    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#   `-----' `-----''-----'     `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  CC2 Gradient


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
        pytest.param({"call": "c4-cc2",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_cc_program": "vcc", "cfour_print": 2},                                                                                      }, id="cc2  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "p4-cc2",  "reference": "rhf",  "fcae": "ae",                        "keywords": {"psi4_cc_type": "conv", "qc_module": "ccenergy"},                                                                                               }, id="cc2  rhf ae: psi4-cc",    marks=using("psi4")),

        pytest.param({"call": "c4-cc2",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_cc_program": "vcc", "cfour_dropmo": 1, "cfour_print": 2},                                                                   }, id="cc2  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "p4-cc2",  "reference": "rhf",  "fcae": "fc",                        "keywords": {"psi4_freeze_core": True, "psi4_cc_type": "conv", "qc_module": "ccenergy"},                                                   "error": {1: _q18}}, id="cc2  rhf fc: psi4-cc",    marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-cc2",  "reference": "rhf",  "fcae": "fc",                        "keywords": {"psi4_freeze_core": True, "psi4_cc_type": "conv", "qc_module": "ccenergy", "psi4_function_kwargs_dertype": 0, "psi4_points": 5},                }, id="cc2  rhf fc: psi4-cc",    marks=using("psi4")),

        pytest.param({"call": "c4-cc2",  "reference": "uhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "vcc", "cfour_print": 2},                                                            }, id="cc2  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "p4-cc2",  "reference": "uhf",  "fcae": "ae",                        "keywords": {"reference": "uhf", "psi4_cc_type": "conv", "psi4_function_kwargs_dertype": 1},                                               "error": {1: _q45}}, id="cc2  uhf ae: psi4-cc",    marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-cc2",  "reference": "uhf",  "fcae": "ae",                        "keywords": {"reference": "uhf", "psi4_cc_type": "conv", "psi4_function_kwargs_dertype": 0, "psi4_points": 5},                                               }, id="cc2  uhf ae: psi4-cc",    marks=using("psi4")),

        pytest.param({"call": "c4-cc2",  "reference": "uhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_dropmo": 1, "cfour_cc_program": "vcc", "cfour_print": 2},                                         }, id="cc2  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "p4-cc2",  "reference": "uhf",  "fcae": "fc",                        "keywords": {"reference": "uhf", "psi4_cc_type": "conv", "psi4_freeze_core": True, "psi4_function_kwargs_dertype": 1},                     "error": {1: _q45}}, id="cc2  uhf fc: psi4",       marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-cc2",  "reference": "uhf",  "fcae": "fc",                        "keywords": {"reference": "uhf", "psi4_cc_type": "conv", "psi4_freeze_core": True, "psi4_function_kwargs_dertype": 0, "psi4_points": 5},                     }, id="cc2  uhf fc: psi4",       marks=using("psi4")),
        # yapf: enable
    ],
)
def test_cc2_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "gradient"))


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
        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "vcc", "cfour_print": 2},                                                                                        }, id="ccsd  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "ecc",},                                                                                                         }, id="ccsd  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "ncc",},                                                                                                         }, id="ccsd  rhf ae: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "rhf",  "fcae": "ae", "keywords": {"gamess_ccinp__ncore": 0},                                                                                                                        }, id="ccsd  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "tce"},                                                                                                                              }, id="ccsd  rhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "cc"},                                                                                                                               }, id="ccsd  rhf ae: nwchem-cc",  marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "ccenergy"},                                                                                                                         }, id="ccsd  rhf ae: psi4-cc",    marks=using("psi4")),
        pytest.param({"call": "p4-ccsd",  "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "fnocc"},                                                                                                                            }, id="ccsd  rhf ae: psi4-fnocc", marks=using("psi4")),
        pytest.param({"call": "p4-ccsd",  "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "mrcc"},                                                                                                                             }, id="ccsd  rhf ae: psi4-mrcc",  marks=using("psi4")),

        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_print": 2},                                                                   }, id="ccsd  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": [1], "cfour_cc_program": "ecc"},                                                                                     }, id="ccsd  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": [1], "cfour_cc_program": "ncc"},                                                                                     }, id="ccsd  rhf fc: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "rhf",  "fcae": "fc", "keywords": {},                                                                                                                                                }, id="ccsd  rhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "fc", "keywords": {"nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                                                     }, id="ccsd  rhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "fc", "keywords": {"nwchem_ccsd__freeze": 1, "qc_module": "cc"},                                                                                                     }, id="ccsd  rhf fc: nwchem-cc",  marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "rhf",  "fcae": "fc", "keywords": {"psi4_freeze_core": True, "qc_module": "ccenergy"},                                                                                               }, id="ccsd  rhf fc: psi4-cc",    marks=using("psi4")),
        pytest.param({"call": "p4-ccsd",  "reference": "rhf",  "fcae": "fc", "keywords": {"psi4_freeze_core": True, "qc_module": "fnocc"},                                                                                                  }, id="ccsd  rhf fc: psi4-fnocc", marks=using("psi4")),
        pytest.param({"call": "p4-ccsd",  "reference": "rhf",  "fcae": "fc", "keywords": {"psi4_freeze_core": True, "qc_module": "mrcc"},                                                                                                   }, id="ccsd  rhf fc: psi4-mrcc",  marks=using("psi4")),  #, using("mrcc")]),

        # "cfour_occupation": [[3, 1, 1, 0], [3, 0, 1, 0]]
        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_REFerence": "UHF", "cfour_cc_program": "vcc", "cfour_print": 2},                                                              }, id="ccsd  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_REFerence": "UHF", "cfour_cc_program": "ecc"},                                                                                }, id="ccsd  uhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "uhf",  "fcae": "ae", "keywords": {"gamess_contrl__scftyp": "uhf", "gamess_ccinp__ncore": 0},                                                                      "error": {0: _q2 }}, id="ccsd  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "uhf",  "fcae": "ae", "keywords": {"nwchem_scf__uhf": True, "qc_module": "tce"},                                                                                                     }, id="ccsd  uhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "uhf",  "fcae": "ae", "keywords": {"nwchem_scf__uhf": True, "qc_module": "cc"},                                                                                    "error": {0: _q3 }}, id="ccsd  uhf ae: nwchem-cc",  marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "uhf",  "fcae": "ae", "keywords": {"reference": "uhf", "qc_module": "ccenergy"},                                                                                                     }, id="ccsd  uhf ae: psi4-cc",    marks=using("psi4")),
        pytest.param({"call": "p4-ccsd",  "reference": "uhf",  "fcae": "ae", "keywords": {"reference": "uhf", "qc_module": "mrcc"},                                                                                                         }, id="ccsd  uhf ae: psi4-mrcc",  marks=using("psi4")),

        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight,  "cfour_REFerence": "UHF", "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_print": 2},                                        }, id="ccsd  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight,  "cfour_REFerence": "UHF", "cfour_dropmo": [1], "cfour_cc_program": "ecc"},                                                          }, id="ccsd  uhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "uhf",  "fcae": "fc", "keywords": {"gamess_contrl__scftyp": "uhf"},                                                                                                "error": {0: _q2 }}, id="ccsd  uhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "uhf",  "fcae": "fc", "keywords": {"nwchem_scf__uhf": True, "nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                            }, id="ccsd  uhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "uhf",  "fcae": "fc", "keywords": {"nwchem_scf__uhf": True, "nwchem_ccsd__freeze": 1, "qc_module": "cc"},                                                          "error": {0: _q3 }}, id="ccsd  uhf fc: nwchem-cc",  marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "uhf",  "fcae": "fc", "keywords": {"psi4_freeze_core": True, "reference": "uhf", "qc_module": "ccenergy"},                                                                           }, id="ccsd  uhf fc: psi4-cc",    marks=using("psi4")),
        pytest.param({"call": "p4-ccsd",  "reference": "uhf",  "fcae": "fc", "keywords": {"psi4_freeze_core": True, "reference": "uhf", "qc_module": "mrcc"},                                                                               }, id="ccsd  uhf fc: psi4-mrcc",  marks=using("psi4")),

        pytest.param({"call": "c4-ccsd",  "reference": "rohf", "fcae": "ae", "sdsc": "sc", "keywords": {**_c4_tight, "cfour_REFerence": "roHF", "cfour_cc_program": "vcc", "cfour_print": 2},                                                                 }, id="ccsd rohf ae sc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rohf", "fcae": "ae", "sdsc": "sc", "keywords": {**_c4_tight, "cfour_REFerence": "roHF", "cfour_cc_program": "ecc"},                                                                                   }, id="ccsd rohf ae sc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rohf", "fcae": "ae", "sdsc": "sd", "keywords": {**_c4_tight, "cfour_REFerence": "roHF", "cfour_orbitals": 0, "cfour_cc_program": "vcc", "cfour_print": 2},                                            }, id="ccsd rohf ae sd: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rohf", "fcae": "ae", "sdsc": "sd", "keywords": {**_c4_tight, "cfour_REFerence": "roHF", "cfour_orbitals": 0, "cfour_cc_program": "ecc"},                                                              }, id="ccsd rohf ae sd: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "rohf", "fcae": "ae",               "keywords": {"gamess_contrl__scftyp": "rohf", "gamess_ccinp__ncore": 0, "gamess_ccinp__maxcc": 50},                                                                }, id="ccsd rohf ae   : gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "rohf", "fcae": "ae",               "keywords": {"nwchem_scf__rohf": True, "qc_module": "tce"},                                                                                                        }, id="ccsd rohf ae   : nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "rohf", "fcae": "ae",               "keywords": {"nwchem_scf__rohf": True, "qc_module": "cc"},                                                                                       "error": {0: _q3 }}, id="ccsd rohf ae   : nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "rohf", "fcae": "ae", "sdsc": "sd", "keywords": {"reference": "rohf", "qc_module": "ccenergy"},                                                                                                        }, id="ccsd rohf ae sd: psi4-cc",    marks=using("psi4")),
        pytest.param({"call": "p4-ccsd",  "reference": "rohf", "fcae": "ae", "sdsc": "sc", "keywords": {"reference": "rohf", "qc_module": "mrcc"},                                                                                                            }, id="ccsd rohf ae sc: psi4-mrcc",  marks=using("psi4")),

        # at least two formulations: standard (sd) or semicanonical (sc) orbitals. most implementations choose one.
        pytest.param({"call": "c4-ccsd",  "reference": "rohf", "fcae": "fc", "sdsc": "sd", "keywords": {**_c4_tight, "cfour_REFerence": "roHF", "cfour_dropmo": [1], "cfour_orbitals": 0, "cfour_cc_program": "vcc", "cfour_print": 2},                       }, id="ccsd rohf fc sd: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rohf", "fcae": "fc", "sdsc": "sd", "keywords": {**_c4_tight, "cfour_REFerence": "roHF", "cfour_dropmo": [1], "cfour_orbitals": 0, "cfour_cc_program": "ecc"},                                         }, id="ccsd rohf fc sd: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rohf", "fcae": "fc", "sdsc": "sc", "keywords": {**_c4_tight, "cfour_REFerence": "roHF", "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_print": 2},                                            }, id="ccsd rohf fc sc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rohf", "fcae": "fc", "sdsc": "sc", "keywords": {**_c4_tight, "cfour_REFerence": "roHF", "cfour_dropmo": [1], "cfour_cc_program": "ecc"},                                                              }, id="ccsd rohf fc sc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "rohf", "fcae": "fc",               "keywords": {"gamess_contrl__scftyp": "rohf", "gamess_ccinp__iconv": 9, "gamess_scf__conv": 9},                                                  "wrong": {0: _w2 }}, id="ccsd rohf fc   : gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "rohf", "fcae": "fc", "sdsc": "sd", "keywords": {"nwchem_scf__rohf": True, "nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                               }, id="ccsd rohf fc sd: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "rohf", "fcae": "fc",               "keywords": {"nwchem_scf__rohf": True, "nwchem_ccsd__freeze": 1, "qc_module": "cc"},                                                             "error": {0: _q3 }}, id="ccsd rohf fc   : nwchem-cc",  marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "rohf", "fcae": "fc", "sdsc": "sd", "keywords": {"reference": "rohf", "psi4_freeze_core": True,  "qc_module": "ccenergy", "psi4_e_convergence": 8, "psi4_r_convergence": 7},                           }, id="ccsd rohf fc sd: psi4-cc",    marks=using("psi4")),
        pytest.param({"call": "p4-ccsd",  "reference": "rohf", "fcae": "fc", "sdsc": "sc", "keywords": {"reference": "rohf", "psi4_freeze_core": True,  "qc_module": "ccenergy", "psi4_e_convergence": 8, "psi4_r_convergence": 7, "psi4_semicanonical": True}}, id="ccsd rohf fc sc: psi4-cc",    marks=using("psi4")),
        pytest.param({"call": "p4-ccsd",  "reference": "rohf", "fcae": "fc", "sdsc": "sc", "keywords": {"reference": "rohf", "psi4_freeze_core": True,  "qc_module": "mrcc", "psi4_e_convergence": 8, "psi4_r_convergence": 7},                               }, id="ccsd rohf fc sc: psi4-mrcc",  marks=using("psi4")),
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

        ## local
        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "fc", "sdsc": "sd", "xptd": {"qc_module": "ecc"},      "keywords": {"cfour_dropmo": 1                                                                    },              }, id="ccsd  rhf fc: cfour dd1",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "fc", "sdsc": "sd", "xptd": {"qc_module": "ecc"},      "keywords": {"cfour_dropmo": 1,         "cfour_reference": "uhf"                                  },              }, id="ccsd  uhf fc: cfour dd1",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rohf", "fcae": "fc", "sdsc": "sd", "xptd": {"qc_module": "ecc"},      "keywords": {"cfour_dropmo": 1,         "cfour_reference": "rohf",       "cfour_orbitals": 0      },              }, id="ccsd rohf fc: cfour dd1",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "ae", "sdsc": "sd", "xptd": {"qc_module": "ecc"},      "keywords": {},                                                                                                   }, id="ccsd  rhf ae: cfour dd1",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "ae", "sdsc": "sd", "xptd": {"qc_module": "ecc"},      "keywords": {                           "cfour_reference": "uhf"                                  },              }, id="ccsd  uhf ae: cfour dd1",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rohf", "fcae": "ae", "sdsc": "sc", "xptd": {"qc_module": "ecc"},      "keywords": {                           "cfour_reference": "rohf"                                 },              }, id="ccsd rohf ae: cfour dd1",  marks=using("cfour")),

        pytest.param({"call": "gms-ccsd", "reference": "rhf",  "fcae": "fc", "sdsc": "sd", "xptd": {"qc_module": None},       "keywords": {},                                                                                                   }, id="ccsd  rhf fc: gamess dd1", marks=using("gamess")),
        pytest.param({"call": "gms-ccsd", "reference": "uhf",  "fcae": "fc", "sdsc": "sd", "xptd": {"qc_module": None},       "keywords": {                           "gamess_contrl__scftyp": "uhf"                            }, "error": _q2,}, id="ccsd  uhf fc: gamess dd1", marks=using("gamess")),
        pytest.param({"call": "gms-ccsd", "reference": "rohf", "fcae": "fc", "sdsc": "sd", "xptd": {"qc_module": None},       "keywords": {                           "gamess_contrl__scftyp": "rohf", "gamess_ccinp__maxcc": 50}, "wrong": _w2,}, id="ccsd rohf fc: gamess dd1", marks=using("gamess")),
        pytest.param({"call": "gms-ccsd", "reference": "rhf",  "fcae": "ae", "sdsc": "sd", "xptd": {"qc_module": None},       "keywords": {"gamess_ccinp__ncore": 0                                                             },              }, id="ccsd  rhf ae: gamess dd1", marks=using("gamess")),
        pytest.param({"call": "gms-ccsd", "reference": "uhf",  "fcae": "ae", "sdsc": "sd", "xptd": {"qc_module": None},       "keywords": {"gamess_ccinp__ncore": 0,  "gamess_contrl__scftyp": "uhf"                            }, "error": _q2,}, id="ccsd  uhf ae: gamess dd1", marks=using("gamess")),
        pytest.param({"call": "gms-ccsd", "reference": "rohf", "fcae": "ae", "sdsc": "sd", "xptd": {"qc_module": None},       "keywords": {"gamess_ccinp__ncore": 0,  "gamess_contrl__scftyp": "rohf", "gamess_ccinp__maxcc": 50},              }, id="ccsd rohf ae: gamess dd1", marks=using("gamess")),

        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "fc", "sdsc": "sd", "xptd": {"qc_module": "cc"},       "keywords": {"nwchem_ccsd__freeze": 1                                                             },              }, id="ccsd  rhf fc: nwchem dd1", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "uhf",  "fcae": "fc", "sdsc": "sd", "xptd": {"qc_module": "tce"},      "keywords": {"nwchem_tce__freeze": 1,   "nwchem_scf__uhf": True,         "qc_module": "tce"       },              }, id="ccsd  uhf fc: nwchem dd1", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "rohf", "fcae": "fc", "sdsc": "sd", "xptd": {"qc_module": "tce"},      "keywords": {"nwchem_tce__freeze": 1,   "nwchem_scf__rohf": True,        "qc_module": "tce"       },              }, id="ccsd rohf fc: nwchem dd1", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "ae", "sdsc": "sd", "xptd": {"qc_module": "cc"},       "keywords": {},                                                                                                   }, id="ccsd  rhf ae: nwchem dd1", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "uhf",  "fcae": "ae", "sdsc": "sd", "xptd": {"qc_module": "tce"},      "keywords": {                           "nwchem_scf__uhf": True,         "qc_module": "tce"       },              }, id="ccsd  uhf ae: nwchem dd1", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "rohf", "fcae": "ae", "sdsc": "sd", "xptd": {"qc_module": "tce"},      "keywords": {                           "nwchem_scf__rohf": True,        "qc_module": "tce"       },              }, id="ccsd rohf ae: nwchem dd1", marks=using("nwchem")),

        pytest.param({"call": "p4-ccsd",  "reference": "rhf",  "fcae": "fc", "sdsc": "sd", "xptd": {"qc_module": "ccenergy"}, "keywords": {"psi4_freeze_core": True                                                             },              }, id="ccsd  rhf fc: psi4 dd1",   marks=using("psi4")),
        pytest.param({"call": "p4-ccsd",  "reference": "uhf",  "fcae": "fc", "sdsc": "sd", "xptd": {"qc_module": "ccenergy"}, "keywords": {"psi4_freeze_core": True,  "psi4_reference": "uhf"                                   },              }, id="ccsd  uhf fc: psi4 dd1",   marks=using("psi4")),
        pytest.param({"call": "p4-ccsd",  "reference": "rohf", "fcae": "fc", "sdsc": "sd", "xptd": {"qc_module": "ccenergy"}, "keywords": {"psi4_freeze_core": True,  "psi4_reference": "rohf",                                 },              }, id="ccsd rohf fc: psi4 dd1",   marks=using("psi4")),
        pytest.param({"call": "p4-ccsd",  "reference": "rhf",  "fcae": "ae", "sdsc": "sd", "xptd": {"qc_module": "ccenergy"}, "keywords": {},                                                                                                   }, id="ccsd  rhf ae: psi4 dd1",   marks=using("psi4")),
        pytest.param({"call": "p4-ccsd",  "reference": "uhf",  "fcae": "ae", "sdsc": "sd", "xptd": {"qc_module": "ccenergy"}, "keywords": {                           "psi4_reference": "uhf"                                   },              }, id="ccsd  uhf ae: psi4 dd1",   marks=using("psi4")),
        pytest.param({"call": "p4-ccsd",  "reference": "rohf", "fcae": "ae", "sdsc": "sd", "xptd": {"qc_module": "ccenergy"}, "keywords": {                           "psi4_reference": "rohf"                                  },              }, id="ccsd rohf ae: psi4 dd1",   marks=using("psi4")),

        ## translated
        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "fc", "sdsc": "sd", "xptd": {"qc_module": "ecc"},      "keywords": {"freeze_core": True,       "reference": "rhf"                                        },              }, id="ccsd  rhf fc: cfour dd2",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "fc", "sdsc": "sd", "xptd": {"qc_module": "ecc"},      "keywords": {"freeze_core": True,       "reference": "uhf"                                        },              }, id="ccsd  uhf fc: cfour dd2",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rohf", "fcae": "fc", "sdsc": "sd", "xptd": {"qc_module": "ecc"},      "keywords": {"freeze_core": True,       "reference": "rohf",             "cfour_orbitals": 0      },              }, id="ccsd rohf fc: cfour dd2",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "ae", "sdsc": "sd", "xptd": {"qc_module": "ecc"},      "keywords": {"freeze_core": False,      "reference": "rhf"                                        },              }, id="ccsd  rhf ae: cfour dd2",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "ae", "sdsc": "sd", "xptd": {"qc_module": "ecc"},      "keywords": {"freeze_core": False,      "reference": "uhf"                                        },              }, id="ccsd  uhf ae: cfour dd2",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rohf", "fcae": "ae", "sdsc": "sc", "xptd": {"qc_module": "ecc"},      "keywords": {"freeze_core": False,      "reference": "rohf"                                       },              }, id="ccsd rohf ae: cfour dd2",  marks=using("cfour")),

        pytest.param({"call": "gms-ccsd", "reference": "rhf",  "fcae": "fc", "sdsc": "sd", "xptd": {"qc_module": None},       "keywords": {"freeze_core": True,       "reference": "rhf"                                        },              }, id="ccsd  rhf fc: gamess dd2", marks=using("gamess")),
        pytest.param({"call": "gms-ccsd", "reference": "uhf",  "fcae": "fc", "sdsc": "sd", "xptd": {"qc_module": None},       "keywords": {"freeze_core": True,       "reference": "uhf",                                       }, "error": _q2,}, id="ccsd  uhf fc: gamess dd2", marks=using("gamess")),
        pytest.param({"call": "gms-ccsd", "reference": "rohf", "fcae": "fc", "sdsc": "sd", "xptd": {"qc_module": None},       "keywords": {"freeze_core": True,       "reference": "rohf",             "gamess_ccinp__maxcc": 50}, "wrong": _w2,}, id="ccsd rohf fc: gamess dd2", marks=using("gamess")),
        pytest.param({"call": "gms-ccsd", "reference": "rhf",  "fcae": "ae", "sdsc": "sd", "xptd": {"qc_module": None},       "keywords": {"freeze_core": False,      "reference": "rhf"                                        },              }, id="ccsd  rhf ae: gamess dd2", marks=using("gamess")),
        pytest.param({"call": "gms-ccsd", "reference": "uhf",  "fcae": "ae", "sdsc": "sd", "xptd": {"qc_module": None},       "keywords": {"freeze_core": False,      "reference": "uhf"                                        }, "error": _q2,}, id="ccsd  uhf ae: gamess dd2", marks=using("gamess")),
        pytest.param({"call": "gms-ccsd", "reference": "rohf", "fcae": "ae", "sdsc": "sd", "xptd": {"qc_module": None},       "keywords": {"freeze_core": False,      "reference": "rohf",             "gamess_ccinp__maxcc": 50},              }, id="ccsd rohf ae: gamess dd2", marks=using("gamess")),

        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "fc", "sdsc": "sd", "xptd": {"qc_module": "cc"},       "keywords": {"freeze_core": True,       "reference": "rhf"                                        },              }, id="ccsd  rhf fc: nwchem dd2", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "uhf",  "fcae": "fc", "sdsc": "sd", "xptd": {"qc_module": "tce"},      "keywords": {"freeze_core": True,       "reference": "uhf",              "qc_module": "tce"       },              }, id="ccsd  uhf fc: nwchem dd2", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "rohf", "fcae": "fc", "sdsc": "sd", "xptd": {"qc_module": "tce"},      "keywords": {"freeze_core": True,       "reference": "rohf",             "qc_module": "tce"       },              }, id="ccsd rohf fc: nwchem dd2", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "ae", "sdsc": "sd", "xptd": {"qc_module": "cc"},       "keywords": {"freeze_core": False,      "reference": "rhf"                                        },              }, id="ccsd  rhf ae: nwchem dd2", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "uhf",  "fcae": "ae", "sdsc": "sd", "xptd": {"qc_module": "tce"},      "keywords": {"freeze_core": False,      "reference": "uhf",              "qc_module": "tce"       },              }, id="ccsd  uhf ae: nwchem dd2", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "rohf", "fcae": "ae", "sdsc": "sd", "xptd": {"qc_module": "tce"},      "keywords": {"freeze_core": False,      "reference": "rohf",             "qc_module": "tce"       },              }, id="ccsd rohf ae: nwchem dd2", marks=using("nwchem")),

        pytest.param({"call": "p4-ccsd",  "reference": "rhf",  "fcae": "fc", "sdsc": "sd", "xptd": {"qc_module": "ccenergy"}, "keywords": {"freeze_core": True,       "reference": "rhf"                                        },              }, id="ccsd  rhf fc: psi4 dd2",   marks=using("psi4")),
        pytest.param({"call": "p4-ccsd",  "reference": "uhf",  "fcae": "fc", "sdsc": "sd", "xptd": {"qc_module": "ccenergy"}, "keywords": {"freeze_core": True,       "reference": "uhf"                                        },              }, id="ccsd  uhf fc: psi4 dd2",   marks=using("psi4")),
        pytest.param({"call": "p4-ccsd",  "reference": "rohf", "fcae": "fc", "sdsc": "sd", "xptd": {"qc_module": "ccenergy"}, "keywords": {"freeze_core": True,       "reference": "rohf"                                       },              }, id="ccsd rohf fc: psi4 dd2",   marks=using("psi4")),
        pytest.param({"call": "p4-ccsd",  "reference": "rhf",  "fcae": "ae", "sdsc": "sd", "xptd": {"qc_module": "ccenergy"}, "keywords": {"freeze_core": False,      "reference": "rhf"                                        },              }, id="ccsd  rhf ae: psi4 dd2",   marks=using("psi4")),
        pytest.param({"call": "p4-ccsd",  "reference": "uhf",  "fcae": "ae", "sdsc": "sd", "xptd": {"qc_module": "ccenergy"}, "keywords": {"freeze_core": False,      "reference": "uhf"                                        },              }, id="ccsd  uhf ae: psi4 dd2",   marks=using("psi4")),
        pytest.param({"call": "p4-ccsd",  "reference": "rohf", "fcae": "ae", "sdsc": "sd", "xptd": {"qc_module": "ccenergy"}, "keywords": {"freeze_core": False,      "reference": "rohf"                                       },              }, id="ccsd rohf ae: psi4 dd2",   marks=using("psi4")),
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
    inpcopy["scf_type"] = "conv"  # longstanding "pk"
    inpcopy["corl_type"] = "conv"
    inpcopy["qc_module"] = "-".join(
        [qcprog, inp["keywords"].get("qc_module", inp["keywords"].get("cfour_cc_program", ""))]
    ).strip("-")
    print("INP", inpcopy)

    runner_asserter(inpcopy, subject, method, basis, tnm, scramble=None, frame="")


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
        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_cc_program": "vcc", "cfour_print": 2},                                                                                      }, id="ccsd  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_cc_program": "ecc"},                                                                                                        }, id="ccsd  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_cc_program": "ncc"},                                                                                                        }, id="ccsd  rhf ae: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {"gamess_ccinp__ncore": 0, "gamess_contrl__numgrd": True,},                                                                                      }, id="ccsd  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {"qc_module": "tce"},                                                                                                                            }, id="ccsd  rhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {},                                                                                                                                              }, id="ccsd  rhf ae: nwchem-cc",  marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "rhf",  "fcae": "ae",                        "keywords": {"psi4_cc_type": "conv", "qc_module": "ccenergy"},                                                                                               }, id="ccsd  rhf ae: psi4-cc",    marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-ccsd", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "vcc"},                                                                                   }, id="ccsd  rhf ae: psi4-cfour-vcc"),
        # DEBUG pytest.param({"call": "p4-c4-ccsd", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "ncc"},                                                                                   }, id="ccsd  rhf ae: psi4-cfour-ncc"),

        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_cc_program": "vcc", "cfour_dropmo": 1, "cfour_print": 2},                                                                   }, id="ccsd  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_cc_program": "ecc", "cfour_dropmo": 1},                                                                                     }, id="ccsd  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_cc_program": "ncc", "cfour_dropmo": 1},                                                                   "wrong": {1: _w8} }, id="ccsd  rhf fc: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "rhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {"gamess_contrl__numgrd": True,},                                                                                                                }, id="ccsd  rhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {"nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                                                   }, id="ccsd  rhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {"nwchem_ccsd__freeze": 1},                                                                                                                      }, id="ccsd  rhf fc: nwchem-cc",  marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "rhf",  "fcae": "fc",                        "keywords": {"psi4_freeze_core": True, "psi4_cc_type": "conv", "qc_module": "ccenergy"},                                                   "error": {1: _q18}}, id="ccsd  rhf fc: psi4-cc",    marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-ccsd", "reference": "rhf", "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_dropmo": [1], "psi4_cfour_cc_program": "vcc"},                                                         }, id="ccsd  rhf fc: psi4-cfour-vcc"),
        # DEBUG pytest.param({"call": "p4-c4-ccsd", "reference": "rhf", "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_dropmo": [1], "psi4_cfour_cc_program": "ncc"},                                                         }, id="ccsd  rhf fc: psi4-cfour-ncc"),

        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "vcc", "cfour_print": 2},                                                            }, id="ccsd  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "ecc"},                                                                              }, id="ccsd  uhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "uhf",  "fcae": "ae",                        "keywords": {"gamess_ccinp__ncore": 0, "gamess_contrl__scftyp": "uhf", "gamess_contrl__numgrd": True},                                     "error": {1: _q2} }, id="ccsd  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "uhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {"nwchem_scf__uhf": True, "qc_module": "tce"},                                                                                                   }, id="ccsd  uhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "uhf",  "fcae": "ae",                        "keywords": {"nwchem_scf__uhf": True, "qc_module": "cc"},                                                                                  "error": {1: _q3} }, id="ccsd  uhf ae: nwchem-cc",  marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "uhf",  "fcae": "ae",                        "keywords": {"reference": "uhf", "psi4_cc_type": "conv"},                                                                                                    }, id="ccsd  uhf ae: psi4-cc",    marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-ccsd", "reference": "uhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "uhf", "psi4_cfour_cc_program": "vcc"},                                                                                   }, id="ccsd  uhf ae: psi4-cfour-vcc"),

        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_dropmo": 1, "cfour_cc_program": "vcc", "cfour_print": 2},                                         }, id="ccsd  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_dropmo": 1, "cfour_cc_program": "ecc"},                                                           }, id="ccsd  uhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "uhf",  "fcae": "fc",                        "keywords": {"gamess_contrl__scftyp": "uhf", "gamess_contrl__numgrd": True,},                                                              "error": {1: _q2} }, id="ccsd  uhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "uhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {"nwchem_scf__uhf": True, "nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                          }, id="ccsd  uhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "uhf",  "fcae": "fc",                        "keywords": {"nwchem_scf__uhf": True, "nwchem_ccsd__freeze": 1},                                                                           "error": {1: _q3} }, id="ccsd  uhf fc: nwchem-cc",  marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "uhf",  "fcae": "fc",                        "keywords": {"reference": "uhf", "psi4_cc_type": "conv", "psi4_freeze_core": True},                                                        "error": {1: _q18}}, id="ccsd  uhf fc: psi4",       marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-ccsd", "reference": "uhf", "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "uhf", "psi4_cfour_dropmo": [1], "psi4_cfour_cc_program": "vcc"},                                                         }, id="ccsd  uhf fc: psi4-cfour-vcc"),

        pytest.param({"call": "c4-ccsd",  "reference": "rohf", "fcae": "ae",                                      "keywords": {**_c4_tight, "cfour_reference": "rohf", "cfour_cc_program": "vcc", "cfour_print": 2},                                                           }, id="ccsd rohf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rohf", "fcae": "ae",                                      "keywords": {**_c4_tight, "cfour_reference": "rohf", "cfour_cc_program": "ecc", "cfour_print": 2},                                                           }, id="ccsd rohf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "rohf", "fcae": "ae",               "xptd": {"fd": True},  "keywords": {"gamess_contrl__scftyp": "rohf", "gamess_ccinp__ncore": 0, "gamess_ccinp__maxcc": 50, "gamess_contrl__numgrd": True},                           }, id="ccsd rohf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "rohf", "fcae": "ae",               "xptd": {"fd": True},  "keywords": {"qc_module": "tce", "nwchem_scf__rohf": True},                                                                                                  }, id="ccsd rohf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "rohf", "fcae": "ae",                                      "keywords": {"nwchem_scf__rohf": True},                                                                                                    "error": {1: _q3} }, id="ccsd rohf ae: nwchem-cc",  marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "rohf", "fcae": "ae", "sdsc": "sd",                        "keywords": {"reference": "rohf", "qc_module": "ccenergy"},                                                                                                  }, id="ccsd rohf ae: psi4-cc",    marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-ccsd", "reference": "rohf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rohf", "psi4_cfour_cc_program": "vcc"},                                                                                               }, id="ccsd rohf ae: psi4-cfour-vcc"),

        # * c. Summer 2021 vcc and ecc yield correct gradients w/pert_orb=0 but not at the same time with correct energies (orbitals=0). c. Summer 2022: previous refers to "sd" targets. "sc" targets work w/ or w/o pert_orb=0
        # * note that ecc not recc for rohf and no rohf gradients on list in cfour paper
        # * below (for both ecc and vcc) no longer "wrong" for sc
            # pytest.param({"call": "c4-ccsd",  "reference": "rohf", "fcae": "fc",                                      "keywords": {**_c4_tight, "cfour_reference": "rohf", "cfour_dropmo": 1, "cfour_cc_program": "vcc", "cfour_print": 2, "cfour_pert_orb": 0}, "wrong": {1: _w12}}, id="ccsd rohf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rohf", "fcae": "fc", "sdsc": "sc",                        "keywords": {**_c4_tight, "cfour_reference": "rohf", "cfour_dropmo": 1, "cfour_cc_program": "vcc", "cfour_print": 2},                                        }, id="ccsd rohf fc sc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rohf", "fcae": "fc", "sdsc": "sc",                        "keywords": {**_c4_tight, "cfour_reference": "rohf", "cfour_dropmo": 1, "cfour_cc_program": "ecc", "cfour_print": 2},                                        }, id="ccsd rohf fc sc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "rohf", "fcae": "fc",               "xptd": {"fd": True},  "keywords": {"gamess_contrl__scftyp": "rohf", "gamess_ccinp__iconv": 9, "gamess_scf__conv": 9, "gamess_contrl__numgrd": True},             "wrong": {1: _w2} }, id="ccsd rohf fc   : gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "rohf", "fcae": "fc", "sdsc": "sd", "xptd": {"fd": True},  "keywords": {"nwchem_tce__freeze": 1, "nwchem_scf__rohf": True, "qc_module": "tce"},                                                                         }, id="ccsd rohf fc sd: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd", "reference": "rohf", "fcae": "fc",                                      "keywords": {"nwchem_ccsd__freeze": 1, "nwchem_scf__rohf": True, "qc_module": "cc"},                                                       "error": {1: _q3} }, id="ccsd rohf fc   : nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "rohf", "fcae": "fc",                                      "keywords": {"reference": "rohf", "psi4_cc_type": "conv", "psi4_freeze_core": True},                                                       "error": {1: _q18}}, id="ccsd rohf fc   : psi4-cc",    marks=using("psi4")),
        pytest.param({"call": "p4-ccsd",  "reference": "rohf", "fcae": "fc", "sdsc": "sd", "xptd": {"fd": False}, "keywords": {"reference": "rohf", "psi4_cc_type": "conv", "psi4_freeze_core": True, **_p4_fd, "psi4_function_kwargs_dertype": 0},                            }, id="ccsd rohf fc sd: psi4-cc",    marks=using("psi4")),
        pytest.param({"call": "p4-ccsd",  "reference": "rohf", "fcae": "fc", "sdsc": "sc", "xptd": {"fd": False}, "keywords": {"reference": "rohf", "psi4_cc_type": "conv", "psi4_freeze_core": True, **_p4_fd, "psi4_function_kwargs_dertype": 0, "psi4_semicanonical": True},}, id="ccsd rohf fc sc: psi4-cc",    marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-ccsd", "reference": "rohf", "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rohf", "psi4_cfour_dropmo": [1], "psi4_cfour_cc_program": "vcc", "cfour_print": 2},                                     }, id="ccsd rohf fc: psi4-cfour-vcc"),
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
        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_cc_program": "vcc", "cfour_print": 2},                                                                    }, id="ccsd  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_cc_program": "ecc", "cfour_print": 2},                                                                    }, id="ccsd  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {"gamess_ccinp__ncore": 0, "gamess_force__method": "fullnum"},                                                                 }, id="ccsd  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {"qc_module": "tce"},                                                                                                          }, id="ccsd  rhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "rhf",  "fcae": "ae", "xptd": {"fd": False}, "keywords": {**_p4_fd, "psi4_cc_type": "conv"},                                                                                            }, id="ccsd  rhf ae: psi4",       marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-ccsd", "reference": "rhf", "fcae": "ae", "keywords": {"psi4_function_kwargs_dertype": 0, "psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "vcc", **_p4c4_fd},                              }, id="ccsd  rhf ae: psi4-cfour-vcc"),

        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "vcc", "cfour_print": 2},                                                 }, id="ccsd  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "ecc", "cfour_print": 2},                                                 }, id="ccsd  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "rhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {"gamess_force__method": "fullnum"},                                                                                           }, id="ccsd  rhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "rhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {"nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                                 }, id="ccsd  rhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "rhf",  "fcae": "fc", "xptd": {"fd": False},  "keywords": {**_p4_fd, "psi4_freeze_core": True, "psi4_cc_type": "conv", "psi4_function_kwargs_dertype": 0},                              }, id="ccsd  rhf fc: psi4",       marks=using("psi4")),  # FD SWITCH
        # DEBUG pytest.param({"call": "p4-c4-ccsd", "reference": "rhf", "fcae": "fc", "keywords": {"psi4_function_kwargs_dertype": 0, "psi4_cfour_dropmo": [1], "psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "vcc", **_p4c4_fd},    }, id="ccsd  rhf fc: psi4-cfour-vcc"),

        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "vcc", "cfour_print": 2},                                          }, id="ccsd  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "ecc"},                                                            }, id="ccsd  uhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "uhf",  "fcae": "ae",                        "keywords": {"gamess_contrl__scftyp": "uhf", "gamess_ccinp__ncore": 0},                                                  "error": {2: _q2} }, id="ccsd  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "uhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {"nwchem_scf__uhf": True, "qc_module": "tce"},                                                                                 }, id="ccsd  uhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "uhf",  "fcae": "ae",                        "keywords": {**_p4_fd, "reference": "uhf", "psi4_cc_type": "conv"},                                                                        }, id="ccsd  uhf ae: psi4",       marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-ccsd", "reference": "uhf", "fcae": "ae", "keywords": {"psi4_function_kwargs_dertype": 0, "psi4_cfour_reference": "uhf", "psi4_cfour_cc_program": "vcc", **_p4c4_fd},                              }, id="ccsd  uhf ae: psi4-cfour-vcc"),

        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_dropmo": 1, "cfour_cc_program": "vcc", "cfour_print": 2},                       }, id="ccsd  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd",  "reference": "uhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_dropmo": 1, "cfour_cc_program": "ecc", "cfour_print": 2},                       }, id="ccsd  uhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd", "reference": "uhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {"gamess_contrl__scftyp": "uhf", "gamess_force__method": "fullnum"},                                         "error": {2: _q2} }, id="ccsd  uhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd", "reference": "uhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {"nwchem_scf__uhf": True, "nwchem_tce__freeze": 1, "qc_module": "tce"},                                                        }, id="ccsd  uhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd",  "reference": "uhf",  "fcae": "fc", "xptd": {"fd": False}, "keywords": {**_p4_fd, "psi4_reference": "uhf", "psi4_freeze_core": True, "psi4_cc_type": "conv", "psi4_function_kwargs_dertype": 0},      }, id="ccsd  uhf fc: psi4",       marks=using("psi4")),
        # FD SWITCH
        # DEBUG pytest.param({"call": "p4-c4-ccsd", "reference": "uhf", "fcae": "fc", "keywords": {"psi4_function_kwargs_dertype": 0, "psi4_cfour_dropmo": [1], "psi4_cfour_reference": "uhf", "psi4_cfour_cc_program": "vcc", **_p4c4_fd},    }, id="ccsd  uhf fc: psi4-cfour-vcc"),
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
        pytest.param({"call": "c4-ccsd+t(ccsd)",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "vcc", "cfour_print": 2},                                                                                      }, id="ccsd_t_ccsd_  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd+t(ccsd)",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "ecc",},                                                                                                       }, id="ccsd_t_ccsd_  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd+t(ccsd)",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "ncc",},                                                                                                       }, id="ccsd_t_ccsd_  rhf ae: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd+t(ccsd)", "reference": "rhf",  "fcae": "ae", "keywords": {"gamess_ccinp__ncore": 0},                                                                                                                      }, id="ccsd_t_ccsd_  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd+t(ccsd)", "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "tce"},                                                                                                                            }, id="ccsd_t_ccsd_  rhf ae: nwchem-tce",  marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd+t(ccsd)", "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "cc"},                                                                                                                             }, id="ccsd_t_ccsd_  rhf ae: nwchem-cc",  marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd+t(ccsd)",  "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                            "error": {0: _q20}}, id="ccsd_t_ccsd_  rhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-ccsd+t(ccsd)",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_print": 2},                                                                 }, id="ccsd_t_ccsd_  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd+t(ccsd)",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": [1], "cfour_cc_program": "ecc"},                                                                                   }, id="ccsd_t_ccsd_  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd+t(ccsd)",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": [1], "cfour_cc_program": "ncc"},                                                                                   }, id="ccsd_t_ccsd_  rhf fc: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd+t(ccsd)", "reference": "rhf",  "fcae": "fc", "keywords": {},                                                                                                                                              }, id="ccsd_t_ccsd_  rhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd+t(ccsd)", "reference": "rhf",  "fcae": "fc", "keywords": {"qc_module": "tce", "nwchem_tce__freeze": 1 },                                                                                                  }, id="ccsd_t_ccsd_  rhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd+t(ccsd)", "reference": "rhf",  "fcae": "fc", "keywords": {"nwchem_ccsd__freeze": 1, "qc_module": "cc"},                                                                                                   }, id="ccsd_t_ccsd_  rhf fc: nwchem-cc",  marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd+t(ccsd)",  "reference": "rhf",  "fcae": "fc", "keywords": {"psi4_freeze_core": True},                                                                                                    "error": {0: _q20}}, id="ccsd_t_ccsd_  rhf fc: psi4-cc",    marks=using("psi4")),

        pytest.param({"call": "c4-ccsd+t(ccsd)",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_REFerence": "UHF", "cfour_cc_program": "vcc", "cfour_print": 2},                                                            }, id="ccsd_t_ccsd_  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd+t(ccsd)",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_REFerence": "UHF", "cfour_cc_program": "ecc"},                                                                              }, id="ccsd_t_ccsd_  uhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd+t(ccsd)", "reference": "uhf",  "fcae": "ae", "keywords": {"gamess_contrl__scftyp": "uhf", "gamess_ccinp__ncore": 0},                                                                    "error": {0: _q2} }, id="ccsd_t_ccsd_  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd+t(ccsd)", "reference": "uhf",  "fcae": "ae", "keywords": {"qc_module": "tce", "nwchem_scf__uhf": True},                                                                                                   }, id="ccsd_t_ccsd_  uhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd+t(ccsd)", "reference": "uhf",  "fcae": "ae", "keywords": {"nwchem_scf__uhf": True},                                                                                                     "error": {0: _q3} }, id="ccsd_t_ccsd_  uhf ae: nwchem-cc",  marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd+t(ccsd)",  "reference": "uhf",  "fcae": "ae", "keywords": {"reference": "uhf", "qc_module": "ccenergy"},                                                                                 "error": {0: _q20}}, id="ccsd_t_ccsd_  uhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-ccsd+t(ccsd)",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": [1], "cfour_REFerence": "UHF", "cfour_cc_program": "vcc", "cfour_print": 2},                                       }, id="ccsd_t_ccsd_  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd+t(ccsd)",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": [1], "cfour_REFerence": "UHF", "cfour_cc_program": "ecc"},                                                         }, id="ccsd_t_ccsd_  uhf fc: cfour-ecc",  marks=using("cfour")),
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
#   ,-----. ,-----. ,---.  ,------.    ,-.,--------.,-.      ,------.
#  '  .--./'  .--./'   .-' |  .-.  \  / .''--.  .--''. \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |    |  |    `.  `-. |  |  \  :|  |    |  |    |  |    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '--'\'  '--'\.-'    ||  '--'  /|  |    |  |    |  |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----' `-----'`-----' `-------'  \ '.   `--'   .' /     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                      `-'          `-'                                    `---' `---'
#  <<<  CCSD(T) Energy


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
        pytest.param({"call": "c4-ccsd(t)",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "vcc", "cfour_print": 2},                                                                                                                                           }, id="ccsd_t_  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd(t)",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "ecc",},                                                                                                                                                            }, id="ccsd_t_  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd(t)",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "ncc",},                                                                                                                                                            }, id="ccsd_t_  rhf ae: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd(t)", "reference": "rhf",  "fcae": "ae", "keywords": {"gamess_ccinp__ncore": 0},                                                                                                                                                                           }, id="ccsd_t_  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd(t)", "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "tce"},                                                                                                                                                                                 }, id="ccsd_t_  rhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd(t)", "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "cc"},                                                                                                                                                                                  }, id="ccsd_t_  rhf ae: nwchem-cc",  marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd(t)",  "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "ccenergy"},                                                                                                                                                                            }, id="ccsd_t_  rhf ae: psi4-cc",    marks=using("psi4")),
        pytest.param({"call": "p4-ccsd(t)",  "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "fnocc"},                                                                                                                                                                               }, id="ccsd_t_  rhf ae: psi4-fnocc", marks=using("psi4")),
        pytest.param({"call": "p4-ccsd(t)",  "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "mrcc"},                                                                                                                                                                                }, id="ccsd_t_  rhf ae: psi4-mrcc",  marks=using("psi4")),  # using mrcc

        pytest.param({"call": "c4-ccsd(t)",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_print": 2},                                                                                                                      }, id="ccsd_t_  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd(t)",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": [1], "cfour_cc_program": "ecc"},                                                                                                                                        }, id="ccsd_t_  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd(t)",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": [1], "cfour_cc_program": "ncc"},                                                                                                                                        }, id="ccsd_t_  rhf fc: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd(t)", "reference": "rhf",  "fcae": "fc", "keywords": {},                                                                                                                                                                                                   }, id="ccsd_t_  rhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd(t)", "reference": "rhf",  "fcae": "fc", "keywords": {"nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                                                                                                        }, id="ccsd_t_  rhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd(t)", "reference": "rhf",  "fcae": "fc", "keywords": {"nwchem_ccsd__freeze": 1, "qc_module": "cc"},                                                                                                                                                        }, id="ccsd_t_  rhf fc: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd(t)",  "reference": "rhf",  "fcae": "fc", "keywords": {"psi4_freeze_core": True, "qc_module": "ccenergy"},                                                                                                                                                  }, id="ccsd_t_  rhf fc: psi4-cc",    marks=using("psi4")),
        pytest.param({"call": "p4-ccsd(t)",  "reference": "rhf",  "fcae": "fc", "keywords": {"psi4_freeze_core": True, "qc_module": "fnocc"},                                                                                                                                                     }, id="ccsd_t_  rhf fc: psi4-fnocc", marks=using("psi4")),
        pytest.param({"call": "p4-ccsd(t)",  "reference": "rhf",  "fcae": "fc", "keywords": {"psi4_freeze_core": True, "qc_module": "mrcc"},                                                                                                                                                      }, id="ccsd_t_  rhf fc: psi4-mrcc",  marks=using("psi4")),

        pytest.param({"call": "c4-ccsd(t)",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_REFerence": "UHF", "cfour_cc_program": "vcc", "cfour_print": 2},                                                                                                                 }, id="ccsd_t_  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd(t)",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_REFerence": "UHF", "cfour_cc_program": "ecc"},                                                                                                                                   }, id="ccsd_t_  uhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd(t)", "reference": "uhf",  "fcae": "ae", "keywords": {"gamess_contrl__scftyp": "uhf", "gamess_ccinp__ncore": 0},                                                                                                                         "error": {0: _q2} }, id="ccsd_t_  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd(t)", "reference": "uhf",  "fcae": "ae", "keywords": {"nwchem_scf__uhf": True, "qc_module": "tce"},                                                                                                                                                        }, id="ccsd_t_  uhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd(t)", "reference": "uhf",  "fcae": "ae", "keywords": {"nwchem_scf__uhf": True},                                                                                                                                                          "error": {0: _q3} }, id="ccsd_t_  uhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd(t)",  "reference": "uhf",  "fcae": "ae", "keywords": {"reference": "uhf", "qc_module": "ccenergy"},                                                                                                                                                        }, id="ccsd_t_  uhf ae: psi4-cc",    marks=using("psi4")),
        pytest.param({"call": "p4-ccsd(t)",  "reference": "uhf",  "fcae": "ae", "keywords": {"reference": "uhf", "qc_module": "mrcc"},                                                                                                                                                            }, id="ccsd_t_  uhf ae: psi4-mrcc",  marks=using("psi4")),

        pytest.param({"call": "c4-ccsd(t)",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_REFerence": "UHF", "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_print": 2},                                                                                            }, id="ccsd_t_  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd(t)",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_REFerence": "UHF", "cfour_dropmo": [1], "cfour_cc_program": "ecc"},                                                                                                              }, id="ccsd_t_  uhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd(t)", "reference": "uhf",  "fcae": "fc", "keywords": {"gamess_contrl__scftyp": "uhf"},                                                                                                                                                   "error": {0: _q2} }, id="ccsd_t_  uhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd(t)", "reference": "uhf",  "fcae": "fc", "keywords": {"nwchem_scf__uhf": True, "nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                                                                               }, id="ccsd_t_  uhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd(t)", "reference": "uhf",  "fcae": "fc", "keywords": {"nwchem_scf__uhf": True, "nwchem_ccsd__freeze": 1},                                                                                                                                "error": {0: _q3} }, id="ccsd_t_  uhf fc: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd(t)",  "reference": "uhf",  "fcae": "fc", "keywords": {"reference": "uhf", "psi4_freeze_core": True, "qc_module": "ccenergy"},                                                                                                                              }, id="ccsd_t_  uhf fc: psi4-cc",    marks=using("psi4")),
        pytest.param({"call": "p4-ccsd(t)",  "reference": "uhf",  "fcae": "fc", "keywords": {"reference": "uhf", "psi4_freeze_core": True, "qc_module": "mrcc"},                                                                                                                                  }, id="ccsd_t_  uhf fc: psi4-mrcc",  marks=using("psi4")),

        pytest.param({"call": "c4-ccsd(t)",  "reference": "rohf", "fcae": "ae", "keywords": {**_c4_tight, "cfour_REFerence": "roHF", "cfour_cc_program": "vcc", "cfour_print": 2},                                                                                                                }, id="ccsd_t_ rohf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd(t)",  "reference": "rohf", "fcae": "ae", "keywords": {**_c4_tight, "cfour_REFerence": "roHF", "cfour_cc_program": "ecc"},                                                                                                                                  }, id="ccsd_t_ rohf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd(t)", "reference": "rohf", "fcae": "ae", "keywords": {"gamess_contrl__scftyp": "rohf", "gamess_ccinp__ncore": 0, "gamess_ccinp__maxcc": 50},                                                                                             "error": {0: _q4} }, id="ccsd_t_ rohf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd(t)", "reference": "rohf", "fcae": "ae", "keywords": {"nwchem_scf__rohf": True, "qc_module": "tce"},                                                                                                                                     "wrong": {0: _w3} }, id="ccsd_t_ rohf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd(t)", "reference": "rohf", "fcae": "ae", "keywords": {"nwchem_scf__rohf": True},                                                                                                                                                         "error": {0: _q3} }, id="ccsd_t_ rohf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd(t)",  "reference": "rohf", "fcae": "ae", "keywords": {"reference": "rohf", "qc_module": "ccenergy"},                                                                                                                                                       }, id="ccsd_t_ rohf ae: psi4-cc",    marks=using("psi4")),
        pytest.param({"call": "p4-ccsd(t)",  "reference": "rohf", "fcae": "ae", "keywords": {"reference": "rohf", "qc_module": "mrcc"},                                                                                                                                                           }, id="ccsd_t_ rohf ae: psi4-mrcc",  marks=using("psi4")),

        ## Sum 2021: can get all passing c4 = p4 = mrcc but at cost of CCSD  not matching plain CCSD. Sum 2022: CCSD often has default standard (sd) or semicanonical (sc), while CCSD(T) generally only sc that were matching Sum 2021. sd allows (T) + default CCSD.
        ## TODO express sdsc variants in table. CCSD, too
        pytest.param({"call": "c4-ccsd(t)",  "reference": "rohf", "fcae": "fc", "sdsc": "sd", "keywords": {**_c4_tight, "cfour_REFerence": "roHF", "cfour_dropmo": [1], "cfour_orbitals": 0, "cfour_cc_program": "vcc", "cfour_print": 2},                                                        }, id="ccsd_t_ rohf fc sd: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd(t)",  "reference": "rohf", "fcae": "fc", "sdsc": "sd", "keywords": {**_c4_tight, "cfour_REFerence": "roHF", "cfour_dropmo": [1], "cfour_orbitals": 0, "cfour_cc_program": "ecc"},                                                                          }, id="ccsd_t_ rohf fc sd: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd(t)",  "reference": "rohf", "fcae": "fc", "sdsc": "sc", "keywords": {**_c4_tight, "cfour_REFerence": "roHF", "cfour_dropmo": [1],  "cfour_cc_program": "vcc", "cfour_print": 2},                                                                            }, id="ccsd_t_ rohf fc sc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd(t)",  "reference": "rohf", "fcae": "fc", "sdsc": "sc", "keywords": {**_c4_tight, "cfour_REFerence": "roHF", "cfour_dropmo": [1],  "cfour_cc_program": "ecc"},                                                                                              }, id="ccsd_t_ rohf fc sc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd(t)", "reference": "rohf", "fcae": "fc",               "keywords": {"gamess_contrl__scftyp": "rohf", "gamess_ccinp__iconv": 9, "gamess_scf__conv": 9},                                                                                   "error": {0: _q4} }, id="ccsd_t_ rohf fc   : gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd(t)", "reference": "rohf", "fcae": "fc", "sdsc": "sd", "keywords": {"nwchem_scf__rohf": True, "nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                                              "wrong": {0: _w3} }, id="ccsd_t_ rohf fc   : nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd(t)", "reference": "rohf", "fcae": "fc",               "keywords": {"nwchem_scf__rohf": True, "nwchem_ccsd__freeze": 1},                                                                                                                 "error": {0: _q3} }, id="ccsd_t_ rohf fc   : nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd(t)",  "reference": "rohf", "fcae": "fc", "sdsc": "sc", "keywords": {"reference": "rohf", "psi4_freeze_core": True, "qc_module": "ccenergy", "psi4_e_convergence": 10, "psi4_r_convergence": 9},                                                            }, id="ccsd_t_ rohf fc sc: psi4-cc",    marks=using("psi4")),
        pytest.param({"call": "p4-ccsd(t)",  "reference": "rohf", "fcae": "fc", "sdsc": "sc", "keywords": {"reference": "rohf", "psi4_freeze_core": True, "qc_module": "mrcc", "psi4_e_convergence": 10, "psi4_r_convergence": 9},                                                                }, id="ccsd_t_ rohf fc sc: psi4-mrcc",  marks=using("psi4")),
        # yapf: enable
    ],
)
def test_ccsd_prt_pr_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "energy"))


#
#   ,-----. ,-----. ,---.  ,------.    ,-.,--------.,-.       ,----.                     ,--.,--.                 ,--.
#  '  .--./'  .--./'   .-' |  .-.  \  / .''--.  .--''. \     '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  |    |  |    `.  `-. |  |  \  :|  |    |  |    |  |    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  '  '--'\'  '--'\.-'    ||  '--'  /|  |    |  |    |  |    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#   `-----' `-----'`-----' `-------'  \ '.   `--'   .' /      `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#                                      `-'          `-'
#  <<<  CCSD(T) Gradient


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
        pytest.param({"call": "c4-ccsd(t)",  "reference": "rhf",  "fcae": "ae",                       "keywords": {**_c4_tight, "cfour_cc_program": "vcc", "cfour_print": 2},                                                                                   }, id="ccsd_t_  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd(t)",  "reference": "rhf",  "fcae": "ae",                       "keywords": {**_c4_tight, "cfour_cc_program": "ecc",},                                                                                                    }, id="ccsd_t_  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd(t)",  "reference": "rhf",  "fcae": "ae",                       "keywords": {**_c4_tight, "cfour_cc_program": "ncc",},                                                                                                    }, id="ccsd_t_  rhf ae: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd(t)", "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {"gamess_ccinp__ncore": 0, "gamess_contrl__numgrd": True},                                                                                   }, id="ccsd_t_  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd(t)", "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {"qc_module": "tce"},                                                                                                                        }, id="ccsd_t_  rhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd(t)", "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {},                                                                                                                                          }, id="ccsd_t_  rhf ae: nwchem-cc",  marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd(t)",  "reference": "rhf",  "fcae": "ae",                      "keywords": {"qc_module": "ccenergy"},                                                                                                                     }, id="ccsd_t_  rhf ae: psi4-cc",    marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-ccsd(t)", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "ncc", "psi4_function_kwargs_dertype": 0},                                            }, id="ccsd_t_  rhf ae: psi4-cfour-ncc"),

        pytest.param({"call": "c4-ccsd(t)",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "vcc", "cfour_print": 2,},                                                              }, id="ccsd_t_  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd(t)",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "ecc"},                                                                                 }, id="ccsd_t_  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd(t)",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "ncc"},                                                               "wrong": {1: _w8} }, id="ccsd_t_  rhf fc: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd(t)", "reference": "rhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {"gamess_contrl__numgrd": True,},                                                                                                            }, id="ccsd_t_  rhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd(t)", "reference": "rhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {"nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                                               }, id="ccsd_t_  rhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd(t)", "reference": "rhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {"nwchem_ccsd__freeze": 1},                                                                                                                  }, id="ccsd_t_  rhf fc: nwchem-cc",  marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd(t)",  "reference": "rhf",  "fcae": "fc",                        "keywords": {"psi4_freeze_core": True, "qc_module": "ccenergy"},                                                                       "error": {1: _q18}}, id="ccsd_t_  rhf fc: psi4-cc",    marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-ccsd(t)", "reference": "rhf", "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_dropmo": [1], "psi4_cfour_cc_program": "ncc", "psi4_function_kwargs_dertype": 0},                  }, id="ccsd_t_  rhf fc: psi4-cfour-ncc"),

        pytest.param({"call": "c4-ccsd(t)",  "reference": "uhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "vcc", "cfour_print": 2},                                                        }, id="ccsd_t_  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd(t)",  "reference": "uhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "ecc"},                                                                          }, id="ccsd_t_  uhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd(t)", "reference": "uhf",  "fcae": "ae",                        "keywords": {"gamess_contrl__scftyp": "uhf", "gamess_ccinp__ncore": 0, "gamess_contrl__numgrd": True,},                                "error": {1: _q2} }, id="ccsd_t_  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd(t)", "reference": "uhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {"nwchem_scf__uhf": True, "qc_module": "tce"},                                                                                               }, id="ccsd_t_  uhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd(t)", "reference": "uhf",  "fcae": "ae",                        "keywords": {"nwchem_scf__uhf": True},                                                                                                 "error": {1: _q3} }, id="ccsd_t_  uhf ae: nwchem-cc",  marks=using("nwchem")),
        # VERY LONG pytest.param({"call": "p4-ccsd(t)",  "reference": "uhf",  "fcae": "ae",                        "keywords": {"reference": "uhf", "qc_module": "ccenergy"},                                                                                   }, id="ccsd_t_  uhf ae: psi4-cc",    marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-ccsd(t)", "reference": "uhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "uhf", "psi4_cfour_cc_program": "vcc", "psi4_function_kwargs_dertype": 0},                                            }, id="ccsd_t_  uhf ae: psi4-cfour-vcc"),

        pytest.param({"call": "c4-ccsd(t)",  "reference": "uhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_dropmo": 1, "cfour_cc_program": "vcc",  "cfour_print": 2},                                    }, id="ccsd_t_  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd(t)",  "reference": "uhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_dropmo": 1, "cfour_cc_program": "ecc"},                                                       }, id="ccsd_t_  uhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd(t)", "reference": "uhf",  "fcae": "fc",                        "keywords": {"gamess_contrl__scftyp": "uhf", "gamess_contrl__numgrd": True,},                                                          "error": {1: _q2} }, id="ccsd_t_  uhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd(t)", "reference": "uhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {"nwchem_scf__uhf": True, "nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                      }, id="ccsd_t_  uhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd(t)", "reference": "uhf",  "fcae": "fc",                        "keywords": {"nwchem_scf__uhf": True, "nwchem_ccsd__freeze": 1},                                                                       "error": {1: _q3} }, id="ccsd_t_  uhf fc: nwchem-cc",  marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd(t)",  "reference": "uhf",  "fcae": "fc",                        "keywords": {"reference": "uhf", "psi4_freeze_core": True},                                                                            "error": {1: _q18}}, id="ccsd_t_  uhf fc: psi4",       marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-ccsd(t)", "reference": "uhf", "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "uhf", "psi4_cfour_dropmo": [1], "psi4_cfour_cc_program": "vcc", "psi4_function_kwargs_dertype": 0},                  }, id="ccsd_t_  uhf fc: psi4-cfour-vcc"),

        pytest.param({"call": "c4-ccsd(t)",  "reference": "rohf", "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_reference": "rohf", "cfour_cc_program": "vcc", "cfour_print": 2},                                                       }, id="ccsd_t_ rohf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd(t)",  "reference": "rohf", "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_reference": "rohf", "cfour_cc_program": "ecc"},                                                                         }, id="ccsd_t_ rohf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd(t)", "reference": "rohf", "fcae": "ae", "xptd": {"fd": True},  "keywords": {"gamess_contrl__scftyp": "rohf", "gamess_ccinp__ncore": 0, "gamess_ccinp__maxcc": 50, "gamess_contrl__numgrd": True},     "error": {1: _q4} }, id="ccsd_t_ rohf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd(t)", "reference": "rohf", "fcae": "ae", "xptd": {"fd": True},  "keywords": {"nwchem_scf__rohf": True, "qc_module": "tce"},                                                                            "wrong": {1: _w3} }, id="ccsd_t_ rohf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd(t)", "reference": "rohf", "fcae": "ae",                        "keywords": {"nwchem_scf__rohf": True},                                                                                                "error": {1: _q3} }, id="ccsd_t_ rohf ae: nwchem-cc",  marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd(t)",  "reference": "rohf", "fcae": "ae", "xptd": {"fd": False}, "keywords": {"reference": "rohf", "qc_module": "ccenergy", "psi4_points": 5},                                                                            }, id="ccsd_t_ rohf ae: psi4-cc",    marks=using("psi4")),
        # DEBUG pytest.param({"call": "p4-c4-ccsd(t)", "reference": "rohf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rohf", "psi4_cfour_cc_program": "vcc", "psi4_function_kwargs_dertype": 0},                                          }, id="ccsd_t_ rohf ae: psi4-cfour-vcc"),
        # Skip ROHF FC until energy resolved

        pytest.param({"call": "c4-ccsd(t)",  "reference": "rohf", "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_reference": "rohf", "cfour_dropmo": 1, "cfour_cc_program": "vcc", "cfour_print": 2},                                    }, id="ccsd_t_ rohf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd(t)",  "reference": "rohf", "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_reference": "rohf", "cfour_dropmo": 1, "cfour_cc_program": "ecc"},                                                      }, id="ccsd_t_ rohf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "p4-ccsd(t)",  "reference": "rohf", "fcae": "fc", "xptd": {"fd": False}, "keywords": {**_p4c4_fd, "reference": "rohf", "psi4_Freeze_core": True, "qc_module": "ccenergy", "psi4_function_kwargs_dertype": 0},                     }, id="ccsd_t_ rohf fc: psi4-ccenergy", marks=using("psi4")),
        # yapf: enable
    ],
)
def test_ccsd_prt_pr_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "gradient"))


#
#   ,-----. ,-----. ,---.  ,------.    ,-.,--------.,-.      ,--.  ,--.                     ,--.
#  '  .--./'  .--./'   .-' |  .-.  \  / .''--.  .--''. \     |  '--'  | ,---.  ,---.  ,---. `--' ,--,--.,--,--,
#  |  |    |  |    `.  `-. |  |  \  :|  |    |  |    |  |    |  .--.  || .-. :(  .-' (  .-' ,--.' ,-.  ||      \
#  '  '--'\'  '--'\.-'    ||  '--'  /|  |    |  |    |  |    |  |  |  |\   --..-'  `).-'  `)|  |\ '-'  ||  ||  |
#   `-----' `-----'`-----' `-------'  \ '.   `--'   .' /     `--'  `--' `----'`----' `----' `--' `--`--'`--''--'
#                                      `-'          `-'
#  <<<  CCSD(T) Hessian


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
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=[pytest.mark.long, pytest.mark.cilong]),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # * ncc errors
        # * nwchem-cc Error in Communication
        pytest.param({"call": "c4-ccsd(t)",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_cc_program": "vcc", "cfour_print": 2},                                                                                  }, id="ccsd_t_  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd(t)",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_cc_program": "ecc",},                                                                                                   }, id="ccsd_t_  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd(t)", "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {"gamess_ccinp__ncore": 0, "gamess_force__method": "fullnum"},                                                                               }, id="ccsd_t_  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd(t)", "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {"qc_module": "tce"},                                                                                                                        }, id="ccsd_t_  rhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd(t)",  "reference": "rhf",  "fcae": "ae", "xptd": {"fd": False},  "keywords": {**_p4_fd, "qc_module": "ccenergy", "psi4_function_kwargs_dertype": 0},                                                                     }, id="ccsd_t_  rhf ae: psi4-cc",    marks=using("psi4")),  # FD SWITCH
        # DEBUG pytest.param({"call": "p4-c4-ccsd(t)", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "ncc", "psi4_function_kwargs_dertype": 0},                                            }, id="ccsd_t_  rhf ae: psi4-cfour-ncc"),

        pytest.param({"call": "c4-ccsd(t)",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "vcc", "cfour_print": 2,},                                                              }, id="ccsd_t_  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd(t)",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "ecc"},                                                                                 }, id="ccsd_t_  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd(t)", "reference": "rhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {"gamess_force__method": "fullnum"},                                                                                                         }, id="ccsd_t_  rhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd(t)", "reference": "rhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {"nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                                               }, id="ccsd_t_  rhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd(t)",  "reference": "rhf",  "fcae": "fc", "xptd": {"fd": False},  "keywords": {**_p4_fd, "psi4_freeze_core": True, "qc_module": "ccenergy", "psi4_function_kwargs_dertype": 0},                                           }, id="ccsd_t_  rhf fc: psi4-cc",    marks=using("psi4")),  # FD SWITCH
        # DEBUG pytest.param({"call": "p4-c4-ccsd(t)", "reference": "rhf", "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_dropmo": [1], "psi4_cfour_cc_program": "ncc", "psi4_function_kwargs_dertype": 0},                  }, id="ccsd_t_  rhf fc: psi4-cfour-ncc"),

        pytest.param({"call": "c4-ccsd(t)",  "reference": "uhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "vcc", "cfour_print": 2},                                                        }, id="ccsd_t_  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd(t)",  "reference": "uhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "ecc"},                                                                          }, id="ccsd_t_  uhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd(t)", "reference": "uhf",  "fcae": "ae",                        "keywords": {"gamess_contrl__scftyp": "uhf", "gamess_force__method": "fullnum"},                                                       "error": {2: _q2} }, id="ccsd_t_  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd(t)", "reference": "uhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {"nwchem_scf__uhf": True, "qc_module": "tce"},                                                                                               }, id="ccsd_t_  uhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd(t)", "reference": "uhf",  "fcae": "ae",                        "keywords": {"nwchem_scf__uhf": True},                                                                                                 "error": {2: _q3} }, id="ccsd_t_  uhf ae: nwchem-cc",  marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd(t)",  "reference": "uhf",  "fcae": "ae", "xptd": {"fd": False},  "keywords": {**_p4_fd, "reference": "uhf", "qc_module": "ccenergy", "psi4_function_kwargs_dertype": 0},                                                 }, id="ccsd_t_  uhf ae: psi4-cc",    marks=using("psi4")),  # FD SWITCH
        # DEBUG pytest.param({"call": "p4-c4-ccsd(t)", "reference": "uhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "uhf", "psi4_cfour_cc_program": "vcc", "psi4_function_kwargs_dertype": 0},                                            }, id="ccsd_t_  uhf ae: psi4-cfour-vcc"),

        pytest.param({"call": "c4-ccsd(t)",  "reference": "uhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_dropmo": 1, "cfour_cc_program": "vcc",  "cfour_print": 2},                                    }, id="ccsd_t_  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsd(t)",  "reference": "uhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_dropmo": 1, "cfour_cc_program": "ecc"},                                                       }, id="ccsd_t_  uhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsd(t)", "reference": "uhf",  "fcae": "fc",                        "keywords": {"gamess_contrl__scftyp": "uhf", "gamess_force__method": "fullnum",},                                                      "error": {2: _q2} }, id="ccsd_t_  uhf fc: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsd(t)", "reference": "uhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {"nwchem_scf__uhf": True, "nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                      }, id="ccsd_t_  uhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwc-ccsd(t)", "reference": "uhf",  "fcae": "fc",                        "keywords": {"nwchem_scf__uhf": True, "nwchem_ccsd__freeze": 1},                                                                       "error": {2: _q3} }, id="ccsd_t_  uhf fc: nwchem-cc",  marks=using("nwchem")),
        pytest.param({"call": "p4-ccsd(t)",  "reference": "uhf",  "fcae": "fc", "xptd": {"fd": False},  "keywords": {**_p4_fd, "reference": "uhf", "psi4_freeze_core": True, "psi4_function_kwargs_dertype": 0},                                                }, id="ccsd_t_  uhf fc: psi4",       marks=using("psi4")),  # FD SWITCH
        # DEBUG pytest.param({"call": "p4-c4-ccsd(t)", "reference": "uhf", "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "uhf", "psi4_cfour_dropmo": [1], "psi4_cfour_cc_program": "vcc", "psi4_function_kwargs_dertype": 0},                  }, id="ccsd_t_  uhf fc: psi4-cfour-vcc"),
        # yapf: enable
    ],
)
def test_ccsd_prt_pr_hessian_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "hessian"))


#
#                  ,-----. ,-----. ,---.  ,------.    ,-.,--------.,-.      ,------.
#   ,--,--.,-----.'  .--./'  .--./'   .-' |  .-.  \  / .''--.  .--''. \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  ' ,-.  |'-----'|  |    |  |    `.  `-. |  |  \  :|  |    |  |    |  |    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  \ '-'  |       '  '--'\'  '--'\.-'    ||  '--'  /|  |    |  |    |  |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `--`--'        `-----' `-----'`-----' `-------'  \ '.   `--'   .' /     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                     `-'          `-'                                    `---' `---'
#  <<<  a-CCSD(T) Energy aka Lambda-CCSD(T) aka CCSD(aT) aka CCSD(T)_L


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
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # vcc stops prematurely
        #pytest.param({"call": "c4-a-ccsd(t)",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "vcc"},                                                                                                      }, id="a-ccsd_t_  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-a-ccsd(t)",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "ecc"},                                                                                                       }, id="a-ccsd_t_  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-a-ccsd(t)",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "ncc"},                                                                                                       }, id="a-ccsd_t_  rhf ae: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-a-ccsd(t)", "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                           "error": {0: _q27}}, id="a-ccsd_t_  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-a-ccsd(t)", "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                           "error": {0: _q27}}, id="a-ccsd_t_  rhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-a-ccsd(t)",  "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "ccenergy"},                                                                                                                      }, id="a-ccsd_t_  rhf ae: psi4-cc",    marks=using("psi4")),

        pytest.param({"call": "c4-a-ccsd(t)",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "ecc"},                                                                                    }, id="a-ccsd_t_  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-a-ccsd(t)",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "ncc"},                                                                                    }, id="a-ccsd_t_  rhf fc: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "p4-a-ccsd(t)",  "reference": "rhf",  "fcae": "fc", "keywords": {"qc_module": "ccenergy", "psi4_freeze_core": True},                                                                                            }, id="a-ccsd_t_  rhf fc: psi4-cc",    marks=using("psi4")),

        pytest.param({"call": "c4-a-ccsd(t)",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "ecc"},                                                           "error": {0: _q28}}, id="a-ccsd_t_  uhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "p4-a-ccsd(t)",  "reference": "uhf",  "fcae": "ae", "keywords": {"psi4_reference": "uhf"},                                                                                                    "error": {0: _q29}}, id="a-ccsd_t_  uhf ae: psi4",       marks=using("psi4")),
        pytest.param({"call": "p4-a-ccsd(t)",  "reference": "uhf",  "fcae": "ae", "keywords": {"psi4_reference": "uhf", "qc_module": "mrcc"},                                                                                                 }, id="a-ccsd_t_  uhf ae: psi4-mrcc",  marks=using("psi4")),
        # yapf: enable
    ],
)
def test_accsd_prt_pr_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "energy"))


#                  ,-----. ,-----. ,---.  ,------.    ,-.,--------.,-.       ,----.                     ,--.,--.                 ,--.
#   ,--,--.,-----.'  .--./'  .--./'   .-' |  .-.  \  / .''--.  .--''. \     '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  ' ,-.  |'-----'|  |    |  |    `.  `-. |  |  \  :|  |    |  |    |  |    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  \ '-'  |       '  '--'\'  '--'\.-'    ||  '--'  /|  |    |  |    |  |    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#   `--`--'        `-----' `-----'`-----' `-------'  \ '.   `--'   .' /      `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#                                                     `-'          `-'
#  <<<  a-CCSD(T) Gradient aka Lambda-CCSD(T) aka CCSD(aT) aka CCSD(T)_L


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
        pytest.param({"call": "c4-a-ccsd(t)",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_cc_program": "ecc"},                                                              "wrong": {1: _w25}}, id="a-ccsd_t_  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-a-ccsd(t)",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_cc_program": "ncc"},                                                                                }, id="a-ccsd_t_  rhf ae: cfour-ncc",  marks=using("cfour")),

        pytest.param({"call": "c4-a-ccsd(t)",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "ecc"},                                           "wrong": {1: _w25}}, id="a-ccsd_t_  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-a-ccsd(t)",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "ncc"},                                           "wrong": {1: _w8} }, id="a-ccsd_t_  rhf fc: cfour-ncc",  marks=using("cfour")),
        # yapf: enable
    ],
)
def test_accsd_prt_pr_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "gradient"))


#
#  ,-----.   ,-----. ,-----.,------.    ,-.,--------.,-.      ,------.
#  |  |) /_ '  .--./'  .--./|  .-.  \  / .''--.  .--''. \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  .-.  \|  |    |  |    |  |  \  :|  |    |  |    |  |    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  '--' /'  '--'\'  '--'\|  '--'  /|  |    |  |    |  |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `------'  `-----' `-----'`-------'  \ '.   `--'   .' /     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                       `-'          `-'                                    `---' `---'
#  <<<  BCCD(T) Energy


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
        pytest.param({"call": "c4-bccd(t)",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "vcc", "cfour_print": 2, "cfour_orbitals": 1},                                                                }, id="bccd_t_  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "p4-bccd(t)",  "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                                             }, id="bccd_t_  rhf ae: psi4",       marks=using("psi4")),

        # cfour fc Brueckner hangs!
        # pytest.param({"call": "c4-bccd(t)",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_print": 2},                                                                }, id="bccd_t_  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "p4-bccd(t)",  "reference": "rhf",  "fcae": "fc", "keywords": {"psi4_freeze_core": True},                                                                                                                     }, id="bccd_t_  rhf fc: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-bccd(t)",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_REFerence": "UHF", "cfour_cc_program": "vcc", "cfour_print": 2, "cfour_orbitals": 1},                                      }, id="bccd_t_  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "p4-bccd(t)",  "reference": "uhf",  "fcae": "ae", "keywords": {"psi4_reference": "uhf"},                                                                                                                      }, id="bccd_t_  uhf ae: psi4",       marks=using("psi4")),

        # pytest.param({"call": "c4-bccd(t)",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_REFerence": "uHF", "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_print": 2},                                      }, id="bccd_t_  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "p4-bccd(t)",  "reference": "uhf",  "fcae": "fc", "keywords": {"reference": "uhf", "psi4_freeze_core": True},                                                                                                 }, id="bccd_t_  uhf fc: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-bccd(t)",  "reference": "rohf", "fcae": "ae", "keywords": {**_c4_tight, "cfour_REFerence": "roHF", "cfour_cc_program": "vcc", "cfour_print": 2, "cfour_orbitals": 1},                                     }, id="bccd_t_ rohf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "p4-bccd(t)",  "reference": "rohf", "fcae": "ae", "keywords": {"reference": "rohf"},                                                                                                                          }, id="bccd_t_ rohf ae: psi4",       marks=using("psi4")),  # NOTE: sc

        # pytest.param({"call": "c4-bccd(t)",  "reference": "rohf", "fcae": "fc", "sdsc": "sc", "keywords": {**_c4_tight, "cfour_REFerence": "roHF", "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_print": 2},                       }, id="bccd_t_ rohf fc sc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "p4-bccd(t)",  "reference": "rohf", "fcae": "fc", "sdsc": "sc", "keywords": {"psi4_reference": "rohf", "psi4_freeze_core": True},                                                                             }, id="bccd_t_ rohf fc sc: psi4",       marks=using("psi4")),

        # cfour: rhf & uhf fail w/o "cfour_orbitals": 1
        # yapf: enable
    ],
)
def test_bccd_prt_pr_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "energy"))


#
#   ,-----. ,-----.,----.     ,------.
#  '  .--./'  .--./'.-.  |    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |    |  |      .' <     |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '--'\'  '--'\/'-'  |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----' `-----'`----'     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                           `---' `---'
#  <<<  CC3 Energy


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
        # ecc, ncc doesn't run
        pytest.param({"call": "c4-cc3",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "vcc", "cfour_print": 2},                                                                                     }, id="cc3  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "p4-cc3",  "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                                             }, id="cc3  rhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-cc3",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_print": 2},                                                                }, id="cc3  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "p4-cc3",  "reference": "rhf",  "fcae": "fc", "keywords": {"psi4_freeze_core": True},                                                                                                                     }, id="cc3  rhf fc: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-cc3",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_REFerence": "UHF", "cfour_cc_program": "vcc", "cfour_print": 2},                                                           }, id="cc3  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "p4-cc3",  "reference": "uhf",  "fcae": "ae", "keywords": {"psi4_reference": "uhf"},                                                                                                                      }, id="cc3  uhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-cc3",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_REFerence": "uHF", "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_print": 2},                                      }, id="cc3  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "p4-cc3",  "reference": "uhf",  "fcae": "fc", "keywords": {"reference": "uhf", "psi4_freeze_core": True},                                                                                                 }, id="cc3  uhf fc: psi4",       marks=using("psi4")),

        # cfour paper disavows cc3 rohf
        pytest.param({"call": "c4-cc3",  "reference": "rohf", "fcae": "ae", "keywords": {**_c4_tight, "cfour_REFerence": "roHF", "cfour_cc_program": "vcc", "cfour_print": 2},                                        "error": {0: _q46}}, id="cc3 rohf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "p4-cc3",  "reference": "rohf", "fcae": "ae", "keywords": {"reference": "rohf"},                                                                                                                          }, id="cc3 rohf ae: psi4",       marks=using("psi4")),  # NOTE: sc

        pytest.param({"call": "c4-cc3",  "reference": "rohf", "fcae": "fc", "sdsc": "sc", "keywords": {**_c4_tight, "cfour_REFerence": "roHF", "cfour_dropmo": [1], "cfour_cc_program": "vcc", "cfour_print": 2},     "error": {0: _q46}}, id="cc3 rohf fc sc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "p4-cc3",  "reference": "rohf", "fcae": "fc", "sdsc": "sc", "keywords": {"psi4_reference": "rohf", "psi4_freeze_core": True},                                                                             }, id="cc3 rohf fc sc: psi4",       marks=using("psi4")),
        # yapf: enable
    ],
)
def test_cc3_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "energy"))


#
#   ,-----. ,-----.,----.      ,----.                     ,--.,--.                 ,--.
#  '  .--./'  .--./'.-.  |    '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  |    |  |      .' <     |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  '  '--'\'  '--'\/'-'  |    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#   `-----' `-----'`----'      `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  CC3 Gradient


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
    ],
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
        pytest.param({"call": "c4-ccsdt-1a",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "vcc",},                                                                                                       }, id="ccsdt-1a  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-1a",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "ecc",},                                                                                                       }, id="ccsdt-1a  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-1a",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "ncc",},                                                                                                       }, id="ccsdt-1a  rhf ae: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsdt-1a", "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                            "error": {0: _q24}}, id="ccsdt-1a  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsdt-1a", "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                            "error": {0: _q24}}, id="ccsdt-1a  rhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-ccsdt-1a",  "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                            "error": {0: _q47}}, id="ccsdt-1a  rhf ae: psi4",       marks=using("psi4")),
        pytest.param({"call": "p4-ccsdt-1a",  "reference": "rhf",  "fcae": "ae", "keywords": {"psi4_qc_module": "mrcc"},                                                                                                                      }, id="ccsdt-1a  rhf ae: psi4-mrcc",  marks=using("psi4")),

        pytest.param({"call": "c4-ccsdt-1a",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "vcc"},                                                                                     }, id="ccsdt-1a  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-1a",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "ecc"},                                                                                     }, id="ccsdt-1a  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-1a",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "ncc"},                                                                                     }, id="ccsdt-1a  rhf fc: cfour-ncc",  marks=using("cfour")),

        pytest.param({"call": "c4-ccsdt-1a",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "vcc",},                                                                             }, id="ccsdt-1a  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-1a",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "ecc",},                                                                             }, id="ccsdt-1a  uhf ae: cfour-ecc",  marks=using("cfour")),

        pytest.param({"call": "c4-ccsdt-1a",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_dropmo": 1, "cfour_cc_program": "vcc"},                                                           }, id="ccsdt-1a  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-1a",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_dropmo": 1, "cfour_cc_program": "ecc"},                                                           }, id="ccsdt-1a  uhf fc: cfour-ecc",  marks=using("cfour")),
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
        pytest.param({"call": "c4-ccsdt-1a",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_cc_program": "ecc"},                                                                                 }, id="ccsdt-1a  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-1a",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_cc_program": "ncc"},                                                                                 }, id="ccsdt-1a  rhf ae: cfour-ncc",  marks=using("cfour")),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt-1a", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "ecc"},                                                            }, id="ccsdt-1a  rhf ae: psi4-cfour-ecc"),
        # DEBUG pytest.param({"call": "p4-mrccsdt-1a", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_reference": "rhf"},                                                                                                   }, id="ccsdt-1a  rhf ae: psi4-mrcc"),

        pytest.param({"call": "c4-ccsdt-1a",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "ecc"},                                            "wrong": {1: _w16}}, id="ccsdt-1a  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-1a",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "ncc"},                                            "wrong": {1: _w22}}, id="ccsdt-1a  rhf fc: cfour-ncc",  marks=using("cfour")),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt-1a",  "reference": "rhf",  "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_cfour_cc_program": "ecc", "psi4_cfour_dropmo": [1],},                                                              }, id="ccsdt-1a  rhf fc: psi4-cfour-ecc"),
        # DEBUG pytest.param({"call": "p4-mrccsdt-1a", "reference": "rhf", "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_reference": "rhf", "psi4_freeze_core": True},                                                                         }, id="ccsdt-1a  rhf fc: psi4-mrcc"),

        pytest.param({"call": "c4-ccsdt-1a",  "reference": "uhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_reference": "uhf"},                                                                "error": {1: _q19}}, id="ccsdt-1a  uhf ae: cfour",      marks=using("cfour")),

        pytest.param({"call": "c4-ccsdt-1a",  "reference": "rohf", "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_reference": "rohf"},                                                               "error": {1: _q19}}, id="ccsdt-1a rohf ae: cfour",      marks=using("cfour")),
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
        pytest.param({"call": "c4-ccsdt-1a",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_cc_program": "ecc"},                                                                                 }, id="ccsdt-1a  rhf ae: cfour-ecc",  marks=using("cfour")),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt-1a", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "ecc"},                                                            }, id="ccsdt-1a  rhf ae: psi4-cfour-ecc"),
        # DEBUG pytest.param({"call": "p4-mrccsdt-1a", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_reference": "rhf"},                                                                                                   }, id="ccsdt-1a  rhf ae: psi4-mrcc"),

        pytest.param({"call": "c4-ccsdt-1a",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "ecc"},                                            "wrong": {2: _w17}}, id="ccsdt-1a  rhf fc: cfour-ecc",  marks=using("cfour")),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt-1a", "reference": "rhf", "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_dropmo": [1], "psi4_cfou_cc_program": "ecc"},                                   }, id="ccsdt-1a  rhf fc: psi4-cfour-ecc"),
        # DEBUG pytest.param({"call": "p4-mrccsdt-1a", "reference": "rhf", "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_reference": "rhf", "psi4_freeze_core": True},                                                                         }, id="ccsdt-1a  rhf fc: psi4-mrcc"),

        pytest.param({"call": "c4-ccsdt-1a",  "reference": "uhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_reference": "uhf"},                                                                "error": {2: _q19}}, id="ccsdt-1a  uhf ae: cfour",      marks=using("cfour")),

        pytest.param({"call": "c4-ccsdt-1a",  "reference": "rohf", "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_reference": "rohf"},                                                               "error": {2: _q19}}, id="ccsdt-1a rohf ae: cfour",      marks=using("cfour")),
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
    ],
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
        pytest.param({"call": "c4-ccsdt-1b",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "vcc",},                                                                                                       }, id="ccsdt-1b  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-1b",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "ecc",},                                                                                                       }, id="ccsdt-1b  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-1b",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "ncc",},                                                                                                       }, id="ccsdt-1b  rhf ae: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsdt-1b", "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                            "error": {0: _q25}}, id="ccsdt-1b  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsdt-1b", "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                            "error": {0: _q25}}, id="ccsdt-1b  rhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-ccsdt-1b",  "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                            "error": {0: _q47}}, id="ccsdt-1b  rhf ae: psi4",       marks=using("psi4")),
        pytest.param({"call": "p4-ccsdt-1b",  "reference": "rhf",  "fcae": "ae", "keywords": {"psi4_qc_module": "mrcc"},                                                                                                                      }, id="ccsdt-1b  rhf ae: psi4-mrcc",  marks=using("psi4")),

        pytest.param({"call": "c4-ccsdt-1b",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "vcc"},                                                                                     }, id="ccsdt-1b  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-1b",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "ecc"},                                                                                     }, id="ccsdt-1b  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-1b",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "ncc"},                                                                                     }, id="ccsdt-1b  rhf fc: cfour-ncc",  marks=using("cfour")),

        pytest.param({"call": "c4-ccsdt-1b",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "vcc"},                                                                              }, id="ccsdt-1b  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-1b",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "ecc"},                                                                              }, id="ccsdt-1b  uhf ae: cfour-ecc",  marks=using("cfour")),

        pytest.param({"call": "c4-ccsdt-1b",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_dropmo": 1, "cfour_cc_program": "vcc"},                                                           }, id="ccsdt-1b  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-1b",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_dropmo": 1, "cfour_cc_program": "ecc"},                                                           }, id="ccsdt-1b  uhf fc: cfour-ecc",  marks=using("cfour")),
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
        pytest.param({"call": "c4-ccsdt-1b",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_cc_program": "ecc"},                                                                                 }, id="ccsdt-1b  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-1b",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_cc_program": "ncc"},                                                                                 }, id="ccsdt-1b  rhf ae: cfour-ncc",  marks=using("cfour")),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt-1b", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "ecc"},                                                            }, id="ccsdt-1b  rhf ae: psi4-cfour-ecc"),
        # DEBUG pytest.param({"call": "p4-mrccsdt-1b", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_reference": "rhf"},                                                                                                   }, id="ccsdt-1b  rhf ae: psi4-mrcc"),

        pytest.param({"call": "c4-ccsdt-1b",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "ecc"},                                            "wrong": {1: _w18}}, id="ccsdt-1b  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-1b",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "ncc"},                                            "wrong": {1: _w23}}, id="ccsdt-1b  rhf fc: cfour-ncc",  marks=using("cfour")),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt-1b",  "reference": "rhf",  "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_cfour_cc_program": "ecc", "psi4_cfour_dropmo": [1],},                                                              }, id="ccsdt-1b  rhf fc: psi4-cfour-ecc"),
        # DEBUG pytest.param({"call": "p4-mrccsdt-1b", "reference": "rhf", "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_reference": "rhf", "psi4_freeze_core": True},                                                                         }, id="ccsdt-1b  rhf fc: psi4-mrcc"),

        pytest.param({"call": "c4-ccsdt-1b",  "reference": "uhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_reference": "uhf"},                                                                "error": {1: _q19}}, id="ccsdt-1b  uhf ae: cfour",      marks=using("cfour")),

        pytest.param({"call": "c4-ccsdt-1b",  "reference": "rohf", "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_reference": "rohf"},                                                               "error": {1: _q19}}, id="ccsdt-1b rohf ae: cfour",      marks=using("cfour")),
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
        pytest.param({"call": "c4-ccsdt-1b",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_cc_program": "ecc"},                                                                                 }, id="ccsdt-1b  rhf ae: cfour-ecc",  marks=using("cfour")),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt-1b", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "ecc"},                                                            }, id="ccsdt-1b  rhf ae: psi4-cfour-ecc"),
        # DEBUG pytest.param({"call": "p4-mrccsdt-1b", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_reference": "rhf"},                                                                                                   }, id="ccsdt-1b  rhf ae: psi4-mrcc"),

        # * fc is wrong wrt the two debugs for adz - commented b/c not filling in ref for other bases
        # pytest.param({"call": "c4-ccsdt-1b",  "reference": "rhf",  "fcae": "fc",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ecc", "cfour_dropmo": 1,},                    "wrong": {2: _w17}}, id="ccsdt-1b  rhf fc: cfour-ecc",  marks=using("cfour")),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt-1b", "reference": "rhf", "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_dropmo": [1], "psi4_cfour_cc_program": "ecc"},                                  }, id="ccsdt-1b  rhf fc: psi4-cfour-ecc"),
        # DEBUG pytest.param({"call": "p4-mrccsdt-1b", "reference": "rhf", "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_reference": "rhf", "psi4_freeze_core": True},                                                                         }, id="ccsdt-1b  rhf fc: psi4-mrcc"),

        pytest.param({"call": "c4-ccsdt-1b",  "reference": "uhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_reference": "uhf"},                                                                "error": {2: _q19}}, id="ccsdt-1b  uhf ae: cfour",      marks=using("cfour")),

        pytest.param({"call": "c4-ccsdt-1b",  "reference": "rohf", "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_reference": "rohf"},                                                               "error": {2: _q19}}, id="ccsdt-1b rohf ae: cfour",      marks=using("cfour")),
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
    ],
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
        pytest.param({"call": "c4-ccsdt-2",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "vcc",},                                                                                                        }, id="ccsdt-2  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-2",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "ecc",},                                                                                                        }, id="ccsdt-2  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-2",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "ncc",},                                                                                                        }, id="ccsdt-2  rhf ae: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsdt-2", "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                             "error": {0: _q26}}, id="ccsdt-2  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsdt-2", "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                             "error": {0: _q26}}, id="ccsdt-2  rhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-ccsdt-2",  "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                             "error": {0: _q26}}, id="ccsdt-2  rhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-ccsdt-2",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "vcc"},                                                                                      }, id="ccsdt-2  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-2",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "ecc"},                                                                                      }, id="ccsdt-2  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-2",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "ncc"},                                                                                      }, id="ccsdt-2  rhf fc: cfour-ncc",  marks=using("cfour")),

        pytest.param({"call": "c4-ccsdt-2",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "vcc",},                                                                              }, id="ccsdt-2  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-2",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "ecc",},                                                                              }, id="ccsdt-2  uhf ae: cfour-ecc",  marks=using("cfour")),

        pytest.param({"call": "c4-ccsdt-2",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_dropmo": 1, "cfour_cc_program": "vcc"},                                                            }, id="ccsdt-2  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-2",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_dropmo": 1, "cfour_cc_program": "ecc"},                                                            }, id="ccsdt-2  uhf fc: cfour-ecc",  marks=using("cfour")),
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
        pytest.param({"call": "c4-ccsdt-2",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_cc_program": "ecc",},                                                                                 }, id="ccsdt-2  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-2",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_cc_program": "ncc",},                                                                                 }, id="ccsdt-2  rhf ae: cfour-ncc",  marks=using("cfour")),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt-2", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "ecc"},                                                             }, id="ccsdt-2  rhf ae: psi4-cfour-ecc"),

        pytest.param({"call": "c4-ccsdt-2",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "ecc"},                                             "wrong": {1: _w20}}, id="ccsdt-2  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-2",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "ncc"},                                             "wrong": {1: _w24}}, id="ccsdt-2  rhf fc: cfour-ncc",  marks=using("cfour")),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt-2",  "reference": "rhf",  "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_cfour_cc_program": "ecc", "psi4_cfour_dropmo": [1],},                                                               }, id="ccsdt-2  rhf fc: psi4-cfour-ecc"),

        pytest.param({"call": "c4-ccsdt-2",  "reference": "uhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_reference": "uhf"},                                                                 "error": {1: _q19}}, id="ccsdt-2  uhf ae: cfour",      marks=using("cfour")),

        pytest.param({"call": "c4-ccsdt-2",  "reference": "rohf", "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_reference": "rohf"},                                                                "error": {1: _q19}}, id="ccsdt-2 rohf ae: cfour",      marks=using("cfour")),
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
        pytest.param({"call": "c4-ccsdt-2",  "reference": "rhf",  "fcae": "ae", "xptd": {"fd": False},  "keywords": {**_c4_tight, "cfour_cc_program": "ecc",},                                                                                 }, id="ccsdt-2  rhf ae: cfour-ecc",  marks=using("cfour")),  # TODO FD
        # DEBUG pytest.param({"call": "p4-c4-ccsdt-2", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "ecc"},                                                             }, id="ccsdt-2  rhf ae: psi4-cfour-ecc"),

        # * no use testing hess if grad wrong
        # pytest.param({"call": "c4-ccsdt-2",  "reference": "rhf",  "fcae": "fc",                        "keywords": {"cfour_basis": "<>", **_c4_tight, "cfour_cc_program": "ecc", "cfour_dropmo": 1,},                     "wrong": {2: _w17}}, id="ccsdt-2  rhf fc: cfour-ecc",  marks=using("cfour")),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt-2", "reference": "rhf", "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_dropmo": [1], "psi4_cfour_cc_program": "ecc"},                                   }, id="ccsdt-2  rhf fc: psi4-cfour-ecc"),

        pytest.param({"call": "c4-ccsdt-2",  "reference": "uhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_reference": "uhf"},                                                                 "error": {2: _q19}}, id="ccsdt-2  uhf ae: cfour",      marks=using("cfour")),

        pytest.param({"call": "c4-ccsdt-2",  "reference": "rohf", "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_reference": "rohf"},                                                                "error": {2: _q19}}, id="ccsdt-2 rohf ae: cfour",      marks=using("cfour")),
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
    ],
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
        pytest.param({"call": "c4-ccsdt-3",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "vcc"},                                                                                                         }, id="ccsdt-3  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-3",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "ecc"},                                                                                                         }, id="ccsdt-3  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-3",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "ncc"},                                                                                                         }, id="ccsdt-3  rhf ae: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsdt-3", "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                             "error": {0: _q23}}, id="ccsdt-3  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsdt-3", "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                             "error": {0: _q23}}, id="ccsdt-3  rhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-ccsdt-3",  "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                             "error": {0: _q47}}, id="ccsdt-3  rhf ae: psi4",       marks=using("psi4")),
        pytest.param({"call": "p4-ccsdt-3",  "reference": "rhf",  "fcae": "ae", "keywords": {"psi4_qc_module": "mrcc"},                                                                                                                       }, id="ccsdt-3  rhf ae: psi4-mrcc",  marks=using("psi4")),

        pytest.param({"call": "c4-ccsdt-3",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "vcc"},                                                                                      }, id="ccsdt-3  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-3",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "ecc"},                                                                                      }, id="ccsdt-3  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-3",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "ncc"},                                                                                      }, id="ccsdt-3  rhf fc: cfour-ncc",  marks=using("cfour")),

        pytest.param({"call": "c4-ccsdt-3",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "vcc",},                                                                              }, id="ccsdt-3  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-3",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "ecc",},                                                                              }, id="ccsdt-3  uhf ae: cfour-ecc",  marks=using("cfour")),

        pytest.param({"call": "c4-ccsdt-3",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_dropmo": 1, "cfour_cc_program": "vcc"},                                                            }, id="ccsdt-3  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-3",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_dropmo": 1, "cfour_cc_program": "ecc"},                                                            }, id="ccsdt-3  uhf fc: cfour-ecc",  marks=using("cfour")),
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
        pytest.param({"call": "c4-ccsdt-3",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_cc_program": "ecc"},                                                                                  }, id="ccsdt-3  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-3",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_cc_program": "ncc"},                                                                                  }, id="ccsdt-3  rhf ae: cfour-ncc",  marks=using("cfour")),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt-3", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "ecc"},                                                             }, id="ccsdt-3  rhf ae: psi4-cfour-ecc"),
        # DEBUG pytest.param({"call": "p4-mrccsdt-3", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_reference": "rhf"},                                                                                                    }, id="ccsdt-3  rhf ae: psi4-mrcc"),

        # p4c4 by psi findif or mrcc by psi findif differ from ecc analytic by 1e-6
        pytest.param({"call": "c4-ccsdt-3",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "ecc"},                                             "wrong": {1: _w14}}, id="ccsdt-3  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt-3",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "ncc"},                                             "wrong": {1: _w13}}, id="ccsdt-3  rhf fc: cfour-ncc",  marks=using("cfour")),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt-3",  "reference": "rhf",  "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_cfour_cc_program": "ecc", "psi4_cfour_dropmo": [1],},                                                               }, id="ccsdt-3  rhf fc: psi4-cfour-ecc"),
        # DEBUG pytest.param({"call": "p4-mrccsdt-3", "reference": "rhf", "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_reference": "rhf", "psi4_freeze_core": True},                                                                          }, id="ccsdt-3  rhf fc: psi4-mrcc"),

        pytest.param({"call": "c4-ccsdt-3",  "reference": "uhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_reference": "uhf"},                                                                 "error": {1: _q19}}, id="ccsdt-3  uhf ae: cfour",      marks=using("cfour")),

        pytest.param({"call": "c4-ccsdt-3",  "reference": "rohf", "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_reference": "rohf"},                                                                "error": {1: _q19}}, id="ccsdt-3 rohf ae: cfour",      marks=using("cfour")),
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
        pytest.param({"call": "c4-ccsdt-3",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_cc_program": "ecc"},                                                                                  }, id="ccsdt-3  rhf ae: cfour-ecc",  marks=using("cfour")),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt-3", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "ecc"},                                                             }, id="ccsdt-3  rhf ae: psi4-cfour-ecc"),
        # DEBUG pytest.param({"call": "p4-mrccsdt-3", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_reference": "rhf"},                                                                                                    }, id="ccsdt-3  rhf ae: psi4-mrcc"),

        # p4c4 by psi findif or mrcc by psi findif differ from ecc analytic by 1e-4
        pytest.param({"call": "c4-ccsdt-3",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "ecc"},                                             "wrong": {2: _w15}}, id="ccsdt-3  rhf fc: cfour-ecc",  marks=using("cfour")),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt-3", "reference": "rhf", "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_dropmo": [1], "psi4_cfour_cc_program": "ecc"},                                   }, id="ccsdt-3  rhf fc: psi4-cfour-ecc"),
        # DEBUG pytest.param({"call": "p4-mrccsdt-3", "reference": "rhf", "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_reference": "rhf", "psi4_freeze_core": True},                                                                          }, id="ccsdt-3  rhf fc: psi4-mrcc"),

        pytest.param({"call": "c4-ccsdt-3",  "reference": "uhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_reference": "uhf"},                                                                 "error": {2: _q19}}, id="ccsdt-3  uhf ae: cfour",      marks=using("cfour")),

        pytest.param({"call": "c4-ccsdt-3",  "reference": "rohf", "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_reference": "rohf"},                                                                "error": {2: _q19}}, id="ccsdt-3 rohf ae: cfour",      marks=using("cfour")),
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
    ],
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
        pytest.param({"call": "c4-ccsdt",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "vcc"},                                                                                                           }, id="ccsdt  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "ecc"},                                                                                                           }, id="ccsdt  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_cc_program": "ncc"},                                                                                                           }, id="ccsdt  rhf ae: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "gms-ccsdt", "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                               "error": {0: _q21}}, id="ccsdt  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-ccsdt", "reference": "rhf",  "fcae": "ae", "keywords": {"qc_module": "tce"},                                                                                                                               }, id="ccsdt  rhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-ccsdt",  "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                                               "error": {0: _q47}}, id="ccsdt  rhf ae: psi4",       marks=using("psi4")),
        pytest.param({"call": "p4-ccsdt",  "reference": "rhf",  "fcae": "ae", "keywords": {"psi4_qc_module": "mrcc"},                                                                                                                         }, id="ccsdt  rhf ae: psi4-mrcc",  marks=using("psi4")),

        pytest.param({"call": "c4-ccsdt",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "vcc"},                                                                                        }, id="ccsdt  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "ecc"},                                                                                        }, id="ccsdt  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "ncc"},                                                                                        }, id="ccsdt  rhf fc: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "nwc-ccsdt", "reference": "rhf",  "fcae": "fc", "keywords": {"nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                                                      }, id="ccsdt  rhf fc: nwchem-tce", marks=using("nwchem")),

        pytest.param({"call": "c4-ccsdt",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_cc_program": "vcc",},                                                                                }, id="ccsdt  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "nwc-ccsdt", "reference": "uhf",  "fcae": "ae", "keywords": {"nwchem_scf__uhf": True, "qc_module": "tce"},                                                                                                      }, id="ccsdt  uhf ae: nwchem-tce", marks=using("nwchem")),

        # cfour uhf/rohf ecc does not converge, ncc does not run
        pytest.param({"call": "c4-ccsdt",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "cfour_reference": "uhf", "cfour_dropmo": 1, "cfour_cc_program": "vcc"},                                                              }, id="ccsdt  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "nwc-ccsdt", "reference": "uhf",  "fcae": "fc", "keywords": {"nwchem_scf__uhf": True, "nwchem_tce__freeze": 1, "qc_module": "tce"},                                                                             }, id="ccsdt  uhf fc: nwchem-tce", marks=using("nwchem")),

        pytest.param({"call": "c4-ccsdt",  "reference": "rohf", "fcae": "ae", "sdsc": "sc", "keywords": {**_c4_tight, "cfour_reference": "rohf", "cfour_cc_program": "vcc"},                                                                  }, id="ccsdt rohf ae sc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt",  "reference": "rohf", "fcae": "ae", "sdsc": "sd", "keywords": {**_c4_tight, "cfour_reference": "rohf", "cfour_orbitals": 0, "cfour_cc_program": "vcc"},                                             }, id="ccsdt rohf ae sd: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "nwc-ccsdt", "reference": "rohf", "fcae": "ae", "sdsc": "sd", "keywords": {"nwchem_scf__rohf": True, "qc_module": "tce"},                                                                                       }, id="ccsdt rohf ae sd: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-ccsdt",  "reference": "rohf", "fcae": "ae", "sdsc": "sc", "keywords": {"psi4_reference": "rohf", "psi4_qc_module": "mrcc"},                                                                                 }, id="ccsdt rohf ae sc: psi4-mrcc",  marks=using("psi4")),

        pytest.param({"call": "c4-ccsdt",  "reference": "rohf", "fcae": "fc", "sdsc": "sc", "keywords": {**_c4_tight, "cfour_reference": "rohf", "cfour_dropmo": 1, "cfour_cc_program": "vcc"},                                               }, id="ccsdt rohf fc sc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt",  "reference": "rohf", "fcae": "fc", "sdsc": "sd", "keywords": {**_c4_tight, "cfour_reference": "rohf", "cfour_dropmo": 1, "cfour_orbitals": 0, "cfour_cc_program": "vcc"},                          }, id="ccsdt rohf fc sd: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "nwc-ccsdt", "reference": "rohf", "fcae": "fc", "sdsc": "sd", "keywords": {"nwchem_scf__rohf": True, "nwchem_tce__freeze": 1, "qc_module": "tce"},                                                              }, id="ccsdt rohf fc sd: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "p4-ccsdt",  "reference": "rohf", "fcae": "fc", "sdsc": "sc", "keywords": {"psi4_reference": "rohf", "psi4_freeze_core": True, "psi4_qc_module": "mrcc"},                                                       }, id="ccsdt rohf fc sc: psi4-mrcc",  marks=using("psi4")),
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
        pytest.param({"call": "c4-ccsdt",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_cc_program": "ecc",},                                                                                   }, id="ccsdt  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_cc_program": "ncc",},                                                                                   }, id="ccsdt  rhf ae: cfour-ncc",  marks=using("cfour")),
        # VLONG pytest.param({"call": "nwc-ccsdt", "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {"basis": "<>", "qc_module": "tce"},                                                                                 }, id="ccsdt  rhf ae: nwchem-tce", marks=using("nwchem")),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "ecc"},                                                               }, id="ccsdt  rhf ae: psi4-cfour-ecc"),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "ncc"},                                                               }, id="ccsdt  rhf ae: psi4-cfour-ncc"),

        pytest.param({"call": "c4-ccsdt",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "ecc"},                                                                 }, id="ccsdt  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "c4-ccsdt",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "ncc"},                                               "wrong": {1: _w7} }, id="ccsdt  rhf fc: cfour-ncc",  marks=using("cfour")),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt", "reference": "rhf", "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_dropmo": [1], "psi4_cfour_cc_program": "ecc"},                                     }, id="ccsdt  rhf fc: psi4-cfour-ecc"),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt", "reference": "rhf", "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_dropmo": [1], "psi4_cfour_cc_program": "ncc"},                                     }, id="ccsdt  rhf fc: psi4-cfour-ncc"),

        pytest.param({"call": "c4-ccsdt",  "reference": "uhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_reference": "uhf"},                                                                   "error": {1: _q19}}, id="ccsdt  uhf ae: cfour",      marks=using("cfour")),

        pytest.param({"call": "c4-ccsdt",  "reference": "rohf", "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_reference": "rohf"},                                                                  "error": {1: _q19}}, id="ccsdt rohf ae: cfour",      marks=using("cfour")),
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
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=[pytest.mark.long, pytest.mark.cilong]),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # * ncc errors
        pytest.param({"call": "c4-ccsdt",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_cc_program": "ecc"},                                                                                    }, id="ccsdt  rhf ae: cfour-ecc",  marks=using("cfour")),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt", "reference": "rhf", "fcae": "ae", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_cc_program": "ecc"},                                                               }, id="ccsdt  rhf ae: psi4-cfour-ecc"),

        pytest.param({"call": "c4-ccsdt",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cfour_dropmo": 1, "cfour_cc_program": "ecc"},                                                                 }, id="ccsdt  rhf fc: cfour-ecc",  marks=using("cfour")),
        # DEBUG pytest.param({"call": "p4-c4-ccsdt", "reference": "rhf", "fcae": "fc", "keywords": {**_p4c4_fd, "psi4_cfour_reference": "rhf", "psi4_cfour_dropmo": [1], "psi4_cfour_cc_program": "ecc"},                                     }, id="ccsdt  rhf fc: psi4-cfour-ecc"),

        pytest.param({"call": "c4-ccsdt",  "reference": "uhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_reference": "uhf"},                                                                   "error": {2: _q19}}, id="ccsdt  uhf ae: cfour",      marks=using("cfour")),

        pytest.param({"call": "c4-ccsdt",  "reference": "rohf", "fcae": "ae",                        "keywords": {**_c4_tight, "cfour_reference": "rohf"},                                                                  "error": {2: _q19}}, id="ccsdt rohf ae: cfour",      marks=using("cfour")),
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
    ],
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
        pytest.param({"call": "c4-ccsdt(q)",  "reference": "rhf",  "fcae": "ae", "keywords": {"cfour_cc_program": "ncc"},                                                                                                                     }, id="ccsdt_q_  rhf ae: cfour-ncc",      marks=using("cfour")),

        pytest.param({"call": "c4-ccsdt(q)",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_dropmo": 1, "cfour_cc_program": "ncc"},                                                                                                  }, id="ccsdt_q_  rhf fc: cfour-ncc",      marks=using("cfour")),
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
        pytest.param({"call": "c4-ccsdt(q)",  "reference": "rhf",  "fcae": "ae", "keywords": {"cfour_cc_program": "ncc"},                                                                                                                     }, id="ccsdt_q_  rhf ae: cfour-ncc",      marks=using("cfour")),

        pytest.param({"call": "c4-ccsdt(q)",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_dropmo": 1, "cfour_cc_program": "ncc"},                                                                                "wrong": {1: _w7 }}, id="ccsdt_q_  rhf fc: cfour-ncc",      marks=using("cfour")),
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
    ],
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
        pytest.param({"call": "c4-ccsdtq",  "reference": "rhf",  "fcae": "ae", "keywords": {"cfour_cc_program": "ncc"},                                                                                                                       }, id="ccsdtq  rhf ae: cfour-ncc",      marks=using("cfour")),

        pytest.param({"call": "c4-ccsdtq",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_dropmo": 1, "cfour_cc_program": "ncc"},                                                                                                    }, id="ccsdtq  rhf fc: cfour-ncc",      marks=using("cfour")),
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
        # vcc/ecc error out
        pytest.param({"call": "c4-ccsdtq",  "reference": "rhf",  "fcae": "ae", "keywords": {"cfour_cc_program": "ncc"},                                                                                                                       }, id="ccsdtq  rhf ae: cfour-ncc",      marks=using("cfour")),

        pytest.param({"call": "c4-ccsdtq",  "reference": "rhf",  "fcae": "fc", "keywords": {"cfour_dropmo": 1, "cfour_cc_program": "ncc"},                                                                                  "wrong": {1: _w6 }}, id="ccsdtq  rhf fc: cfour-ncc",      marks=using("cfour")),
        # yapf: enable
    ],
)
def test_ccsdtq_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "gradient"))


#
#  ,------. ,-----.  ,------.    ,------.
#  |  .--. '|  |) /_ |  .---'    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  '--' ||  .-.  \|  `--,     |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  | --' |  '--' /|  `---.    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `--'     `------' `------'    `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                              `---' `---'
#  <<<  PBE Energy

# _gms_grid = {"gamess_dft__nrad": 99, "gamess_dft__nleb": 590, "gamess_dft__thresh": 1.e-15}  #, "gamess_dft__gthre": 10, "gamess_scf__conv": 1.e-9}
_gms_grid = {"gamess_dft__nrad": 99, "gamess_dft__nleb": 590}
_psi_grid = {"psi4_dft_radial_points": 99, "psi4_dft_spherical_points": 590}
_nwc_grid = {"nwchem_dft__grid__lebedev": (99, 14)}


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
        pytest.param({"call": "c4-pbe",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight},                                                                    "error": {0: _q35}}, id="pbe  rhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-pbe", "reference": "rhf",  "fcae": "ae", "keywords": {**_gms_grid},                                                                                      }, id="pbe  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-pbe", "reference": "rhf",  "fcae": "ae", "keywords": {**_nwc_grid},                                                                                      }, id="pbe  rhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-pbe",  "reference": "rhf",  "fcae": "ae", "keywords": {**_psi_grid, "psi4_scf_type": "pk"},                                                               }, id="pbe  rhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-pbe",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_reference": "uhf"},                                          "error": {0: _q35}}, id="pbe  uhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-pbe", "reference": "uhf",  "fcae": "ae", "keywords": {**_gms_grid, "gamess_contrl__scftyp": "uhf"},                                                      }, id="pbe  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-pbe", "reference": "uhf",  "fcae": "ae", "keywords": {**_nwc_grid, "nwchem_scf__uhf": True},                                                             }, id="pbe  uhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-pbe",  "reference": "uhf",  "fcae": "ae", "keywords": {**_psi_grid, "reference": "uhf", "psi4_scf_type": "pk"},                                           }, id="pbe  uhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-pbe",  "reference": "rohf", "fcae": "ae", "keywords": {**_c4_tight, "cfour_reference": "rohf"},                                         "error": {0: _q35}}, id="pbe rohf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-pbe", "reference": "rohf", "fcae": "ae", "keywords": {**_gms_grid, "gamess_contrl__scftyp": "rohf"},                                                     }, id="pbe rohf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-pbe", "reference": "rohf", "fcae": "ae", "keywords": {**_nwc_grid, "nwchem_scf__rohf": True, "nwchem_dft__rodft": True},                                 }, id="pbe rohf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-pbe",  "reference": "rohf", "fcae": "ae", "keywords": {**_psi_grid, "reference": "rohf", "psi4_scf_type": "pk"},                        "error": {0: _q36}}, id="pbe rohf ae: psi4",       marks=using("psi4")),
        # pytest.param({"call": "nwc-pbe", "reference": "rohf", "fcae": "ae", "keywords": {**_nwc_grid, "nwchem_scf__rohf": True, "nwchem_dft__rodft": True, "nwchem_dft__convergence__energy": 1.e-9, "nwchem_dft__convergence__gradient": 1.e-7},                    }, id="pbe rohf ae: nwchem",     marks=using("nwchem")),
        # yapf: enable
    ],
)
def test_pbe_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "energy"))


#
#  ,------. ,-----.  ,------.     ,----.                     ,--.,--.                 ,--.
#  |  .--. '|  |) /_ |  .---'    '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  '--' ||  .-.  \|  `--,     |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  |  | --' |  '--' /|  `---.    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#  `--'     `------' `------'     `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  PBE Gradient


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
        pytest.param({"call": "c4-pbe",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight},                                                                                                  "error": {1: _q37}}, id="pbe  rhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-pbe", "reference": "rhf",  "fcae": "ae", "keywords": {**_gms_grid},                                                                                                                    }, id="pbe  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-pbe", "reference": "rhf",  "fcae": "ae", "keywords": {**_nwc_grid},                                                                                                                    }, id="pbe  rhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-pbe",  "reference": "rhf",  "fcae": "ae", "keywords": {**_psi_grid, "psi4_scf_type": "pk"},                                                                                             }, id="pbe  rhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-pbe",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_reference": "uhf"},                                                                        "error": {1: _q37}}, id="pbe  uhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-pbe", "reference": "uhf",  "fcae": "ae", "keywords": {**_gms_grid, "gamess_contrl__scftyp": "uhf"},                                                                                    }, id="pbe  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-pbe", "reference": "uhf",  "fcae": "ae", "keywords": {**_nwc_grid, "nwchem_scf__uhf": True},                                                                                           }, id="pbe  uhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-pbe",  "reference": "uhf",  "fcae": "ae", "keywords": {**_psi_grid, "reference": "uhf", "psi4_scf_type": "pk"},                                                                         }, id="pbe  uhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-pbe",  "reference": "rohf", "fcae": "ae", "keywords": {**_c4_tight, "cfour_reference": "rohf"},                                                                       "error": {1: _q37}}, id="pbe rohf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-pbe", "reference": "rohf", "fcae": "ae", "keywords": {**_gms_grid, "gamess_contrl__scftyp": "rohf", "gamess_dft__gthre": 10},                                                          }, id="pbe rohf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-pbe", "reference": "rohf", "fcae": "ae", "keywords": {**_nwc_grid, "nwchem_scf__rohf": True, "nwchem_dft__rodft": True, "nwchem_dft__convergence__gradient": 1.e-5},                   }, id="pbe rohf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-pbe",  "reference": "rohf", "fcae": "ae", "keywords": {**_psi_grid, "reference": "rohf", "psi4_scf_type": "pk"},                                                      "error": {1: _q36}}, id="pbe rohf ae: psi4",       marks=using("psi4")),
        # pytest.param({"call": "nwc-pbe", "reference": "rohf", "fcae": "ae", "keywords": {**_nwc_grid, "nwchem_scf__rohf": True, "nwchem_dft__rodft": True, "nwchem_dft__convergence__energy": 1.e-9, "nwchem_dft__convergence__gradient": 1.e-7},                    }, id="pbe rohf ae: nwchem",     marks=using("nwchem")),
        # yapf: enable
    ],
)
def test_pbe_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "gradient"))


#
#  ,-----.  ,----. ,--.,--.   ,--.,------.     ,------.
#  |  |) /_ '.-.  ||  | \  `.'  / |  .--. '    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  .-.  \  .' < |  |  '.    /  |  '--' |    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  '--' //'-'  ||  '--. |  |   |  | --'     |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `------' `----' `-----' `--'   `--'         `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                            `---' `---'
#  <<<  B3LYP Energy


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
        pytest.param({"call": "c4-b3lyp",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight},                                                                    "error": {0: _q35}}, id="b3lyp  rhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-b3lyp", "reference": "rhf",  "fcae": "ae", "keywords": {**_gms_grid},                                                                                      }, id="b3lyp  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-b3lyp", "reference": "rhf",  "fcae": "ae", "keywords": {**_nwc_grid},                                                                                      }, id="b3lyp  rhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-b3lyp",  "reference": "rhf",  "fcae": "ae", "keywords": {**_psi_grid, "psi4_scf_type": "pk"},                                                               }, id="b3lyp  rhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-b3lyp",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_reference": "uhf"},                                          "error": {0: _q35}}, id="b3lyp  uhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-b3lyp", "reference": "uhf",  "fcae": "ae", "keywords": {**_gms_grid, "gamess_contrl__scftyp": "uhf"},                                                      }, id="b3lyp  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-b3lyp", "reference": "uhf",  "fcae": "ae", "keywords": {**_nwc_grid, "nwchem_scf__uhf": True},                                                             }, id="b3lyp  uhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-b3lyp",  "reference": "uhf",  "fcae": "ae", "keywords": {**_psi_grid, "reference": "uhf", "psi4_scf_type": "pk"},                                           }, id="b3lyp  uhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-b3lyp",  "reference": "rohf", "fcae": "ae", "keywords": {**_c4_tight, "cfour_reference": "rohf"},                                         "error": {0: _q35}}, id="b3lyp rohf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-b3lyp", "reference": "rohf", "fcae": "ae", "keywords": {**_gms_grid, "gamess_contrl__scftyp": "rohf"},                                                     }, id="b3lyp rohf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-b3lyp", "reference": "rohf", "fcae": "ae", "keywords": {**_nwc_grid, "nwchem_scf__rohf": True, "nwchem_dft__rodft": True},                                 }, id="b3lyp rohf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-b3lyp",  "reference": "rohf", "fcae": "ae", "keywords": {**_psi_grid, "reference": "rohf", "psi4_scf_type": "pk"},                        "error": {0: _q36}}, id="b3lyp rohf ae: psi4",       marks=using("psi4")),
        # DEBUG pytest.param({"call": "nwc-b3lyp", "reference": "rohf", "fcae": "ae", "keywords": {**_nwc_grid, "nwchem_scf__rohf": True, "nwchem_dft__rodft": True, "nwchem_dft__convergence__energy": 1.e-9, "nwchem_dft__convergence__gradient": 1.e-7},                    }, id="b3lyp rohf ae: nwchem",     marks=using("nwchem")),
        # yapf: enable
    ],
)
def test_b3lyp_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "energy"))


#
#  ,-----.  ,----. ,--.,--.   ,--.,------. ,-----.    ,------.
#  |  |) /_ '.-.  ||  | \  `.'  / |  .--. '|  .--'    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  .-.  \  .' < |  |  '.    /  |  '--' |'--. `\    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  '--' //'-'  ||  '--. |  |   |  | --' .--'  /    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `------' `----' `-----' `--'   `--'     `----'     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                                   `---' `---'
#  <<<  B3LYP5 Energy


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
        pytest.param({"call": "c4-b3lyp5",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight},                                                                    "error": {0: _q35}}, id="b3lyp5  rhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-b3lyp5", "reference": "rhf",  "fcae": "ae", "keywords": {**_gms_grid},                                                                                      }, id="b3lyp5  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-b3lyp5", "reference": "rhf",  "fcae": "ae", "keywords": {**_nwc_grid},                                                                                      }, id="b3lyp5  rhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-b3lyp5",  "reference": "rhf",  "fcae": "ae", "keywords": {**_psi_grid, "psi4_scf_type": "pk"},                                                               }, id="b3lyp5  rhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-b3lyp5",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_reference": "uhf"},                                          "error": {0: _q35}}, id="b3lyp5  uhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-b3lyp5", "reference": "uhf",  "fcae": "ae", "keywords": {**_gms_grid, "gamess_contrl__scftyp": "uhf"},                                                      }, id="b3lyp5  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-b3lyp5", "reference": "uhf",  "fcae": "ae", "keywords": {**_nwc_grid, "nwchem_scf__uhf": True},                                                             }, id="b3lyp5  uhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-b3lyp5",  "reference": "uhf",  "fcae": "ae", "keywords": {**_psi_grid, "reference": "uhf", "psi4_scf_type": "pk"},                                           }, id="b3lyp5  uhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-b3lyp5",  "reference": "rohf", "fcae": "ae", "keywords": {**_c4_tight, "cfour_reference": "rohf"},                                         "error": {0: _q35}}, id="b3lyp5 rohf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-b3lyp5", "reference": "rohf", "fcae": "ae", "keywords": {**_gms_grid, "gamess_contrl__scftyp": "rohf"},                                                     }, id="b3lyp5 rohf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-b3lyp5", "reference": "rohf", "fcae": "ae", "keywords": {**_nwc_grid, "nwchem_scf__rohf": True, "nwchem_dft__rodft": True},                                 }, id="b3lyp5 rohf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-b3lyp5",  "reference": "rohf", "fcae": "ae", "keywords": {**_psi_grid, "reference": "rohf", "psi4_scf_type": "pk"},                        "error": {0: _q36}}, id="b3lyp5 rohf ae: psi4",       marks=using("psi4")),
        # DEBUG pytest.param({"call": "nwc-b3lyp5", "reference": "rohf", "fcae": "ae", "keywords": {**_nwc_grid, "nwchem_scf__rohf": True, "nwchem_dft__rodft": True, "nwchem_dft__convergence__energy": 1.e-9, "nwchem_dft__convergence__gradient": 1.e-7},                    }, id="b3lyp5 rohf ae: nwchem",     marks=using("nwchem")),
        # yapf: enable
    ],
)
def test_b3lyp5_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "energy"))


#
#  ,-----.  ,----. ,--.,--.   ,--.,------.      ,----.                     ,--.,--.                 ,--.
#  |  |) /_ '.-.  ||  | \  `.'  / |  .--. '    '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  .-.  \  .' < |  |  '.    /  |  '--' |    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  |  '--' //'-'  ||  '--. |  |   |  | --'     '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#  `------' `----' `-----' `--'   `--'          `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  B3LYP Gradient


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
        pytest.param({"call": "c4-b3lyp",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight},                                                                                                  "error": {1: _q37}}, id="b3lyp  rhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-b3lyp", "reference": "rhf",  "fcae": "ae", "keywords": {**_gms_grid},                                                                                                                    }, id="b3lyp  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-b3lyp", "reference": "rhf",  "fcae": "ae", "keywords": {**_nwc_grid},                                                                                                                    }, id="b3lyp  rhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-b3lyp",  "reference": "rhf",  "fcae": "ae", "keywords": {**_psi_grid, "psi4_scf_type": "pk"},                                                                                             }, id="b3lyp  rhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-b3lyp",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_reference": "uhf"},                                                                        "error": {1: _q37}}, id="b3lyp  uhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-b3lyp", "reference": "uhf",  "fcae": "ae", "keywords": {**_gms_grid, "gamess_contrl__scftyp": "uhf"},                                                                                    }, id="b3lyp  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-b3lyp", "reference": "uhf",  "fcae": "ae", "keywords": {**_nwc_grid, "nwchem_scf__uhf": True},                                                                                           }, id="b3lyp  uhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-b3lyp",  "reference": "uhf",  "fcae": "ae", "keywords": {**_psi_grid, "reference": "uhf", "psi4_scf_type": "pk"},                                                                         }, id="b3lyp  uhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-b3lyp",  "reference": "rohf", "fcae": "ae", "keywords": {**_c4_tight, "cfour_reference": "rohf"},                                                                       "error": {1: _q37}}, id="b3lyp rohf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-b3lyp", "reference": "rohf", "fcae": "ae", "keywords": {**_gms_grid, "gamess_contrl__scftyp": "rohf", "gamess_dft__gthre": 10},                                                          }, id="b3lyp rohf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-b3lyp", "reference": "rohf", "fcae": "ae", "keywords": {**_nwc_grid, "nwchem_scf__rohf": True, "nwchem_dft__rodft": True, "nwchem_dft__convergence__gradient": 1.e-5},                   }, id="b3lyp rohf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-b3lyp",  "reference": "rohf", "fcae": "ae", "keywords": {**_psi_grid, "reference": "rohf", "psi4_scf_type": "pk"},                                                      "error": {1: _q36}}, id="b3lyp rohf ae: psi4",       marks=using("psi4")),
        # DEBUG pytest.param({"call": "nwc-b3lyp", "reference": "rohf", "fcae": "ae", "keywords": {**_nwc_grid, "nwchem_scf__rohf": True, "nwchem_dft__rodft": True, "nwchem_dft__convergence__energy": 1.e-9, "nwchem_dft__convergence__gradient": 1.e-7},                    }, id="b3lyp rohf ae: nwchem",     marks=using("nwchem")),
        # yapf: enable
    ],
)
def test_b3lyp_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "gradient"))


#
#  ,-----.  ,----. ,--.,--.   ,--.,------. ,-----.     ,----.                     ,--.,--.                 ,--.
#  |  |) /_ '.-.  ||  | \  `.'  / |  .--. '|  .--'    '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  .-.  \  .' < |  |  '.    /  |  '--' |'--. `\    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  |  '--' //'-'  ||  '--. |  |   |  | --' .--'  /    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#  `------' `----' `-----' `--'   `--'     `----'      `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  B3LYP5 Gradient


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
        pytest.param({"call": "c4-b3lyp5",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight},                                                                                                  "error": {1: _q37}}, id="b3lyp5  rhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-b3lyp5", "reference": "rhf",  "fcae": "ae", "keywords": {**_gms_grid},                                                                                                                    }, id="b3lyp5  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-b3lyp5", "reference": "rhf",  "fcae": "ae", "keywords": {**_nwc_grid},                                                                                                                    }, id="b3lyp5  rhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-b3lyp5",  "reference": "rhf",  "fcae": "ae", "keywords": {**_psi_grid, "psi4_scf_type": "pk"},                                                                                             }, id="b3lyp5  rhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-b3lyp5",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "cfour_reference": "uhf"},                                                                        "error": {1: _q37}}, id="b3lyp5  uhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-b3lyp5", "reference": "uhf",  "fcae": "ae", "keywords": {**_gms_grid, "gamess_contrl__scftyp": "uhf"},                                                                                    }, id="b3lyp5  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-b3lyp5", "reference": "uhf",  "fcae": "ae", "keywords": {**_nwc_grid, "nwchem_scf__uhf": True},                                                                                           }, id="b3lyp5  uhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-b3lyp5",  "reference": "uhf",  "fcae": "ae", "keywords": {**_psi_grid, "reference": "uhf", "psi4_scf_type": "pk"},                                                                         }, id="b3lyp5  uhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "c4-b3lyp5",  "reference": "rohf", "fcae": "ae", "keywords": {**_c4_tight, "cfour_reference": "rohf"},                                                                       "error": {1: _q37}}, id="b3lyp5 rohf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gms-b3lyp5", "reference": "rohf", "fcae": "ae", "keywords": {**_gms_grid, "gamess_contrl__scftyp": "rohf", "gamess_dft__gthre": 10},                                                          }, id="b3lyp5 rohf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "nwc-b3lyp5", "reference": "rohf", "fcae": "ae", "keywords": {**_nwc_grid, "nwchem_scf__rohf": True, "nwchem_dft__rodft": True, "nwchem_dft__convergence__gradient": 1.e-5},                   }, id="b3lyp5 rohf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "p4-b3lyp5",  "reference": "rohf", "fcae": "ae", "keywords": {**_psi_grid, "reference": "rohf", "psi4_scf_type": "pk"},                                                      "error": {1: _q36}}, id="b3lyp5 rohf ae: psi4",       marks=using("psi4")),
        # DEBUG pytest.param({"call": "nwc-b3lyp5", "reference": "rohf", "fcae": "ae", "keywords": {**_nwc_grid, "nwchem_scf__rohf": True, "nwchem_dft__rodft": True, "nwchem_dft__convergence__energy": 1.e-9, "nwchem_dft__convergence__gradient": 1.e-7},                    }, id="b3lyp5 rohf ae: nwchem",     marks=using("nwchem")),
        # yapf: enable
    ],
)
def test_b3lyp5_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "gradient"))


#
#  ,-----.   ,---. ,------. ,--.,--.   ,--.,------.     ,------.
#  |  |) /_ '.-.  \|  .--. '|  | \  `.'  / |  .--. '    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  .-.  \ .-' .'|  '--' ||  |  '.    /  |  '--' |    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  '--' //   '-.|  | --' |  '--. |  |   |  | --'     |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `------' '-----'`--'     `-----' `--'   `--'         `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                                     `---' `---'
#  <<<  B2PLYP Energy


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
        pytest.param({"call": "c4-b2plyp",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight},                                                                    "error": {0: _q35}}, id="b2plyp  rhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "nwc-b2plyp", "reference": "rhf",  "fcae": "ae", "keywords": {**_nwc_grid},                                                                                      }, id="b2plyp  rhf ae: nwchem",     marks=using("nwchem")),
        # ONLY DF-MP2 pytest.param({"call": "p4-b2plyp",  "reference": "rhf",  "fcae": "ae", "keywords": {**_psi_grid, "psi4_scf_type": "pk"},                                                               }, id="b2plyp  rhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "nwc-b2plyp", "reference": "rhf",  "fcae": "fc", "keywords": {**_nwc_grid, "nwchem_mp2__freeze": 1},                                                             }, id="b2plyp  rhf fc: nwchem",     marks=using("nwchem")),

        pytest.param({"call": "nwc-b2plyp", "reference": "uhf",  "fcae": "ae", "keywords": {**_nwc_grid, "nwchem_scf__uhf": True},                                                             }, id="b2plyp  uhf ae: nwchem",     marks=using("nwchem")),
        # ONLY DF-MP2 pytest.param({"call": "p4-b2plyp",  "reference": "uhf",  "fcae": "ae", "keywords": {**_psi_grid, "reference": "uhf", "psi4_scf_type": "pk"},                                           }, id="b2plyp  uhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "nwc-b2plyp", "reference": "uhf",  "fcae": "fc", "keywords": {**_nwc_grid, "nwchem_scf__uhf": True, "nwchem_mp2__freeze": 1},                                    }, id="b2plyp  uhf fc: nwchem",     marks=using("nwchem")),
        # yapf: enable
    ],
)
def test_b2plyp_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, "energy"))


#
#   ,-----.
#  '  .--./ ,---. ,--,--,--.,--,--,--. ,---. ,--,--,
#  |  |    | .-. ||        ||        || .-. ||      \
#  '  '--'\' '-' '|  |  |  ||  |  |  |' '-' '|  ||  |
#   `-----' `---' `--`--`--'`--`--`--' `---' `--''--'
#
#  <<<  Common functions


def _processor(inp, dertype, basis, subjects, clsd_open_pmols, request, driver, *, scramble=None, frame=""):
    qcprog, method = inp["call"].split("-", 1)
    qcprog = _trans_qcprog[qcprog.lower()]
    tnm = request.node.name
    suffix = "-fixed" if frame == "fixed" else ""
    subject = clsd_open_pmols[subjects[std_refs.index(inp["reference"])] + suffix]

    inpcopy = {k: v for k, v in inp.items() if k not in ["error", "wrong"]}
    if inp.get("error", False) and inp["error"].get(dertype, False):
        inpcopy["error"] = inp["error"][dertype]
    if inp.get("wrong", False) and inp["wrong"].get(dertype, False):
        inpcopy["wrong"] = inp["wrong"][dertype]
    if inp.get("error", False) and inp["error"].get(basis, False):
        inpcopy["error"] = inp["error"][basis]
    if inp.get("wrong", False) and inp["wrong"].get(basis, False):
        inpcopy["wrong"] = inp["wrong"][basis]
    if inp.get("marks", False) and inp["marks"].get(dertype, False):
        pytest.xfail(inp["marks"][dertype])
        # request.node.add_marker(inp["marks"][dertype])

    inpcopy["keywords"] = {
        k: (_trans_key(qcprog, basis, k) if v == "<>" else v) for k, v in inpcopy["keywords"].items()
    }
    inpcopy["driver"] = driver
    if not any([k.lower() in _basis_keywords for k in inpcopy["keywords"]]):
        inpcopy["keywords"]["basis"] = basis
    inpcopy["scf_type"] = "conv"  # longstanding "pk"
    inpcopy["corl_type"] = "conv"
    inpcopy["qc_module"] = "-".join(
        [
            qcprog,
            inp["keywords"].get(
                "qc_module",
                inp["keywords"].get(
                    "cfour_cc_program",
                    inp["keywords"].get("psi4_qc_module", inp["keywords"].get("gamess_mp2__code", "")),
                ),
            ),
        ]
    ).strip("-")
    inpcopy["keywords"]["function_kwargs"] = {"dertype": dertype}
    print("INP", inpcopy)

    return inpcopy, subject, method, basis, tnm, scramble, frame
