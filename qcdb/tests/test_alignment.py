import pytest
from qcengine.programs.tests.standard_suite_ref import std_molecules, std_refs
from qcengine.testing import using

import qcdb

from .standard_suite_runner import runner_asserter
from .test_standard_suite import _processor


@pytest.fixture
def clsd_open_pmols():
    frame_not_important = {
        name[:-4]: qcdb.Molecule.from_string(smol, name=name[:-4])
        for name, smol in std_molecules.items()
        if name.endswith("-xyz")
    }
    frame_part_of_spec = {
        name[:-4] + "-fixed": qcdb.Molecule.from_string(smol + "\nno_com\nno_reorient\n", name=name[:-4])
        for name, smol in std_molecules.items()
        if name.endswith("-xyz")
    }

    return {**frame_not_important, **frame_part_of_spec}


#                                                                                          
#  ,--.  ,--.,------.      ,---.  ,--.,--.                                          ,--.   
#  |  '--'  ||  .---'     /  O  \ |  |`--' ,---. ,--,--, ,--,--,--. ,---. ,--,--, ,-'  '-. 
#  |  .--.  ||  `--,     |  .-.  ||  |,--.| .-. ||      \|        || .-. :|      \'-.  .-' 
#  |  |  |  ||  |`       |  | |  ||  ||  |' '-' '|  ||  ||  |  |  |\   --.|  ||  |  |  |   
#  `--'  `--'`--'        `--' `--'`--'`--'.`-  / `--''--'`--`--`--' `----'`--''--'  `--'   
#                                         `---'                                           
#  <<<  HF Alignment

@pytest.mark.parametrize(
    "scramble",
    # * this parameter alters the input molecule by Cartesian frame and/or atom ordering to imitate real-world inputs.
    # * scramble dictionary are arguments to qcdb.Molecule.scramble() or qcel.models.Molecule.scramble() that computes a shifted/rotated/atom-mapped input molecule from `subjects` below and the transformations to be applied to the reference data.
    # * arguments do_shift=T/F and other boolean values generate random perturbations, so multiple lines or runs makes a fuller test.
    # * specific shifts etc or problematic runs can be reproduced by specifying arrays like the commented example below. note that b/c do_resort is natom-dependent, may need to exclude subjects while debugging.
    # * `id`s below are numbered since `pytest -k` doesn't recognize case and so duplicate entries can be added. the 0, 1, 2, 3, 9 progression is handy for debugging since perturbations are added systematically
    [
        # pytest.param({"do_shift": [ 2.660760432055,  1.477336796939, -2.098045335573], "do_rotate": [[ 0.321861140022,  0.445246880671,  0.835560064745], [-0.447874157747,  0.849136107328, -0.279958229125], [-0.834154749052, -0.28411808546 ,  0.472718487209]], "do_resort": [2, 0, 1]}, id="mTTT3"),
        # pytest.param({"do_shift": False, "do_rotate": False, "do_resort": False}, id="srm0"),
        # pytest.param({"do_shift": True, "do_rotate": False, "do_resort": False}, id="Srm1"),
        # pytest.param({"do_shift": True, "do_rotate": True, "do_resort": False}, id="SRm2"),
        # pytest.param({"do_shift": False, "do_rotate": False, "do_resort": True}, id="srM3"),
        # pytest.param({"do_shift": True, "do_rotate": True, "do_resort": True}, id="SRM8"),
        pytest.param({"do_shift": True, "do_rotate": True, "do_resort": True}, id="SRM9"),
    ],
)
@pytest.mark.parametrize(
    "frame",
    # * this parameter alters the input molecule by fix_com and fix_orientation to imitate user signalling frame matters or not.
    [
        pytest.param("fixed"),  # fix_=True (no_com/no_reorient); atres.mol.geom = atin.mol.geom aka scrambled
        pytest.param("free"),   # fix_=False (def)              ; atres.mol.geom = qcdb (aka psi4) interna orientation
    ],
)
@pytest.mark.parametrize(
    "driver,dertype",
    [
        pytest.param("energy", 0, id="ene0"),
        pytest.param("gradient", 1, id="grd1"),
        pytest.param("hessian", 2, id="hes2"),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    # this parameter, along with rhf/uhf below covers four different molecules, with nat=2-4.
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
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
        # yapf: enable
    ],
)
def test_hf_alignment(inp, scramble, frame, driver, dertype, basis, subjects, clsd_open_pmols, request):
    runner_asserter(*_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, driver, scramble=scramble, frame=frame))
