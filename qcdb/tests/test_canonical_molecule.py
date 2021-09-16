# """
# Tests the DQM compute dispatch module
# """
import pprint

# import re
# import tempfile
# from pathlib import Path

import pytest

# import qcelemental as qcel
from qcelemental.models import AtomicInput, Molecule

import qcengine as qcng
from qcengine.testing import has_program, using
from qcengine.programs.tests.standard_suite_ref import std_molecules, std_suite

import qcdb
from .utils import compare_values

_canonical_methods = [
    ("cfour", {"method": "hf", "basis": "cc-pvdz"}, {"reference": "uhf"}),
    ("gamess", {"method": "hf", "basis": "cc-pvdz"}, {"reference": "uhf"}),
    ("nwchem", {"method": "hf", "basis": "cc-pvdz"}, {"reference": "uhf"}),
    ("psi4", {"method": "hf", "basis": "cc-pvdz"}, {"reference": "uhf", "scf_type": "pk"}),
    #    ("cfour", {"method": "hf", "basis": "cc-pvdz"}, {}),  # program defaults to uhf  # "cfour_reference": "uhf"
    #    ("gamess", {"method": "hf", "basis": "cc-pvdz"}, {"gamess_contrl__scftyp": "uhf"}),
    #    ("nwchem", {"method": "hf", "basis": "cc-pvdz"}, {"nwchem_scf__uhf": True}),
    #    ("psi4", {"method": "hf", "basis": "cc-pvdz"}, {"scf_type": "pk"}),  # harness defaults to uhf  # "psi4_reference": "uhf"
]


# @pytest.fixture
# def clsd_open_pmols():
#    return {name: qcdb.Molecule.from_string(smol, name=name) for name, smol in std_molecules.items()}


@pytest.fixture
def clsd_open_pmols():
    return {
        name[:-4]: Molecule.from_data(smol, name=name[:-4])
        for name, smol in std_molecules.items()
        if name.endswith("-xyz")
    }


@pytest.mark.parametrize(
    "mol_trickery",
    [
        pytest.param({}, id="none"),
        pytest.param(
            # native keywords consistent with molecule BH3+ below
            {
                "cfour": {"cfour_charge": 1},
                "gamess": {"gamess_contrl__icharg": 1},
                "nwchem": {"nwchem_charge": 1},
                "psi4": None,  # no charge keyword in psi
            },
            id="dsl-chg",
        ),
        pytest.param(
            {
                "cfour": {"cfour_multiplicity": 2},
                "gamess": {"gamess_contrl__mult": 2},
                "nwchem": {"nwchem_scf__nopen": 1},
                "psi4": None,  # no multiplicity keyword in psi
            },
            id="dsl-mult",
        ),
        pytest.param(
            # native keywords that CONTRADICT molecule BH3+ below
            {
                "cfour": {"cfour_charge": 0},
                "gamess": {"gamess_contrl__icharg": 0},
                "nwchem": {"nwchem_charge": 0},
                "psi4": None,
            },
            id="dsl-chg-contra",
        ),
        pytest.param(
            {
                "cfour": {"cfour_multiplicity": 4},
                "gamess": {"gamess_contrl__mult": 4},
                "nwchem": {"nwchem_scf__nopen": 3},
                "psi4": None,
            },
            id="dsl-mult-contra",
        ),
    ],
)
@pytest.mark.parametrize("program, model, keywords", _canonical_methods)
def test_keyword_reconciliation_molecule(program, model, keywords, mol_trickery, clsd_open_pmols, request):
    """Ensure molecule handling implemented in harness.

        For available harnesses, run small UHF calc on non-neutral singlet, both through molecule alone
          and with clashing (and non-QCEngine-like) keyword spec.

    #    New Harness Instructions
    #    ------------------------
    #    * Make sure minimal calc is in _canonical_methods above.
    #    * If ``managed_memory=True`` in harness, add regex to ``stdout_ref`` below to check that memory
    #      is specifiable.
    #    * If this test doesn't work, implement or adjust ``config.memory`` in your harness.

    """
    if not has_program(program):
        pytest.skip(f"Program '{program}' not found.")

    harness = qcng.get_program(program)
    molecule = clsd_open_pmols["bh3p"]

    addl_keywords = mol_trickery.get(program, mol_trickery)
    if addl_keywords is None:
        pytest.skip(f"Nothing to test for '{program}' id '{request.node.name}'")
    use_keywords = {**keywords, **addl_keywords}

    #    #  <<  Config
    #
    #    config = qcng.config.get_config(
    #        hostname="something",
    #        local_options={
    #            "ncores": 1,
    #            "nnodes": 1,
    #            "memory": 1.555,
    #        },
    #    )

    #  <<  Run

    inp = AtomicInput(molecule=molecule, driver="energy", model=model, keywords=use_keywords)
    ret = qcdb.compute(inp, program)  # , local_options=config.dict())
    pprint.pprint(ret.dict(), width=200)

    if "contra" in request.node.name:
        assert ret.success is False
        assert "KeywordReconciliationError" in ret.error.error_message
        return

    assert ret.success is True

    #  <<  Reference

    ref_e = std_suite["bh3p_cc-pvdz_pk_uhf_ae_conv"]["HF TOTAL ENERGY"]
    #    stdout_ref = {  # 1.555 GiB = 208708567 quad-words
    #        "cfour": "Allocated    1592 MB of main memory",
    #        "gamess": "208000000 WORDS OF MEMORY AVAILABLE",
    #        "nwchem": r"total    =  2087085\d\d doubles =   1592.3 Mbytes",  # doubles is quad-words. Mbytes is MiB
    #        "psi4": "1592 MiB Core",
    #    }

    #  <<  Test

    assert compare_values(ref_e, ret.return_result, atol=1.0e-6, label="uhf ene")


#    assert config.ncores == 1
#    assert pytest.approx(config.memory, 0.1) == 1.555
#
#    if harness._defaults["managed_memory"] is True:
#        assert re.search(stdout_ref[program], ret.stdout), f"Memory pattern not found: {stdout_ref[program]}"
# assert 0


@pytest.mark.parametrize("program, model, keywords", _canonical_methods)
def hide_test_mode(program, model, keywords, clsd_open_pmols):
    #    """Ensure scratch handling implemented in harness (if applicable).
    #
    #    For available harnesses, run minimal calc at specific scratch directory name (randomly generated
    #      during test) and skip scratch clean-up. Check scratch settings show up in ``TaskConfig``.
    #    For ``scratch``-active harnesses, check that an expected file is written to and left behind in
    #      scratch directory. Check any scratch-related printing in output.
    #
    #    New Harness Instructions
    #    ------------------------
    #    * Make sure minimal calc is in _canonical_methods above.
    #    * If ``scratch=True`` in harness, add single file (preferrably output) glob to ``scratch_sample``
    #      below to check that program scratch is directable.
    #    * If ``scratch=True`` in harness, if scratch directory mentioned in output, add regex to
    #      ``stdout_ref`` below to check that program scratch is directable. Otherwise, add an
    #      always-passing regex.
    #    * If this test doesn't work, implement or adjust ``config.scratch_directory`` and
    #      ``config.scratch_messy`` in your harness.
    #
    #    """
    #    if not has_program(program):
    #        pytest.skip(f"Program '{program}' not found.")

    harness = qcng.get_program(program)
    molecule = clsd_open_pmols["hf"]

    #    #  <<  Config
    #
    #    scratch_directory = tempfile.mkdtemp(suffix="_" + program)
    #
    #    config = qcng.config.get_config(
    #        hostname="something",
    #        local_options={
    #            "scratch_directory": scratch_directory,
    #            "scratch_messy": True,
    #        },
    #    )
    #
    #    #  <<  Run
    #
    #    inp = AtomicInput(molecule=molecule, driver="energy", model=model, keywords=keywords)
    #    ret = qcdb.compute(inp, program, local_options=config.dict())
    #    pprint.pprint(ret.dict(), width=200)
    #    assert ret.success is True
    #
    #    #  <<  Reference
    #
    #    stdout_ref = {
    #        "cfour": "University of Florida",  # freebie
    #        "dftd3": "Grimme",  # freebie
    #        "gamess": "IOWA STATE UNIVERSITY",  # freebie
    #        "gcp": "Grimme",  # freebie
    #        "mp2d": "Beran",  # freebie
    #        "nwchem": "E. Apra",  # freebie
    #        "psi4": rf"Scratch directory: {scratch_directory}/tmp\w+_psi_scratch/",
    #    }
    #
    #    # a scratch file (preferrably output) expected after job if scratch not cleaned up
    #    scratch_sample = {
    #        "cfour": "*/NEWFOCK",
    #        "dftd3": "*/dftd3_geometry.xyz",  # no outfiles
    #        "gamess": "*/gamess.dat",
    #        "gcp": "*/gcp_geometry.xyz",  # no outfiles
    #        "mp2d": "*/mp2d_geometry",  # no outfiles
    #        "nwchem": "*/nwchem.db",
    #        "psi4": "*/psi.*.35",
    #    }
    #
    #    #  <<  Test
    #
    #    assert config.scratch_directory.endswith(program)
    #
    #    if harness._defaults["scratch"] is True:
    #        sample_file = list(Path(scratch_directory).glob(scratch_sample[program]))
    #        assert len(sample_file) == 1, f"Scratch sample not found: {scratch_sample[program]} in {scratch_directory}"
    #
    #        assert re.search(stdout_ref[program], ret.stdout), f"Scratch pattern not found: {stdout_ref[program]}"
    assert 0


# @pytest.mark.parametrize("ncores", [1, 3])
# @pytest.mark.parametrize("program, model, keywords", _canonical_methods)
# def test_local_options_ncores(program, model, keywords, ncores):
#    """Ensure multithreading implemented in harness (if applicable) or multithreaded runs don't
#       break harness (if inapplicable).
#
#    For available harnesses, run minimal calc with single and multiple cores; check ncores count
#      shows up in ``TaskConfig``.
#    For ``thread_parallel``-active harnesses, check ncores count registers in output.
#
#    New Harness Instructions
#    ------------------------
#    * Make sure minimal calc is in _canonical_methods above.
#    * If ``thread_parallel=True`` in harness, add regex to ``stdout_ref`` below to check ncores the
#      program sees.
#    * If this test doesn't work, implement or adjust ``config.ncores`` in your harness.
#
#    """
#    if not has_program(program):
#        pytest.skip(f"Program '{program}' not found.")
#
#    harness = qcng.get_program(program)
#    molecule = _get_molecule(program, model["method"])
#
#    #  <<  Config
#
#    config = qcng.config.get_config(
#        hostname="something",
#        local_options={
#            "ncores": ncores,
#            "nnodes": 1,
#            # "memory": 0.04
#        },
#    )
#
#    #  <<  Run
#
#    inp = AtomicInput(molecule=molecule, driver="energy", model=model, keywords=keywords)
#    ret = qcdb.compute(inp, program, local_options=config.dict())
#    pprint.pprint(ret.dict(), width=200)
#    assert ret.success is True
#
#    #  <<  Reference
#
#    stdout_ref = {
#        "cfour": rf"Running with {ncores} threads/proc",
#        "gamess": rf"MEMDDI DISTRIBUTED OVER\s+{ncores} PROCESSORS",
#        # "gamess": rf"PARALLEL VERSION RUNNING ON\s+{ncores} PROCESSORS IN\s+1 NODES",  # no line for serial
#        # nwchem is node_parallel only
#        "psi4": rf"Threads:\s+{ncores}",
#    }
#
#    #  <<  Test
#
#    assert config.ncores == ncores
#    assert config.nnodes == 1
#
#    if harness._defaults["thread_parallel"] is True:
#        assert re.search(stdout_ref[program], ret.stdout), f"Thread pattern not found: {stdout_ref[program]}"
