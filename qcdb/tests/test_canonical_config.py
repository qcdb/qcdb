"""
Tests the DQM compute dispatch module
"""
import pprint
import re
import tempfile
from pathlib import Path

import pytest
import qcengine as qcng
from qcelemental.models import AtomicInput
from qcengine.testing import has_program, using

import qcdb

_canonical_methods = [
    ("cfour", {"method": "hf", "basis": "6-31G"}, {}),
    #    ("dftd3", {"method": "b3lyp-d3"}, {}),
    ("gamess", {"method": "hf", "basis": "6-31g"}, {}),
    ("gamess", {"method": "mp2", "basis": "aug-cc-pvdz"}, {}),
    ("gamess", {"method": "ccsd", "basis": "aug-cc-pvdz"}, {}),
    ("gamess", {"method": "ccsd(t)", "basis": "aug-cc-pvtz"}, {}),
    #    ("gamess", {"method": "hf", "basis": "n31"}, {"basis__NGAUSS": 6}),
    #    ("gamess", {"method": "mp2", "basis": "accd"}, {"contrl__ispher": 1}),
    #    ("gamess", {"method": "ccsd", "basis": "accd"}, {"contrl__ispher": 1}),
    #    ("gamess", {"method": "ccsd(t)", "basis": "acct"}, {"contrl__ispher": 1}),
    ("nwchem", {"method": "hf", "basis": "6-31G"}, {}),
    ("psi4", {"method": "hf", "basis": "6-31G"}, {}),
]


def _get_molecule(program, method):
    if program in ["openmm", "terachem_pbs"]:
        return qcng.get_molecule("water")
    elif program == "gamess" and method == "ccsd(t)":
        return qcng.get_molecule("water")
    else:
        return qcng.get_molecule("hydrogen")


@pytest.mark.parametrize(
    "memory_trickery",
    [
        pytest.param({}, id="none"),
        pytest.param({"memory": "1.555 gib"}, id="qcdb"),
        pytest.param(
            {
                "cfour": {"memory": "5 gb"},
                "gamess": {"memory": "5 gb"},
                "nwchem": {"memory": "5 gb"},
                "psi4": None,
            },
            id="qcdb-contra",
        ),
        pytest.param(
            # native keywords consistent with config.memory below
            {
                "cfour": {"cfour_memory_size": 208708567},
                "gamess": {"gamess_system__mwords": 208},
                "nwchem": {"nwchem_memory": 1669668536},
                "psi4": None,  # no memory keyword in psi
            },
            id="dsl",
        ),
        pytest.param(
            # native keywords that CONTRADICT config.memory below
            {
                "cfour": {"cfour_memory_size": "5000"},
                "gamess": {"gamess_system__mwords": 500},
                "nwchem": {"nwchem_memory": 500000000},
                "psi4": None,  # no contradictory memory keyword in psi
            },
            id="dsl-contra",
        ),
    ],
)
@pytest.mark.parametrize("program, model, keywords", _canonical_methods)
def test_local_options_memory_gib(program, model, keywords, memory_trickery, request):
    """Ensure memory handling implemented in harness (if applicable).

    For available harnesses, run minimal calc at specific total node memory, both through runtime
      config alone and with clashing (and non-QCEngine-like) keyword spec. Check memory quantity
      shows up in ``TaskConfig``.
    For ``managed-memory``-active harnesses, check that memory registers in output.

    New Harness Instructions
    ------------------------
    * Make sure minimal calc is in _canonical_methods above.
    * If ``managed_memory=True`` in harness, add regex to ``stdout_ref`` below to check that memory
      is specifiable.
    * If this test doesn't work, implement or adjust ``config.memory`` in your harness.

    """
    if not has_program(program):
        pytest.skip(f"Program '{program}' not found.")

    harness = qcng.get_program(program)
    molecule = _get_molecule(program, model["method"])

    addl_keywords = memory_trickery.get(program, memory_trickery)
    if addl_keywords is None:
        pytest.skip(f"Nothing to test for '{program}' id '{request.node.name}'")
    use_keywords = {**keywords, **addl_keywords}

    #  <<  Config

    config = qcng.config.get_config(
        hostname="something",
        local_options={
            "ncores": 1,
            "nnodes": 1,
            "memory": 1.555,
        },
    )

    #  <<  Run

    inp = AtomicInput(molecule=molecule, driver="energy", model=model, keywords=use_keywords)
    ret = qcdb.compute(inp, program, local_options=config.dict())
    pprint.pprint(ret.dict(), width=200)

    if "contra" in request.node.name:
        assert ret.success is False
        assert "KeywordReconciliationError" in ret.error.error_message
        return

    assert ret.success is True

    #  <<  Reference

    stdout_ref = {  # 1.555 GiB = 208708567 quad-words
        "cfour": "Allocated    1592 MB of main memory",
        "gamess": "208000000 WORDS OF MEMORY AVAILABLE",
        "nwchem": r"total    =  2087085\d\d doubles =   1592.3 Mbytes",  # doubles is quad-words. Mbytes is MiB
        "psi4": "1592 MiB Core",
    }

    #  <<  Test

    assert config.ncores == 1
    assert pytest.approx(config.memory, 0.1) == 1.555

    if harness._defaults["managed_memory"] is True:
        assert re.search(stdout_ref[program], ret.stdout), f"Memory pattern not found: {stdout_ref[program]}"
    # assert 0


@pytest.mark.parametrize("program, model, keywords", _canonical_methods)
def test_local_options_scratch(program, model, keywords):
    """Ensure scratch handling implemented in harness (if applicable).

    For available harnesses, run minimal calc at specific scratch directory name (randomly generated
      during test) and skip scratch clean-up. Check scratch settings show up in ``TaskConfig``.
    For ``scratch``-active harnesses, check that an expected file is written to and left behind in
      scratch directory. Check any scratch-related printing in output.

    New Harness Instructions
    ------------------------
    * Make sure minimal calc is in _canonical_methods above.
    * If ``scratch=True`` in harness, add single file (preferrably output) glob to ``scratch_sample``
      below to check that program scratch is directable.
    * If ``scratch=True`` in harness, if scratch directory mentioned in output, add regex to
      ``stdout_ref`` below to check that program scratch is directable. Otherwise, add an
      always-passing regex.
    * If this test doesn't work, implement or adjust ``config.scratch_directory`` and
      ``config.scratch_messy`` in your harness.

    """
    if not has_program(program):
        pytest.skip(f"Program '{program}' not found.")

    harness = qcng.get_program(program)
    molecule = _get_molecule(program, model["method"])

    #  <<  Config

    scratch_directory = tempfile.mkdtemp(suffix="_" + program)

    config = qcng.config.get_config(
        hostname="something",
        local_options={
            "scratch_directory": scratch_directory,
            "scratch_messy": True,
        },
    )

    #  <<  Run

    inp = AtomicInput(molecule=molecule, driver="energy", model=model, keywords=keywords)
    ret = qcdb.compute(inp, program, local_options=config.dict())
    pprint.pprint(ret.dict(), width=200)
    assert ret.success is True

    #  <<  Reference

    stdout_ref = {
        "cfour": "University of Florida",  # freebie
        "dftd3": "Grimme",  # freebie
        "gamess": "IOWA STATE UNIVERSITY",  # freebie
        "gcp": "Grimme",  # freebie
        "mp2d": "Beran",  # freebie
        "nwchem": "E. Apra",  # freebie
        "psi4": rf"Scratch directory: {scratch_directory}/tmp\w+_psi_scratch/",
    }

    # a scratch file (preferrably output) expected after job if scratch not cleaned up
    scratch_sample = {
        "cfour": "*/NEWFOCK",
        "dftd3": "*/dftd3_geometry.xyz",  # no outfiles
        "gamess": "*/gamess.dat",
        "gcp": "*/gcp_geometry.xyz",  # no outfiles
        "mp2d": "*/mp2d_geometry",  # no outfiles
        "nwchem": "*/nwchem.db",
        "psi4": "*/psi.*.35",
    }

    #  <<  Test

    assert config.scratch_directory.endswith(program)

    if harness._defaults["scratch"] is True:
        sample_file = list(Path(scratch_directory).glob(scratch_sample[program]))
        assert len(sample_file) == 1, f"Scratch sample not found: {scratch_sample[program]} in {scratch_directory}"

        assert re.search(stdout_ref[program], ret.stdout), f"Scratch pattern not found: {stdout_ref[program]}"


@pytest.mark.parametrize("ncores", [1, 3])
@pytest.mark.parametrize("program, model, keywords", _canonical_methods)
def test_local_options_ncores(program, model, keywords, ncores):
    """Ensure multithreading implemented in harness (if applicable) or multithreaded runs don't
       break harness (if inapplicable).

    For available harnesses, run minimal calc with single and multiple cores; check ncores count
      shows up in ``TaskConfig``.
    For ``thread_parallel``-active harnesses, check ncores count registers in output.

    New Harness Instructions
    ------------------------
    * Make sure minimal calc is in _canonical_methods above.
    * If ``thread_parallel=True`` in harness, add regex to ``stdout_ref`` below to check ncores the
      program sees.
    * If this test doesn't work, implement or adjust ``config.ncores`` in your harness.

    """
    if not has_program(program):
        pytest.skip(f"Program '{program}' not found.")

    harness = qcng.get_program(program)
    molecule = _get_molecule(program, model["method"])

    #  <<  Config

    config = qcng.config.get_config(
        hostname="something",
        local_options={
            "ncores": ncores,
            "nnodes": 1,
            # "memory": 0.04
        },
    )

    #  <<  Run

    inp = AtomicInput(molecule=molecule, driver="energy", model=model, keywords=keywords)
    ret = qcdb.compute(inp, program, local_options=config.dict())
    pprint.pprint(ret.dict(), width=200)
    assert ret.success is True

    #  <<  Reference

    stdout_ref = {
        "cfour": rf"Running with {ncores} threads/proc",
        "gamess": rf"MEMDDI DISTRIBUTED OVER\s+{ncores} PROCESSORS",
        # "gamess": rf"PARALLEL VERSION RUNNING ON\s+{ncores} PROCESSORS IN\s+1 NODES",  # no line for serial
        # nwchem is node_parallel only
        "psi4": rf"Threads:\s+{ncores}",
    }

    #  <<  Test

    assert config.ncores == ncores
    assert config.nnodes == 1

    if harness._defaults["thread_parallel"] is True:
        assert re.search(stdout_ref[program], ret.stdout), f"Thread pattern not found: {stdout_ref[program]}"
