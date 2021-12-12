import pprint
import re

import pytest
import qcengine
from qcelemental import constants
from qcengine.programs.tests.test_dftd3_mp2d import eneyne_ne_qcdbmols
from qcengine.programs.tests.test_ghost import bimol_ref
from qcengine.testing import using

import qcdb

from .utils import *


@pytest.mark.parametrize(
    "qcp",
    [
        pytest.param("p4-", marks=using("psi4")),
        pytest.param("c4-", marks=using("cfour")),
        pytest.param("nwc-", marks=using("nwchem")),
        pytest.param("gms-", marks=using("gamess")),
    ],
)
def test_simple_ghost(qcp):

    dimer = qcdb.set_molecule(
        f"""
      He
    --
      Ne 1 R
    """
    )
    dimer.R = 2.5
    monomer = dimer.extract_fragments(1, 2)

    qcdb.set_options(
        {
            "basis": "aug-cc-pVDZ",
            "freeze_core": True,
            "scf_type": "conv",  # longstanding "pk"
            "d_convergence": 8,
        }
    )

    ene, wfn = qcdb.energy(qcp + "hf", return_wfn=True, molecule=monomer)
    pprint.pprint(wfn, width=200)

    atol = 1.0e-6
    assert compare_values(0.0, qcdb.variable("nuclear repulsion energy"), atol=atol, label="nre")
    assert compare_values(0.0, wfn["properties"]["nuclear_repulsion_energy"], atol=atol, label="nre")
    assert compare(32, wfn["properties"]["calcinfo_nbasis"], label="nbas")
    assert compare(32, wfn["properties"]["calcinfo_nmo"], label="nmo")
    assert compare_values(-2.8557143339397539, ene, atol=atol, label="ene")


@pytest.mark.parametrize("subject", ["dimer", "mA", "mB", "mAgB", "gAmB"])
@pytest.mark.parametrize(
    "qcprog, keywords",
    [
        pytest.param("cfour", {}, id="cfour", marks=using("cfour")),
        pytest.param("gamess", {"gamess_mp2__nacore": 0}, id="gamess", marks=using("gamess")),
        pytest.param("nwchem", {}, id="nwchem", marks=using("nwchem")),
        pytest.param(
            "psi4",
            {
                "psi4_scf_type": "pk",
                "psi4_mp2_type": "conv",
            },
            id="psi4",
            marks=using("psi4"),
        ),
    ],
)
def test_tricky_ghost(qcprog, subject, keywords):
    qmol = eneyne_ne_qcdbmols()["eneyne"][subject]
    ref = bimol_ref["eneyne"]

    assert qmol.natom() == ref["natom"][subject]
    assert sum([int(qmol.Z(at) > 0) for at in range(qmol.natom())]) == ref["nreal"][subject]

    ene, wfn = qcdb.energy(
        qcdb.util.program_prefix(qcprog) + "mp2/6-31g*", return_wfn=True, molecule=qmol, options=keywords
    )
    pprint.pprint(wfn, width=200)

    assert compare_values(
        ref["nre"][subject], wfn["properties"]["nuclear_repulsion_energy"], atol=1.0e-4, label="nre"
    ), f'nre: {wfn["properties"]["nuclear_repulsion_energy"]} != {ref["nre"][subject]}'
    assert compare(
        ref["nbasis"][subject], wfn["properties"]["calcinfo_nbasis"], label="nbasis"
    ), f'nbasis: {wfn["properties"]["calcinfo_nbasis"]} != {ref["nbasis"][subject]}'
    assert compare(
        ref["nmo"][subject], wfn["properties"]["calcinfo_nmo"], label="nmo"
    ), f'nmo: {wfn["properties"]["calcinfo_nmo"]} != {ref["nmo"][subject]}'
    assert compare_values(ref["mp2"][subject], ene, atol=1.0e-6, label="ene"), f'ene: {ene} != {ref["mp2"][subject]}'

    pgline = {
        "cfour": r"Computational point group: (?P<pg>\w+)",
        "gamess": r"THE POINT GROUP IS (?P<pg>[\w\s,=]+)",
        "nwchem": r"Group name\s+(?P<pg>\w+)",
        "psi4": r"Running in (?P<pg>\w+) symmetry.",
    }
    mobj = re.search(pgline[qcprog], wfn["stdout"])
    if mobj:
        pg = mobj.group("pg").strip()
        if pg == "CNV, NAXIS= 2, ORDER= 4":
            pg = "C2v"
        elif pg == "C1 , NAXIS= 0, ORDER= 1":
            pg = "C1"
        pg = pg.capitalize()

    if qcprog == "gamess" and subject in ["mAgB", "gAmB"]:
        # don't know how to get master frame w/ghosts in gamess, so C1 forced
        assert pg == "C1", f"pg: {pg} != C1"
    else:
        assert pg in ref["pg"][subject], f'pg: {pg} != {ref["pg"][subject]}'
