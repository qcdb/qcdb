import pprint

import pytest
import qcengine
from qcelemental import constants
from qcengine.testing import using

import qcdb

from .utils import *

tu6_ie_scan = {2.5: 0.757717, 3.0: 0.015685, 4.0: -0.016266}


@pytest.mark.parametrize(
    "qcp",
    [
        pytest.param("p4-", marks=using("psi4")),
        pytest.param("c4-", marks=using("cfour")),
        pytest.param("nwc-", marks=using("nwchem")),
        pytest.param("gms-", marks=using("gamess")),
    ],
)
def test_ghost(qcp):

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
            "scf_type": "pk",
            "d_convergence": 8,
        }
    )

    if qcp == "gms-":
        with pytest.raises(qcengine.exceptions.InputError) as e:
            ene, wfn = qcdb.energy(qcp + "hf", return_wfn=True, molecule=monomer)
        return None

    ene, wfn = qcdb.energy(qcp + "hf", return_wfn=True, molecule=monomer)
    pprint.pprint(wfn, width=200)

    assert compare_values(0.0, qcdb.variable("nuclear repulsion energy"), atol=1.0e-6, label="nre")
    assert compare_values(0.0, wfn["properties"]["nuclear_repulsion_energy"], atol=1.0e-6, label="nre")
    assert compare(32, wfn["properties"]["calcinfo_nbasis"], label="nbas")
    assert compare(32, wfn["properties"]["calcinfo_nmo"], label="nmo")
    assert compare_values(-2.8557143339397539, ene, atol=1.0e-6, label="ene")
