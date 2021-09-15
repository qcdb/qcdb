"""
from https://github.com/psi4/psi4/blob/master/tests/tu2-ch2-energy/input.dat
Sample UHF/6-31G** CH2 computation

"""
import pprint
import pytest

import qcdb

from .utils import *

tu2_scf_ene = -38.9250886434
tu2_scf_ene_df = -38.9253346246


@using("psi4")
def test_tu2_uhf_psi4():
    ch2 = qcdb.set_molecule(
        """
        0 3
        C
        H 1 R
        H 1 R 2 A

        R = 2.05
        A = 133.93
        units au
    """
    )

    qcdb.set_keywords(
        {
            "basis": "6-31g**",
            "reference": "uhf",
            "scf_type": "pk",
        }
    )
    print(ch2)
    print(qcdb.get_active_options().print_changed())

    ene, wfn = qcdb.energy("p4-scf", return_wfn=True)
    pprint.pprint(wfn, width=200)  # debug printing

    assert compare_values(tu2_scf_ene, qcdb.variable("hf total energy"), 6, "energy")
    assert compare("Psi4", wfn["provenance"]["creator"], "harness")


@using("cfour")
def test_tu2_uhf_cfour():
    ch2 = qcdb.set_molecule(
        """
        0 3
        C
        H 1 R
        H 1 R 2 A

        R = 2.05
        A = 133.93
        units au
    """
    )

    qcdb.set_keywords(
        {
            "basis": "6-31g**",
            "reference": "uhf",
            "puream": "cart",
        }
    )
    print(ch2)
    print(qcdb.get_active_options().print_changed())

    ene, wfn = qcdb.energy("c4-scf", return_wfn=True)
    pprint.pprint(wfn, width=200)  # debug printing

    assert compare_values(tu2_scf_ene, qcdb.variable("hf total energy"), 6, "energy")
    assert compare("CFOUR", wfn["provenance"]["creator"], "harness")


@using("nwchem")
def test_tu2_uhf_nwchem():
    ch2 = qcdb.set_molecule(
        """
        0 3
        C
        H 1 R
        H 1 R 2 A

        R = 2.05
        A = 133.93
        units au
    """
    )

    qcdb.set_keywords(
        {
            "basis": "6-31g**",
            "reference": "uhf",
        }
    )
    print(ch2)
    print(qcdb.get_active_options().print_changed())

    ene, wfn = qcdb.energy("nwc-scf", return_wfn=True)
    pprint.pprint(wfn, width=200)  # debug printing

    assert compare_values(tu2_scf_ene, qcdb.variable("hf total energy"), 6, "energy")
    assert compare("NWChem", wfn["provenance"]["creator"], "harness")


@using("gamess")
def test_tu2_uhf_gamess():
    ch2 = qcdb.set_molecule(
        """
        0 3
        C
        H 1 R
        H 1 R 2 A

        R = 2.05
        A = 133.93
        units au
    """
    )

    qcdb.set_keywords(
        {
            "basis": "6-31g**",
            "reference": "uhf",
        }
    )
    print(ch2)
    print(qcdb.get_active_options().print_changed())

    ene, wfn = qcdb.energy("gms-scf", return_wfn=True)
    pprint.pprint(wfn, width=200)  # debug printing

    assert compare_values(tu2_scf_ene, qcdb.variable("hf total energy"), 6, "energy")
    assert compare("GAMESS", wfn["provenance"]["creator"], "harness")
