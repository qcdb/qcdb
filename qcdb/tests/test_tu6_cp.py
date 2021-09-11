"""
from https://github.com/psi4/psi4/blob/master/tests/tu6-cp-ne2/input.dat
Example potential energy surface scan and CP-correction for Ne2

"""
import pprint
import pytest

import qcengine
from qcengine.testing import using
import qcdb
from qcelemental import constants

from .addons import *
from .utils import *

tu6_ie_scan = {2.5: 0.757717, 3.0: 0.015685, 4.0: -0.016266}




@using("psi4")
def test_tu6_cp_psi4():
    dimer = qcdb.set_molecule(
        """
      Ne
    --
      Ne 1 R
    """
    )

    qcdb.set_options(
        {
            "basis": "aug-cc-pVDZ",
            "freeze_core": True,
        }
    )

    Rvals = [2.5, 3.0, 4.0]
    ecp = {}
    for R in Rvals:
        dimer.R = R
        ecp[R], wfn = qcdb.energy("p4-ccsd(t)", bsse_type="cp", return_wfn=True)
        pprint.pprint(wfn, width=200)
        assert compare_values(
            tu6_ie_scan[R], ecp[R] * constants.hartree2kcalmol, atol=1.0e-4, label=f"CP-CCSD(T) [{R:3.1f}]"
        )

    # assert compare("Psi4", wfn["provenance"]["creator"], "harness")

    print("\nCP-corrected CCSD(T)/aug-cc-pVDZ interaction energies\n\n")
    print("        R [Ang]         E_int [kcal/mol]             \n")
    print("-----------------------------------------------------\n")
    for R in Rvals:
        e = ecp[R] * constants.hartree2kcalmol
        print(f"        {R:3.1f}            {e:10.6f}\n")




@using("cfour")
def test_tu6_cp_cfour():
    dimer = qcdb.set_molecule(
        """
      Ne
    --
      Ne 1 R
    """
    )

    qcdb.set_options(
        {
            "basis": "aug-cc-pVDZ",
            "freeze_core": True,
        }
    )

    Rvals = [2.5, 3.0, 4.0]
    ecp = {}
    for R in Rvals:
        dimer.R = R
        ecp[R], wfn = qcdb.energy("c4-ccsd(t)", bsse_type="cp", return_wfn=True)
        pprint.pprint(wfn, width=200)
        assert compare_values(
            tu6_ie_scan[R], ecp[R] * constants.hartree2kcalmol, atol=1.0e-4, label=f"CP-CCSD(T) [{R:3.1f}]"
        )

    # assert compare("CFOUR", wfn["provenance"]["creator"], "harness")

    print("\nCP-corrected CCSD(T)/aug-cc-pVDZ interaction energies\n\n")
    print("        R [Ang]         E_int [kcal/mol]             \n")
    print("-----------------------------------------------------\n")
    for R in Rvals:
        e = ecp[R] * constants.hartree2kcalmol
        print(f"        {R:3.1f}            {e:10.6f}\n")




@using("nwchem")
def test_tu6_cp_nwchem():
    dimer = qcdb.set_molecule(
        """
      Ne
    --
      Ne 1 R
    """
    )

    qcdb.set_options(
        {
            "basis": "aug-cc-pVDZ",
            "freeze_core": True,
        }
    )

    Rvals = [2.5, 3.0, 4.0]
    ecp = {}
    for R in Rvals:
        dimer.R = R
        ecp[R], wfn = qcdb.energy("nwc-ccsd(t)", bsse_type="cp", return_wfn=True)
        pprint.pprint(wfn, width=200)
        assert compare_values(
            tu6_ie_scan[R], ecp[R] * constants.hartree2kcalmol, atol=1.0e-4, label=f"CP-CCSD(T) [{R:3.1f}]"
        )

    # assert compare("NWChem", wfn["provenance"]["creator"], "harness")

    print("\nCP-corrected CCSD(T)/aug-cc-pVDZ interaction energies\n\n")
    print("        R [Ang]         E_int [kcal/mol]             \n")
    print("-----------------------------------------------------\n")
    for R in Rvals:
        e = ecp[R] * constants.hartree2kcalmol
        print(f"        {R:3.1f}            {e:10.6f}\n")


@using("gamess")
def test_tu6_cp_gamess():
    dimer = qcdb.set_molecule(
        """
      Ne
    --
      Ne 1 R
    """
    )

    qcdb.set_options(
        {
            "basis": "aug-cc-pVDZ",
            "freeze_core": True,
        }
    )

    Rvals = [2.5, 3.0, 4.0]
    ecp = {}
    for R in Rvals:
        dimer.R = R
        with pytest.raises(qcengine.exceptions.InputError) as e:
            ecp[R], wfn = qcdb.energy("gms-ccsd(t)", bsse_type="cp", return_wfn=True)
        return 0
        pprint.pprint(wfn, width=200)
        assert compare_values(
            tu6_ie_scan[R], ecp[R] * constants.hartree2kcalmol, atol=1.0e-4, label=f"CP-CCSD(T) [{R:3.1f}]"
        )

    # assert compare("GAMESS", wfn["provenance"]["creator"], "harness")

    print("\nCP-corrected CCSD(T)/aug-cc-pVDZ interaction energies\n\n")
    print("        R [Ang]         E_int [kcal/mol]             \n")
    print("-----------------------------------------------------\n")
    for R in Rvals:
        e = ecp[R] * constants.hartree2kcalmol
        print(f"        {R:3.1f}            {e:10.6f}\n")
