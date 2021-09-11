"""
from https://github.com/psi4/psi4/blob/master/tests/tu4-h2o-freq/input.dat
Optimization followed by frequencies H2O HF/cc-pVDZ

"""
import pprint
import pytest

import qcdb

from .addons import *
from .utils import *

tu3_nre_start = 9.16819
tu3_nre_opt = 9.30079
tu3_scf_ene = -76.0270535127645
tu4_scf_freqs = [1775.33, 4113.46, 4212.16]
tu4_scf_freqs_df = [1775.65, 4113.38, 4212.18]


@using_psi4
def test_tu4_freq_psi4():
    h2o = qcdb.set_molecule(
        """
      O
      H 1 0.96
      H 1 0.96 2 104.5
    """
    )
    h2o.update_geometry()  # only needed for pre-calc NRE check next line
    assert compare_values(tu3_nre_start, h2o.nuclear_repulsion_energy(), 3, "Nuclear repulsion energy")

    qcdb.optking("p4-scf/cc-pvdz")
    ene, wfn = qcdb.frequency("p4-scf/cc-pvdz", return_wfn=True)
    freqs = wfn["frequency_analysis"]["omega"].data
    pprint.pprint(wfn, width=200)  # debug printing

    assert compare_values(tu3_nre_opt, h2o.nuclear_repulsion_energy(), 3, "Nuclear repulsion energy")
    assert compare_values(tu3_scf_ene, qcdb.variable("CURRENT ENERGY"), 4, "opt energy")  # loose check since DF in psi
    assert compare_values(tu4_scf_freqs, freqs[-3:], atol=0.4, label="freq omegas")  # loose check since DF in psi
    assert compare("Psi4", wfn["provenance"]["creator"], "harness")


@using_cfour
def test_tu4_freq_cfour():
    h2o = qcdb.set_molecule(
        """
      O
      H 1 0.96
      H 1 0.96 2 104.5
    """
    )
    h2o.update_geometry()  # only needed for pre-calc NRE check next line
    assert compare_values(tu3_nre_start, h2o.nuclear_repulsion_energy(), 3, "Nuclear repulsion energy")

    qcdb.optking("c4-scf/cc-pvdz")
    ene, wfn = qcdb.frequency("c4-scf/cc-pvdz", return_wfn=True)
    freqs = wfn["frequency_analysis"]["omega"].data
    pprint.pprint(wfn, width=200)  # debug printing

    assert compare_values(tu3_nre_opt, h2o.nuclear_repulsion_energy(), 3, "Nuclear repulsion energy")
    assert compare_values(tu3_scf_ene, qcdb.variable("CURRENT ENERGY"), 6, "opt energy")
    assert compare_values(tu4_scf_freqs, freqs[-3:], atol=0.1, label="freq omegas")
    assert compare("CFOUR", wfn["provenance"]["creator"], "harness")


@using_nwchem
def test_tu4_freq_nwchem():
    h2o = qcdb.set_molecule(
        """
      O
      H 1 0.96
      H 1 0.96 2 104.5
    """
    )
    h2o.update_geometry()  # only needed for pre-calc NRE check next line
    assert compare_values(tu3_nre_start, h2o.nuclear_repulsion_energy(), 3, "Nuclear repulsion energy")

    qcdb.optking("nwc-scf/cc-pvdz")
    ene, wfn = qcdb.frequency("nwc-scf/cc-pvdz", return_wfn=True)
    freqs = wfn["frequency_analysis"]["omega"].data
    pprint.pprint(wfn, width=200)  # debug printing

    assert compare_values(tu3_nre_opt, h2o.nuclear_repulsion_energy(), 3, "Nuclear repulsion energy")
    assert compare_values(tu3_scf_ene, qcdb.variable("CURRENT ENERGY"), 6, "opt energy")
    assert compare_values(tu4_scf_freqs, freqs[-3:], atol=0.1, label="freq omegas")
    assert compare("NWChem", wfn["provenance"]["creator"], "harness")


@using_gamess
def test_tu4_freq_gamess():
    h2o = qcdb.set_molecule(
        """
      O
      H 1 0.96
      H 1 0.96 2 104.5
    """
    )
    h2o.update_geometry()  # only needed for pre-calc NRE check next line
    assert compare_values(tu3_nre_start, h2o.nuclear_repulsion_energy(), 3, "Nuclear repulsion energy")

    qcdb.optking("gms-scf/cc-pvdz")
    ene, wfn = qcdb.frequency("gms-scf/cc-pvdz", return_wfn=True)
    freqs = wfn["frequency_analysis"]["omega"].data
    pprint.pprint(wfn, width=200)  # debug printing

    assert compare_values(tu3_nre_opt, h2o.nuclear_repulsion_energy(), 3, "Nuclear repulsion energy")
    assert compare_values(tu3_scf_ene, qcdb.variable("CURRENT ENERGY"), 6, "opt energy")
    assert compare_values(tu4_scf_freqs, freqs[-3:], atol=0.1, label="freq omegas")
    assert compare("GAMESS", wfn["provenance"]["creator"], "harness")
