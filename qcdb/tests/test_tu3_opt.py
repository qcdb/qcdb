"""
from https://github.com/psi4/psi4/blob/master/tests/tu3-h2o-opt/input.dat
Optimize H2O HF/cc-pVDZ

"""
import pytest

import qcdb

from .utils import *

tu3_nre_start = 9.16819
tu3_nre_opt = 9.30079
tu3_scf_ene = -76.0270535127645


@using("psi4")
def test_tu3_opt_psi4():
    water1 = qcdb.set_molecule(
        """
      O
      H 1 0.96
      H 1 0.96 2 104.5
    """
    )
    water1.update_geometry()  # only needed for pre-calc NRE check next line
    assert compare_values(tu3_nre_start, water1.nuclear_repulsion_energy(), 3, "Nuclear repulsion energy")

    qcdb.optimize("p4-scf/cc-pvdz")

    assert compare_values(tu3_nre_opt, water1.nuclear_repulsion_energy(), 3, "Nuclear repulsion energy")
    assert compare_values(tu3_scf_ene, qcdb.variable("CURRENT ENERGY"), 4, "opt energy")  # loose check since DF in psi


@using("cfour")
def test_tu3_opt_cfour():
    water2 = qcdb.set_molecule(
        """
      O
      H 1 0.96
      H 1 0.96 2 104.5
    """
    )
    water2.update_geometry()  # only needed for pre-calc NRE check next line
    assert compare_values(tu3_nre_start, water2.nuclear_repulsion_energy(), 3, "Nuclear repulsion energy")

    qcdb.optimize("c4-scf/cc-pvdz")

    assert compare_values(tu3_nre_opt, water2.nuclear_repulsion_energy(), 3, "Nuclear repulsion energy")
    assert compare_values(tu3_scf_ene, qcdb.variable("CURRENT ENERGY"), 6, "opt energy")


@using("nwchem")
def test_tu3_opt_nwchem():
    water3 = qcdb.set_molecule(
        """
      O
      H 1 0.96
      H 1 0.96 2 104.5
    """
    )
    water3.update_geometry()  # only needed for pre-calc NRE check next line
    assert compare_values(tu3_nre_start, water3.nuclear_repulsion_energy(), 3, "Nuclear repulsion energy")

    qcdb.optimize("nwc-scf/cc-pvdz")

    assert compare_values(tu3_nre_opt, water3.nuclear_repulsion_energy(), 3, "Nuclear repulsion energy")
    assert compare_values(tu3_scf_ene, qcdb.variable("CURRENT ENERGY"), 6, "opt energy")


@using("gamess")
def test_tu3_opt_gamess():
    water4 = qcdb.set_molecule(
        """
      O
      H 1 0.96
      H 1 0.96 2 104.5
    """
    )
    water4.update_geometry()  # only needed for pre-calc NRE check next line
    assert compare_values(tu3_nre_start, water4.nuclear_repulsion_energy(), 3, "Nuclear repulsion energy")

    qcdb.optimize("gms-scf/cc-pvdz")

    assert compare_values(tu3_nre_opt, water4.nuclear_repulsion_energy(), 3, "Nuclear repulsion energy")
    assert compare_values(tu3_scf_ene, qcdb.variable("CURRENT ENERGY"), 6, "opt energy")
