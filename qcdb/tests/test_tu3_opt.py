"""
from https://github.com/psi4/psi4/blob/master/tests/tu3-h2o-opt/input.dat
Optimize H2O HF/cc-pVDZ

"""
import pytest

import qcdb

from .addons import *
from .utils import *

tu3_nre_start = 9.16819
tu3_nre_opt = 9.30079
tu3_scf_ene = -76.0270535127645


@using_psi4
def test_tu3_opt_psi4():
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

    assert compare_values(tu3_nre_opt, h2o.nuclear_repulsion_energy(), 3, "Nuclear repulsion energy")
    assert compare_values(tu3_scf_ene, qcdb.variable("CURRENT ENERGY"), 4, "opt energy")  # loose check since DF in psi


@using_cfour
def test_tu3_opt_cfour():
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

    assert compare_values(tu3_nre_opt, h2o.nuclear_repulsion_energy(), 3, "Nuclear repulsion energy")
    assert compare_values(tu3_scf_ene, qcdb.variable("CURRENT ENERGY"), 6, "opt energy")


@using_nwchem
def test_tu3_opt_nwchem():
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

    assert compare_values(tu3_nre_opt, h2o.nuclear_repulsion_energy(), 3, "Nuclear repulsion energy")
    assert compare_values(tu3_scf_ene, qcdb.variable("CURRENT ENERGY"), 6, "opt energy")


@using_gamess
def test_tu3_opt_gamess():
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

    assert compare_values(tu3_nre_opt, h2o.nuclear_repulsion_energy(), 3, "Nuclear repulsion energy")
    assert compare_values(tu3_scf_ene, qcdb.variable("CURRENT ENERGY"), 6, "opt energy")
