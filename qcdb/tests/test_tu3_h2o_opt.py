import pytest

import qcdb

from .addons import *
from .utils import *

nucenergy =   9.3007948234239
refenergy = -76.0270535127645

@pytest.mark.xfail(True, reason='Old Driver, Spring 2019', run=True)
@using_psi4
def test_1a():
    h2o = qcdb.set_molecule("""
      O
      H 1 0.96
      H 1 0.96 2 104.5
    """)

    qcdb.optking('scf/cc-pvdz')

    assert compare_values(nucenergy, h2o.nuclear_repulsion_energy(), 3, "Nuclear repulsion energy")
    assert compare_values(refenergy, qcdb.variable("CURRENT ENERGY"), 4, "Reference energy")

@pytest.mark.xfail(True, reason='Old Driver, Spring 2019', run=True)
@using_psi4
def test_1b():
    h2o = qcdb.set_molecule("""
      O
      H 1 0.96
      H 1 0.96 2 104.5
    """)

    qcdb.optking('c4-scf/cc-pvdz')

    assert compare_values(nucenergy, h2o.nuclear_repulsion_energy(), 3, "Nuclear repulsion energy")
    assert compare_values(refenergy, qcdb.variable("CURRENT ENERGY"), 4, "Reference energy")

@using_geometric
def test_2a():
    h2o = qcdb.set_molecule("""
      O
      H 1 0.96
      H 1 0.96 2 104.5
    """)

    qcdb.geometric('scf/cc-pvdz')

    assert compare_values(nucenergy, h2o.nuclear_repulsion_energy(), 3, "Nuclear repulsion energy")
    assert compare_values(refenergy, qcdb.variable("CURRENT ENERGY"), 4, "Reference energy")

@using_geometric
def test_2b():
    h2o = qcdb.set_molecule("""
      O
      H 1 0.96
      H 1 0.96 2 104.5
    """)

    qcdb.geometric('c4-scf/cc-pvdz')

    assert compare_values(nucenergy, h2o.nuclear_repulsion_energy(), 3, "Nuclear repulsion energy")
    assert compare_values(refenergy, qcdb.variable("CURRENT ENERGY"), 4, "Reference energy")


@using_geometric
def test_2c():
    h2o = qcdb.set_molecule("""
O 0.          0.         -0.12
H 0.         -1.52  1.04
H 0.          1.52  1.04
units au
    """)

    qcdb.geometric('scf/cc-pvdz')

    assert compare_values(nucenergy, h2o.nuclear_repulsion_energy(), 3, "Nuclear repulsion energy")
    assert compare_values(refenergy, qcdb.variable("CURRENT ENERGY"), 4, "Reference energy")

@using_geometric
def test_2d():
    h2o = qcdb.set_molecule("""
O 0.          0.         -0.12
H 0.         -1.52  1.04
H 0.          1.52  1.04
units au
    """)

    qcdb.geometric('c4-scf/cc-pvdz')

    assert compare_values(nucenergy, h2o.nuclear_repulsion_energy(), 3, "Nuclear repulsion energy")
    assert compare_values(refenergy, qcdb.variable("CURRENT ENERGY"), 4, "Reference energy")
