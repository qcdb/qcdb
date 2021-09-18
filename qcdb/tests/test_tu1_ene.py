"""
from https://github.com/psi4/psi4/blob/master/tests/tu1-h2o-energy/input.dat
Sample HF/cc-pVDZ H2O computation

"""
import pprint

import pytest

import qcdb

from .utils import *

tu1_scf_ene = -76.02665366
tu1_scf_ene_df = -76.0266327341


@using("psi4")
def test_tu1_ene_psi4():

    # memory 600 mb
    #
    # molecule h2o {
    #  O
    #  H 1 0.96
    #  H 1 0.96 2 104.5
    # }
    #
    # set basis cc-pVDZ
    # energy('scf')

    # compare_values(-76.0266327341067125, variable('SCF TOTAL ENERGY'), 6, 'SCF energy')  #TEST

    h2o = qcdb.set_molecule(
        """
  O 
  H 1 0.96
  H 1 0.96 2 104.5
"""
    )

    qcdb.set_keywords({"scf_type": "pk"})

    ene = qcdb.energy("p4-hf/cc-pVDZ")

    assert compare_values(tu1_scf_ene, ene, 6, "energy")


@using("cfour")
def test_tu1_ene_cfour():
    h2o = qcdb.set_molecule(
        """
  O 
  H 1 0.96
  H 1 0.96 2 104.5
"""
    )

    ene = qcdb.energy("c4-HF/cc-pVDZ")

    assert compare_values(tu1_scf_ene, ene, 6, "energy")


@using("nwchem")
def test_tu1_ene_nwchem():
    h2o = qcdb.set_molecule(
        """
  O 
  H 1 0.96
  H 1 0.96 2 104.5
"""
    )

    ene = qcdb.energy("nwc-HF/cc-pVDZ")

    assert compare_values(tu1_scf_ene, ene, 6, "energy")


@using("gamess")
def test_tu1_ene_gamess():
    h2o = qcdb.set_molecule(
        """
  O 
  H 1 0.96
  H 1 0.96 2 104.5
"""
    )

    ene = qcdb.energy("gms-HF/cc-pVDZ")

    assert compare_values(tu1_scf_ene, ene, 6, "energy")
