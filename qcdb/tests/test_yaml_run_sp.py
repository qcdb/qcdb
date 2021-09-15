"""
from https://github.com/psi4/psi4/blob/master/tests/tu2-ch2-energy/input.dat

"""
import pprint
import pytest

import qcdb

from .utils import *



@using("psi4")
def test_uhf_yaml_psi4():
    yamlin = """
molecule: |
  0 3
  C
  H 1 R
  H 1 R 2 A

  R = 1.075
  A = 133.93

driver: !!python/name:qcdb.energy
method: p4-hf
options:
  basis: '6-31g**'
  reference: uhf
  puream: cart
  scf_type: pk
"""

    import yaml

    asdf = yaml.load(yamlin, Loader=yaml.FullLoader)

    ene = asdf["driver"](asdf["method"], options=asdf["options"], molecule=qcdb.Molecule(asdf["molecule"]))

    assert compare_values(-38.9253416082900827, ene, 6, "calc from yaml str")


@using("cfour")
def test_uhf_yaml_cfour():

    yamlin = """
molecule: |
  0 3
  C
  H 1 R
  H 1 R 2 A

  R = 1.075
  A = 133.93

driver: !!python/name:qcdb.energy
method: c4-hf
options:
  basis: '6-31g**'
  reference: uhf
  puream: cart
  scf_type: pk
"""

    import yaml
    asdf = yaml.load(yamlin, Loader=yaml.FullLoader)

    ene = asdf['driver'](asdf['method'],
                         options=asdf['options'],
                         molecule=qcdb.Molecule(asdf['molecule']))

    assert compare_values(-38.9253416082900827, ene, 6, 'calc from yaml str')


@using("nwchem")
def test_uhf_yaml_nwchem():

    yamlin = """
molecule: |
  0 3
  C
  H 1 R
  H 1 R 2 A

  R = 1.075
  A = 133.93

driver: !!python/name:qcdb.energy
method: nwc-hf
options:
  basis: '6-31g**'
  reference: uhf
  puream: cart
  scf_type: pk
"""

    import yaml
    asdf = yaml.load(yamlin, Loader=yaml.FullLoader)

    ene = asdf['driver'](asdf['method'],
                         options=asdf['options'],
                         molecule=qcdb.Molecule(asdf['molecule']))

    assert compare_values(-38.9253416082900827, ene, 6, 'calc from yaml str')


@using("gamess")
def test_uhf_yaml_gamess():

    yamlin = """
molecule: |
  0 3
  C
  H 1 R
  H 1 R 2 A

  R = 1.075
  A = 133.93

driver: !!python/name:qcdb.energy
method: gms-hf
options:
  basis: '6-31g**'
  reference: uhf
  puream: cart
  scf_type: pk
"""

    import yaml
    asdf = yaml.load(yamlin, Loader=yaml.FullLoader)

    ene = asdf['driver'](asdf['method'],
                         options=asdf['options'],
                         molecule=qcdb.Molecule(asdf['molecule']))

    assert compare_values(-38.9253416082900827, ene, 6, 'calc from yaml str')

