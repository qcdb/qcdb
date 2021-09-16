import os
import sys

import qcdb

from ..utils import *

import pytest
import numpy as np

# import pprinter as pp


def check_dipole(return_value):

    dip_x = -0.0000000000
    dip_y = -0.0000000000
    dip_z = -2.2247923027
    dip_grad = np.array([-0.0000000, -0.000000, -2.2247923027])

    assert compare_values(dip_x, qcdb.variable("CURRENT DIPOLE X"), 5, "dip x")
    assert compare_values(dip_y, qcdb.variable("CURRENT DIPOLE Y"), 5, "dip y")
    assert compare_values(dip_z, qcdb.variable("CURRENT DIPOLE Z"), 5, "dip z")
    assert compare_values(dip_grad, qcdb.variable("CURRENT DIPOLE GRADIENT"), 5, "dip array")


@using("nwchem")
@pytest.mark.xfail(True, reason="properties() NYI", run=True)
def test_1_dipole():
    qcdb.set_molecule(
        """ 
     O 0       0        0
     H 0       1.430   -1.107
     H 0      -1.430   -1.107
     units au
     """
    )

    qcdb.set_options(
        {"basis": "6-31g*", "nwchem_scf__rohf": True, "nwchem_scf__singlet": True, "nwchem_property__dipole": True}
    )

    val = qcdb.properties("nwc-scf")
    check_dipole(val)
    print(qcdb.get_active_options().print_changed())
