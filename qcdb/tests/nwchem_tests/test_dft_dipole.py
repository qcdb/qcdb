import os
import sys

import pytest

import qcdb

from ..utils import *


def check_dft(return_value):
    dft = -76.420359078925
    x = 0.0000000000
    y = 0.0000000000
    z = -1.9104758249

    assert compare_values(dft, qcdb.variable("DFT TOTAL ENERGY"), 5, "dft ref")
    assert compare_values(x, qcdb.variable("CURRENT DIPOLE X"), 5, "dipole x")
    assert compare_values(y, qcdb.variable("CURRENT DIPOLE Y"), 5, "dipole y")
    assert compare_values(z, qcdb.variable("CURRENT DIPOLE Z"), 5, "dipole z")


@pytest.mark.xfail(True, reason="Property not implement | harvester error", run=True)
def test_1_dft():
    h2o = qcdb.set_molecule(
        """
                O      0.00000000     0.00000000     0.11726921
                H      0.75698224     0.00000000    -0.46907685
                H     -0.75698224     0.00000000    -0.46907685
                """
    )

    qcdb.set_options(
        {
            "basis": "cc-pvdz",
            "nwchem_dft__xc": "b3lyp",
            "nwchem_property__dipole": True,
        }
    )

    val = qcdb.properties("nwc-b3lyp")
    check_dft(val)
