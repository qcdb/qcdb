import os
import sys

from ..utils import *

import qcdb

def mcscf_ch2(return_value):

    hf  =   -38.857160053058
    mcscf   =   -38.913118315963

    assert compare_values(mcscf, qcdb.variable("MCSCF TOTAL ENERGY"), 5, "mcscf 1a1")

@using("nwchem")
def test_1_1a1():
    h2o = qcdb.set_molecule('''
            C    0   0     0
            H    0  1.87  -0.82
            H    0  -1.87  -0.82
            symmetry c2v
            units au''')

    qcdb.set_options({
        'basis' :   '6-31g**',
        'nwchem_mcscf__active'  :   6,
        'nwchem_mcscf__actelec' :   6,
        'nwchem_mcscf__state'   :   '1a1',
        })
    val = qcdb.energy('nwc-mcscf')
    mcscf_ch2(val)
