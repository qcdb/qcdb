import os
import sys

import qcdb

from ..addons import *
from ..utils import *

#import pprinter as pp

def check_dipole(return_value):

    dip_x = -0.0000000000
    dip_y = -0.0000000000
    dip_z = -2.2247923027

    assert compare_values(dip_x, qcdb.variable('CURRENT DIPOLE X'), 5, 'dip x')
    assert compare_values(dip_y, qcdb.variable('CURRENT DIPOLE Y'), 5, 'dip y')
    assert compare_values(dip_z, qcdb.variable('CURRENT DIPOLE Z'), 5, 'dip z')

@using_nwchem
def test_1_dipole():
    qcdb.set_molecule(''' 
     O 0       0        0
     H 0       1.430   -1.107
     H 0      -1.430   -1.107
     units au
     ''')

    qcdb.set_options({
        'basis' : '6-31g*',
        'nwchem_scf__rohf': True,
        'nwchem_scf__singlet': True,
        'nwchem_property__dipole': True
        })

