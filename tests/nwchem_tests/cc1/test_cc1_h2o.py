import os
import sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))
from utils import *
from addons import *

@using_nwchem
def test_cc1():
    """ cc1-h20-energy/input.dat
    global testing
    """
    h2o = qcdb.set_molecule("""
    O
    H 1 0.97
    H 1 0.97 2 103.0
    """)
    print(h2o)
    print(qcdb.get_active_options.print_changed())
