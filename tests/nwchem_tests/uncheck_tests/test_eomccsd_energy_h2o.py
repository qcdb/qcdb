#LR-CCSD on h2o

import os
import sys
from addons import *
from utils import *
import qcdb

h2o = qcdb.set_molecule('''
        O      0.000000000000     0.000000000000    -0.123909374404
        H      0.000000000000     1.429936611037     0.983265845431
        H      0.000000000000    -1.429936611037     0.983265845431
        ''')

print(h2o)

def test_lr_ccsd(return_value, is_df):
    if is_df:
        ref     =
