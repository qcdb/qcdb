#Xenon ZORA calcualtion

import os
import sys
import qcdb
from addons import *
from utils import *

def check_zora(return_value):
    dft =   -7499.170349126238 #dzvp orbital basis

    assert compare_values(dft, qcdb.get_variable('DFT TOTAL ENERGY'), 5, 'dft total')

def test_1_xe():
    xe = qcdb.set_molecule('Xe 0 0 0')

    qcdb.set_options({
        'basis' :   'dzvp (dft orbital)',
        'nwchem_dft__iterations' : '100',
        'nwchem_relativistic__zora': True,
        })
    val = qcdb.energy('nwc-dft')
    check_zora(val)
    
