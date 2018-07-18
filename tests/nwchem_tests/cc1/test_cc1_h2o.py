import os
import sys
from utils import *
from addons import *

print(" cc1-h20-energy/input.dat global testing ")
h2o = qcdb.set_molecule("""
    O
    H 1 0.97
    H 1 0.97 2 103.0
    """)
print(h2o)
print(qcdb.get_active_options.print_changed())

qcdb.set_options({basis: '6-31G**'})

qcdb.energy('ccsd')

refnuc   =   9.1654609427539  #TEST
refscf   = -76.0229427274435  #TEST
refccsd  = -0.20823570806196  #TEST
reftotal = -76.2311784355056  #TEST

assert compare_values(refnuc,   h2o.nuclear_repulsion_energy(),          3, "Nuclear repulsion energy") #TEST
assert compare_values(refscf,   qcdb.get_variable("SCF total energy"),        5, "SCF energy")               #TEST
assert compare_values(refccsd,  qcdb.get_variable("CCSD correlation energy"), 4, "CCSD contribution")        #TEST
assert compare_values(reftotal, qcdb.get_variable("Current energy"),          7, "Total energy")             #TEST

