#631g** H2O CCSD optimization by energies with z-matrix input

import os
import sys
import utils
import addons
import qcdb

print("cc2-h2o-energy/input.dat global testing")
h2o= qcdb.set_molecule("""
        O
        H 1 0.97
        H 1 0.97 2 103.0
        """)
print(h2o)
print(qcdb.get_active_options())

qcdb.set_options({'basis':"6-31G**"})
qcdb.energy('ccsd')


refnuc   =  9.166137300098260  #TEST
refscf   = -76.02293995043653  #TEST
refccsd  = -0.208238532935749  #TEST
reftotal = -76.2311784833722   #TEST

assert compare_values(refnuc, h2o.nuclear_repulsion_energy(),3, "Nuclear repulsion energy") #TEST
assert compare_values(refscf, qcdb.get_variable("SCF total energy"),5, "SCF energy")       #TEST
assert compare_values(refccsd,qcdb.get_variable("CCSD correlation energy"), 4, "CCSD contribution")        #TEST
assert compare_values(reftotal,qcdb.get_variable("Current energy"), 7, "Total energy")   #TEST

