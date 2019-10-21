#EOM CCSD
import os
import sys

import qcdb

from ..addons import *
from ..utils import *


def check_eomccsd(return_value):
    ref     =       -75.633713836043
    ccsd_tot=       -75.942608958748735
    ccsd_corl=      -0.308895122705897
    root1_exc=      0.317644190387170
    root1_tot=      -75.89014254
    root2_exc=      0.399118315871657
    root2_tot=      -75.80866841
    root3_exc=      0.419648919776505
    root3_tot=      -75.78813781
    root4_exc=      0.505153089073363
    root4_tot=      -75.70263364
    root5_exc=      0.569331468468014
    root5_tot=      -75.63845526
    root6_exc=      0.695937098625459
    root6_tot=      -75.51184963
    root7_exc=      0.995506121301140
    root7_tot=      -75.21228061
    root8_exc=      1.018576341192903
    root8_tot=      -75.18921039
    root9_exc=      1.048599044061786
    root9_tot=      -75.15918769
    root10_exc=     1.072014344987968
    root10_tot=     -75.13577239
    root11_exc=     1.123688743818949
    root11_tot=     -75.08409799
    root12_exc=     1.128736537578306
    root12_tot=     -75.07905019

#add all excited states being tested
    assert compare_values(ref, qcdb.get_variable('HF TOTAL ENERGY'), 6, 'hf ref')
    assert compare_values(ccsd_tot, qcdb.get_variable('CCSD TOTAL ENERGY'), 6, 'ccsd tot')
    assert compare_values(ccsd_corl, qcdb.get_variable('CCSD CORRELATION ENERGY'), 6, 'ccsd corl')
    assert compare_values(root1_exc, qcdb.get_variable('EOM-CCSD ROOT 0 -> ROOT 1 EXCITATION ENERGY - B1 SYMMETRY'), 6, 'EOM-CCSD 0 -> 1 ext')
    assert compare_values(root1_tot, qcdb.get_variable('EOM-CCSD ROOT 0 -> ROOT 1 TOTAL ENERGY - B1 SYMMETRY'), 6, 'EOM-CCSD 0 -> 1 total')
    assert compare_values(root2_exc, qcdb.get_variable('EOM-CCSD ROOT 0 -> ROOT 2 EXCITATION ENERGY - A2 SYMMETRY'), 6, 'EOM-CCSD 0 -> 2 ext')
    assert compare_values(root2_tot, qcdb.get_variable('EOM-CCSD ROOT 0 -> ROOT 2 TOTAL ENERGY - A2 SYMMETRY'), 6, 'EOM-CCSD 0 -> 2 total')
    assert compare_values(root3_exc, qcdb.get_variable('EOM-CCSD ROOT 0 -> ROOT 3 EXCITATION ENERGY - A1 SYMMETRY'), 6, 'EOM-CCSD 0 -> 3 ext')
    assert compare_values(root3_tot, qcdb.get_variable('EOM-CCSD ROOT 0 -> ROOT 3 TOTAL ENERGY - A1 SYMMETRY'), 6, 'EOM-CCSD 0 -> 3 total')
    assert compare_values(root4_exc, qcdb.get_variable('EOM-CCSD ROOT 0 -> ROOT 4 EXCITATION ENERGY - B2 SYMMETRY'), 6, 'EOM-CCSD 0 -> 4 ext')
    assert compare_values(root4_tot, qcdb.get_variable('EOM-CCSD ROOT 0 -> ROOT 4 TOTAL ENERGY - B2 SYMMETRY'), 6, 'EOM-CCSD 0 -> 4 total')
    assert compare_values(root5_exc, qcdb.get_variable('EOM-CCSD ROOT 0 -> ROOT 5 EXCITATION ENERGY - B2 SYMMETRY'), 6, 'EOM-CCSD 0 -> 5 ext')
    assert compare_values(root5_tot, qcdb.get_variable('EOM-CCSD ROOT 0 -> ROOT 5 TOTAL ENERGY - B2 SYMMETRY'), 6, 'EOM-CCSD 0 -> 5 total')
    assert compare_values(root6_exc, qcdb.get_variable('EOM-CCSD ROOT 0 -> ROOT 6 EXCITATION ENERGY - A1 SYMMETRY'), 6, 'EOM-CCSD 0 -> 6 ext')
    assert compare_values(root6_tot, qcdb.get_variable('EOM-CCSD ROOT 0 -> ROOT 6 TOTAL ENERGY - A1 SYMMETRY'), 6, 'EOM-CCSD 0 -> 6 total')
    assert compare_values(root7_exc, qcdb.get_variable('EOM-CCSD ROOT 0 -> ROOT 7 EXCITATION ENERGY - A1 SYMMETRY'), 6, 'EOM-CCSD 0 -> 7 ext')
    assert compare_values(root7_tot, qcdb.get_variable('EOM-CCSD ROOT 0 -> ROOT 7 TOTAL ENERGY - A1 SYMMETRY'), 6, 'EOM-CCSD 0 ->7 total')
    assert compare_values(root8_exc, qcdb.get_variable('EOM-CCSD ROOT 0 -> ROOT 8 EXCITATION ENERGY - A2 SYMMETRY'), 6, 'EOM-CCSD 0 -> 8 ext')
    assert compare_values(root8_tot, qcdb.get_variable('EOM-CCSD ROOT 0 -> ROOT 8 TOTAL ENERGY - A2 SYMMETRY'), 6, 'EOM-CCSD 0 -> 8 total')
    assert compare_values(root9_exc, qcdb.get_variable('EOM-CCSD ROOT 0 -> ROOT 9 EXCITATION ENERGY - B1 SYMMETRY'), 6, 'EOM-CCSD 0 -> 9 ext')
    assert compare_values(root9_tot, qcdb.get_variable('EOM-CCSD ROOT 0 -> ROOT 9 TOTAL ENERGY - B1 SYMMETRY'), 6, 'EOM-CCSD 0 -> 9 total')
    assert compare_values(root10_exc, qcdb.get_variable('EOM-CCSD ROOT 0 -> ROOT 10 EXCITATION ENERGY - A2 SYMMETRY'), 6, 'EOM-CCSD 0 -> 10 ext')
    assert compare_values(root10_tot, qcdb.get_variable('EOM-CCSD ROOT 0 -> ROOT 10 TOTAL ENERGY - A2 SYMMETRY'), 6, 'EOM-CCSD 0 -> 10 total')
    assert compare_values(root11_exc, qcdb.get_variable('EOM-CCSD ROOT 0 -> ROOT 11 EXCITATION ENERGY - B1 SYMMETRY'), 6, 'EOM-CCSD 0 -> 11 ext')
    assert compare_values(root11_tot, qcdb.get_variable('EOM-CCSD ROOT 0 -> ROOT 11 TOTAL ENERGY - B1 SYMMETRY'), 6, 'EOM-CCSD 0 -> 11 total')
    assert compare_values(root12_exc, qcdb.get_variable('EOM-CCSD ROOT 0 -> ROOT 12 EXCITATION ENERGY - B2 SYMMETRY'), 6, 'EOM-CCSD 0 -> 12 ext')
    assert compare_values(root12_tot, qcdb.get_variable('EOM-CCSD ROOT 0 -> ROOT 12 TOTAL ENERGY - B2 SYMMETRY'), 6, 'EOM-CCSD 0 -> 12 total')

@using_nwchem
def test_1_eomccsd():
    h2o = qcdb.set_molecule('''
        O      0.000000000000     0.000000000000    -0.123909374404
        H      0.000000000000     1.429936611037     0.983265845431
        H      0.000000000000    -1.429936611037     0.983265845431
        ''')

    qcdb.set_options({
        'basis'     :       '6-31g*',
        'memory'    :       '15000 mb',
        'scf__e_convergence'    : 1.0e-10,
        #'nwchem_memory' :   '1500 mb',
        #'nwchem_memory' :   '[total, 1500, stack, 400, heap, 400, global, 700, mb]', #The way nwchem speak for memory may need to change
        #'nwchem_stack_memory'  :   '400 mb',
        #'nwchem_heap_memory'   :   '400 mb',
        #'nwchem_global_memory' :   '700 mb',
        'nwchem_scf__thresh'     :   1.0e-10,
        'nwchem_scf__tol2e'      :   1.0e-10,
        'nwchem_scf__rhf'       :   True,
        'qc_module'             :   'tce',
        'nwchem_tce__ccsd'      :   True,
        'nwchem_tce__nroots'     :   12
        })
    print('Testing EOM-CCSD...')
    val = qcdb.energy('nwc-eom-ccsd')
    check_eomccsd(val)
