import os
import sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))
from .utils import *
from .addons import *

import numpy as np

from qcengine.programs.util.hessparse import load_hessian

import qcdb

#! Various gradients for a strained helium dimer and water molecule

ref_e_scf_dz = -5.705182806503631
ref_scf_dz = np.array(
               [[ 0.0,  0.0,   0.01233095],
                [ 0.0,  0.0,  -0.01233095]])
ref_hess_scf_dz = load_hessian("""
    2    6
       -0.0036251364        0.0000000000        0.0000000000
        0.0036251364        0.0000000000        0.0000000000
        0.0000000000       -0.0036251364        0.0000000000
        0.0000000000        0.0036251364        0.0000000000
        0.0000000000        0.0000000000        0.0283246495
        0.0000000000        0.0000000000       -0.0283246495
        0.0036251364        0.0000000000        0.0000000000
       -0.0036251364        0.0000000000        0.0000000000
        0.0000000000        0.0036251364        0.0000000000
        0.0000000000       -0.0036251364        0.0000000000
        0.0000000000        0.0000000000       -0.0283246495
        0.0000000000        0.0000000000        0.0283246495
""", dtype='fcmfinal')
assert np.count_nonzero(ref_hess_scf_dz) == 12

ref_e_scf_tz = -5.716880280818783
ref_grad_scf_tz = np.array( 
               [[ 0.0,  0.0,   0.01246097],
                [ 0.0,  0.0,  -0.01246097]])
ref_hess_scf_tz = load_hessian("""
    2    6
       -0.0036633603        0.0000000000        0.0000000000
        0.0036633603        0.0000000000        0.0000000000
        0.0000000000       -0.0036633603        0.0000000000
        0.0000000000        0.0036633603        0.0000000000
        0.0000000000        0.0000000000        0.0279412505
        0.0000000000        0.0000000000       -0.0279412505
        0.0036633603        0.0000000000        0.0000000000
       -0.0036633603        0.0000000000        0.0000000000
        0.0000000000        0.0036633603        0.0000000000
        0.0000000000       -0.0036633603        0.0000000000
        0.0000000000        0.0000000000       -0.0279412505
        0.0000000000        0.0000000000        0.0279412505
""", dtype='fcmfinal')

ref_e_scf_dtz = -5.71973063
ref_grad_scf_dtz = np.array( 
               [[ 0.0,  0.0,   0.01249265],
                [ 0.0,  0.0,  -0.01249265]])
ref_hess_scf_dtz = np.array(
               [[-0.00367267,  0.        ,  0.        ,  0.00367267,  0.        ,  0.        ],
                [ 0.        , -0.00367267,  0.        ,  0.        ,  0.00367267,  0.        ],
                [ 0.        ,  0.        ,  0.02784783,  0.        ,  0.        , -0.02784783],
                [ 0.00367267,  0.        ,  0.        , -0.00367267,  0.        ,  0.        ],
                [ 0.        ,  0.00367267,  0.        ,  0.        , -0.00367267,  0.        ],
                [ 0.        ,  0.        , -0.02784783,  0.        ,  0.        ,  0.02784783]])

ref_e_scf_dtqz = -5.7176447
ref_grad_scf_dtqz = np.array( 
               [[ 0.0,  0.0,   0.01244412],
                [ 0.0,  0.0,  -0.01244412]])
ref_hess_scf_dtqz = np.array(
               [[-0.00365852,  0.        ,  0.        ,  0.00365852,  0.        ,  0.        ],
                [ 0.        , -0.00365852,  0.        ,  0.        ,  0.00365852,  0.        ],
                [ 0.        ,  0.        ,  0.02798016,  0.        ,  0.        , -0.02798016],
                [ 0.00365852,  0.        ,  0.        , -0.00365852,  0.        ,  0.        ],
                [ 0.        ,  0.00365852,  0.        ,  0.        , -0.00365852,  0.        ],
                [ 0.        ,  0.        , -0.02798016,  0.        ,  0.        ,  0.02798016]])

ref_e_mp2_dtz = -5.7898452
ref_grad_mp2_dtz = np.array( 
               [[ 0.0,  0.0,   0.01155124],
                [ 0.0,  0.0,  -0.01155124]])
ref_hess_mp2_dtz = np.array(
               [[-0.00339591,  0.        ,  0.        ,  0.00339591,  0.        ,  0.        ],
                [ 0.        , -0.00339591,  0.        ,  0.        ,  0.00339591,  0.        ],
                [ 0.        ,  0.        ,  0.02630844,  0.        ,  0.        , -0.02630844],
                [ 0.00339591,  0.        ,  0.        , -0.00339591,  0.        ,  0.        ],
                [ 0.        ,  0.00339591,  0.        ,  0.        , -0.00339591,  0.        ],
                [ 0.        ,  0.        , -0.02630844,  0.        ,  0.        ,  0.02630844]])

# y-axis, fixed
permuted_indices = [ 3, 5, 4, 0, 2, 1]
permuted_indices_col = [ 0, 2, 1]
permuted_indices_row = [ 1, 0]

ref_scf_dz_y = ref_scf_dz[:, permuted_indices_col][permuted_indices_row, :]
ref_hess_scf_dz_y = ref_hess_scf_dz[:, permuted_indices][permuted_indices, :]
ref_mp2_dtz_y = np.array( 
               [[ 0.0, -0.01155124,   0.0],
                [ 0.0,  0.01155124,   0.0]])
#ref_scf_dz_yf = np.array(
#               [[ 0.0,  0.02466190,   0.0],
#                [ 0.0,  0.0,          0.0]])

nucenergy_ref = 1.17594935242


def system1():
    he_dimer = qcdb.set_molecule("""
        He 0 0 0
        He 0 0 1.8
    """)
    
    # Get a reasonable guess, to save some iterations
    qcdb.set_options({
        'scf_type': 'pk',
        'mp2_type': 'conv',
        'reference': 'rhf',
    })
    
    he_dimer.update_geometry()
    compare_values(nucenergy_ref, he_dimer.nuclear_repulsion_energy(), 8, "Nuclear repulsion energy") #TEST


def system2():
    he_dimer = qcdb.set_molecule("""
        He 0 1.8 0
        He 0 0 0 
        no_reorient
        no_com
    """)
    
    # Get a reasonable guess, to save some iterations
    qcdb.set_options({
        'scf_type': 'pk',
        'mp2_type': 'conv',
        'reference': 'rhf',
    })
    
    he_dimer.update_geometry()
    compare_values(nucenergy_ref, he_dimer.nuclear_repulsion_energy(), 8, "Nuclear repulsion energy") #TEST


# SCF TESTS

@using_psi4
def test_1a():
    system1()
    lbl = "[1a] SCF/cc-pVDZ, Psi4"

    scf_dz, jrec = qcdb.hessian('SCF/cc-pVDZ', return_wfn=True)
    assert compare_arrays(ref_hess_scf_dz, scf_dz, 6, lbl)
    assert compare_arrays(ref_hess_scf_dz, qcdb.get_variable('CURRENT HESSIAN'), 6, lbl)
    assert compare_arrays(ref_hess_scf_dz, jrec['qcvars']['CURRENT HESSIAN'].data, 6, lbl)
    assert compare_arrays(ref_scf_dz, jrec['qcvars']['CURRENT GRADIENT'].data, 6, lbl)
    assert compare_arrays(ref_scf_dz, qcdb.get_variable('CURRENT GRADIENT'), 6, lbl)
    assert compare_values(ref_e_scf_dz, jrec['qcvars']['CURRENT ENERGY'].data, 6, lbl)
    assert compare_values(ref_e_scf_dz, qcdb.get_variable('CURRENT ENERGY'), 6, lbl)
    assert compare_arrays(ref_hess_scf_dz, jrec['qcvars']['HF TOTAL HESSIAN'].data, 6, lbl)
#    assert compare_arrays(ref_scf_dz, jrec['qcvars']['HF TOTAL GRADIENT'].data, 6, lbl)
    assert compare_values(ref_e_scf_dz, jrec['qcvars']['HF TOTAL ENERGY'].data, 6, lbl)
    #assert compare_arrays(ref_hess_scf_dz, jrec['qcvars']['HF/CC-PVDZ TOTAL HESSIAN'].data, 6, lbl)
    #assert compare_arrays(ref_scf_dz, jrec['qcvars']['HF/CC-PVDZ TOTAL GRADIENT'].data, 6, lbl)
    #assert compare_values(ref_e_scf_dz, jrec['qcvars']['HF/CC-PVDZ TOTAL ENERGY'].data, 6, lbl)
    # TODO provenance kill list
    assert ['QCDB', 'Psi4'] == [d['creator'] for d in jrec['provenance']], "[1a] prov"
    print(jrec['provenance'])


@using_cfour
def test_1b():
    system1()
    lbl = "[1b] SCF/cc-pVDZ, Cfour"

    scf_dz, jrec = qcdb.hessian('c4-SCF/cc-pVDZ', return_wfn=True)
    assert compare_arrays(ref_hess_scf_dz, scf_dz, 6, lbl)
    assert compare_arrays(ref_hess_scf_dz, qcdb.get_variable('CURRENT HESSIAN'), 6, lbl)
    assert compare_arrays(ref_hess_scf_dz, jrec['qcvars']['CURRENT HESSIAN'].data, 6, lbl)
    assert compare_arrays(ref_scf_dz, jrec['qcvars']['CURRENT GRADIENT'].data, 6,lbl)
    assert compare_arrays(ref_scf_dz, qcdb.get_variable('CURRENT GRADIENT'), 6, lbl)
#    assert compare_values(ref_e_scf_dz, jrec['qcvars']['CURRENT ENERGY'].data, 6, lbl)
#    assert compare_values(ref_e_scf_dz, qcdb.get_variable('CURRENT ENERGY'), 6, lbl)
#    assert compare_arrays(ref_hess_scf_dz, jrec['qcvars']['HF TOTAL HESSIAN'].data, 6, lbl)
#    assert compare_arrays(ref_scf_dz, jrec['qcvars']['HF TOTAL GRADIENT'].data, 6, lbl)
    assert compare_values(ref_e_scf_dz, jrec['qcvars']['HF TOTAL ENERGY'].data, 6, lbl)
    #assert compare_arrays(ref_hess_scf_dz, jrec['qcvars']['HF/CC-PVDZ TOTAL HESSIAN'].data, 6, lbl)
    #assert compare_arrays(ref_scf_dz, jrec['qcvars']['HF/CC-PVDZ TOTAL GRADIENT'].data, 6, lbl)
    #assert compare_values(ref_e_scf_dz, jrec['qcvars']['HF/CC-PVDZ TOTAL ENERGY'].data, 6, lbl)
    assert 'CFOUR' == jrec['provenance']['creator'], "[1b] prov"


@using_psi4
def test_1c():
    system2()
    lbl = "[1c] SCF/cc-pVDZ, Psi4"

    scf_dz, jrec = qcdb.hessian('SCF/cc-pVDZ', return_wfn=True)
    assert compare_arrays(ref_hess_scf_dz_y, scf_dz, 6, lbl)
    assert compare_arrays(ref_hess_scf_dz_y, qcdb.get_variable('CURRENT HESSIAN'), 6, lbl)
    assert compare_arrays(ref_hess_scf_dz_y, jrec['qcvars']['CURRENT HESSIAN'].data, 6, lbl)
    assert compare_arrays(ref_scf_dz_y, jrec['qcvars']['CURRENT GRADIENT'].data, 6, lbl)
    assert compare_arrays(ref_scf_dz_y, qcdb.get_variable('CURRENT GRADIENT'), 6, lbl)
    #assert compare_arrays(ref_hess_scf_dz_y, jrec['qcvars']['HF/CC-PVDZ TOTAL HESSIAN'].data, 6, lbl)
    #assert compare_arrays(ref_scf_dz_y, jrec['qcvars']['HF/CC-PVDZ TOTAL GRADIENT'].data, 6, lbl)
    # TODO provenance kill list
    assert ['QCDB', 'Psi4'] == [d['creator'] for d in jrec['provenance']], "[1c] prov"
    print(jrec['provenance'])


@using_cfour
def test_1d():
    system2()
    lbl = "[1d] SCF/cc-pVDZ, Cfour"

    scf_dz, jrec = qcdb.hessian('c4-SCF/cc-pVDZ', return_wfn=True)
    assert compare_arrays(ref_hess_scf_dz_y, scf_dz, 6, lbl)
    assert compare_arrays(ref_hess_scf_dz_y, qcdb.get_variable('CURRENT HESSIAN'), 6, lbl)
    assert compare_arrays(ref_hess_scf_dz_y, jrec['qcvars']['CURRENT HESSIAN'].data, 6, lbl)
    assert compare_arrays(ref_scf_dz_y, jrec['qcvars']['CURRENT GRADIENT'].data, 6, lbl)
    assert compare_arrays(ref_scf_dz_y, qcdb.get_variable('CURRENT GRADIENT'), 6, lbl)
    #assert compare_arrays(ref_hess_scf_dz_y, jrec['qcvars']['HF/CC-PVDZ TOTAL HESSIAN'].data, 6, lbl)
    #assert compare_arrays(ref_scf_dz_y, jrec['qcvars']['HF/CC-PVDZ TOTAL GRADIENT'].data, 6, lbl)
    assert 'CFOUR' == jrec['provenance']['creator'], "[1d] prov"
    print(jrec['provenance'])


@using_psi4
def test_2a():
    system1()
    lbl = "[2a] SCF/cc-pVDZ, Psi4, dertype=1"

    scf_dz, jrec = qcdb.hessian('SCF/cc-pVDZ', return_wfn=True, hess_dertype=1)
    print('HESS OUT')
    print(scf_dz)
    assert compare_arrays(ref_hess_scf_dz, scf_dz, 6, lbl)
    assert compare_arrays(ref_hess_scf_dz, qcdb.get_variable('CURRENT HESSIAN'), 6, lbl)
    assert compare_arrays(ref_hess_scf_dz, jrec['qcvars']['CURRENT HESSIAN'].data, 6, lbl)
    assert compare_arrays(ref_scf_dz, jrec['qcvars']['CURRENT GRADIENT'].data, 6, lbl)
    assert compare_arrays(ref_scf_dz, qcdb.get_variable('CURRENT GRADIENT'), 6, lbl)
    assert compare_values(ref_e_scf_dz, jrec['qcvars']['CURRENT ENERGY'].data, 6, lbl)
    assert compare_values(ref_e_scf_dz, qcdb.get_variable('CURRENT ENERGY'), 6, lbl)
#    assert compare_arrays(ref_hess_scf_dz, jrec['qcvars']['HF TOTAL HESSIAN'].data, 6, lbl)
#    assert compare_arrays(ref_scf_dz, jrec['qcvars']['HF TOTAL GRADIENT'].data, 6, lbl)
    assert compare_values(ref_e_scf_dz, jrec['qcvars']['HF TOTAL ENERGY'].data, 6, lbl)
    #assert compare_arrays(ref_hess_scf_dz, jrec['qcvars']['HF/CC-PVDZ TOTAL HESSIAN'].data, 6, lbl)
    #assert compare_arrays(ref_scf_dz, jrec['qcvars']['HF/CC-PVDZ TOTAL GRADIENT'].data, 6, lbl)
    #assert compare_values(ref_e_scf_dz, jrec['qcvars']['HF/CC-PVDZ TOTAL ENERGY'].data, 6, lbl)
    # TODO provenance kill list
    assert ['QCDB', 'Psi4'] == [d['creator'] for d in jrec['provenance']], "[1a] prov"
    print(jrec['provenance'])


@using_cfour
def test_2b():
    system1()
    lbl = "[2b] SCF/cc-pVDZ, Cfour, dertype=1"

    scf_dz, jrec = qcdb.hessian('c4-SCF/cc-pVDZ', return_wfn=True, hess_dertype=1)
    assert compare_arrays(ref_hess_scf_dz, scf_dz, 6, lbl)
    assert compare_arrays(ref_hess_scf_dz, qcdb.get_variable('CURRENT HESSIAN'), 6, lbl)
    assert compare_arrays(ref_hess_scf_dz, jrec['qcvars']['CURRENT HESSIAN'].data, 6, lbl)
    assert compare_arrays(ref_scf_dz, jrec['qcvars']['CURRENT GRADIENT'].data, 6,lbl)
    assert compare_arrays(ref_scf_dz, qcdb.get_variable('CURRENT GRADIENT'), 6, lbl)
#    assert compare_values(ref_e_scf_dz, jrec['qcvars']['CURRENT ENERGY'].data, 6, lbl)
#    assert compare_values(ref_e_scf_dz, qcdb.get_variable('CURRENT ENERGY'), 6, lbl)
#    assert compare_arrays(ref_hess_scf_dz, jrec['qcvars']['HF TOTAL HESSIAN'].data, 6, lbl)
#    assert compare_arrays(ref_scf_dz, jrec['qcvars']['HF TOTAL GRADIENT'].data, 6, lbl)
    assert compare_values(ref_e_scf_dz, jrec['qcvars']['HF TOTAL ENERGY'].data, 6, lbl)
    #assert compare_arrays(ref_hess_scf_dz, jrec['qcvars']['HF/CC-PVDZ TOTAL HESSIAN'].data, 6, lbl)
    #assert compare_arrays(ref_scf_dz, jrec['qcvars']['HF/CC-PVDZ TOTAL GRADIENT'].data, 6, lbl)
    #assert compare_values(ref_e_scf_dz, jrec['qcvars']['HF/CC-PVDZ TOTAL ENERGY'].data, 6, lbl)
    assert 'CFOUR' == jrec['provenance']['creator'], "[1b] prov"


@using_psi4
def test_2c():
    system2()
    lbl = "[2c] SCF/cc-pVDZ, Psi4, dertype=1"

    scf_dz, jrec = qcdb.hessian('SCF/cc-pVDZ', return_wfn=True, dertype=1)
    assert compare_arrays(ref_hess_scf_dz_y, scf_dz, 6, lbl)
    assert compare_arrays(ref_hess_scf_dz_y, qcdb.get_variable('CURRENT HESSIAN'), 6, lbl)
    assert compare_arrays(ref_hess_scf_dz_y, jrec['qcvars']['CURRENT HESSIAN'].data, 6, lbl)
    assert compare_arrays(ref_scf_dz_y, jrec['qcvars']['CURRENT GRADIENT'].data, 6, lbl)
    assert compare_arrays(ref_scf_dz_y, qcdb.get_variable('CURRENT GRADIENT'), 6, lbl)
    #assert compare_arrays(ref_hess_scf_dz_y, jrec['qcvars']['HF/CC-PVDZ TOTAL HESSIAN'].data, 6, lbl)
    #assert compare_arrays(ref_scf_dz_y, jrec['qcvars']['HF/CC-PVDZ TOTAL GRADIENT'].data, 6, lbl)
    # TODO provenance kill list
    assert ['QCDB', 'Psi4'] == [d['creator'] for d in jrec['provenance']], "[1c] prov"
    print(jrec['provenance'])


@using_cfour
def test_2d():
    system2()
    lbl = "[2d] SCF/cc-pVDZ, Cfour, dertype=1"

    scf_dz, jrec = qcdb.hessian('c4-SCF/cc-pVDZ', return_wfn=True, dertype=1)
    assert compare_arrays(ref_hess_scf_dz_y, scf_dz, 6, lbl)
    assert compare_arrays(ref_hess_scf_dz_y, qcdb.get_variable('CURRENT HESSIAN'), 6, lbl)
    assert compare_arrays(ref_hess_scf_dz_y, jrec['qcvars']['CURRENT HESSIAN'].data, 6, lbl)
    assert compare_arrays(ref_scf_dz_y, jrec['qcvars']['CURRENT GRADIENT'].data, 6, lbl)
    assert compare_arrays(ref_scf_dz_y, qcdb.get_variable('CURRENT GRADIENT'), 6, lbl)
    #assert compare_arrays(ref_hess_scf_dz_y, jrec['qcvars']['HF/CC-PVDZ TOTAL HESSIAN'].data, 6, lbl)
    #assert compare_arrays(ref_scf_dz_y, jrec['qcvars']['HF/CC-PVDZ TOTAL GRADIENT'].data, 6, lbl)
    assert 'CFOUR' == jrec['provenance']['creator'], "[1d] prov"
    print(jrec['provenance'])


##def hide_test_2():
##    system1()
##
##    scf_tz = qcdb.hessian('SCF/cc-pVTZ', dertype=0)
##    assert compare_arrays(ref_grad_scf_tz, scf_tz, 6, "[2] SCF/cc-pVTZ Gradient, dertype=0")


##def hide_test_3():
##    system1()
##
##    scf_dtz = qcdb.hessian('SCF/cc-pV[23]Z', dertype=0)
##    assert compare_arrays(ref_grad_scf_dtz, scf_dtz, 6, "[3] SCF/cc-pV[DT]Z Gradient, dertype=0")


@using_psi4
def test_4a():
    system1()
    lbl = "[4a] SCF/cc-pV[DT]Z Gradient, dertype=1"

    scf_dtz, jrec = qcdb.hessian('HF/cc-pV[23]Z', return_wfn=True)
    assert compare_arrays(ref_hess_scf_dtz, scf_dtz, 6, lbl)
    assert compare_arrays(ref_hess_scf_dtz, qcdb.get_variable('CURRENT HESSIAN'), 6, lbl)
    assert compare_arrays(ref_hess_scf_dtz, jrec['qcvars']['CURRENT HESSIAN'].data, 6, lbl)
    assert compare_arrays(ref_grad_scf_dtz, jrec['qcvars']['CURRENT GRADIENT'].data, 6, lbl)
    assert compare_values(ref_e_scf_dtz, jrec['qcvars']['CURRENT ENERGY'].data, 6, lbl)
    assert compare_arrays(ref_grad_scf_dtz, qcdb.get_variable('CURRENT GRADIENT'), 6, lbl)
    assert compare_values(ref_e_scf_dtz, qcdb.get_variable('CURRENT ENERGY'), 6, lbl)
    #assert compare_arrays(ref_hess_scf_dtz, jrec['qcvars']['HF/CC-PV[DTQ]Z TOTAL HESSIAN'].data, 6, lbl)
    #assert compare_arrays(ref_grad_scf_dtz, jrec['qcvars']['HF/CC-PV[DTQ]Z TOTAL GRADIENT'].data, 6, lbl)
    #assert compare_arrays(ref_e_scf_dtz, jrec['qcvars']['HF/CC-PV[DTQ]Z TOTAL ENNERGY'].data, 6, lbl)


@using_cfour
def test_4b():
    system1()
    lbl = "[4b] SCF/cc-pV[DT]Z Cfour, dertype=1"

    scf_dtz, jrec = qcdb.hessian('c4-HF/cc-pV[23]Z', return_wfn=True)
    assert compare_arrays(ref_hess_scf_dtz, scf_dtz, 6, lbl)
    assert compare_arrays(ref_hess_scf_dtz, qcdb.get_variable('CURRENT HESSIAN'), 6, lbl)
    assert compare_arrays(ref_hess_scf_dtz, jrec['qcvars']['CURRENT HESSIAN'].data, 6, lbl)
    assert compare_arrays(ref_grad_scf_dtz, jrec['qcvars']['CURRENT GRADIENT'].data, 6, lbl)
#    assert compare_values(ref_e_scf_dtz, jrec['qcvars']['CURRENT ENERGY'].data, 6, lbl)
    assert compare_arrays(ref_grad_scf_dtz, qcdb.get_variable('CURRENT GRADIENT'), 6, lbl)
#    assert compare_values(ref_e_scf_dtz, qcdb.get_variable('CURRENT ENERGY'), 6, lbl)
    #assert compare_arrays(ref_hess_scf_dtz, jrec['qcvars']['HF/CC-PV[DTQ]Z TOTAL HESSIAN'].data, 6, lbl)
    #assert compare_arrays(ref_grad_scf_dtz, jrec['qcvars']['HF/CC-PV[DTQ]Z TOTAL GRADIENT'].data, 6, lbl)
    #assert compare_arrays(ref_e_scf_dtz, jrec['qcvars']['HF/CC-PV[DTQ]Z TOTAL ENNERGY'].data, 6, lbl)


@using_psi4
def test_5a():
    system1()
    lbl = "[5a] HF/cc-pV[DTQ]Z, Psi4"

    scf_dtqz, jrec = qcdb.hessian('HF/cc-pV[DTQ]Z', return_wfn=True)
    assert compare_arrays(ref_hess_scf_dtqz, scf_dtqz, 6, lbl)
    assert compare_arrays(ref_hess_scf_dtqz, qcdb.get_variable('CURRENT HESSIAN'), 6, lbl)
    assert compare_arrays(ref_hess_scf_dtqz, jrec['qcvars']['CURRENT HESSIAN'].data, 6, lbl)
    assert compare_arrays(ref_grad_scf_dtqz, jrec['qcvars']['CURRENT GRADIENT'].data, 6, lbl)
    assert compare_values(ref_e_scf_dtqz, jrec['qcvars']['CURRENT ENERGY'].data, 6, lbl)
    assert compare_arrays(ref_grad_scf_dtqz, qcdb.get_variable('CURRENT GRADIENT'), 6, lbl)
    assert compare_values(ref_e_scf_dtqz, qcdb.get_variable('CURRENT ENERGY'), 6, lbl)
    #assert compare_arrays(ref_hess_scf_dtqz, jrec['qcvars']['HF/CC-PV[DTQ]Z TOTAL HESSIAN'].data, 6, lbl)
    #assert compare_arrays(ref_grad_scf_dtqz, jrec['qcvars']['HF/CC-PV[DTQ]Z TOTAL GRADIENT'].data, 6, lbl)
    #assert compare_arrays(ref_e_scf_dtqz, jrec['qcvars']['HF/CC-PV[DTQ]Z TOTAL ENNERGY'].data, 6, lbl)


@using_cfour
def test_5b():
    system1()
    lbl = "[5b] SCF/cc-pV[DTQ]Z, Cfour"

    scf_dtqz, jrec = qcdb.hessian('c4-HF/cc-pV[DTQ]Z', return_wfn=True)
    assert compare_arrays(ref_hess_scf_dtqz, scf_dtqz, 6, lbl)
    assert compare_arrays(ref_hess_scf_dtqz, qcdb.get_variable('CURRENT HESSIAN'), 6, lbl)
    assert compare_arrays(ref_hess_scf_dtqz, jrec['qcvars']['CURRENT HESSIAN'].data, 6, lbl)
    assert compare_arrays(ref_grad_scf_dtqz, jrec['qcvars']['CURRENT GRADIENT'].data, 6, lbl)
#    assert compare_values(ref_e_scf_dtqz, jrec['qcvars']['CURRENT ENERGY'].data, 6, lbl)
    assert compare_arrays(ref_grad_scf_dtqz, qcdb.get_variable('CURRENT GRADIENT'), 6, lbl)
#    assert compare_values(ref_e_scf_dtqz, qcdb.get_variable('CURRENT ENERGY'), 6, lbl)
    #assert compare_arrays(ref_hess_scf_dtqz, jrec['qcvars']['HF/CC-PV[DTQ]Z TOTAL HESSIAN'].data, 6, lbl)
    #assert compare_arrays(ref_grad_scf_dtqz, jrec['qcvars']['HF/CC-PV[DTQ]Z TOTAL GRADIENT'].data, 6, lbl)
    #assert compare_arrays(ref_e_scf_dtqz, jrec['qcvars']['HF/CC-PV[DTQ]Z TOTAL ENNERGY'].data, 6, lbl)


# MP2 TESTS

@using_psi4
def test_6a():
    system1()
    lbl = '[6a] MP2/cc-pV[DT]Z, Psi4'

    # No MP2 freqs
    mp2_dtz, jrec = qcdb.hessian('MP2/cc-pV[DT]Z', return_wfn=True)
    assert compare_arrays(ref_hess_mp2_dtz, mp2_dtz, 6, lbl)
    assert compare_arrays(ref_hess_mp2_dtz, qcdb.get_variable('CURRENT HESSIAN'), 6, lbl)
    assert compare_arrays(ref_hess_mp2_dtz, jrec['qcvars']['CURRENT HESSIAN'].data, 6, lbl)
    assert compare_arrays(ref_grad_mp2_dtz, jrec['qcvars']['CURRENT GRADIENT'].data, 6, lbl)
    assert compare_arrays(ref_grad_mp2_dtz, qcdb.get_variable('CURRENT GRADIENT'), 6, lbl)
    assert compare_values(ref_e_mp2_dtz, jrec['qcvars']['CURRENT ENERGY'].data, 6, lbl)
    assert compare_values(ref_e_mp2_dtz, qcdb.get_variable('CURRENT ENERGY'), 6, lbl)
    #assert compare_arrays(ref_hess_mp2_dtz, jrec['qcvars']['MP2/CC-PV[DT]Z TOTAL HESSIAN'].data, 6, lbl)
    #assert compare_arrays(ref_grad_mp2_dtz, jrec['qcvars']['MP2/CC-PV[DT]Z TOTAL GRADIENT'].data, 6, lbl)
    #assert compare_values(ref_e_mp2_dtz, jrec['qcvars']['MP2/CC-PV[DT]Z TOTAL ENERGY'].data, 6, lbl)
    #assert ['QCDB', 'Psi4'] == [d['creator'] for d in jrec['provenance']], "[1a] prov"


@using_cfour
def test_6b():
    system1()
    lbl = '[6b] MP2/cc-pV[DT]Z, Cfour'

    mp2_dtz, jrec = qcdb.hessian('c4-MP2/cc-pV[DT]Z', return_wfn=True)
    assert compare_arrays(ref_hess_mp2_dtz, mp2_dtz, 6, lbl)
    assert compare_arrays(ref_hess_mp2_dtz, qcdb.get_variable('CURRENT HESSIAN'), 6, lbl)
    assert compare_arrays(ref_hess_mp2_dtz, jrec['qcvars']['CURRENT HESSIAN'].data, 6, lbl)
    assert compare_arrays(ref_grad_mp2_dtz, jrec['qcvars']['CURRENT GRADIENT'].data, 6, lbl)
    assert compare_arrays(ref_grad_mp2_dtz, qcdb.get_variable('CURRENT GRADIENT'), 6, lbl)
# cfour parse is finding MP2 as currE instead of HF
#    assert compare_values(ref_e_mp2_dtz, jrec['qcvars']['CURRENT ENERGY'].data, 6, lbl)
#    assert compare_values(ref_e_mp2_dtz, qcdb.get_variable('CURRENT ENERGY'), 6, lbl)
    #assert compare_arrays(ref_hess_mp2_dtz, jrec['qcvars']['MP2/CC-PV[DT]Z TOTAL HESSIAN'].data, 6, lbl)
    #assert compare_arrays(ref_grad_mp2_dtz, jrec['qcvars']['MP2/CC-PV[DT]Z TOTAL GRADIENT'].data, 6, lbl)
    #assert compare_values(ref_e_mp2_dtz, jrec['qcvars']['MP2/CC-PV[DT]Z TOTAL ENERGY'].data, 6, lbl)


#def test_6c():
#    system2()
#
#    mp2_dtz = qcdb.hessian('MP2/cc-pV[DT]Z')
#    assert compare_arrays(ref_mp2_dtz_y, mp2_dtz, 6, "[6] MP2/cc-pV[DT]Z Gradient")
#
#
#def test_6d():
#    system2()
#
#    mp2_dtz = qcdb.hessian('c4-MP2/cc-pV[DT]Z')
#    assert compare_arrays(ref_mp2_dtz_y, mp2_dtz, 6, "[6] MP2/cc-pV[DT]Z Gradient")
#
#
##def hide_test_7():
##    system1()
##
##    mp2_dtz = qcdb.hessian('MP2/cc-pV[DT]Z', dertype='energy')
##    assert compare_arrays(ref_mp2_dtz, mp2_dtz, 6, "[7] MP2/cc-pV[DT]Z Gradient, dertype=0")
#
