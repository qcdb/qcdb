import os
import sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))
from .utils import *
from .addons import *

import numpy as np

import qcdb

#! Various gradients for a strained helium dimer and water molecule

ref_scf_dz = np.array(
               [[ 0.0,  0.0,   0.01233095],
                [ 0.0,  0.0,  -0.01233095]])
ref_scf_tz = np.array( 
               [[ 0.0,  0.0,   0.01246097],
                [ 0.0,  0.0,  -0.01246097]])
ref_scf_dtz = np.array( 
               [[ 0.0,  0.0,   0.01249265],
                [ 0.0,  0.0,  -0.01249265]])
ref_scf_dtqz = np.array( 
               [[ 0.0,  0.0,   0.01244412],
                [ 0.0,  0.0,  -0.01244412]])
ref_mp2_dtz = np.array( 
               [[ 0.0,  0.0,   0.01155124],
                [ 0.0,  0.0,  -0.01155124]])

# y-axis, exchanged
permuted_indices_col = [ 0, 2, 1]
permuted_indices_row = [ 1, 0]
ref_scf_dz_y = ref_scf_dz[:, permuted_indices_col][permuted_indices_row, :]
ref_mp2_dtz_y = ref_mp2_dtz[:, permuted_indices_col][permuted_indices_row, :]

# y-axis, exchanged, fixed
#ref_scf_dz_yf = np.array(
#               [[ 0.0,  0.02466190,   0.0],
#                [ 0.0,  0.0,          0.0]])

nucenergy_ref = 1.17594935242 * a2a


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
    assert compare_values(nucenergy_ref, he_dimer.nuclear_repulsion_energy(), 9, "Nuclear repulsion energy")


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
    assert compare_values(nucenergy_ref, he_dimer.nuclear_repulsion_energy(), 9, "Nuclear repulsion energy")


# SCF TESTS

@using_psi4
def test_1a():
    system1()

    scf_dz, jrec = qcdb.gradient('SCF/cc-pVDZ', return_wfn=True)
    assert compare_arrays(ref_scf_dz, scf_dz, 6, "[1a] SCF/cc-pVDZ Gradient")
    assert compare_arrays(ref_scf_dz, jrec['qcvars']['CURRENT GRADIENT'].data, 6, "[1a] SCF/cc-pVDZ Hessian")
    assert compare_arrays(ref_scf_dz, qcdb.get_variable('CURRENT GRADIENT'), 6, "[1a] SCF/cc-pVDZ Hessian")
    #assert compare_arrays(ref_scf_dz, jrec['qcvars']['HF/CC-PVDZ TOTAL GRADIENT'].data, 6, "[1a] SCF/cc-pVDZ Hessian")
    # TODO provenance kill list
    assert ['QCDB', 'Psi4'] == [d['creator'] for d in jrec['provenance']], "[1a] prov"

    print(jrec['provenance'])


@using_cfour
def test_1b():
    system1()

    scf_dz, jrec = qcdb.gradient('c4-SCF/cc-pVDZ', return_wfn=True)
    assert compare_arrays(ref_scf_dz, scf_dz, 6, "[1b] SCF/cc-pVDZ Gradient")
    assert compare_arrays(ref_scf_dz, jrec['qcvars']['CURRENT GRADIENT'].data, 6, "[1b] SCF/cc-pVDZ Hessian")
    assert compare_arrays(ref_scf_dz, qcdb.get_variable('CURRENT GRADIENT'), 6, "[1b] SCF/cc-pVDZ Hessian")
    #assert compare_arrays(ref_scf_dz, jrec['qcvars']['HF/CC-PVDZ TOTAL GRADIENT'].data, 6, "[1b] SCF/cc-pVDZ Hessian")
    assert 'CFOUR' == jrec['provenance']['creator'], "[1b] prov"


@using_psi4
def test_1c():
    system2()

    scf_dz, jrec = qcdb.gradient('SCF/cc-pVDZ', return_wfn=True)
    assert compare_arrays(ref_scf_dz_y, scf_dz, 6, "[1c] SCF/cc-pVDZ Gradient")
    assert compare_arrays(ref_scf_dz_y, jrec['qcvars']['CURRENT GRADIENT'].data, 6, "[1c] SCF/cc-pVDZ Hessian")
    assert compare_arrays(ref_scf_dz_y, qcdb.get_variable('CURRENT GRADIENT'), 6, "[1c] SCF/cc-pVDZ Hessian")
    #assert compare_arrays(ref_scf_dz_y, jrec['qcvars']['HF/CC-PVDZ TOTAL GRADIENT'].data, 6, "[1c] SCF/cc-pVDZ Hessian")
    # TODO provenance kill list
    assert ['QCDB', 'Psi4'] == [d['creator'] for d in jrec['provenance']], "[1c] prov"
    print(jrec['provenance'])


@using_cfour
def test_1d():
    system2()

    scf_dz, jrec = qcdb.gradient('c4-SCF/cc-pVDZ', return_wfn=True)
    assert compare_arrays(ref_scf_dz_y, scf_dz, 6, "[1b] SCF/cc-pVDZ Gradient")
    assert compare_arrays(ref_scf_dz_y, jrec['qcvars']['CURRENT GRADIENT'].data, 6, "[1d] SCF/cc-pVDZ Hessian")
    assert compare_arrays(ref_scf_dz_y, qcdb.get_variable('CURRENT GRADIENT'), 6, "[1d] SCF/cc-pVDZ Hessian")
    #assert compare_arrays(ref_scf_dz_y, jrec['qcvars']['HF/CC-PVDZ TOTAL GRADIENT'].data, 6, "[1d] SCF/cc-pVDZ Hessian")
    assert 'CFOUR' == jrec['provenance']['creator'], "[1d] prov"


#def hide_test_2():
#    system1()
#
#    scf_tz = qcdb.gradient('SCF/cc-pVTZ', dertype=0)
#    assert compare_arrays(ref_scf_tz, scf_tz, 6, "[2] SCF/cc-pVTZ Gradient, dertype=0")
#
#
#def hide_test_3():
#    system1()
#
#    scf_dtz = qcdb.gradient('SCF/cc-pV[23]Z', dertype=0)
#    assert compare_arrays(ref_scf_dtz, scf_dtz, 6, "[3] SCF/cc-pV[DT]Z Gradient, dertype=0")


@using_psi4
def test_4a():
    system1()

    scf_dtz, jrec = qcdb.gradient('HF/cc-pV[23]Z', return_wfn=True)
    assert compare_arrays(ref_scf_dtz, scf_dtz, 6, "[4] SCF/cc-pV[DT]Z Gradient, dertype=1")
    pp.pprint(jrec)


@using_cfour
def test_4b():
    system1()

    scf_dtz, jrec = qcdb.gradient('c4-HF/cc-pV[23]Z', return_wfn=True)
    assert compare_arrays(ref_scf_dtz, scf_dtz, 6, "[4] SCF/cc-pV[DT]Z Gradient, dertype=1")
    pp.pprint(jrec)


@using_psi4
def test_5a():
    system1()

    scf_dtqz = qcdb.gradient('HF/cc-pV[DTQ]Z')
    assert compare_arrays(ref_scf_dtqz, scf_dtqz, 6, "[5] SCF/cc-pV[DTQ]Z Gradient")


@using_cfour
def test_5b():
    system1()

    scf_dtqz = qcdb.gradient('c4-HF/cc-pV[DTQ]Z')
    assert compare_arrays(ref_scf_dtqz, scf_dtqz, 6, "[5] SCF/cc-pV[DTQ]Z Gradient")


@using_psi4
def test_6a():
    system1()

    mp2_dtz = qcdb.gradient('MP2/cc-pV[DT]Z')
    assert compare_arrays(ref_mp2_dtz, mp2_dtz, 6, "[6] MP2/cc-pV[DT]Z Gradient")


@using_cfour
def test_6b():
    system1()

    mp2_dtz = qcdb.gradient('c4-MP2/cc-pV[DT]Z')
    assert compare_arrays(ref_mp2_dtz, mp2_dtz, 6, "[6] MP2/cc-pV[DT]Z Gradient")


@using_psi4
def test_6c():
    system2()

    mp2_dtz = qcdb.gradient('MP2/cc-pV[DT]Z')
    assert compare_arrays(ref_mp2_dtz_y, mp2_dtz, 6, "[6] MP2/cc-pV[DT]Z Gradient")


@using_cfour
def test_6d():
    system2()

    mp2_dtz = qcdb.gradient('c4-MP2/cc-pV[DT]Z')
    assert compare_arrays(ref_mp2_dtz_y, mp2_dtz, 6, "[6] MP2/cc-pV[DT]Z Gradient")


#def hide_test_7():
#    system1()
#
#    mp2_dtz = qcdb.gradient('MP2/cc-pV[DT]Z', dertype='energy')
#    assert compare_arrays(ref_mp2_dtz, mp2_dtz, 6, "[7] MP2/cc-pV[DT]Z Gradient, dertype=0")

