import os

import numpy as np

import qcdb

from .addons import *
from .utils import *

##! Various gradients for a strained helium dimer and water molecule
#
#ref_scf_dz = np.array(
#               [[ 0.0,  0.0,   0.01233095],
#                [ 0.0,  0.0,  -0.01233095]])
#ref_scf_tz = np.array( 
#               [[ 0.0,  0.0,   0.01246097],
#                [ 0.0,  0.0,  -0.01246097]])
#ref_scf_dtz = np.array( 
#               [[ 0.0,  0.0,   0.01249265],
#                [ 0.0,  0.0,  -0.01249265]])
#ref_scf_dtqz = np.array( 
#               [[ 0.0,  0.0,   0.01244412],
#                [ 0.0,  0.0,  -0.01244412]])
#ref_mp2_dtz = np.array( 
#               [[ 0.0,  0.0,   0.01155124],
#                [ 0.0,  0.0,  -0.01155124]])
#
## y-axis, exchanged
#permuted_indices_col = [ 0, 2, 1]
#permuted_indices_row = [ 1, 0]
#ref_scf_dz_y = ref_scf_dz[:, permuted_indices_col][permuted_indices_row, :]
#ref_mp2_dtz_y = ref_mp2_dtz[:, permuted_indices_col][permuted_indices_row, :]
#
## y-axis, exchanged, fixed
##ref_scf_dz_yf = np.array(
##               [[ 0.0,  0.02466190,   0.0],
##                [ 0.0,  0.0,          0.0]])



def system1():
    h2 = qcdb.set_molecule("""
    H
    H 1 R
    R = 1
    """)
    
    qcdb.set_options({
        'scf_type': 'pk',
        'mp2_type': 'conv',
        'g_convergence': 'GAU_VERYTIGHT',
        'e_convergence': 1.e-10,
    })
    
    h2.update_geometry()
    assert compare_values(0.529177208590000 * a2a, h2.nuclear_repulsion_energy(), 9, "Nuclear repulsion energy") #TEST

    return h2





test1_ene = -1.128746115958
test1_R = 0.747953788665

test2_ene = -1.1340651615708
test2_R = 0.730953222371

test3_ene = -1.1335602448947
test3_R = 0.7337084536

test4_ene = -1.1552182822114
test4_R = 0.754362336677

test5_ene = -1.1668917766245
test5_R = 0.735403792164

test6_ene = -1.1676372240400
test6_R = 0.735872194986

#test7_ene =
test7_R = 0.740686885481


@using_psi4
def test_1a():
    h2 = system1()
    refene = test1_ene
    lbl = tnm() + ' [1] SCF/cc-pVDZ Optimized R'
    assert compare_values(0.529177208590000 * a2a, h2.nuclear_repulsion_energy(), 9, "Nuclear repulsion energy")

    ene, jrec = qcdb.optking('SCF/cc-pVDZ', return_wfn=True, molecule=h2)
    assert compare_values(refene, ene, 6, lbl)
    assert compare_values(refene, jrec['qcvars']['CURRENT ENERGY'].data, 4, lbl)
    assert compare_values(refene, qcdb.variable('CURRENT ENERGY'), 4, lbl)
    assert compare_values(0.0, np.sqrt(np.mean(np.square(jrec['qcvars']['CURRENT GRADIENT'].data))), 4, lbl)
    assert compare_values(0.0, np.sqrt(np.mean(np.square(qcdb.variable('CURRENT GRADIENT')))), 4, lbl)
    assert compare_values(test1_R, h2.R, 4, lbl)
    #print(jrec['provenance'])


@using_cfour
def test_1b():
    h2 = system1()
    refene = test1_ene
    lbl = tnm() + ' [1] SCF/cc-pVDZ Optimized R'
    assert compare_values(0.529177208590000 * a2a, h2.nuclear_repulsion_energy(), 9, "Nuclear repulsion energy")

    ene, jrec = qcdb.optking('c4-SCF/cc-pVDZ', return_wfn=True, molecule=h2)
    assert compare_values(refene, ene, 6, lbl)
    assert compare_values(refene, jrec['qcvars']['CURRENT ENERGY'].data, 4, lbl)
    assert compare_values(refene, qcdb.variable('CURRENT ENERGY'), 4, lbl)
    assert compare_values(0.0, np.sqrt(np.mean(np.square(jrec['qcvars']['CURRENT GRADIENT'].data))), 4, lbl)
    assert compare_values(0.0, np.sqrt(np.mean(np.square(qcdb.variable('CURRENT GRADIENT')))), 4, lbl)
    assert compare_values(test1_R, h2.R, 4, lbl)
    #print(jrec['provenance'])


@using_psi4
def test_2a():
    h2 = system1()
    refene = test2_ene
    lbl = tnm() + ' [2] SCF/cc-pV[DT]Z Optimized R'
    assert compare_values(0.529177208590000 * a2a, h2.nuclear_repulsion_energy(), 9, "Nuclear repulsion energy")

    ene, jrec = qcdb.optking('SCF/cc-pV[DT]Z', return_wfn=True, molecule=h2)
    assert compare_values(refene, ene, 6, lbl)
    assert compare_values(refene, jrec['qcvars']['CURRENT ENERGY'].data, 4, lbl)
    assert compare_values(refene, qcdb.variable('CURRENT ENERGY'), 4, lbl)
    assert compare_values(0.0, np.sqrt(np.mean(np.square(jrec['qcvars']['CURRENT GRADIENT'].data))), 4, lbl)
    assert compare_values(0.0, np.sqrt(np.mean(np.square(qcdb.variable('CURRENT GRADIENT')))), 4, lbl)
    assert compare_values(test2_R, h2.R, 4, lbl)
    #print(jrec['provenance'])


@using_cfour
def test_2b():
    h2 = system1()
    refene = test2_ene
    lbl = tnm() + ' [2] SCF/cc-pV[DT]Z Optimized R'
    assert compare_values(0.529177208590000 * a2a, h2.nuclear_repulsion_energy(), 9, "Nuclear repulsion energy")

    ene, jrec = qcdb.optking('c4-SCF/cc-pV[DT]Z', return_wfn=True, molecule=h2)
    assert compare_values(refene, ene, 6, lbl)
    assert compare_values(refene, jrec['qcvars']['CURRENT ENERGY'].data, 4, lbl)
    assert compare_values(refene, qcdb.variable('CURRENT ENERGY'), 4, lbl)
    assert compare_values(0.0, np.sqrt(np.mean(np.square(jrec['qcvars']['CURRENT GRADIENT'].data))), 4, lbl)
    assert compare_values(0.0, np.sqrt(np.mean(np.square(qcdb.variable('CURRENT GRADIENT')))), 4, lbl)
    assert compare_values(test2_R, h2.R, 4, lbl)
    #print(jrec['provenance'])


@using_psi4
def test_3a():
    h2 = system1()
    refene = test3_ene
    lbl = tnm() + " [3] SCF/cc-pV[DTQ]Z Optimized R"
    assert compare_values(0.529177208590000 * a2a, h2.nuclear_repulsion_energy(), 9, "Nuclear repulsion energy")

    ene, jrec = qcdb.optking('SCF/cc-pV[DTQ]Z', return_wfn=True, molecule=h2)
    assert compare_values(refene, ene, 6, lbl)
    assert compare_values(refene, jrec['qcvars']['CURRENT ENERGY'].data, 4, lbl)
    assert compare_values(refene, qcdb.variable('CURRENT ENERGY'), 4, lbl)
    assert compare_values(0.0, np.sqrt(np.mean(np.square(jrec['qcvars']['CURRENT GRADIENT'].data))), 4, lbl)
    assert compare_values(0.0, np.sqrt(np.mean(np.square(qcdb.variable('CURRENT GRADIENT')))), 4, lbl)
    assert compare_values(test3_R, h2.R, 4, lbl)
    #print(jrec['provenance'])


@using_cfour
def test_3b():
    h2 = system1()
    refene = test3_ene
    lbl = tnm() + " [3] SCF/cc-pV[DTQ]Z Optimized R"
    assert compare_values(0.529177208590000 * a2a, h2.nuclear_repulsion_energy(), 9, "Nuclear repulsion energy")

    ene, jrec = qcdb.optking('c4-SCF/cc-pV[DTQ]Z', return_wfn=True, molecule=h2)
    assert compare_values(refene, ene, 6, lbl)
    assert compare_values(refene, jrec['qcvars']['CURRENT ENERGY'].data, 4, lbl)
    assert compare_values(refene, qcdb.variable('CURRENT ENERGY'), 4, lbl)
    assert compare_values(0.0, np.sqrt(np.mean(np.square(jrec['qcvars']['CURRENT GRADIENT'].data))), 4, lbl)
    assert compare_values(0.0, np.sqrt(np.mean(np.square(qcdb.variable('CURRENT GRADIENT')))), 4, lbl)
    assert compare_values(test3_R, h2.R, 4, lbl)
    #print(jrec['provenance'])

@using_psi4
def test_4a():
    h2 = system1()
    refene = test4_ene
    lbl = tnm() + " [4] MP2/cc-pVDZ Optimized R"
    assert compare_values(0.529177208590000 * a2a, h2.nuclear_repulsion_energy(), 9, "Nuclear repulsion energy")

    ene, jrec = qcdb.optking('MP2/cc-pVDZ', return_wfn=True, molecule=h2)
    assert compare_values(refene, ene, 6, lbl)
    assert compare_values(refene, jrec['qcvars']['CURRENT ENERGY'].data, 4, lbl)
    assert compare_values(refene, qcdb.variable('CURRENT ENERGY'), 4, lbl)
    assert compare_values(0.0, np.sqrt(np.mean(np.square(jrec['qcvars']['CURRENT GRADIENT'].data))), 4, lbl)
    assert compare_values(0.0, np.sqrt(np.mean(np.square(qcdb.variable('CURRENT GRADIENT')))), 4, lbl)
    assert compare_values(test4_R, h2.R, 4, lbl)
    #print(jrec['provenance'])


@using_cfour
def test_4b():
    h2 = system1()
    refene = test4_ene
    lbl = tnm() + " [4] MP2/cc-pVDZ Optimized R"
    assert compare_values(0.529177208590000 * a2a, h2.nuclear_repulsion_energy(), 9, "Nuclear repulsion energy")

    ene, jrec = qcdb.optking('c4-MP2/cc-pVDZ', return_wfn=True, molecule=h2)
    assert compare_values(refene, ene, 6, lbl)
    assert compare_values(refene, jrec['qcvars']['CURRENT ENERGY'].data, 4, lbl)
    assert compare_values(refene, qcdb.variable('CURRENT ENERGY'), 4, lbl)
    assert compare_values(0.0, np.sqrt(np.mean(np.square(jrec['qcvars']['CURRENT GRADIENT'].data))), 4, lbl)
    assert compare_values(0.0, np.sqrt(np.mean(np.square(qcdb.variable('CURRENT GRADIENT')))), 4, lbl)
    assert compare_values(test4_R, h2.R, 4, lbl)
    #print(jrec['provenance'])

@using_psi4
def test_5a():
    h2 = system1()
    refene = test5_ene
    lbl = tnm() + " [5] MP2/cc-pV[DT]Z Optimized R"
    assert compare_values(0.529177208590000 * a2a, h2.nuclear_repulsion_energy(), 9, "Nuclear repulsion energy")

    ene, jrec = qcdb.optking('MP2/cc-pV[DT]Z', return_wfn=True, molecule=h2)
    assert compare_values(refene, ene, 6, lbl)
    assert compare_values(refene, jrec['qcvars']['CURRENT ENERGY'].data, 4, lbl)
    assert compare_values(refene, qcdb.variable('CURRENT ENERGY'), 4, lbl)
    assert compare_values(0.0, np.sqrt(np.mean(np.square(jrec['qcvars']['CURRENT GRADIENT'].data))), 4, lbl)
    assert compare_values(0.0, np.sqrt(np.mean(np.square(qcdb.variable('CURRENT GRADIENT')))), 4, lbl)
    assert compare_values(test5_R, h2.R, 4, lbl)
    #print(jrec['provenance'])


@using_cfour
def test_5b():
    h2 = system1()
    refene = test5_ene
    lbl = tnm() + " [5] MP2/cc-pV[DT]Z Optimized R"
    assert compare_values(0.529177208590000 * a2a, h2.nuclear_repulsion_energy(), 9, "Nuclear repulsion energy")

    ene, jrec = qcdb.optking('c4-MP2/cc-pV[DT]Z', return_wfn=True, molecule=h2)
    assert compare_values(refene, ene, 6, lbl)
    assert compare_values(refene, jrec['qcvars']['CURRENT ENERGY'].data, 4, lbl)
    assert compare_values(refene, qcdb.variable('CURRENT ENERGY'), 4, lbl)
    assert compare_values(0.0, np.sqrt(np.mean(np.square(jrec['qcvars']['CURRENT GRADIENT'].data))), 4, lbl)
    assert compare_values(0.0, np.sqrt(np.mean(np.square(qcdb.variable('CURRENT GRADIENT')))), 4, lbl)
    assert compare_values(test5_R, h2.R, 4, lbl)
    #print(jrec['provenance'])

@using_psi4
def test_6a():
    h2 = system1()
    refene = test6_ene
    lbl = tnm() + " [6] MP2/cc-pV[TQ]Z Optimized R"
    assert compare_values(0.529177208590000 * a2a, h2.nuclear_repulsion_energy(), 9, "Nuclear repulsion energy")

    ene, jrec = qcdb.optking('MP2/cc-pV[TQ]Z', return_wfn=True, molecule=h2)
    assert compare_values(refene, ene, 6, lbl)
    assert compare_values(refene, jrec['qcvars']['CURRENT ENERGY'].data, 4, lbl)
    assert compare_values(refene, qcdb.variable('CURRENT ENERGY'), 4, lbl)
    assert compare_values(0.0, np.sqrt(np.mean(np.square(jrec['qcvars']['CURRENT GRADIENT'].data))), 4, lbl)
    assert compare_values(0.0, np.sqrt(np.mean(np.square(qcdb.variable('CURRENT GRADIENT')))), 4, lbl)
    assert compare_values(test6_R, h2.R, 4, lbl)
    #print(jrec['provenance'])


@using_cfour
def test_6b():
    h2 = system1()
    refene = test6_ene
    lbl = tnm() + " [6] MP2/cc-pV[TQ]Z Optimized R"
    assert compare_values(0.529177208590000 * a2a, h2.nuclear_repulsion_energy(), 9, "Nuclear repulsion energy")

    ene, jrec = qcdb.optking('c4-MP2/cc-pV[TQ]Z', return_wfn=True, molecule=h2)
    assert compare_values(refene, ene, 6, lbl)
    assert compare_values(refene, jrec['qcvars']['CURRENT ENERGY'].data, 4, lbl)
    assert compare_values(refene, qcdb.variable('CURRENT ENERGY'), 4, lbl)
    assert compare_values(0.0, np.sqrt(np.mean(np.square(jrec['qcvars']['CURRENT GRADIENT'].data))), 4, lbl)
    assert compare_values(0.0, np.sqrt(np.mean(np.square(qcdb.variable('CURRENT GRADIENT')))), 4, lbl)
    assert compare_values(test6_R, h2.R, 4, lbl)
    #print(jrec['provenance'])


    lbl = tnm() + " [7] CI2/cc-pV[DT]Z Optimized R"
##optimize('ci2/cc-pv[dt]z')
