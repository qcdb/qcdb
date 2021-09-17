import os
import sys

import numpy as np
import pytest
import qcengine as qcng

import qcdb

from .addons import *
from .utils import *


@pytest.fixture
def h2o():
    return """
    O
    H 1 R
    H 1 R 2 A

    R=0.958
    A=104.5
    symmetry c1
"""

@pytest.fixture
def nh2():
    return """
0 2
N
H 1 R
H 1 R 2 A

R=1.008
A=105.0
"""

rhf_scf = np.array(
    [[ 0.00000000,  0.00000000, -0.02642604],
     [ 0.00000000, -0.01394333,  0.01321302],
     [ 0.00000000,  0.01394333,  0.01321302]])
uhf_scf = np.array(
    [[ 0.00000000,  0.00000000,  0.00155083],
     [ 0.00000000,  0.00009240, -0.00077542],
     [ 0.00000000, -0.00009240, -0.00077542]])
rohf_scf = np.array(
    [[ 0.00000000,  0.00000000,  0.00068377],
     [ 0.00000000, -0.00052602, -0.00034188],
     [ 0.00000000,  0.00052602, -0.00034188]])


@pytest.mark.parametrize('mtd,opts', [
    pytest.param('c4-hf', {'cfour_basis': 'qz2p', 'cfour_scf_conv': 12}, marks=using_cfour),
    pytest.param('c4-hf', {'basis': 'cfour-qz2p', 'cfour_scf_conv': 12}, marks=using_cfour),
    pytest.param('gms-hf', {'basis': 'cfour-qz2p'}, marks=using_gamess),
    pytest.param('nwc-hf', {'basis': 'cfour-qz2p', 'qc_module': 'tce', 'nwchem_scf__thresh': 12}, marks=using_nwchem), # tce doesn't matter here
    pytest.param('nwc-hf', {'basis': 'cfour-qz2p', 'nwchem_scf__thresh': 12}, marks=using_nwchem),
    pytest.param('p4-hf', {'basis': 'cfour-qz2p', 'psi4_scf_type' : 'direct'}, marks=using_psi4),
])
def test_sp_hf_rhf(mtd, opts, h2o):
    """cfour/???/input.dat
    #! single point HF/qz2p on water

    """
    h2o = qcdb.set_molecule(h2o)
    qcdb.set_options(opts)

    g, jrec = qcdb.gradient(mtd, return_wfn=True, molecule=h2o)

    print(g)
    assert compare_arrays(rhf_scf, g, atol=1.0e-5)
    # from cfour
    #scf_tot = -76.0627484601
    #atol = 1.e-6

    #assert compare_values(scf_tot, e, tnm() + ' SCF', atol=atol)
    #assert compare_values(scf_tot, qcdb.variable('hf total energy'), tnm() + ' SCF', atol=atol)
    #assert compare_values(scf_tot, qcdb.variable('scf total energy'), tnm() + ' SCF', atol=atol)
    #assert compare_values(scf_tot, qcdb.variable('current energy'), tnm() + ' SCF', atol=atol)


@pytest.mark.parametrize('mtd,opts', [
    pytest.param('c4-hf', {'cfour_basis': 'qz2p', 'cfour_reference': 'uhf', 'cfour_occupation': [[3,1,1,0],[3,0,1,0]], 'cfour_scf_conv': 12}, marks=using_cfour),
    pytest.param('c4-hf', {'basis': 'cfour-qz2p', 'cfour_reference': 'uhf'}, marks=using_cfour),
    pytest.param('gms-hf', {'basis': 'cfour-qz2p', 'gamess_contrl__scftyp': 'uhf'}, marks=using_gamess),
    pytest.param('nwc-hf', {'basis': 'cfour-qz2p', 'qc_module': 'tce', 'nwchem_scf__uhf': True}, marks=using_nwchem),
    pytest.param('nwc-hf', {'basis': 'cfour-qz2p', 'nwchem_scf__uhf': True}, marks=using_nwchem),
    pytest.param('p4-hf', {'basis': 'cfour-qz2p', 'reference': 'uhf', 'psi4_scf_type' : 'direct'}, marks=using_psi4),
])
def test_sp_hf_uhf(mtd, opts, nh2):
    """cfour/???/input.dat
    #! single-point HF/qz2p on NH2

    """
    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options(opts)

    g, jrec = qcdb.gradient(mtd, return_wfn=True, molecule=nh2)

    print(g)
    assert compare_arrays(uhf_scf, g, atol=1.0e-5)

    ## from cfour
    #scf_tot = -55.5893469688
    #atol = 1.e-6

    #assert compare_values(scf_tot, e, tnm() + ' SCF', atol=atol)
    #assert compare_values(scf_tot, qcdb.variable('hf total energy'), tnm() + ' SCF', atol=atol)
    #assert compare_values(scf_tot, qcdb.variable('scf total energy'), tnm() + ' SCF', atol=atol)
    #assert compare_values(scf_tot, qcdb.variable('current energy'), tnm() + ' SCF', atol=atol)


@pytest.mark.parametrize('mtd,opts', [
    pytest.param('c4-hf', {'cfour_basis': 'qz2p', 'cfour_reference': 'rohf', 'cfour_occupation': [[3,1,1,0],[3,0,1,0]], 'cfour_scf_conv': 12}, marks=using_cfour),
    pytest.param('c4-hf', {'basis': 'cfour-qz2p', 'cfour_reference': 'rohf'}, marks=using_cfour),
    pytest.param('gms-hf', {'basis': 'cfour-qz2p', 'gamess_contrl__scftyp': 'rohf'}, marks=using_gamess),
    pytest.param('nwc-hf', {'basis': 'cfour-qz2p', 'qc_module': 'tce', 'nwchem_scf__rohf': True}, marks=using_nwchem),
    pytest.param('nwc-hf', {'basis': 'cfour-qz2p', 'nwchem_scf__rohf': True}, marks=using_nwchem),
    pytest.param('p4-hf', {'basis': 'cfour-qz2p', 'reference': 'rohf', 'psi4_scf_type' : 'direct'}, marks=using_psi4),
])
def test_sp_hf_rohf(mtd, opts, nh2):
    """cfour/???/input.dat
    #! single-point HF/qz2p on NH2

    """
    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options(opts)

    g, jrec = qcdb.gradient(mtd, return_wfn=True, molecule=nh2)

    print(g)
    assert compare_arrays(rohf_scf, g, atol=1.0e-5)

    ## from cfour
    #scf_tot = -55.5847372601
    #atol = 1.e-6

    #assert compare_values(scf_tot, e, tnm() + ' SCF', atol=atol)
    #assert compare_values(scf_tot, qcdb.variable('hf total energy'), tnm() + ' SCF', atol=atol)
    #assert compare_values(scf_tot, qcdb.variable('scf total energy'), tnm() + ' SCF', atol=atol)
    #assert compare_values(scf_tot, qcdb.variable('current energy'), tnm() + ' SCF', atol=atol)
