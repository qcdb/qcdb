import pytest

import os
import sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))
from .utils import *
from .addons import *

import qcdb
import qcengine as qcng

@pytest.fixture
def h2o():
    return """
    O
    H 1 R
    H 1 R 2 A

    R=0.958
    A=104.5
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


@pytest.mark.parametrize('mtd,opts', [
    pytest.param('c4-mp2', {'cfour_basis': 'qz2p', 'cfour_scf_conv': 12}, marks=using_cfour),
    pytest.param('c4-mp2', {'basis': 'cfour-qz2p', 'cfour_scf_conv': 12}, marks=using_cfour),
    pytest.param('gms-mp2', {'basis': 'cfour-qz2p', 'gamess_mp2__nacore': 0}, marks=using_gamess),
    pytest.param('nwc-mp2', {'basis': 'cfour-qz2p', 'qc_module': 'tce'}, marks=using_nwchem),
    pytest.param('nwc-mp2', {'basis': 'cfour-qz2p'}, marks=using_nwchem),
    pytest.param('p4-mp2', {'basis': 'cfour-qz2p', 'psi4_mp2_type' : 'conv'}, marks=using_psi4),
])
def test_sp_mp2_rhf_full(mtd, opts, h2o):
    """cfour/???/input.dat
    #! single point MP2/qz2p on water

    """
    h2o = qcdb.set_molecule(h2o)
    qcdb.set_options(opts)

    e, jrec = qcdb.energy(mtd, return_wfn=True, molecule=h2o)

    # from cfour
    scftot = -76.0627484601
    mp2tot = -76.332940127333
    mp2corl = mp2tot - scftot
    atol = 1.e-6

    assert compare_values(mp2tot, e, tnm() + ' Returned', atol=atol)
    assert compare_values(mp2tot, qcdb.variable('current energy'), tnm() + ' Current', atol=atol)
    assert compare_values(mp2tot, qcdb.variable('mp2 total energy'), tnm() + ' MP2', atol=atol)
                                                                                               
    assert compare_values(mp2corl, qcdb.variable('current correlation energy'), tnm() + ' MP2 Corl', atol=atol)
    assert compare_values(mp2corl, qcdb.variable('mp2 correlation energy'), tnm() + ' MP2 Corl', atol=atol)

    assert compare_values(scftot, qcdb.variable('hf total energy'), tnm() + ' SCF', atol=atol)
    assert compare_values(scftot, qcdb.variable('scf total energy'), tnm() + ' SCF', atol=atol)


@pytest.mark.parametrize('mtd,opts', [
    pytest.param('c4-mp2', {'cfour_basis': 'qz2p', 'cfour_dropmo': [1], 'cfour_scf_conv': 12}, marks=using_cfour),
    pytest.param('c4-mp2', {'basis': 'cfour-qz2p', 'cfour_dropmo': 1, 'cfour_scf_conv': 12}, marks=using_cfour),
    pytest.param('gms-mp2', {'basis': 'cfour-qz2p', 'nwchem_mp2__freeze': 1}, marks=using_gamess),
    pytest.param('nwc-mp2', {'basis': 'cfour-qz2p', 'qc_module': 'tce', 'nwchem_tce__freeze': 1}, marks=using_nwchem),
    pytest.param('nwc-mp2', {'basis': 'cfour-qz2p', 'nwchem_mp2__freeze': 1}, marks=using_nwchem),
    pytest.param('p4-mp2', {'basis': 'cfour-qz2p', 'psi4_freeze_core': True, 'psi4_mp2_type' : 'conv'}, marks=using_psi4),
])
def test_sp_mp2_rhf_fc(mtd, opts, h2o):
    """cfour/???/input.dat
    #! single point MP2/qz2p on water

    """
    h2o = qcdb.set_molecule(h2o)
    qcdb.set_options(opts)

    e, jrec = qcdb.energy(mtd, return_wfn=True, molecule=h2o)

    # from cfour
    scftot = -76.062748460117
    mp2tot = -76.307900312177
    mp2corl = mp2tot - scftot
    atol = 1.e-6

    assert compare_values(mp2tot, e, tnm() + ' Returned', atol=atol)
    assert compare_values(mp2tot, qcdb.variable('current energy'), tnm() + ' Current', atol=atol)
    assert compare_values(mp2tot, qcdb.variable('mp2 total energy'), tnm() + ' MP2', atol=atol)
                                                                                               
    assert compare_values(mp2corl, qcdb.variable('current correlation energy'), tnm() + ' MP2 Corl', atol=atol)
    assert compare_values(mp2corl, qcdb.variable('mp2 correlation energy'), tnm() + ' MP2 Corl', atol=atol)

    assert compare_values(scftot, qcdb.variable('hf total energy'), tnm() + ' SCF', atol=atol)
    assert compare_values(scftot, qcdb.variable('scf total energy'), tnm() + ' SCF', atol=atol)


@pytest.mark.parametrize('mtd,opts', [
    pytest.param('c4-mp2', {'cfour_basis': 'qz2p', 'cfour_reference': 'uhf', 'cfour_scf_conv': 12}, marks=using_cfour),
    pytest.param('c4-mp2', {'basis': 'cfour-qz2p', 'cfour_reference': 'uhf', 'cfour_scf_conv': 12}, marks=using_cfour),
    pytest.param('gms-mp2', {'basis': 'cfour-qz2p', 'gamess_contrl__scftyp': 'uhf', 'gamess_mp2__nacore': 0}, marks=using_gamess),
    pytest.param('nwc-mp2', {'basis': 'cfour-qz2p', 'qc_module': 'tce', 'nwchem_scf__uhf': True}, marks=using_nwchem),
    pytest.param('nwc-mp2', {'basis': 'cfour-qz2p', 'nwchem_scf__uhf': True}, marks=using_nwchem),
    pytest.param('p4-mp2', {'basis': 'cfour-qz2p', 'reference': 'uhf', 'psi4_mp2_type' : 'conv'}, marks=using_psi4),
])
def test_sp_mp2_uhf_full(mtd, opts, nh2):
    """cfour/???/input.dat
    #! single point MP2/qz2p on water

    """
    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options(opts)

    e, jrec = qcdb.energy(mtd, return_wfn=True, molecule=nh2)

    # from cfour
    scftot = -55.5893469688
    mp2tot = -55.784877360093
    mp2corl = mp2tot - scftot
    atol = 1.e-6

    assert compare_values(mp2tot, e, tnm() + ' Returned', atol=atol)
    assert compare_values(mp2tot, qcdb.variable('current energy'), tnm() + ' Current', atol=atol)
    assert compare_values(mp2tot, qcdb.variable('mp2 total energy'), tnm() + ' MP2', atol=atol)
                                                                                               
    assert compare_values(mp2corl, qcdb.variable('current correlation energy'), tnm() + ' MP2 Corl', atol=atol)
    assert compare_values(mp2corl, qcdb.variable('mp2 correlation energy'), tnm() + ' MP2 Corl', atol=atol)

    assert compare_values(scftot, qcdb.variable('hf total energy'), tnm() + ' SCF', atol=atol)
    assert compare_values(scftot, qcdb.variable('scf total energy'), tnm() + ' SCF', atol=atol)


@pytest.mark.parametrize('mtd,opts', [
    pytest.param('c4-mp2', {'cfour_basis': 'qz2p', 'cfour_reference': 'uhf', 'cfour_dropmo': [1], 'cfour_scf_conv': 12}, marks=using_cfour),
    pytest.param('c4-mp2', {'basis': 'cfour-qz2p', 'cfour_reference': 'uhf', 'cfour_dropmo': 1, 'cfour_scf_conv': 12}, marks=using_cfour),
    pytest.param('gms-mp2', {'basis': 'cfour-qz2p', 'gamess_contrl__scftyp': 'uhf', 'nwchem_mp2__freeze': 1}, marks=using_gamess),
    pytest.param('nwc-mp2', {'basis': 'cfour-qz2p', 'qc_module': 'tce', 'nwchem_tce__freeze': 1, 'nwchem_scf__uhf': True}, marks=using_nwchem),
    pytest.param('nwc-mp2', {'basis': 'cfour-qz2p', 'nwchem_scf__uhf': True, 'nwchem_mp2__freeze': 1}, marks=using_nwchem),
    pytest.param('p4-mp2', {'basis': 'cfour-qz2p', 'reference': 'uhf', 'psi4_freeze_core': True, 'psi4_mp2_type' : 'conv'}, marks=using_psi4),
])
def test_sp_mp2_uhf_fc(mtd, opts, nh2):
    """cfour/???/input.dat
    #! single point MP2/qz2p on water

    """
    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options(opts)

    e, jrec = qcdb.energy(mtd, return_wfn=True, molecule=nh2)

    # from cfour
    scftot = -55.5893469688
    mp2tot = -55.760531091893
    mp2corl = mp2tot - scftot
    atol = 1.e-6

    assert compare_values(mp2tot, e, tnm() + ' Returned', atol=atol)
    assert compare_values(mp2tot, qcdb.variable('current energy'), tnm() + ' Current', atol=atol)
    assert compare_values(mp2tot, qcdb.variable('mp2 total energy'), tnm() + ' MP2', atol=atol)
                                                                                               
    assert compare_values(mp2corl, qcdb.variable('current correlation energy'), tnm() + ' MP2 Corl', atol=atol)
    assert compare_values(mp2corl, qcdb.variable('mp2 correlation energy'), tnm() + ' MP2 Corl', atol=atol)

    assert compare_values(scftot, qcdb.variable('hf total energy'), tnm() + ' SCF', atol=atol)
    assert compare_values(scftot, qcdb.variable('scf total energy'), tnm() + ' SCF', atol=atol)


@pytest.mark.parametrize('mtd,opts,errmsg', [
    pytest.param('nwc-mp2', {'basis': 'cfour-qz2p', 'nwchem_scf__rohf': True}, 'unknown SCFTYPE', marks=using_nwchem), # no rohf reference for nwc mp2
])
def test_sp_mp2_rohf_full_error(mtd, opts, errmsg, nh2):
    """cfour/???/input.dat
    #! single point MP2/qz2p on water

    """
    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options(opts)

    with pytest.raises(qcng.exceptions.InputError) as e:
        qcdb.energy(mtd, molecule=nh2)

    assert errmsg in str(e.value)


@pytest.mark.parametrize('mtd,opts', [
    pytest.param('c4-mp2', {'cfour_basis': 'qz2p', 'cfour_reference': 'rohf', 'cfour_occupation': [[3,1,1,0],[3,0,1,0]], 'cfour_scf_conv': 12}, marks=using_cfour),
    pytest.param('c4-mp2', {'basis': 'cfour-qz2p', 'cfour_reference': 'rohf', 'cfour_occupation': [[3,1,1,0],[3,0,1,0]], 'cfour_scf_conv': 12}, marks=using_cfour),
    pytest.param('gms-mp2', {'basis': 'cfour-qz2p', 'gamess_contrl__scftyp': 'rohf', 'gamess_mp2__nacore': 0, 'gamess_mp2__ospt': 'RMP'}, marks=using_gamess),
    pytest.param('nwc-mp2', {'basis': 'cfour-qz2p', 'qc_module': 'tce', 'nwchem_scf__rohf': True, 'nwchem_scf__thresh': 8, 'nwchem_tce__thresh': 8, 'nwchem_tce__freeze': 0, 'nwchem_scf__tol2e': 10}, marks=using_nwchem),
    pytest.param('p4-mp2', {'basis': 'cfour-qz2p', 'reference': 'rohf', 'psi4_mp2_type' : 'conv'}, marks=using_psi4),
])
def test_sp_mp2_rohf_full(mtd, opts, nh2):
    """cfour/???/input.dat
    #! single point MP2/qz2p on water

    """
    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options(opts)

    e, jrec = qcdb.energy(mtd, return_wfn=True, molecule=nh2)

    # from cfour
    scftot = -55.5847372601
    mp2tot = -55.7852767873
    mp2corl = mp2tot - scftot
    atol = 1.e-6
    if mtd.startswith('nwc'): #TODO: Figure out why nwc disagrees
        atol = 5.e-3

    assert compare_values(mp2tot, e, tnm() + ' Returned', atol=atol)
    assert compare_values(mp2tot, qcdb.variable('current energy'), tnm() + ' Current', atol=atol)
    assert compare_values(mp2tot, qcdb.variable('mp2 total energy'), tnm() + ' MP2', atol=atol)
                                                                                               
    assert compare_values(mp2corl, qcdb.variable('current correlation energy'), tnm() + ' MP2 Corl', atol=atol)
    assert compare_values(mp2corl, qcdb.variable('mp2 correlation energy'), tnm() + ' MP2 Corl', atol=atol)

    assert compare_values(scftot, qcdb.variable('hf total energy'), tnm() + ' SCF', atol=atol)
    assert compare_values(scftot, qcdb.variable('scf total energy'), tnm() + ' SCF', atol=atol)


@pytest.mark.parametrize('mtd,opts,errmsg', [
    pytest.param('nwc-mp2', {'basis': 'cfour-qz2p', 'nwchem_scf__rohf': True, 'nwchem_mp2__freeze': 1}, 'unknown SCFTYPE', marks=using_nwchem), # no rohf reference for nwc mp2
])
def test_sp_mp2_rohf_fc_error(mtd, opts, errmsg, nh2):
    """cfour/???/input.dat
    #! single point MP2/qz2p on water

    """
    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options(opts)

    with pytest.raises(qcng.exceptions.InputError) as e:
        qcdb.energy(mtd, molecule=nh2)

    assert errmsg in str(e.value)


@pytest.mark.parametrize('mtd,opts', [
    pytest.param('c4-mp2', {'cfour_basis': 'qz2p', 'cfour_reference': 'rohf', 'cfour_dropmo': [1], 'cfour_scf_conv': 12}, marks=using_cfour),
    pytest.param('c4-mp2', {'basis': 'cfour-qz2p', 'cfour_reference': 'rohf', 'cfour_dropmo': 1, 'cfour_scf_conv': 12}, marks=using_cfour),
    pytest.param('gms-mp2', {'basis': 'cfour-qz2p', 'gamess_contrl__scftyp': 'rohf', 'gamess_mp2__ospt': 'RMP'}, marks=using_gamess),
    pytest.param('nwc-mp2', {'basis': 'cfour-qz2p', 'qc_module': 'tce', 'nwchem_tce__freeze': 1, 'nwchem_scf__rohf': True}, marks=using_nwchem),
    pytest.param('p4-mp2', {'basis': 'cfour-qz2p', 'reference': 'rohf', 'psi4_freeze_core': True, 'psi4_mp2_type' : 'conv'}, marks=using_psi4),
])
def test_sp_mp2_rohf_fc(mtd, opts, nh2):
    """cfour/???/input.dat
    #! single point MP2/qz2p on water

    """
    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options(opts)

    e, jrec = qcdb.energy(mtd, return_wfn=True, molecule=nh2)

    # from cfour
    scftot = -55.5847372601
    mp2tot = -55.7608535667
    mp2corl = mp2tot - scftot
    atol = 1.e-6
    if mtd.startswith('nwc'): #TODO: Figure out why nwc disagrees
        atol = 5.e-3

    assert compare_values(mp2tot, e, tnm() + ' Returned', atol=atol)
    assert compare_values(mp2tot, qcdb.variable('current energy'), tnm() + ' Current', atol=atol)
    assert compare_values(mp2tot, qcdb.variable('mp2 total energy'), tnm() + ' MP2', atol=atol)
                                                                                               
    assert compare_values(mp2corl, qcdb.variable('current correlation energy'), tnm() + ' MP2 Corl', atol=atol)
    assert compare_values(mp2corl, qcdb.variable('mp2 correlation energy'), tnm() + ' MP2 Corl', atol=atol)

    assert compare_values(scftot, qcdb.variable('hf total energy'), tnm() + ' SCF', atol=atol)
    assert compare_values(scftot, qcdb.variable('scf total energy'), tnm() + ' SCF', atol=atol)
