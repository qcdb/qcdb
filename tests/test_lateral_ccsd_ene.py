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
    pytest.param('c4-ccsd', {'cfour_BASIS': 'qz2p', 'cfour_SCF_CONV': 12, 'cfour_CC_CONV': 12}, marks=using_cfour),
    pytest.param('c4-ccsd', {'BASIS': 'cfour-qz2p'}, marks=using_cfour),
    pytest.param('gms-ccsd', {'basis': 'cfour-qz2p', 'gamess_ccinp__ncore': 0}, marks=using_gamess),
    pytest.param('nwc-ccsd', {'basis': 'cfour-qz2p'}, marks=using_nwchem),
    pytest.param('nwc-ccsd', {'basis': 'cfour-qz2p', 'qc_module': 'tce'}, marks=using_nwchem),
    pytest.param('p4-ccsd', {'basis': 'cfour-qz2p'}, marks=using_psi4),
])
def test_sp_ccsd_rhf_ae(mtd, opts, h2o):
    """cfour/sp-rhf-ccsd/input.dat
    #! single point CCSD/qz2p on water

    """
    h2o = qcdb.set_molecule(h2o)
    qcdb.set_options(opts)

    e, jrec = qcdb.energy(mtd, return_wfn=True, molecule=h2o)

    scf_tot = -76.062748460117
    mp2_tot = -76.332940127333
    ccsd_tot = -76.338453951890

    mp2_corl = mp2_tot - scf_tot
    ccsd_corl = ccsd_tot - scf_tot

    atol = 1.e-6

    # cc terms
    assert compare_values(ccsd_tot, e, tnm() + ' Returned', atol=atol)
    assert compare_values(ccsd_tot, qcdb.variable('current energy'), tnm() + ' Current', atol=atol)
    assert compare_values(ccsd_tot, qcdb.variable('ccsd total energy'), tnm() + ' CCSD', atol=atol)
    assert compare_values(ccsd_corl, qcdb.variable('current correlation energy'), tnm() + ' CCSD', atol=atol)
    assert compare_values(ccsd_corl, qcdb.variable('ccsd correlation energy'), tnm() + ' CCSD', atol=atol)

    # mp2 terms (not printed for nwc tce)
    if not (mtd.startswith('nwc') and opts.get('qc_module', 'nein').lower() == 'tce'):
        assert compare_values(mp2_tot, qcdb.variable('mp2 total energy'), tnm() + ' MP2', atol=atol)
        assert compare_values(mp2_corl, qcdb.variable('mp2 correlation energy'), tnm() + ' MP2', atol=atol)

    # hf terms
    assert compare_values(scf_tot, qcdb.variable('hf total energy'), tnm() + ' SCF', atol=atol)
    assert compare_values(scf_tot, qcdb.variable('scf total energy'), tnm() + ' SCF', atol=atol)


@pytest.mark.parametrize('mtd,opts', [
    pytest.param('c4-ccsd', {'cfour_BASIS': 'qz2p', 'cfour_dropmo': [1], 'cfour_SCF_CONV': 12, 'cfour_CC_CONV': 12}, marks=using_cfour),
    pytest.param('c4-ccsd', {'BASIS': 'cfour-qz2p', 'cfour_dropmo': 1}, marks=using_cfour),
    pytest.param('gms-ccsd', {'basis': 'cfour-qz2p'}, marks=using_gamess),
    pytest.param('nwc-ccsd', {'basis': 'cfour-qz2p', 'nwchem_ccsd__freeze': 1}, marks=using_nwchem),
    pytest.param('nwc-ccsd', {'basis': 'cfour-qz2p', 'qc_module': 'tce', 'nwchem_tce__freeze': 1}, marks=using_nwchem),
    pytest.param('p4-ccsd', {'basis': 'cfour-qz2p', 'psi4_freeze_core': True}, marks=using_psi4),
])
def test_sp_ccsd_rhf_fc(mtd, opts, h2o):
    """cfour/sp-rhf-ccsd/input.dat
    #! single point CCSD/qz2p on water

    """
    h2o = qcdb.set_molecule(h2o)
    qcdb.set_options(opts)

    e, jrec = qcdb.energy(mtd, return_wfn=True, molecule=h2o)

    scf_tot = -76.062748460117
    mp2_tot = -76.307900312177
    ccsd_tot = -76.313079023615

    mp2_corl = mp2_tot - scf_tot
    ccsd_corl = ccsd_tot - scf_tot

    atol = 1.e-6

    # cc terms
    assert compare_values(ccsd_tot, e, tnm() + ' Returned', atol=atol)
    assert compare_values(ccsd_tot, qcdb.variable('current energy'), tnm() + ' Current', atol=atol)
    assert compare_values(ccsd_tot, qcdb.variable('ccsd total energy'), tnm() + ' CCSD', atol=atol)
    assert compare_values(ccsd_corl, qcdb.variable('current correlation energy'), tnm() + ' CCSD', atol=atol)
    assert compare_values(ccsd_corl, qcdb.variable('ccsd correlation energy'), tnm() + ' CCSD', atol=atol)

    # mp2 terms (not printed for nwc tce)
    if not (mtd.startswith('nwc') and opts.get('qc_module', 'nein').lower() == 'tce'):
        assert compare_values(mp2_tot, qcdb.variable('mp2 total energy'), tnm() + ' MP2', atol=atol)
        assert compare_values(mp2_corl, qcdb.variable('mp2 correlation energy'), tnm() + ' MP2', atol=atol)

    # hf terms
    assert compare_values(scf_tot, qcdb.variable('hf total energy'), tnm() + ' SCF', atol=atol)
    assert compare_values(scf_tot, qcdb.variable('scf total energy'), tnm() + ' SCF', atol=atol)


@pytest.mark.parametrize('mtd,opts,errmsg', [
    pytest.param('nwc-ccsd', {'basis': 'cfour-qz2p', 'nwchem_scf__uhf': True}, 'ccsd: nopen is not zero', marks=using_nwchem),
    pytest.param('gms-ccsd', {'basis': 'cfour-qz2p', 'gamess_contrl__scftyp': 'uhf', 'gamess_ccinp__ncore': 0}, 'CCTYP IS PROGRAMMED ONLY FOR SCFTYP=RHF OR ROHF', marks=using_gamess),
])
def test_sp_ccsd_uhf_ae_error(mtd, opts, nh2, errmsg):
    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options(opts)

    with pytest.raises(qcng.exceptions.InputError) as e:
        qcdb.energy(mtd, molecule=nh2)

    assert errmsg in str(e.value)


@pytest.mark.parametrize('mtd,opts', [
    pytest.param('c4-ccsd', {'cfour_BASIS': 'qz2p', 'cfour_REFerence': 'UHF', 'cfour_occupation': [[3,1,1,0],[3,0,1,0]], 'cfour_SCF_CONV': 12, 'cfour_CC_CONV': 12}, marks=using_cfour),
    pytest.param('c4-ccsd', {'BASIS': 'cfour-qz2p', 'cfour_reference': 'uhf'}, marks=using_cfour),
    pytest.param('nwc-ccsd', {'basis': 'cfour-qz2p', 'qc_module': 'tce', 'nwchem_scf__uhf': True}, marks=using_nwchem),
    pytest.param('p4-ccsd', {'basis': 'cfour-qz2p', 'reference': 'uhf'}, marks=using_psi4),
])
def test_sp_ccsd_uhf_ae(mtd, opts, nh2):
    """cfour/sp-uhf-ccsd/input.dat
    #! single-point CCSD/qz2p on NH2

    """
    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options(opts)

    e = qcdb.energy(mtd, molecule=nh2)

    scf_tot = -55.5893469688
    mp2_tot = -55.784877360093
    ccsd_tot = -55.802645023972

    mp2_corl = mp2_tot - scf_tot
    ccsd_corl = ccsd_tot - scf_tot

    atol = 1.e-6

    # cc terms
    assert compare_values(ccsd_tot, e, tnm() + ' Returned', atol=atol)
    assert compare_values(ccsd_tot, qcdb.variable('current energy'), tnm() + ' Current', atol=atol)
    assert compare_values(ccsd_tot, qcdb.variable('ccsd total energy'), tnm() + ' CCSD', atol=atol)
    assert compare_values(ccsd_corl, qcdb.variable('current correlation energy'), tnm() + ' CCSD', atol=atol)
    assert compare_values(ccsd_corl, qcdb.variable('ccsd correlation energy'), tnm() + ' CCSD', atol=atol)

    # mp2 terms (not printed for nwc tce)
    if not (mtd.startswith('nwc') and opts.get('qc_module', 'nein').lower() == 'tce'):
        assert compare_values(mp2_tot, qcdb.variable('mp2 total energy'), tnm() + ' MP2', atol=atol)
        assert compare_values(mp2_corl, qcdb.variable('mp2 correlation energy'), tnm() + ' MP2', atol=atol)

    # hf terms
    assert compare_values(scf_tot, qcdb.variable('hf total energy'), tnm() + ' SCF', atol=atol)
    assert compare_values(scf_tot, qcdb.variable('scf total energy'), tnm() + ' SCF', atol=atol)


@pytest.mark.parametrize('mtd,opts,errmsg', [
    pytest.param('nwc-ccsd', {'basis': 'cfour-qz2p', 'nwchem_ccsd__freeze': 1, 'nwchem_scf__uhf': True}, 'ccsd: nopen is not zero', marks=using_nwchem),
    pytest.param('gms-ccsd', {'basis': 'cfour-qz2p', 'gamess_contrl__scftyp': 'uhf'}, 'CCTYP IS PROGRAMMED ONLY FOR SCFTYP=RHF OR ROHF', marks=using_gamess),
])
def test_sp_ccsd_uhf_fc_error(mtd, opts, nh2, errmsg):
    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options(opts)

    with pytest.raises(qcng.exceptions.InputError) as e:
        qcdb.energy(mtd, molecule=nh2)

    assert errmsg in str(e.value)


@pytest.mark.parametrize('mtd,opts', [
    pytest.param('c4-ccsd', {'cfour_BASIS': 'qz2p', 'cfour_dropmo': [1], 'cfour_REFerence': 'UHF', 'cfour_occupation': [[3,1,1,0],[3,0,1,0]], 'cfour_SCF_CONV': 12, 'cfour_CC_CONV': 12}, marks=using_cfour),
    pytest.param('c4-ccsd', {'BASIS': 'cfour-qz2p', 'cfour_dropmo': 1, 'cfour_reference': 'uhf'}, marks=using_cfour),
    pytest.param('nwc-ccsd', {'basis': 'cfour-qz2p', 'nwchem_tce__freeze': 1, 'qc_module': 'tce', 'nwchem_scf__uhf': True}, marks=using_nwchem),
    pytest.param('p4-ccsd', {'basis': 'cfour-qz2p', 'psi4_freeze_core': True, 'reference': 'uhf'}, marks=using_psi4),
])
def test_sp_ccsd_uhf_fc(mtd, opts, nh2):
    """cfour/sp-uhf-ccsd/input.dat
    #! single-point CCSD/qz2p on NH2

    """
    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options(opts)

    e = qcdb.energy(mtd, molecule=nh2)

    scf_tot = -55.5893469688
    mp2_tot = -55.760531091893
    ccsd_tot = -55.777664191533

    mp2_corl = mp2_tot - scf_tot
    ccsd_corl = ccsd_tot - scf_tot

    atol = 1.e-6

    # cc terms
    assert compare_values(ccsd_tot, e, tnm() + ' Returned', atol=atol)
    assert compare_values(ccsd_tot, qcdb.variable('current energy'), tnm() + ' Current', atol=atol)
    assert compare_values(ccsd_tot, qcdb.variable('ccsd total energy'), tnm() + ' CCSD', atol=atol)
    assert compare_values(ccsd_corl, qcdb.variable('current correlation energy'), tnm() + ' CCSD', atol=atol)
    assert compare_values(ccsd_corl, qcdb.variable('ccsd correlation energy'), tnm() + ' CCSD', atol=atol)

    # mp2 terms (not printed for nwc tce)
    if not (mtd.startswith('nwc') and opts.get('qc_module', 'nein').lower() == 'tce'):
        assert compare_values(mp2_tot, qcdb.variable('mp2 total energy'), tnm() + ' MP2', atol=atol)
        assert compare_values(mp2_corl, qcdb.variable('mp2 correlation energy'), tnm() + ' MP2', atol=atol)

    # hf terms
    assert compare_values(scf_tot, qcdb.variable('hf total energy'), tnm() + ' SCF', atol=atol)
    assert compare_values(scf_tot, qcdb.variable('scf total energy'), tnm() + ' SCF', atol=atol)


@pytest.mark.parametrize('mtd,opts,errmsg', [
    pytest.param('nwc-ccsd', {'basis': 'cfour-qz2p', 'nwchem_scf__rohf': True}, 'ccsd: nopen is not zero', marks=using_nwchem),
])
def test_sp_ccsd_rohf_ae_error(mtd, opts, nh2, errmsg):
    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options(opts)

    with pytest.raises(qcng.exceptions.InputError) as e:
        qcdb.energy(mtd, molecule=nh2)

    assert errmsg in str(e.value)


@pytest.mark.parametrize('mtd,opts', [
    pytest.param('c4-ccsd', {'cfour_BASIS': 'qz2p', 'cfour_REFerence': 'roHF', 'cfour_occupation': [[3,1,1,0],[3,0,1,0]], 'cfour_SCF_CONV': 12, 'cfour_CC_CONV': 12}, marks=using_cfour),
    pytest.param('c4-ccsd', {'BASIS': 'cfour-qz2p', 'cfour_reference': 'rohf'}, marks=using_cfour),
    pytest.param('gms-ccsd', {'basis': 'cfour-qz2p', 'gamess_contrl__scftyp': 'rohf', 'gamess_ccinp__ncore': 0}, marks=using_gamess),
    pytest.param('nwc-ccsd', {'basis': 'cfour-qz2p', 'qc_module': 'tce', 'nwchem_scf__rohf': True}, marks=using_nwchem),
    pytest.param('p4-ccsd', {'basis': 'cfour-qz2p', 'reference': 'rohf'}, marks=using_psi4),
])
def test_sp_ccsd_rohf_ae(mtd, opts, nh2):
    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options(opts)

    e, jrec = qcdb.energy(mtd, return_wfn = True, molecule=nh2)
                         
    scf_tot = -55.5847372601
    mp2_tot = -55.7852767873
    ccsd_tot = -55.802586766392

    mp2_corl = mp2_tot - scf_tot
    ccsd_corl = ccsd_tot - scf_tot

    atol = 1.e-6

    # cc terms
    assert compare_values(ccsd_tot, e, tnm() + ' Returned', atol=atol)
    assert compare_values(ccsd_tot, qcdb.variable('current energy'), tnm() + ' Current', atol=atol)
    assert compare_values(ccsd_tot, qcdb.variable('ccsd total energy'), tnm() + ' CCSD', atol=atol)
    assert compare_values(ccsd_corl, qcdb.variable('current correlation energy'), tnm() + ' CCSD', atol=atol)
    assert compare_values(ccsd_corl, qcdb.variable('ccsd correlation energy'), tnm() + ' CCSD', atol=atol)

    # mp2 terms (only printed for c4)
    #if not (mtd.startswith('nwc') and opts.get('qc_module', 'nein').lower() == 'tce'):
    #    assert compare_values(mp2_tot, qcdb.variable('mp2 total energy'), tnm() + ' MP2', atol=atol)
    #    assert compare_values(mp2_corl, qcdb.variable('mp2 correlation energy'), tnm() + ' MP2', atol=atol)

    # hf terms
    assert compare_values(scf_tot, qcdb.variable('hf total energy'), tnm() + ' SCF', atol=atol)
    assert compare_values(scf_tot, qcdb.variable('scf total energy'), tnm() + ' SCF', atol=atol)

#    cfour isn't splitting out the singles from OS. and psi4 isn't incl the singles in either so SS + OS != corl
#    and maybe singles moved from SS to OS in Cfour btwn 2010 and 2014 versions (change in ref)
#    if not (mtd in ['gms-ccsd', 'nwc-ccsd']):
#        assert compare_values(osccsdcorl, qcdb.variable('ccsd opposite-spin correlation energy'), tnm() + ' CCSD OS corl', atol=atol)
#        assert compare_values(ssccsdcorl, qcdb.variable('ccsd same-spin correlation energy'), tnm() + ' CCSD SS corl', atol=atol)


@pytest.mark.parametrize('mtd,opts,errmsg', [
    pytest.param('nwc-ccsd', {'basis': 'cfour-qz2p', 'nwchem_ccsd__freeze': 1, 'nwchem_scf__rohf': True}, 'ccsd: nopen is not zero', marks=using_nwchem), #, pytest.mark.xfail(True, reason='NWChem has no hand-coded ROHF CC. Sensible error msg.', run=True)]),
])
def test_sp_ccsd_rohf_fc_error(mtd, opts, nh2, errmsg):
    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options(opts)

    with pytest.raises(qcng.exceptions.InputError) as e:
        qcdb.energy(mtd, molecule=nh2)

    assert errmsg in str(e.value)
    # check input file, too?


@pytest.mark.parametrize('mtd,opts', [
    pytest.param('c4-ccsd', {'cfour_BASIS': 'qz2p', 'cfour_dropmo': [1], 'cfour_print': 2, 'cfour_REFerence': 'roHF', 'cfour_occupation': [[3,1,1,0],[3,0,1,0]], 'cfour_SCF_CONV': 12, 'cfour_CC_CONV': 12, 'cfour_orbitals': 0}, marks=using_cfour),
    pytest.param('c4-ccsd', {'BASIS': 'cfour-qz2p', 'cfour_dropmo': 1, 'cfour_reference': 'rohf', 'cfour_orbitals': 0}, marks=using_cfour),
    pytest.param('gms-ccsd', {'basis': 'cfour-qz2p', 'gamess_contrl__scftyp': 'rohf', 'gamess_ccinp__iconv': 9, 'gamess_scf__conv': 9}, marks=using_gamess),
    pytest.param('nwc-ccsd', {'basis': 'cfour-qz2p', 'nwchem_tce__freeze': 1, 'qc_module': 'tce', 'nwchem_scf__rohf': True}, marks=using_nwchem),
    pytest.param('p4-ccsd', {'basis': 'cfour-qz2p', 'psi4_e_convergence': 8, 'psi4_r_convergence': 7, 'psi4_freeze_core': True, 'reference': 'rohf'}, marks=using_psi4),
])
def test_sp_ccsd_rohf_fc(mtd, opts, nh2):
    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options(opts)

    e, jrec = qcdb.energy(mtd, return_wfn = True, molecule=nh2)

    # from Cfour
    scf_tot = -55.5847372600528120
    mp2_tot = -55.760613403041812
    ccsd_tot = -55.7775634749542241

    mp2_corl = mp2_tot - scf_tot
    ccsd_corl = ccsd_tot - scf_tot

    atol = 1.e-6
    # gms CCSD correlation disagrees with Cfour, Psi4, and NWChem by ~2.e-4
    # hf is in agreement across all programs
    if mtd.startswith('gms'):
        atol = 2.e-5

    # cc terms
    assert compare_values(ccsd_tot, e, tnm() + ' Returned', atol=atol)
    assert compare_values(ccsd_tot, qcdb.variable('current energy'), tnm() + ' Current', atol=atol)
    assert compare_values(ccsd_tot, qcdb.variable('ccsd total energy'), tnm() + ' CCSD', atol=atol)
    assert compare_values(ccsd_corl, qcdb.variable('current correlation energy'), tnm() + ' CCSD', atol=atol)
    assert compare_values(ccsd_corl, qcdb.variable('ccsd correlation energy'), tnm() + ' CCSD', atol=atol)

    # mp2 terms (only printed for c4)
    #if not (mtd.startswith('nwc') and opts.get('qc_module', 'nein').lower() == 'tce'):
    #    assert compare_values(mp2_tot, qcdb.variable('mp2 total energy'), tnm() + ' MP2', atol=atol)
    #    assert compare_values(mp2_corl, qcdb.variable('mp2 correlation energy'), tnm() + ' MP2', atol=atol)

    # hf terms
    assert compare_values(scf_tot, qcdb.variable('hf total energy'), tnm() + ' SCF', atol=atol)
    assert compare_values(scf_tot, qcdb.variable('scf total energy'), tnm() + ' SCF', atol=atol)
