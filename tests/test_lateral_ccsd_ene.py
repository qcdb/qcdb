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

@pytest.mark.parametrize('basis', [
    'cfour-qz2p',
    'aug-cc-pvdz',
])
@pytest.mark.parametrize('method,keywords', [
    pytest.param('c4-ccsd', {'cfour_SCF_CONV': 12, 'cfour_CC_CONV': 12}, marks=using_cfour),
    pytest.param('c4-ccsd', {}, marks=using_cfour),
    pytest.param('nwc-ccsd', {}, marks=using_nwchem),
    pytest.param('nwc-ccsd', {'qc_module': 'tce'}, marks=using_nwchem),
    pytest.param('p4-ccsd', {}, marks=using_psi4),
    pytest.param('gms-ccsd', {'gamess_ccinp__ncore': 0}, marks=using_gamess),
])
def test_sp_ccsd_rhf_ae(method, basis, keywords, h2o):
    """cfour/sp-rhf-ccsd/input.dat
    #! single point CCSD/qz2p on water

    """
    h2o = qcdb.set_molecule(h2o)
    keywords['basis'] = basis
    qcdb.set_options(keywords)

    e, jrec = qcdb.energy(method, return_wfn=True, molecule=h2o)

    if basis == 'cfour-qz2p':
        scf_tot = -76.062748460117
        mp2_tot = -76.332940127333
        #ccsd_corl = -0.275705491773
        ccsd_tot = -76.338453951890
    elif basis == 'aug-cc-pvdz':
        scf_tot = -76.0413815332
        mp2_tot = -76.2632792578
        #mp2_corl = -0.2218977246
        #ccsd_corl = -0.2294105794
        ccsd_tot = -76.2707921127

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
    if not (method.startswith('nwc') and keywords.get('qc_module', 'nein').lower() == 'tce'):
        assert compare_values(mp2_tot, qcdb.variable('mp2 total energy'), tnm() + ' MP2', atol=atol)
        assert compare_values(mp2_corl, qcdb.variable('mp2 correlation energy'), tnm() + ' MP2', atol=atol)

    # hf terms
    assert compare_values(scf_tot, qcdb.variable('hf total energy'), tnm() + ' SCF', atol=atol)
    assert compare_values(scf_tot, qcdb.variable('scf total energy'), tnm() + ' SCF', atol=atol)


@pytest.mark.parametrize('basis', [
    'cfour-qz2p',
    'aug-cc-pvdz',
])
@pytest.mark.parametrize('method,keywords', [
    pytest.param('c4-ccsd', {'cfour_dropmo': [1], 'cfour_SCF_CONV': 12, 'cfour_CC_CONV': 12}, marks=using_cfour),
    pytest.param('c4-ccsd', {'cfour_dropmo': 1}, marks=using_cfour),
    pytest.param('nwc-ccsd', {'nwchem_ccsd__freeze': 1}, marks=using_nwchem),
    pytest.param('nwc-ccsd', {'qc_module': 'tce', 'nwchem_tce__freeze': 1}, marks=using_nwchem),
    pytest.param('p4-ccsd', {'psi4_freeze_core': True}, marks=using_psi4),
    pytest.param('gms-ccsd', {}, marks=using_gamess),
])
def test_sp_ccsd_rhf_fc(method, keywords, basis, h2o):
    """cfour/sp-rhf-ccsd/input.dat
    #! single point CCSD/qz2p on water

    """
    h2o = qcdb.set_molecule(h2o)
    keywords['basis'] = basis
    qcdb.set_options(keywords)

    e, jrec = qcdb.energy(method, return_wfn=True, molecule=h2o)

    if basis == 'cfour-qz2p':
        scf_tot = -76.062748460117
        #mp2_corl = -0.245151837406
        mp2_tot = -76.307900312177
        #ccsd_corl = -0.250330548844
        ccsd_tot = -76.313079023615
    elif basis == 'aug-cc-pvdz':
        scf_tot = -76.0413815332
        #mp2_corl = -0.2194081478
        mp2_tot = -76.2607896810
        #ccsd_corl = -0.2271733460
        ccsd_tot = -76.2685548793

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
    if not (method.startswith('nwc') and keywords.get('qc_module', 'nein').lower() == 'tce'):
        assert compare_values(mp2_tot, qcdb.variable('mp2 total energy'), tnm() + ' MP2', atol=atol)
        assert compare_values(mp2_corl, qcdb.variable('mp2 correlation energy'), tnm() + ' MP2', atol=atol)

    # hf terms
    assert compare_values(scf_tot, qcdb.variable('hf total energy'), tnm() + ' SCF', atol=atol)
    assert compare_values(scf_tot, qcdb.variable('scf total energy'), tnm() + ' SCF', atol=atol)


@pytest.mark.parametrize('method,keywords,errmsg', [
    pytest.param('nwc-ccsd', {'basis': 'cfour-qz2p', 'nwchem_scf__uhf': True}, 'ccsd: nopen is not zero', marks=using_nwchem),
    pytest.param('gms-ccsd', {'basis': 'cfour-qz2p', 'gamess_contrl__scftyp': 'uhf', 'gamess_ccinp__ncore': 0}, 'CCTYP IS PROGRAMMED ONLY FOR SCFTYP=RHF OR ROHF', marks=using_gamess),
])
def test_sp_ccsd_uhf_ae_error(method, keywords, nh2, errmsg):
    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options(keywords)

    with pytest.raises(qcng.exceptions.InputError) as e:
        qcdb.energy(method, molecule=nh2)

    assert errmsg in str(e.value)


@pytest.mark.parametrize('method,keywords', [
    pytest.param('c4-ccsd', {'cfour_BASIS': 'qz2p', 'cfour_REFerence': 'UHF', 'cfour_occupation': [[3,1,1,0],[3,0,1,0]], 'cfour_SCF_CONV': 12, 'cfour_CC_CONV': 12}, marks=using_cfour),
    pytest.param('c4-ccsd', {'BASIS': 'cfour-qz2p', 'cfour_reference': 'uhf'}, marks=using_cfour),
    pytest.param('nwc-ccsd', {'basis': 'cfour-qz2p', 'qc_module': 'tce', 'nwchem_scf__uhf': True}, marks=using_nwchem),
    pytest.param('p4-ccsd', {'basis': 'cfour-qz2p', 'reference': 'uhf'}, marks=using_psi4),
])
def test_sp_ccsd_uhf_ae(method, keywords, nh2):
    """cfour/sp-uhf-ccsd/input.dat
    #! single-point CCSD/qz2p on NH2

    """
    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options(keywords)

    e = qcdb.energy(method, molecule=nh2)

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
    if not (method.startswith('nwc') and keywords.get('qc_module', 'nein').lower() == 'tce'):
        assert compare_values(mp2_tot, qcdb.variable('mp2 total energy'), tnm() + ' MP2', atol=atol)
        assert compare_values(mp2_corl, qcdb.variable('mp2 correlation energy'), tnm() + ' MP2', atol=atol)

    # hf terms
    assert compare_values(scf_tot, qcdb.variable('hf total energy'), tnm() + ' SCF', atol=atol)
    assert compare_values(scf_tot, qcdb.variable('scf total energy'), tnm() + ' SCF', atol=atol)


@pytest.mark.parametrize('method,keywords,errmsg', [
    pytest.param('nwc-ccsd', {'basis': 'cfour-qz2p', 'nwchem_ccsd__freeze': 1, 'nwchem_scf__uhf': True}, 'ccsd: nopen is not zero', marks=using_nwchem),
    pytest.param('gms-ccsd', {'basis': 'cfour-qz2p', 'gamess_contrl__scftyp': 'uhf'}, 'CCTYP IS PROGRAMMED ONLY FOR SCFTYP=RHF OR ROHF', marks=using_gamess),
])
def test_sp_ccsd_uhf_fc_error(method, keywords, nh2, errmsg):
    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options(keywords)

    with pytest.raises(qcng.exceptions.InputError) as e:
        qcdb.energy(method, molecule=nh2)

    assert errmsg in str(e.value)


@pytest.mark.parametrize('method,keywords', [
    pytest.param('c4-ccsd', {'cfour_BASIS': 'qz2p', 'cfour_dropmo': [1], 'cfour_REFerence': 'UHF', 'cfour_occupation': [[3,1,1,0],[3,0,1,0]], 'cfour_SCF_CONV': 12, 'cfour_CC_CONV': 12}, marks=using_cfour),
    pytest.param('c4-ccsd', {'BASIS': 'cfour-qz2p', 'cfour_dropmo': 1, 'cfour_reference': 'uhf'}, marks=using_cfour),
    pytest.param('nwc-ccsd', {'basis': 'cfour-qz2p', 'nwchem_tce__freeze': 1, 'qc_module': 'tce', 'nwchem_scf__uhf': True}, marks=using_nwchem),
    pytest.param('p4-ccsd', {'basis': 'cfour-qz2p', 'psi4_freeze_core': True, 'reference': 'uhf'}, marks=using_psi4),
])
def test_sp_ccsd_uhf_fc(method, keywords, nh2):
    """cfour/sp-uhf-ccsd/input.dat
    #! single-point CCSD/qz2p on NH2

    """
    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options(keywords)

    e = qcdb.energy(method, molecule=nh2)

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
    if not (method.startswith('nwc') and keywords.get('qc_module', 'nein').lower() == 'tce'):
        assert compare_values(mp2_tot, qcdb.variable('mp2 total energy'), tnm() + ' MP2', atol=atol)
        assert compare_values(mp2_corl, qcdb.variable('mp2 correlation energy'), tnm() + ' MP2', atol=atol)

    # hf terms
    assert compare_values(scf_tot, qcdb.variable('hf total energy'), tnm() + ' SCF', atol=atol)
    assert compare_values(scf_tot, qcdb.variable('scf total energy'), tnm() + ' SCF', atol=atol)


@pytest.mark.parametrize('method,keywords,errmsg', [
    pytest.param('nwc-ccsd', {'basis': 'cfour-qz2p', 'nwchem_scf__rohf': True}, 'ccsd: nopen is not zero', marks=using_nwchem),
])
def test_sp_ccsd_rohf_ae_error(method, keywords, nh2, errmsg):
    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options(keywords)

    with pytest.raises(qcng.exceptions.InputError) as e:
        qcdb.energy(method, molecule=nh2)

    assert errmsg in str(e.value)


@pytest.mark.parametrize('basis', [
    'cfour-qz2p',
    'aug-cc-pvdz',
])
@pytest.mark.parametrize('method,keywords', [
    pytest.param('c4-ccsd', {'cfour_REFerence': 'roHF', 'cfour_occupation': [[3,1,1,0],[3,0,1,0]], 'cfour_SCF_CONV': 12, 'cfour_CC_CONV': 12}, marks=using_cfour),
    pytest.param('c4-ccsd', {'cfour_reference': 'rohf'}, marks=using_cfour),
    pytest.param('gms-ccsd', {'gamess_contrl__scftyp': 'rohf', 'gamess_ccinp__ncore': 0}, marks=using_gamess),
    pytest.param('nwc-ccsd', {'qc_module': 'tce', 'nwchem_scf__rohf': True}, marks=using_nwchem),
    pytest.param('p4-ccsd', {'reference': 'rohf'}, marks=using_psi4),
])
def test_sp_ccsd_rohf_ae(method, basis, keywords, nh2):
    nh2 = qcdb.set_molecule(nh2)
    keywords['basis'] = basis
    qcdb.set_options(keywords)

    e, jrec = qcdb.energy(method, return_wfn=True, molecule=nh2)

    if basis == 'cfour-qz2p':
        scf_tot = -55.5847372601
        mp2_tot = -55.7852767873
        ssccsd_corl = -0.0398915
        osccsd_corl = -0.1779580
        #osccsd_corl = -0.1745752
        #ssccsd_corl = -0.0432743
        ccsd_corl = -0.217849506326
        ccsd_tot = -55.802586766392
    elif basis == 'aug-cc-pvdz':
        scf_tot = -55.570724348574
        mp2_tot = 0.0  #TODO
        ssccsd_corl = -0.0339827
        osccsd_corl = -0.1442533
        ccsd_corl = -0.178236032911
        ccsd_tot = -55.748960381485

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
    #if not (method.startswith('nwc') and keywords.get('qc_module', 'nein').lower() == 'tce'):
    #    assert compare_values(mp2_tot, qcdb.variable('mp2 total energy'), tnm() + ' MP2', atol=atol)
    #    assert compare_values(mp2_corl, qcdb.variable('mp2 correlation energy'), tnm() + ' MP2', atol=atol)

    # hf terms
    assert compare_values(scf_tot, qcdb.variable('hf total energy'), tnm() + ' SCF', atol=atol)
    assert compare_values(scf_tot, qcdb.variable('scf total energy'), tnm() + ' SCF', atol=atol)

#    cfour isn't splitting out the singles from OS. and psi4 isn't incl the singles in either so SS + OS != corl
#    and maybe singles moved from SS to OS in Cfour btwn 2010 and 2014 versions (change in ref)
#    if not (method in ['gms-ccsd', 'nwc-ccsd']):
#        'cfour_PRINT': 2
#        assert compare_values(osccsd_corl, qcdb.variable('ccsd opposite-spin correlation energy'), tnm() + ' CCSD OS corl', atol=atol)
#        assert compare_values(ssccsd_corl, qcdb.variable('ccsd same-spin correlation energy'), tnm() + ' CCSD SS corl', atol=atol)


@pytest.mark.parametrize('method,keywords,errmsg', [
    pytest.param('nwc-ccsd', {'basis': 'cfour-qz2p', 'nwchem_ccsd__freeze': 1, 'nwchem_scf__rohf': True}, 'ccsd: nopen is not zero', marks=using_nwchem), #, pytest.mark.xfail(True, reason='NWChem has no hand-coded ROHF CC. Sensible error msg.', run=True)]),
])
def test_sp_ccsd_rohf_fc_error(method, keywords, nh2, errmsg):
    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options(keywords)

    with pytest.raises(qcng.exceptions.InputError) as e:
        qcdb.energy(method, molecule=nh2)

    assert errmsg in str(e.value)
    # check input file, too?


@pytest.mark.parametrize('method,keywords', [
    pytest.param('c4-ccsd', {'cfour_BASIS': 'qz2p', 'cfour_dropmo': [1], 'cfour_print': 2, 'cfour_REFerence': 'roHF', 'cfour_occupation': [[3,1,1,0],[3,0,1,0]], 'cfour_SCF_CONV': 12, 'cfour_CC_CONV': 12, 'cfour_orbitals': 0}, marks=using_cfour),
    pytest.param('c4-ccsd', {'BASIS': 'cfour-qz2p', 'cfour_dropmo': 1, 'cfour_reference': 'rohf', 'cfour_orbitals': 0}, marks=using_cfour),
    pytest.param('gms-ccsd', {'basis': 'cfour-qz2p', 'gamess_contrl__scftyp': 'rohf', 'gamess_ccinp__iconv': 9, 'gamess_scf__conv': 9}, marks=using_gamess),
    pytest.param('nwc-ccsd', {'basis': 'cfour-qz2p', 'nwchem_tce__freeze': 1, 'qc_module': 'tce', 'nwchem_scf__rohf': True}, marks=using_nwchem),
    pytest.param('p4-ccsd', {'basis': 'cfour-qz2p', 'psi4_e_convergence': 8, 'psi4_r_convergence': 7, 'psi4_freeze_core': True, 'reference': 'rohf'}, marks=using_psi4),
])
def test_sp_ccsd_rohf_fc(method, keywords, nh2):
    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options(keywords)

    e, jrec = qcdb.energy(method, return_wfn=True, molecule=nh2)

    # differing std/semicanonical
    ## from Cfour
    #osccsdcorl = -0.1563275
    #ssccsdcorl = -0.0365048
    ##ssccsdcorl = -0.036502859383024

    ## from Psi4. psi & nwchem agree to 6
    osccsd_corl = -0.152968752464362
    ssccsd_corl = -0.036502859383024

    ## from gamess.
    #ccsdcorl = -0.1928447371
    #ccsdtot = -55.7775819971
    # from Cfour

    scf_tot = -55.5847372600528120
    mp2_tot = -55.760613403041812
    ccsd_tot = -55.7775634749542241

    mp2_corl = mp2_tot - scf_tot
    ccsd_corl = ccsd_tot - scf_tot

    atol = 1.e-6
    # gms CCSD correlation disagrees with Cfour, Psi4, and NWChem by ~2.e-4
    # hf is in agreement across all programs
    if method.startswith('gms'):
        atol = 2.e-5

    # cc terms
    assert compare_values(ccsd_tot, e, tnm() + ' Returned', atol=atol)
    assert compare_values(ccsd_tot, qcdb.variable('current energy'), tnm() + ' Current', atol=atol)
    assert compare_values(ccsd_tot, qcdb.variable('ccsd total energy'), tnm() + ' CCSD', atol=atol)
    assert compare_values(ccsd_corl, qcdb.variable('current correlation energy'), tnm() + ' CCSD', atol=atol)
    assert compare_values(ccsd_corl, qcdb.variable('ccsd correlation energy'), tnm() + ' CCSD', atol=atol)

    # mp2 terms (only printed for c4)
    #if not (method.startswith('nwc') and keywords.get('qc_module', 'nein').lower() == 'tce'):
    #    assert compare_values(mp2_tot, qcdb.variable('mp2 total energy'), tnm() + ' MP2', atol=atol)
    #    assert compare_values(mp2_corl, qcdb.variable('mp2 correlation energy'), tnm() + ' MP2', atol=atol)

    # hf terms
    assert compare_values(scf_tot, qcdb.variable('hf total energy'), tnm() + ' SCF', atol=atol)
    assert compare_values(scf_tot, qcdb.variable('scf total energy'), tnm() + ' SCF', atol=atol)

    if not (method in ['gms-ccsd', 'nwc-ccsd']):
#        assert compare_values(osccsd_corl, qcdb.variable('ccsd opposite-spin correlation energy'), tnm() + ' CCSD OS corl', atol=atol)
        assert compare_values(ssccsd_corl, qcdb.variable('ccsd same-spin correlation energy'), tnm() + ' CCSD SS corl', atol=atol)
