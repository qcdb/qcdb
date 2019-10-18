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


@pytest.mark.parametrize('method,keywords', [
    pytest.param('c4-ccsd(t)', {'cfour_basis': 'qz2p', 'cfour_scf_conv': 12, 'cfour_cc_conv': 12}, marks=using_cfour),
    pytest.param('c4-ccsd(t)', {'basis': 'cfour-qz2p', 'cfour_scf_conv': 12, 'cfour_cc_conv': 12}, marks=using_cfour),
    pytest.param('gms-ccsd(t)', {'basis': 'cfour-qz2p', 'gamess_ccinp__ncore': 0}, marks=using_gamess),
    pytest.param('nwc-ccsd(t)', {'basis': 'cfour-qz2p', 'qc_module': 'tce'}, marks=using_nwchem),
    pytest.param('nwc-ccsd(t)', {'basis': 'cfour-qz2p'}, marks=using_nwchem),
    pytest.param('p4-ccsd(t)', {'basis': 'cfour-qz2p'}, marks=using_psi4),
])
def test_sp_ccsd_t_rhf_ae(method, keywords, h2o):
    """cfour/???/input.dat
    #! single point CCSD(T)/qz2p on water

    """
    h2o = qcdb.set_molecule(h2o)
    qcdb.set_options(keywords)

    e, jrec = qcdb.energy(method, return_wfn=True, molecule=h2o)

    # from cfour
    scf_tot = -76.0627484601
    mp2_tot = -76.3329401286
    ccsd_tot = -76.3384539547
    ccsd_t_tot = -76.3457175511

    mp2_corl = mp2_tot - scf_tot
    ccsd_corl = ccsd_tot - scf_tot
    ccsd_t_corl = ccsd_t_tot - scf_tot
    t_corl = ccsd_t_tot - ccsd_tot

    atol = 1.e-6

    # cc terms
    assert compare_values(ccsd_t_tot, e, tnm() + ' Returned', atol=atol)
    assert compare_values(ccsd_t_tot, qcdb.variable('current energy'), tnm() + ' Current', atol=atol)
    assert compare_values(ccsd_t_tot, qcdb.variable('ccsd(t) total energy'), tnm() + ' CCSD(T)', atol=atol)

    assert compare_values(ccsd_t_corl, qcdb.variable('current correlation energy'), tnm() + ' CCSD(T)', atol=atol)
    assert compare_values(ccsd_t_corl, qcdb.variable('ccsd(t) correlation energy'), tnm() + ' CCSD(T)', atol=atol)

    assert compare_values(ccsd_tot, qcdb.variable('ccsd total energy'), tnm() + ' CCSD', atol=atol)
    assert compare_values(ccsd_corl, qcdb.variable('ccsd correlation energy'), tnm() + ' CCSD', atol=atol)

    assert compare_values(t_corl, qcdb.variable('(t) correction energy'), tnm() + ' (T)', atol=atol)

    # mp2 terms (not printed for nwc tce)
    if not (method.startswith('nwc') and keywords.get('qc_module', 'nein').lower() == 'tce'):
        assert compare_values(mp2_tot, qcdb.variable('mp2 total energy'), tnm() + ' MP2', atol=atol)
        assert compare_values(mp2_corl, qcdb.variable('mp2 correlation energy'), tnm() + ' MP2', atol=atol)

    # hf terms
    assert compare_values(scf_tot, qcdb.variable('hf total energy'), tnm() + ' SCF', atol=atol)
    assert compare_values(scf_tot, qcdb.variable('scf total energy'), tnm() + ' SCF', atol=atol)


@pytest.mark.parametrize('method,keywords', [
    pytest.param('c4-ccsd(t)', {'cfour_basis': 'qz2p', 'cfour_dropmo': [1], 'cfour_scf_conv': 12, 'cfour_cc_conv': 12}, marks=using_cfour),
    pytest.param('c4-ccsd(t)', {'basis': 'cfour-qz2p', 'cfour_dropmo': 1, 'cfour_scf_conv': 12, 'cfour_cc_conv': 12}, marks=using_cfour),
    pytest.param('gms-ccsd(t)', {'basis': 'cfour-qz2p'}, marks=using_gamess),
    pytest.param('nwc-ccsd(t)', {'basis': 'cfour-qz2p', 'qc_module': 'tce', 'nwchem_tce__freeze': 1}, marks=using_nwchem),
    pytest.param('nwc-ccsd(t)', {'basis': 'cfour-qz2p', 'nwchem_ccsd__freeze': 1}, marks=using_nwchem),
    pytest.param('p4-ccsd(t)', {'basis': 'cfour-qz2p', 'psi4_freeze_core': True}, marks=using_psi4),
])
def test_sp_ccsd_t_rhf_fc(method, keywords, h2o):
    """cfour/???/input.dat
    #! single point CCSD(T)/qz2p on water

    """
    h2o = qcdb.set_molecule(h2o)
    qcdb.set_options(keywords)

    e, jrec = qcdb.energy(method, return_wfn=True, molecule=h2o)

    # from cfour
    scf_tot = -76.0627484747
    mp2_tot = -76.3079003122
    ccsd_tot = -76.3130790236
    ccsd_t_tot = -76.3201755322

    mp2_corl = mp2_tot - scf_tot
    ccsd_corl = ccsd_tot - scf_tot
    ccsd_t_corl = ccsd_t_tot - scf_tot
    t_corl = ccsd_t_tot - ccsd_tot

    atol = 1.e-6

    # cc terms
    assert compare_values(ccsd_t_tot, e, tnm() + ' Returned', atol=atol)
    assert compare_values(ccsd_t_tot, qcdb.variable('current energy'), tnm() + ' Current', atol=atol)
    assert compare_values(ccsd_t_tot, qcdb.variable('ccsd(t) total energy'), tnm() + ' CCSD(T)', atol=atol)

    assert compare_values(ccsd_t_corl, qcdb.variable('current correlation energy'), tnm() + ' CCSD(T)', atol=atol)
    assert compare_values(ccsd_t_corl, qcdb.variable('ccsd(t) correlation energy'), tnm() + ' CCSD(T)', atol=atol)

    assert compare_values(ccsd_tot, qcdb.variable('ccsd total energy'), tnm() + ' CCSD', atol=atol)
    assert compare_values(ccsd_corl, qcdb.variable('ccsd correlation energy'), tnm() + ' CCSD', atol=atol)

    assert compare_values(t_corl, qcdb.variable('(t) correction energy'), tnm() + ' (T)', atol=atol)

    # mp2 terms (not printed for nwc tce)
    if not (method.startswith('nwc') and keywords.get('qc_module', 'nein').lower() == 'tce'):
        assert compare_values(mp2_tot, qcdb.variable('mp2 total energy'), tnm() + ' MP2', atol=atol)
        assert compare_values(mp2_corl, qcdb.variable('mp2 correlation energy'), tnm() + ' MP2', atol=atol)

    # hf terms
    assert compare_values(scf_tot, qcdb.variable('hf total energy'), tnm() + ' SCF', atol=atol)
    assert compare_values(scf_tot, qcdb.variable('scf total energy'), tnm() + ' SCF', atol=atol)


@pytest.mark.parametrize('method,keywords,errmsg', [
    pytest.param('gms-ccsd(t)', {'basis': 'cfour-qz2p', 'gamess_contrl__scftyp': 'uhf', 'gamess_ccinp__ncore': 0}, 'CCTYP IS PROGRAMMED ONLY FOR SCFTYP=RHF OR ROHF', marks=using_gamess),
    pytest.param('nwc-ccsd(t)', {'basis': 'cfour-qz2p', 'nwchem_scf__uhf': True}, 'ccsd: nopen is not zero', marks=using_nwchem),
])
def test_sp_ccsd_t_uhf_ae_error(method, keywords, nh2, errmsg):
    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options(keywords)

    with pytest.raises(qcng.exceptions.InputError) as e:
        qcdb.energy(method, molecule=nh2)

    assert errmsg in str(e.value)


@pytest.mark.parametrize('method,keywords', [
    pytest.param('c4-ccsd(t)', {'cfour_basis': 'qz2p', 'cfour_reference': 'uhf', 'cfour_scf_conv': 12, 'cfour_cc_conv': 12}, marks=using_cfour),
    pytest.param('c4-ccsd(t)', {'basis': 'cfour-qz2p', 'cfour_reference': 'uhf', 'cfour_scf_conv': 12, 'cfour_cc_conv': 12}, marks=using_cfour),
    pytest.param('nwc-ccsd(t)', {'basis': 'cfour-qz2p', 'qc_module': 'tce', 'nwchem_scf__uhf': True}, marks=using_nwchem),
    pytest.param('p4-ccsd(t)', {'basis': 'cfour-qz2p', 'reference': 'uhf'}, marks=using_psi4),
])
def test_sp_ccsd_t_uhf_ae(method, keywords, nh2):
    """cfour/???/input.dat
    #! single point CCSD(T)/qz2p on water

    """
    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options(keywords)

    e, jrec = qcdb.energy(method, return_wfn=True, molecule=nh2)

    # from cfour
    scf_tot = -55.5893469682
    mp2_tot = -55.7848773537
    ccsd_tot = -55.8026450157
    ccsd_t_tot = -55.8078116024

    mp2_corl = mp2_tot - scf_tot
    ccsd_corl = ccsd_tot - scf_tot
    ccsd_t_corl = ccsd_t_tot - scf_tot
    t_corl = ccsd_t_tot - ccsd_tot

    atol = 1.e-6

    # cc terms
    assert compare_values(ccsd_t_tot, e, tnm() + ' Returned', atol=atol)
    assert compare_values(ccsd_t_tot, qcdb.variable('current energy'), tnm() + ' Current', atol=atol)
    assert compare_values(ccsd_t_tot, qcdb.variable('ccsd(t) total energy'), tnm() + ' CCSD(T)', atol=atol)

    assert compare_values(ccsd_t_corl, qcdb.variable('current correlation energy'), tnm() + ' CCSD(T)', atol=atol)
    assert compare_values(ccsd_t_corl, qcdb.variable('ccsd(t) correlation energy'), tnm() + ' CCSD(T)', atol=atol)

    assert compare_values(ccsd_tot, qcdb.variable('ccsd total energy'), tnm() + ' CCSD', atol=atol)
    assert compare_values(ccsd_corl, qcdb.variable('ccsd correlation energy'), tnm() + ' CCSD', atol=atol)

    assert compare_values(t_corl, qcdb.variable('(t) correction energy'), tnm() + ' (T)', atol=atol)

    # mp2 terms (not printed for nwc tce)
    if not (method.startswith('nwc') and keywords.get('qc_module', 'nein').lower() == 'tce'):
        assert compare_values(mp2_tot, qcdb.variable('mp2 total energy'), tnm() + ' MP2', atol=atol)
        assert compare_values(mp2_corl, qcdb.variable('mp2 correlation energy'), tnm() + ' MP2', atol=atol)

    # hf terms
    assert compare_values(scf_tot, qcdb.variable('hf total energy'), tnm() + ' SCF', atol=atol)
    assert compare_values(scf_tot, qcdb.variable('scf total energy'), tnm() + ' SCF', atol=atol)


@pytest.mark.parametrize('method,keywords,errmsg', [
    pytest.param('gms-ccsd(t)', {'basis': 'cfour-qz2p', 'gamess_contrl__scftyp': 'uhf'}, 'CCTYP IS PROGRAMMED ONLY FOR SCFTYP=RHF OR ROHF', marks=using_gamess),
    pytest.param('nwc-ccsd(t)', {'basis': 'cfour-qz2p', 'nwchem_scf__uhf': True, 'nwchem_ccsd__freeze': 1}, 'ccsd: nopen is not zero', marks=using_nwchem),
])
def test_sp_ccsd_t_uhf_fc_error(method, keywords, nh2, errmsg):
    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options(keywords)

    with pytest.raises(qcng.exceptions.InputError) as e:
        qcdb.energy(method, molecule=nh2)

    assert errmsg in str(e.value)


@pytest.mark.parametrize('method,keywords', [
    pytest.param('c4-ccsd(t)', {'cfour_basis': 'qz2p', 'cfour_reference': 'uhf', 'cfour_dropmo': [1], 'cfour_scf_conv': 12, 'cfour_cc_conv': 12}, marks=using_cfour),
    pytest.param('c4-ccsd(t)', {'basis': 'cfour-qz2p', 'cfour_reference': 'uhf', 'cfour_dropmo': 1, 'cfour_scf_conv': 12, 'cfour_cc_conv': 12}, marks=using_cfour),
    pytest.param('nwc-ccsd(t)', {'basis': 'cfour-qz2p', 'qc_module': 'tce', 'nwchem_tce__freeze': 1, 'nwchem_scf__uhf': True}, marks=using_nwchem),
    pytest.param('p4-ccsd(t)', {'basis': 'cfour-qz2p', 'reference': 'uhf', 'psi4_freeze_core': True}, marks=using_psi4),
])
def test_sp_ccsd_t_uhf_fc(method, keywords, nh2):
    """cfour/???/input.dat
    #! single point CCSD(T)/qz2p on water

    """
    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options(keywords)

    e, jrec = qcdb.energy(method, return_wfn=True, molecule=nh2)

    # from cfour
    scf_tot = -55.5893469681738990
    mp2_tot = -55.760531091806899
    ccsd_tot = -55.7776641914257188
    ccsd_t_tot = -55.7826468408968990


    mp2_corl = mp2_tot - scf_tot
    ccsd_corl = ccsd_tot - scf_tot
    ccsd_t_corl = ccsd_t_tot - scf_tot
    t_corl = ccsd_t_tot - ccsd_tot

    atol = 1.e-6

    # cc terms
    assert compare_values(ccsd_t_tot, e, tnm() + ' Returned', atol=atol)
    assert compare_values(ccsd_t_tot, qcdb.variable('current energy'), tnm() + ' Current', atol=atol)
    assert compare_values(ccsd_t_tot, qcdb.variable('ccsd(t) total energy'), tnm() + ' CCSD(T)', atol=atol)

    assert compare_values(ccsd_t_corl, qcdb.variable('current correlation energy'), tnm() + ' CCSD(T)', atol=atol)
    assert compare_values(ccsd_t_corl, qcdb.variable('ccsd(t) correlation energy'), tnm() + ' CCSD(T)', atol=atol)

    assert compare_values(ccsd_tot, qcdb.variable('ccsd total energy'), tnm() + ' CCSD', atol=atol)
    assert compare_values(ccsd_corl, qcdb.variable('ccsd correlation energy'), tnm() + ' CCSD', atol=atol)

    assert compare_values(t_corl, qcdb.variable('(t) correction energy'), tnm() + ' (T)', atol=atol)

    # mp2 terms (not printed for nwc tce)
    if not (method.startswith('nwc') and keywords.get('qc_module', 'nein').lower() == 'tce'):
        assert compare_values(mp2_tot, qcdb.variable('mp2 total energy'), tnm() + ' MP2', atol=atol)
        assert compare_values(mp2_corl, qcdb.variable('mp2 correlation energy'), tnm() + ' MP2', atol=atol)

    # hf terms
    assert compare_values(scf_tot, qcdb.variable('hf total energy'), tnm() + ' SCF', atol=atol)
    assert compare_values(scf_tot, qcdb.variable('scf total energy'), tnm() + ' SCF', atol=atol)


@pytest.mark.parametrize('method,keywords,errmsg', [
    pytest.param('gms-ccsd(t)', {'basis': 'cfour-qz2p', 'gamess_contrl__scftyp': 'rohf', 'gamess_ccinp__ncore': 0}, 'ROHF\'S CCTYP MUST BE CCSD OR CR-CCL', marks=using_gamess), # No (T) w/ rohf in gms
    pytest.param('nwc-ccsd(t)', {'basis': 'cfour-qz2p', 'nwchem_scf__rohf': True}, 'ccsd: nopen is not zero', marks=using_nwchem),
])
def test_sp_ccsd_t_rohf_ae_error(method, keywords, errmsg, nh2):
    """cfour/???/input.dat
    #! single point CCSD(T)/qz2p on water

    """
    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options(keywords)

    with pytest.raises(qcng.exceptions.InputError) as e:
        qcdb.energy(method, molecule=nh2)

    assert errmsg in str(e.value)


@pytest.mark.parametrize('method,keywords', [
    pytest.param('c4-ccsd(t)', {'cfour_basis': 'qz2p', 'cfour_reference': 'rohf', 'cfour_occupation': [[3,1,1,0],[3,0,1,0]], 'cfour_scf_conv': 12, 'cfour_cc_conv': 12}, marks=using_cfour),
    pytest.param('c4-ccsd(t)', {'basis': 'cfour-qz2p', 'cfour_reference': 'rohf', 'cfour_occupation': [[3,1,1,0],[3,0,1,0]], 'cfour_scf_conv': 12, 'cfour_cc_conv': 12}, marks=using_cfour),
    pytest.param('nwc-ccsd(t)', {'basis': 'cfour-qz2p', 'qc_module': 'tce', 'nwchem_scf__rohf': True, 'nwchem_scf__thresh': 8, 'nwchem_tce__thresh': 8, 'nwchem_tce__freeze': 0, 'nwchem_scf__tol2e': 10}, marks=using_nwchem),
    pytest.param('p4-ccsd(t)', {'basis': 'cfour-qz2p', 'reference': 'rohf'}, marks=using_psi4),

])
def test_sp_ccsd_t_rohf_ae(method, keywords, nh2):
    """cfour/???/input.dat
    #! single point CCSD(T)/qz2p on water

    """
    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options(keywords)

    e, jrec = qcdb.energy(method, return_wfn=True, molecule=nh2)

    # from Cfour
    scf_tot = -55.5847372601
    mp2_tot = -55.7852767809
    ccsd_tot = -55.8025867548
    ccsd_t_tot = -55.8078206933

    mp2_corl = mp2_tot - scf_tot
    ccsd_corl = ccsd_tot - scf_tot
    ccsd_t_corl = ccsd_t_tot - scf_tot
    t_corl = ccsd_t_tot - ccsd_tot

    atol = 1.e-6

    # nwc CCSD(T) correlation disagrees with Cfour and Psi4 by ~2.e-4
    # hf and ccsd correlation are in agreement
    if method.startswith('nwc'):
        atol = 2.e-4

    # cc terms
    assert compare_values(ccsd_t_tot, e, tnm() + ' Returned', atol=atol)
    assert compare_values(ccsd_t_tot, qcdb.variable('current energy'), tnm() + ' Current', atol=atol)
    assert compare_values(ccsd_t_tot, qcdb.variable('ccsd(t) total energy'), tnm() + ' CCSD(T)', atol=atol)

    assert compare_values(ccsd_t_corl, qcdb.variable('current correlation energy'), tnm() + ' CCSD(T)', atol=atol)
    assert compare_values(ccsd_t_corl, qcdb.variable('ccsd(t) correlation energy'), tnm() + ' CCSD(T)', atol=atol)

    assert compare_values(ccsd_tot, qcdb.variable('ccsd total energy'), tnm() + ' CCSD', atol=atol)
    assert compare_values(ccsd_corl, qcdb.variable('ccsd correlation energy'), tnm() + ' CCSD', atol=atol)

    assert compare_values(t_corl, qcdb.variable('(t) correction energy'), tnm() + ' (T)', atol=atol)

    # mp2 terms (not printed for nwc tce)
    if not (method.startswith('nwc') and keywords.get('qc_module', 'nein').lower() == 'tce'):
        assert compare_values(mp2_tot, qcdb.variable('mp2 total energy'), tnm() + ' MP2', atol=atol)
        assert compare_values(mp2_corl, qcdb.variable('mp2 correlation energy'), tnm() + ' MP2', atol=atol)

    # hf terms
    assert compare_values(scf_tot, qcdb.variable('hf total energy'), tnm() + ' SCF', atol=atol)
    assert compare_values(scf_tot, qcdb.variable('scf total energy'), tnm() + ' SCF', atol=atol)


@pytest.mark.parametrize('method,keywords,errmsg', [
    pytest.param('gms-ccsd(t)', {'basis': 'cfour-qz2p', 'gamess_contrl__scftyp': 'rohf'}, 'ROHF\'S CCTYP MUST BE CCSD OR CR-CCL', marks=using_gamess), # No (T) w/ rohf in gms
    pytest.param('nwc-ccsd(t)', {'basis': 'cfour-qz2p', 'nwchem_scf__rohf': True, 'nwchem_ccsd__freeze': 1}, 'ccsd: nopen is not zero', marks=using_nwchem),
])
def test_sp_ccsd_t_rohf_fc_error(method, keywords, errmsg, nh2):
    """cfour/???/input.dat
    #! single point CCSD(T)/qz2p on water

    """
    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options(keywords)

    with pytest.raises(qcng.exceptions.InputError) as e:
        qcdb.energy(method, molecule=nh2)

    assert errmsg in str(e.value)


@pytest.mark.parametrize('method,keywords', [
    pytest.param('c4-ccsd(t)', {'cfour_basis': 'qz2p', 'cfour_reference': 'rohf', 'cfour_dropmo': [1], 'cfour_scf_conv': 12, 'cfour_cc_conv': 12}, marks=using_cfour),
    pytest.param('c4-ccsd(t)', {'basis': 'cfour-qz2p', 'cfour_reference': 'rohf', 'cfour_dropmo': 1, 'cfour_scf_conv': 12, 'cfour_cc_conv': 12}, marks=using_cfour),
    pytest.param('nwc-ccsd(t)', {'basis': 'cfour-qz2p', 'qc_module': 'tce', 'nwchem_tce__freeze': 1, 'nwchem_scf__rohf': True}, marks=using_nwchem),
    pytest.param('p4-ccsd(t)', {'basis': 'cfour-qz2p', 'reference': 'rohf', 'psi4_freeze_core': True}, marks=using_psi4),
])
def test_sp_ccsd_t_rohf_fc(method, keywords, nh2):
    """cfour/???/input.dat
    #! single point CCSD(T)/qz2p on water

    """
    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options(keywords)

    e, jrec = qcdb.energy(method, return_wfn=True, molecule=nh2)

    # from Cfour
    scf_tot = -55.58473726
    mp2_tot = -55.76085357
    ccsd_tot = -55.77756954
    ccsd_t_tot = -55.78261742

    ## from NWChem
    #scf_tot = -55.58473726
    #mp2_tot =   ???
    #ccsd_tot = -55.77756347
    #ccsd_t_tot = -55.78275606

    ## from Psi4
    #scf_tot = -55.58473726
    #mp2_tot = -55.76085357
    #ccsd_tot = -55.77756954
    #ccsd_t_tot = -55.78261742

    mp2_corl = mp2_tot - scf_tot
    ccsd_corl = ccsd_tot - scf_tot
    ccsd_t_corl = ccsd_t_tot - scf_tot
    t_corl = ccsd_t_tot - ccsd_tot

    atol = 1.e-6

    # nwc CCSD(T) correlation disagrees with Cfour and Psi4 by ~2.e-4
    # nwc CCSD correlation disagrees with Cfour and Psi4 by ~7.e-6
    # hf is in agreement
    if method.startswith('nwc'):
        atol = 2.e-4

    # cc terms
    assert compare_values(ccsd_t_tot, e, tnm() + ' Returned', atol=atol)
    assert compare_values(ccsd_t_tot, qcdb.variable('current energy'), tnm() + ' Current', atol=atol)
    assert compare_values(ccsd_t_tot, qcdb.variable('ccsd(t) total energy'), tnm() + ' CCSD(T)', atol=atol)

    assert compare_values(ccsd_t_corl, qcdb.variable('current correlation energy'), tnm() + ' CCSD(T)', atol=atol)
    assert compare_values(ccsd_t_corl, qcdb.variable('ccsd(t) correlation energy'), tnm() + ' CCSD(T)', atol=atol)

    assert compare_values(ccsd_tot, qcdb.variable('ccsd total energy'), tnm() + ' CCSD', atol=atol)
    assert compare_values(ccsd_corl, qcdb.variable('ccsd correlation energy'), tnm() + ' CCSD', atol=atol)

    assert compare_values(t_corl, qcdb.variable('(t) correction energy'), tnm() + ' (T)', atol=atol)

    # mp2 terms (not printed for nwc tce)
    if not (method.startswith('nwc') and keywords.get('qc_module', 'nein').lower() == 'tce'):
        assert compare_values(mp2_tot, qcdb.variable('mp2 total energy'), tnm() + ' MP2', atol=atol)
        assert compare_values(mp2_corl, qcdb.variable('mp2 correlation energy'), tnm() + ' MP2', atol=atol)

    # hf terms
    assert compare_values(scf_tot, qcdb.variable('hf total energy'), tnm() + ' SCF', atol=atol)
    assert compare_values(scf_tot, qcdb.variable('scf total energy'), tnm() + ' SCF', atol=atol)
