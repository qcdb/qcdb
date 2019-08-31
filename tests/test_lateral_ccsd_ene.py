import os

import pytest
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
    pytest.param('nwc-ccsd', {'basis': 'cfour-qz2p'}, marks=using_nwchem),
    pytest.param('nwc-ccsd', {'basis': 'cfour-qz2p', 'qc_module': 'tce'}, marks=using_nwchem),
    pytest.param('p4-ccsd', {'basis': 'cfour-qz2p'}, marks=using_psi4),
    pytest.param('gms-ccsd', {'basis': 'cfour-qz2p', 'gamess_ccinp__ncore': 0}, marks=using_gamess),
])
def test_sp_ccsd_rhf_full(mtd, opts, h2o):
    """cfour/sp-rhf-ccsd/input.dat
    #! single point CCSD/qz2p on water

    """
    h2o = qcdb.set_molecule(h2o)
    qcdb.set_options(opts)

    e, jrec = qcdb.energy(mtd, return_wfn=True, molecule=h2o)

    scftot = -76.062748460117
    mp2tot = -76.332940127333
    ccsdcorl = -0.275705491773
    ccsdtot = -76.338453951890

    atol = 1.e-6
    assert compare_values(scftot, qcdb.variable('scf total energy'), tnm() + ' SCF', atol=atol)
    if not (mtd == 'nwc-ccsd' and opts.get('qc_module', 'nein').lower() == 'tce'):
        assert compare_values(mp2tot, qcdb.variable('mp2 total energy'), tnm() + ' MP2', atol=atol)
    assert compare_values(ccsdcorl, qcdb.variable('ccsd correlation energy'), tnm() + ' CCSD corl', atol=atol)
    assert compare_values(ccsdtot, qcdb.variable('ccsd total energy'), tnm() + ' CCSD', atol=atol)


@pytest.mark.parametrize('mtd,opts', [
    pytest.param('c4-ccsd', {'cfour_BASIS': 'qz2p', 'cfour_dropmo': [1], 'cfour_SCF_CONV': 12, 'cfour_CC_CONV': 12}, marks=using_cfour),
    pytest.param('c4-ccsd', {'BASIS': 'cfour-qz2p', 'cfour_dropmo': 1}, marks=using_cfour),
    pytest.param('nwc-ccsd', {'basis': 'cfour-qz2p', 'nwchem_ccsd__freeze': 1}, marks=using_nwchem),
    pytest.param('nwc-ccsd', {'basis': 'cfour-qz2p', 'qc_module': 'tce', 'nwchem_tce__freeze': 1}, marks=using_nwchem),
    pytest.param('p4-ccsd', {'basis': 'cfour-qz2p', 'psi4_freeze_core': True}, marks=using_psi4),
    pytest.param('gms-ccsd', {'basis': 'cfour-qz2p'}, marks=using_gamess),
])
def test_sp_ccsd_rhf_fc(mtd, opts, h2o):
    """cfour/sp-rhf-ccsd/input.dat
    #! single point CCSD/qz2p on water

    """
    h2o = qcdb.set_molecule(h2o)
    qcdb.set_options(opts)

    e, jrec = qcdb.energy(mtd, return_wfn=True, molecule=h2o)

    scftot = -76.062748460117
    mp2corl = -0.245151837406
    mp2tot = -76.307900312177
    ccsdcorl = -0.250330548844
    ccsdtot = -76.313079023615

    atol = 1.e-6
    assert compare_values(scftot, qcdb.variable('scf total energy'), tnm() + ' SCF', atol=atol)
    if not (mtd == 'nwc-ccsd' and opts.get('qc_module', 'nein').lower() == 'tce'):
        assert compare_values(mp2corl, qcdb.variable('mp2 correlation energy'), tnm() + ' MP2 corl', atol=atol)
        assert compare_values(mp2tot, qcdb.variable('mp2 total energy'), tnm() + ' MP2', atol=atol)
    assert compare_values(ccsdcorl, qcdb.variable('ccsd correlation energy'), tnm() + ' CCSD corl', atol=atol)
    assert compare_values(ccsdtot, qcdb.variable('ccsd total energy'), tnm() + ' CCSD', atol=atol)


@pytest.mark.parametrize('mtd,opts,errmsg', [
    pytest.param('nwc-ccsd', {'basis': 'cfour-qz2p', 'nwchem_scf__uhf': True}, 'ccsd: nopen is not zero', marks=using_nwchem),
    pytest.param('gms-ccsd', {'basis': 'cfour-qz2p', 'gamess_contrl__scftyp': 'uhf', 'gamess_ccinp__ncore': 0}, 'CCTYP IS PROGRAMMED ONLY FOR SCFTYP=RHF OR ROHF', marks=using_gamess),
])
def test_sp_ccsd_uhf_full_error(mtd, opts, nh2, errmsg):
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
def test_sp_ccsd_uhf_full(mtd, opts, nh2):
    """cfour/sp-uhf-ccsd/input.dat
    #! single-point CCSD/qz2p on NH2

    """
    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options(opts)

    e = qcdb.energy(mtd, molecule=nh2)

    scftot = -55.5893469688
    mp2tot = -55.784877360093
    ccsdcorl = -0.213298055172

    atol = 1.e-6
    assert compare_values(scftot, qcdb.variable('scf total energy'), tnm() + ' SCF', atol=atol)
    if not (mtd == 'nwc-ccsd' and opts.get('qc_module', 'nein').lower() == 'tce'):
        assert compare_values(mp2tot, qcdb.variable('mp2 total energy'), tnm() + ' MP2', atol=atol)
    assert compare_values(ccsdcorl, qcdb.variable('ccsd correlation energy'), tnm() + ' CCSD', atol=atol)
    assert compare_values(ccsdcorl+scftot, e, atol=atol)


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

    scftot = -55.5893469688
    mp2corl = -0.171184123712
    mp2tot = -55.760531091893
    ccsdcorl = -0.188317223352
    ccsdtot = -55.777664191533

    atol = 1.e-6
    assert compare_values(scftot, qcdb.variable('scf total energy'), tnm() + ' SCF', atol=atol)
    if not (mtd == 'nwc-ccsd' and opts.get('qc_module', 'nein').lower() == 'tce'):
        assert compare_values(mp2tot, qcdb.variable('mp2 total energy'), tnm() + ' MP2', atol=atol)
    assert compare_values(ccsdcorl, qcdb.variable('ccsd correlation energy'), tnm() + ' CCSD', atol=atol)
    assert compare_values(ccsdcorl+scftot, e, atol=atol)


@pytest.mark.parametrize('mtd,opts,errmsg', [
    pytest.param('nwc-ccsd', {'basis': 'cfour-qz2p', 'nwchem_scf__rohf': True}, 'ccsd: nopen is not zero', marks=using_nwchem),
])
def test_sp_ccsd_rohf_full_error(mtd, opts, nh2, errmsg):
    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options(opts)

    with pytest.raises(qcng.exceptions.InputError) as e:
        qcdb.energy(mtd, molecule=nh2)

    assert errmsg in str(e.value)


@pytest.mark.parametrize('mtd,opts', [
    pytest.param('c4-ccsd', {'cfour_BASIS': 'qz2p', 'cfour_REFerence': 'roHF', 'cfour_occupation': [[3,1,1,0],[3,0,1,0]], 'cfour_SCF_CONV': 12, 'cfour_CC_CONV': 12}, marks=using_cfour),
    pytest.param('c4-ccsd', {'BASIS': 'cfour-qz2p', 'cfour_reference': 'rohf'}, marks=using_cfour),
    pytest.param('nwc-ccsd', {'basis': 'cfour-qz2p', 'qc_module': 'tce', 'nwchem_scf__rohf': True}, marks=using_nwchem),
    pytest.param('p4-ccsd', {'basis': 'cfour-qz2p', 'reference': 'rohf'}, marks=using_psi4),
    pytest.param('gms-ccsd', {'basis': 'cfour-qz2p', 'gamess_contrl__scftyp': 'rohf', 'gamess_ccinp__ncore': 0}, marks=using_gamess),
])
def test_sp_ccsd_rohf_full(mtd, opts, nh2):
    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options(opts)
#        'cfour_PRINT': 2

    e = qcdb.energy(mtd, molecule=nh2)
    scftot = -55.5847372601
    ssccsdcorl = -0.0398915
    osccsdcorl = -0.1779580
    #osccsdcorl = -0.1745752
    #ssccsdcorl = -0.0432743
    ccsdcorl = -0.217849506326
    ccsdtot = -55.802586766392

    atol = 1.e-6
    assert compare_values(scftot, qcdb.variable('scf total energy'), tnm() + 'SCF', atol=atol)
    # not printed assert compare_values(smp2corl, qcdb.variable('mp2 singles energy'), tnm() + ' MP2 singles', atol=atol)
    # not printed assert compare_values(osmp2corl, qcdb.variable('mp2 opposite-spin correlation energy'), tnm() + ' MP2 OS corl', atol=atol)
    # not printed assert compare_values(ssmp2corl, qcdb.variable('mp2 same-spin correlation energy'), tnm() + ' MP2 SS corl', atol=atol)
    # not printed assert compare_values(mp2corl, qcdb.variable('mp2 correlation energy'), tnm() + ' MP2 corl', atol=atol)
    # not printed assert compare_values(mp2tot, qcdb.variable('mp2 total energy'), tnm() + ' MP2', atol=atol)
    if not (mtd in ['gms-ccsd', 'nwc-ccsd']):
# cfour isn't splitting out the singles from OS. and psi4 isn't incl the singles in either so SS + OS != corl
# and maybe singles moved from SS to OS in Cfour btwn 2010 and 2014 versions (change in ref)
#        assert compare_values(osccsdcorl, qcdb.variable('ccsd opposite-spin correlation energy'), tnm() + ' CCSD OS corl', atol=atol)
        assert compare_values(ssccsdcorl, qcdb.variable('ccsd same-spin correlation energy'), tnm() + ' CCSD SS corl', atol=atol)
    assert compare_values(ccsdcorl, qcdb.variable('ccsd correlation energy'), tnm() + ' CCSD corl', atol=atol)
    assert compare_values(ccsdtot, qcdb.variable('ccsd total energy'), tnm() + ' CCSD', atol=atol)
    assert compare_values(ccsdcorl, qcdb.variable('current correlation energy'), tnm() + ' Current corl', atol=atol)
    assert compare_values(ccsdtot, qcdb.variable('current energy'), tnm() + ' Current', atol=atol)


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
    pytest.param('c4-ccsd', {'cfour_BASIS': 'qz2p', 'cfour_dropmo': [1], 'cfour_print': 2, 'cfour_REFerence': 'roHF', 'cfour_occupation': [[3,1,1,0],[3,0,1,0]], 'cfour_SCF_CONV': 12, 'cfour_CC_CONV': 12}, marks=using_cfour),
    pytest.param('c4-ccsd', {'BASIS': 'cfour-qz2p', 'cfour_dropmo': 1, 'cfour_reference': 'rohf'}, marks=using_cfour),
    pytest.param('nwc-ccsd', {'basis': 'cfour-qz2p', 'nwchem_tce__freeze': 1, 'qc_module': 'tce', 'nwchem_scf__rohf': True}, marks=using_nwchem),
    pytest.param('p4-ccsd', {'basis': 'cfour-qz2p', 'psi4_e_convergence': 8, 'psi4_r_convergence': 7, 'psi4_freeze_core': True, 'reference': 'rohf'}, marks=using_psi4),
    pytest.param('gms-ccsd', {'basis': 'cfour-qz2p', 'gamess_contrl__scftyp': 'rohf', 'gamess_ccinp__iconv': 9, 'gamess_scf__conv': 9}, marks=using_gamess),
])
def test_sp_ccsd_rohf_fc(mtd, opts, nh2):
    nh2 = qcdb.set_molecule(nh2)
    qcdb.set_options(opts)

    e = qcdb.energy(mtd, molecule=nh2)

    # from Cfour
    scftot = -55.5847372601
    osccsdcorl = -0.1563275
    ssccsdcorl = -0.0365048
#   ssccsdcorl = -0.036502859383024
    ccsdcorl = -0.192832282505
    ccsdtot = -55.777569542558

    # from Psi4. psi & nwchem agree to 6
    osccsdcorl = -0.152968752464362
    ssccsdcorl = -0.036502859383024
    ccsdcorl = -0.192826214666648
    ccsdtot = -55.777563474720807

    # from gamess.
    ccsdcorl = -0.1928447371
    ccsdtot = -55.7775819971

    # this VERY BAD!
    weakatol = 2.e-5

    assert compare_values(scftot, qcdb.variable('scf total energy'), tnm() + 'SCF', atol=1.e-6)
    if not (mtd in ['gms-ccsd', 'nwc-ccsd']):
#        assert compare_values(osccsdcorl, qcdb.variable('ccsd opposite-spin correlation energy'), tnm() + ' CCSD OS corl', atol=1.e-6)
        assert compare_values(ssccsdcorl, qcdb.variable('ccsd same-spin correlation energy'), tnm() + ' CCSD SS corl', atol=weakatol)
    assert compare_values(ccsdcorl, qcdb.variable('ccsd correlation energy'), tnm() + ' CCSD corl', atol=weakatol)
    assert compare_values(ccsdtot, qcdb.variable('ccsd total energy'), tnm() + ' CCSD', atol=weakatol)
    assert compare_values(ccsdcorl, qcdb.variable('current correlation energy'), tnm() + ' Current corl', atol=weakatol)
    assert compare_values(ccsdtot, qcdb.variable('current energy'), tnm() + ' Current', atol=weakatol)
