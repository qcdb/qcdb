import uuid
from typing import Dict

from ...util import conv_float2negexp


def muster_inherited_options(ropts: 'Keywords', verbose: int = 1) -> None:
    accession = uuid.uuid4()

    kwgs = {'accession': accession, 'verbose': verbose}
    do_translate = ropts.scroll['QCDB']['TRANSLATE_QCDB'].value

    # qcdb/memory [B] --> gamess/system__mwords [MW]
    qopt = ropts.scroll['QCDB']['MEMORY']
    if do_translate or qopt.is_required():
        mem = int(0.000001 * qopt.value / 8.0)
        print('\n\nMEMORY', mem, '\n\n')
        ropts.suggest('GAMESS', 'system__mwords', mem, **kwgs)

    # qcdb/reference --> gamess/contrl__scftyp
    # TODO ref or scf__ref?
    qref = ropts.scroll['QCDB']['SCF__REFERENCE'].value
    if qref in ['RHF', 'UHF', 'ROHF']:
        #ref = {'RHF': 'RHF',
        #       'UHF': 'UHF',
        #       'ROHF': 'ROHF'}[ropts.scroll['QCDB']['REFERENCE'].value]
        ropts.suggest('GAMESS', 'contrl__scftyp', qref, **kwgs)

    # qcdb/scf__d_convergence --> gamess/scf__conv
    qopt = ropts.scroll['QCDB']['SCF__D_CONVERGENCE']
    if qopt.disputed():
        conv = conv_float2negexp(qopt.value)
        ropts.suggest('GAMESS', 'scf__conv', conv, **kwgs)

    # --> gamess/ccinp__iconv
    qopt = ropts.scroll['QCDB']['E_CONVERGENCE']
    if qopt.disputed():
        conv = conv_float2negexp(qopt.value)
        ropts.suggest('GAMESS', 'ccinp__iconv', conv, **kwgs)


def muster_modelchem(name: str, dertype: int, ropts: 'Keywords', sysinfo: Dict, verbose: int = 1) -> None:
    lowername = name.lower()
    accession = uuid.uuid4()

    #runtyp = {'energy': 'energy',
    #          'gradient': 'gradient',
    #          'hessian': 'hessian',
    #          'properties': 'prop',
    #         }[driver]

    runtyp = {
        0: 'energy',
        1: 'gradient',
        2: 'hessian',
        #'properties': 'prop',
    }[dertype]

    ropts.require('GAMESS', 'contrl__runtyp', runtyp, accession=accession, verbose=verbose)

    if lowername == 'gms-gamess':
        pass

    elif lowername in ['gms-scf', 'gms-hf']:
        ropts.require('GAMESS', 'contrl__mplevl', 0, accession=accession, verbose=verbose)
        ropts.require('GAMESS', 'contrl__cityp', 'none', accession=accession, verbose=verbose)
        ropts.require('GAMESS', 'contrl__cctyp', 'none', accession=accession, verbose=verbose)

    elif lowername == 'gms-mp2':
        ropts.require('GAMESS', 'contrl__mplevl', 2, accession=accession, verbose=verbose)
        ropts.require('GAMESS', 'contrl__cityp', 'none', accession=accession, verbose=verbose)
        ropts.require('GAMESS', 'contrl__cctyp', 'none', accession=accession, verbose=verbose)

    elif lowername == 'gms-ccsd':
        ropts.require('GAMESS', 'contrl__mplevl', 0, accession=accession, verbose=verbose)
        ropts.require('GAMESS', 'contrl__cityp', 'none', accession=accession, verbose=verbose)
        ropts.require('GAMESS', 'contrl__cctyp', 'ccsd', accession=accession, verbose=verbose)

    elif lowername == 'gms-ccsd(t)':
        ropts.require('GAMESS', 'contrl__mplevl', 0, accession=accession, verbose=verbose)
        ropts.require('GAMESS', 'contrl__cityp', 'none', accession=accession, verbose=verbose)
        ropts.require('GAMESS', 'contrl__cctyp', 'ccsd(t)', accession=accession, verbose=verbose)

    elif lowername == 'gms-cr-ccl':
        ropts.require('GAMESS', 'contrl__mplevl', 0, accession=accession, verbose=verbose)
        ropts.require('GAMESS', 'contrl__cityp', 'none', accession=accession, verbose=verbose)
        ropts.require('GAMESS', 'contrl__cctyp', 'cr-ccl', accession=accession, verbose=verbose)

    elif lowername == 'gms-ccsd(tq)':
        ropts.require('GAMESS', 'contrl__mplevl', 0, accession=accession, verbose=verbose)
        ropts.require('GAMESS', 'contrl__cityp', 'none', accession=accession, verbose=verbose)
        ropts.require('GAMESS', 'contrl__cctyp', 'ccsd(tq)', accession=accession, verbose=verbose)

    elif lowername == 'gms-fci':
        ropts.require('GAMESS', 'contrl__mplevl', 0, accession=accession, verbose=verbose)
        ropts.require('GAMESS', 'contrl__cityp', 'aldet', accession=accession, verbose=verbose)
        ropts.require('GAMESS', 'contrl__cctyp', 'none', accession=accession, verbose=verbose)

        ropts.suggest('GAMESS', 'cidet__ncore', sysinfo['ncore'], accession=accession, verbose=verbose)
        ropts.suggest('GAMESS', 'cidet__nact', sysinfo['nact'], accession=accession, verbose=verbose)
        ropts.suggest('GAMESS', 'cidet__nels', sysinfo['nels'], accession=accession, verbose=verbose)

    # unused from Nuwan


#    elif lowername == 'gms-dft':
#        if dertype == 0:
#            options ['GAMESS']['GAMESS_CONTRL_RUNTYP']['value'] = 'energy'
#            options ['GAMESS']['GAMESS_CONTRL_DFTTYP']['value'] = 'b3lyp'
#        elif dertype == 1:
#            options ['GAMESS']['GAMESS_CONTRL_RUNTYP']['value'] = 'gradient'
#            options ['GAMESS']['GAMESS_CONTRL_DFTTYP']['value'] = 'b3lyp'
#        elif dertype == 2:
#            options ['GAMESS']['GAMESS_CONTRL_RUNTYP']['value'] = 'hessian'
#            options ['GAMESS']['GAMESS_CONTRL_DFTTYP']['value'] = 'b3lyp'
#
#    elif lowername == 'gms-eom-ccsd':
#        if dertype == 0:
#            options ['GAMESS']['GAMESS_CONTRL_RUNTYP']['value'] = 'energy'
#            options ['GAMESS']['GAMESS_CONTRL_CCTYP']['value']  = 'eom-ccsd'
#        elif dertype == 1:
#            options ['GAMESS']['GAMESS_CONTRL_RUNTYP']['value'] = 'gradient'
#            options ['GAMESS']['GAMESS_CONTRL_CCTYP']['value']  = 'eom-ccsd'
#        elif dertype == 2:
#            options ['GAMESS']['GAMESS_CONTRL_RUNTYP']['value'] = 'hessian'
#            options ['GAMESS']['GAMESS_CONTRL_CCTYP']['value']  = 'eom-ccsd'
#
#    elif lowername == 'gms-cis':
#        if dertype == 0:
#            options ['GAMESS']['GAMESS_CONTRL_RUNTYP']['value'] = 'energy'
#            options ['GAMESS']['GAMESS_CONTRL_CITYP']['value']  = 'cis'
#        elif dertype == 1:
#            options ['GAMESS']['GAMESS_CONTRL_RUNTYP']['value'] = 'gradient'
#            options ['GAMESS']['GAMESS_CONTRL_CITYP']['value']  = 'cis'
#        elif dertype == 2:
#            options ['GAMESS']['GAMESS_CONTRL_RUNTYP']['value'] = 'hessian'
#            options ['GAMESS']['GAMESS_CONTRL_CITYP']['value']  = 'cis'
#
#    elif lowername == 'gms-efp':
#        if dertype == 0:
#            options ['GAMESS']['GAMESS_CONTRL_RUNTYP']['value'] = 'energy'
#            options ['GAMESS']['GAMESS_CONTRL_COORD']['value']  = 'fragonly'
#        elif dertype == 1:
#            options ['GAMESS']['GAMESS_CONTRL_COORD']['value']  = 'gradient'
#            options ['GAMESS']['GAMESS_CONTRL_COORD']['value']  = 'fragonly'
#        elif dertype == 2:
#            options ['GAMESS']['GAMESS_CONTRL_RUNTYP']['value'] = 'hessian'
#            options ['GAMESS']['GAMESS_CONTRL_COORD']['value']  = 'fragonly'

    else:
        raise ValidationError("""Requested GAMESS computational methods %d is not available.""" % (lowername))
