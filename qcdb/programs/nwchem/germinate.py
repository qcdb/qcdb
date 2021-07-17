import sys
import uuid
from typing import Dict

import qcelemental as qcel

from ...exceptions import ValidationError
from ...util import conv_float2negexp

def muster_molecule(molrec: Dict, ropts: 'Keywords', verbose: int = 1) -> str:
    kwgs = {'accession': uuid.uuid4(), 'verbose': verbose}

    molcmd, moldata = qcel.molparse.to_string(molrec, dtype='nwchem', units='Bohr', return_data=True)

    for key, val in moldata['keywords'].items():
        ropts.require('NWCHEM', key, val, **kwgs)

    return molcmd


def muster_basisset(molrec: Dict, ropts: 'Keywords', qbs: 'BasisSet', verbose: int = 1) -> str:

    # this is bad b/c user can't reset puream. adjust after figuring out anonymous options better
    native_puream = qbs.has_puream()
    nwc_puream = {True: 'spherical', False: 'cartesian'}[native_puream]
    #    nwc_puream = 'cartesian'

    bascmd = f"""basis "ao basis" {nwc_puream} print\n"""  # nwc wants role, not basis name, I guess: f"""basis "{qbs.name}" {nwc_puream}\n"""
    bascmd += qbs.print_detail_nwchem()  # every unique printed. need labels in geometry, too?
    bascmd += "\nend\n"

    #ropts.require('NWCHEM', 'basis__puream', {True: 'spherical', False: 'cartesian'}[native_puream], accession=accession, verbose=verbose)

    return bascmd


def muster_modelchem(name: str, dertype: int, ropts: 'Keywords', verbose: int = 1) -> str:
    accession = sys._getframe().f_code.co_name + '_' + str(uuid.uuid4())
    kwgs = {'accession': accession, 'verbose': verbose}

    lowername = name.lower()

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

    if lowername == 'nwc-nwchem':
        pass

    elif lowername in ['nwc-scf', 'nwc-hf']:
        #ropts.require('NWCHEM', 'task__scf', runtyp, **kwgs)
        mdccmd = f'task scf {runtyp}\n\n'
    
    elif lowername == 'nwc-mcscf':
        ropts.suggest('NWCHEM', 'mcscf__active', 1, **kwgs)
        ropts.suggest('NWCHEM', 'mcscf__actelec', 1, **kwgs)
        mdccmd = f'task mcscf {runtyp} \n\n'

    #MPn options
    elif lowername == 'nwc-mp2':
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            mdccmd = f'task tce {runtyp}\n\n'
            ropts.require('NWCHEM', 'tce__mp2', True, **kwgs)
        elif ropts.scroll["QCDB"]["QC_MODULE"].value == "directmp2":
            mdccmd = f"task direct_mp2 {runtyp} \n\n"
        else:
            mdccmd = f'task mp2 {runtyp} \n\n'

    elif lowername == 'nwc-mp3':
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            mdccmd = f'task tce {runtyp}\n\n'
            ropts.require('NWCHEM', 'tce__mp3', True, **kwgs)

    elif lowername == 'nwc-mp4':
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            mdccmd = f'task tce {runtyp}\n\n'
            ropts.require('NWCHEM', 'tce__mp4', True, **kwgs)

    elif lowername == 'nwc-rimp2':
        #rimp2 requires fitting basis meaning topline must be <basis "ri-mp2 basis">
        mdccmd = f'task rimp2 {runtyp} \n\n'

    #CC options
    elif lowername == 'nwc-ccd':
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            mdccmd = f'task tce {runtyp} \n\n'
            ropts.require('NWCHEM', 'tce__ccd', True, **kwgs)

    elif lowername == 'nwc-ccsd':
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            mdccmd = f'task tce {runtyp}\n\n'
            ropts.require('NWCHEM', 'tce__ccsd', True, **kwgs)
        else:
            mdccmd = f'task ccsd {runtyp}\n\n'

            #if (ropts.scroll["NWCHEM"]["SCF__UHF"].value is True) or \
            #   (ropts.scroll["NWCHEM"]["SCF__ROHF"].value is True):
            #    ropts.suggest('NWCHEM', 'tce__ccsd', True, **kwgs)
            #    mdccmd = f'task tce {runtyp}\n\n'

    elif lowername == 'nwc-ccsd+t(ccsd)':
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            pass
            mdccmd = f'task tce {runtyp} \n\n'
            ropts.require('NWCHEM', 'tce__ccsd(t)', True, **kwgs)
        else:
            mdccmd = f'task ccsd+t(ccsd) {runtyp}\n\n'

    elif lowername == 'nwc-ccsdt':
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            mdccmd = f'task tce {runtyp} \n\n'
            ropts.require('NWCHEM', 'tce__ccsdt', True, **kwgs)
        else:
            mdccmd = f'task ccsdt {runtyp}\n\n'

    elif lowername == 'nwc-ccsd(t)':
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            mdccmd = f'task tce {runtyp} \n\n'
            ropts.require('NWCHEM', 'tce__ccsd(t)', True, **kwgs)
        else:
            mdccmd = f'task ccsd(t) {runtyp}\n\n'

    elif lowername == 'nwc-ccsdtq':
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            mdccmd = f'task tce {runtyp} \n\n'
            ropts.require('NWCHEM', 'tce__ccsdtq', True, **kwgs)
    elif lowername == 'nwc-eom-ccsd':
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            mdccmd = f'task tce {runtyp} \n\n'
            ropts.require('NWCHEM', 'tce__ccsd', True, **kwgs)
            ropts.suggest('NWCHEM', 'tce__nroots', 4, **kwgs)
    elif lowername == 'nwc-eom-ccsdt':
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            ropts.require('NWCHEM', 'tce__ccsdt', True, **kwgs)
            ropts.suggest('NWCHEM', 'tce__nroots', 4, **kwgs)


# elif lowername == 'nwc-ccsd_act':

#TCE only opts
    elif lowername == 'nwc-ccsd[t]':
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            mdccmd = f'task tce {runtyp} \n\n'
            ropts.require('NWCHEM', 'tce__ccsd[t]', True, **kwgs)

    elif lowername == 'nwc-cisd':
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            mdccmd = f'task tce {runtyp} \n\n'
            ropts.require('NWCHEM', 'tce__cisd', True, **kwgs)
    elif lowername == 'nwc-cisdt':
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            mdccmd = f'task tce {runtyp} \n\n'
            ropts.require('NWCHEM', 'tce__cisdt', True, **kwgs)
    elif lowername == 'nwc-cisdtq':
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            mdccmd = f'task tce {runtyp} \n\n'
            ropts.require('NWCHEM', 'tce__cisdtq', True, **kwgs)
    elif lowername == 'nwc-qcisd':
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            mdccmd = f'task tce {runtyp} \n\n'
            ropts.require('NWCHEM', 'tce__qcisd', True, **kwgs)

    elif lowername == 'nwc-lccd':
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            mdccmd = f'task tce {runtyp} \n\n'
            ropts.require('NWCHEM', 'tce__lccd', True, **kwgs)
    elif lowername == 'nwc-lccsd':
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            mdccmd = f'task tce {runtyp} \n\n'
            ropts.require('NWCHEM', 'tce__lccsd', True, **kwgs)
    elif lowername == 'nwc-cc2':
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            mdccmd = f'task tce {runtyp} \n\n'
            ropts.require('NWCHEM', 'tce__cc2', True, **kwgs)
    elif lowername == 'nwc-lr-ccsd':
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            mdccmd = f'task tce {runtyp} \n\n'
            ropts.require('NWCHEM', 'tce__lr-ccsd', True, **kwgs)
    #add more opts that are necessary for ccsdta
    elif lowername == 'nwc-ccsdta':
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            mdccmd = f'task tce {runtyp} \n\n'
            ropts.require('NWCHEM', 'tce__ccsdta', True, **kwgs)

    elif lowername == 'nwc-eaccsd':
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            mdccmd = f'task tce {runtyp} \n\n'
            ropts.require('NWCHEM', 'tce__eaccsd', True, **kwgs)
            ropts.suggest('NWCHEM', 'tce__nroots', 4, **kwgs)

    elif lowername == 'nwc-sodft':
        ropts.suggest('NWCHEM', 'xc', 'b3lyp', **kwgs)
        mdccmd = f'task sodft {runtyp} \n\n'

    elif lowername == 'nwc-dft':
        mdccmd = f'task dft {runtyp} \n\n'

    #DFT xc functionals
    elif lowername == 'nwc-pbe0':
        ropts.require('NWCHEM', 'xc', 'pbe0', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-pbeop':
        ropts.require('NWCHEM', 'xc', 'pbeop', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-acm':
        ropts.require('NWCHEM', 'xc', 'acm', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-bhlyp':
        ropts.require('NWCHEM', 'xc', 'bhlyp', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'

    elif lowername == 'nwc-b3lyp':
        ropts.require('NWCHEM', 'dft__xc', 'b3lyp', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'

    elif lowername == 'nwc-b1b95':
        ropts.require('NWCHEM', 'xc', 'b1b95', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-b97-1':
        ropts.require('NWCHEM', 'xc', 'becke97-1', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-b97-2':
        ropts.require('NWCHEM', 'xc', 'becke97-2', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-b97-d':
        ropts.require('NWCHEM', 'xc', 'becke97-d', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-b97-gga1':
        ropts.require('NWCHEM', 'xc', 'becke97gga1', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-b98':
        ropts.require('NWCHEM', 'xc', 'becke98', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-bb1k':
        ropts.require('NWCHEM', 'xc', 'bb1k', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-bhandh':
        ropts.require('NWCHEM', 'xc', 'beckehandh', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-bop':
        ropts.require('NWCHEM', 'xc', 'bop', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-cft97':
        ropts.require('NWCHEM', 'xc', 'cft97', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-dldf':
        ropts.require('NWCHEM', 'xc', 'dldf', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-ft97':
        ropts.require('NWCHEM', 'xc', 'ft97', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-hcth':
        ropts.require('NWCHEM', 'xc', 'hcth', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-hcth120':
        ropts.require('NWCHEM', 'xc', 'hcth120', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-hcth147':
        ropts.require('NWCHEM', 'xc', 'hcth147', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-hcth407':
        ropts.require('NWCHEM', 'xc', 'hcth407', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-hcth407p':
        ropts.require('NWCHEM', 'xc', 'hcth407p', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-hcthp14':
        ropts.require('NWCHEM', 'xc', 'hcthp14', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-m05':
        ropts.require('NWCHEM', 'xc', 'm05', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-m05-2x':
        ropts.require('NWCHEM', 'xc', 'm05-2x', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-m06':
        ropts.require('NWCHEM', 'xc', 'm06', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-m06-2x':
        ropts.require('NWCHEM', 'xc', 'm06-2x', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-m06-hf':
        ropts.require('NWCHEM', 'xc', 'm06-hf', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-m06-l':
        ropts.require('NWCHEM', 'xc', 'm06-l', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-m08-hx':
        ropts.require('NWCHEM', 'xc', 'm08-hx', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-m08-so':
        ropts.require('NWCHEM', 'xc', 'm08-so', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-m11':
        ropts.require('NWCHEM', 'xc', 'm11', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-m11-l':
        ropts.require('NWCHEM', 'xc', 'm11-l', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-mpw1b95':
        ropts.require('NWCHEM', 'xc', 'mpw1b95', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-mpw1k':
        ropts.require('NWCHEM', 'xc', 'mpw1k', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-mpwb1k':
        ropts.require('NWCHEM', 'xc', 'mpwb1k', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-pw6b95':
        ropts.require('NWCHEM', 'xc', 'pw6b95', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-pwb6k':
        ropts.require('NWCHEM', 'xc', 'pwb6k', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-op':
        ropts.require('NWCHEM', 'xc', 'op', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-optx':
        ropts.require('NWCHEM', 'xc', 'optx', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-tpssh':
        ropts.require('NWCHEM', 'xc', 'xctpssh', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-vs98':
        ropts.require('NWCHEM', 'xc', 'vs98', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-xft97':
        ropts.require('NWCHEM', 'xc', 'xft97', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-xtpss03':
        ropts.require('NWCHEM', 'xc', 'xtpss03', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    #DFT functionals potential issues - multiple options or req conditions, need to be on one line or format as:
    #dft [block start]
    #xc functional required_numeric
    #func2
    #Current fix: dft__xc is taking full string
    #end #TODO
    elif lowername == 'nwc-pw91':
        ropts.suggest('NWCHEM', 'xc', 'xperdew91 perdew91', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-pbe96':
        ropts.suggest('NWCHEM', 'xc', 'xpbe96 perdew91', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-bp91':
        ropts.suggest('NWCHEM', 'xc', 'becke88 perdew91', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-blyp':
        ropts.suggest('NWCHEM', 'xc', 'becke88 lyp', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-b97':
        ropts.suggest('NWCHEM', 'xc', 'becke97 hfexch 0.1943', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-mpw1pw':
        ropts.suggest('NWCHEM', 'xc', 'mpw91 0.75 hfexch 0.25 perdew91', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-mpwlyp1m':
        ropts.suggest('NWCHEM', 'xc', 'mpw91 0.95 hfexch 0.05 lyp', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-mpwlyp1w':
        ropts.suggest('NWCHEM', 'xc', 'mpw91 vwn_5 0.12 lyp 0.88', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-b1lyp':
        ropts.suggest('NWCHEM', 'xc', 'hfexch 0.25 becke88 0.75 lyp', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-b1pw91':
        ropts.suggest('NWCHEM', 'xc', 'hfexch 0.25 becke88 0.75 perdew91', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-b2plyp':
        ropts.suggest('NWCHEM', 'xc', 'hfexch 0.53 becke88 0.47 lyp 0.73', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'

    elif lowername == 'nwc-b3lyp5':
        # ropts.suggest('NWCHEM', 'dft__xc', 'hfexch 0.2 slater 0.8 becke88 0.72 vwn_5 0.190 lyp 0.81', **kwgs)  # ATL
        ropts.suggest('NWCHEM', 'dft__xc', 'hfexch 0.2 slater 0.8 becke88 nonlocal 0.72 vwn_5 0.190 lyp 0.81', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n' 

    elif lowername == 'nwc-b86bpbe':
        ropts.suggest('NWCHEM', 'xc', 'becke86 cpbe96', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-b97-1p':
        ropts.suggest('NWCHEM', 'xc', 'becke97gga1 hfexch', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-b97-0':
        ropts.suggest('NWCHEM', 'xc', 'becke97 hfexch', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername in ['nwc-bhandhlyp', 'nwc-bhhlyp']:
        ropts.suggest('NWCHEM', 'xc', 'becke88 0.50 hfexch 0.50 lyp', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-bp86':
        ropts.suggest('NWCHEM', 'xc', 'becke88 perdew86', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-cam-b3lyp':
        ropts.suggest('NWCHEM', 'xc', 'xcamb88 1.00 lyp 0.81 vwn_5 0.19 hfexch 1.0 cam 0.33 cam_alpha 0.19 cam_beta 0.46', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-lrc-wpbe':
        ropts.suggest('NWCHEM', 'xc', 'xwpbe 1.0 cpbe 1.0 hfexch 1.0 cam 0.3 cam_alpha 0 cam_beta 1.0', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-lrc-wpbeh':
        ropts.suggest('NWCHEM', 'xc', 'xwpbe 0.8 cpbe96 1.0 hfexch 1.0 cam 0.2 cam_alpha 0.20 cam_beta 0.80', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'

    elif lowername == 'nwc-pbe':
        ropts.require('NWCHEM', 'dft__xc', 'xpbe96 cpbe96', **kwgs)  
        mdccmd = f'task dft {runtyp} \n\n'

    elif lowername == 'nwc-pbe0-13':
        ropts.suggest('NWCHEM', 'xc', 'mpw91 0.95 HFexch 0.05 lyp', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-pbeh':
        ropts.suggest('NWCHEM', 'xc', 'xpbe96 0.75 hfexch 0.25 cpbe96' **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-pbelyp1w':
        ropts.suggest('NWCHEM', 'xc', 'xpbe96 vwn_5 0.26 lyp 0.74', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-pw86pbe':
        ropts.suggest('NWCHEM', 'xc', 'xperdew86 cpbe96' **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-tpsslyp1w':
        ropts.suggest('NWCHEM', 'xc', 'xtpss03 vwn_5 0.26 lyp 0.74', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-wpbe':
        ropts.suggest('NWCHEM', 'xc', 'xwpbe 1.0 cpbe96 1.0 hfexch 1.0 cam 0.40 cam_alpha 0 cam_beta 1.0', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-xlyp':
        ropts.suggest('NWCHEM', 'xc', 'slater -0.0690 becke88 0.722 xperdew91 0.347 lyp', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    
    # TODO dft__xc suggest or require

    elif lowername == 'nwc-tddft':
        mdccmd = f'task tddft {runtyp} \n\n'

    else:
        raise ValidationError(f"""Requested NWChem computational method {lowername} is not available.""")

    return mdccmd


def muster_inherited_keywords(ropts: 'Keywords', verbose: int = 1) -> None:
    accession = sys._getframe().f_code.co_name + '_' + str(uuid.uuid4())
    kwgs = {'accession': accession, 'verbose': verbose}

    do_translate = ropts.scroll['QCDB']['TRANSLATE_QCDB'].value

    # qcdb/memory [B] --> nwchem/total_memory [B]
    qopt = ropts.scroll['QCDB']['MEMORY']
    if do_translate or qopt.is_required():
        mem = qopt.value
        ropts.suggest('NWCHEM', 'memory', mem, **kwgs)

    # qcdb/reference --> nwchem/scf_[r|u|ro]hf
    # TODO ref or scf__ref?
    qref = ropts.scroll['QCDB']['SCF__REFERENCE'].value
    print('<<<< QREF {} >>>'.format(qref))
    if qref in ['RHF', 'UHF', 'ROHF']:
        rhf = (qref == 'RHF')
        uhf = (qref == 'UHF')
        rohf = (qref == 'ROHF')
        ropts.suggest('NWCHEM', 'scf__rhf', rhf, **kwgs)
        ropts.suggest('NWCHEM', 'scf__uhf', uhf, **kwgs)
        ropts.suggest('NWCHEM', 'scf__rohf', rohf, **kwgs)

    # qcdb/scf__d_convergence --> nwchem/scf__thresh
    # TODO should be e_conv?
    qopt = ropts.scroll['QCDB']['SCF__D_CONVERGENCE']
    if qopt.disputed():
        conv = conv_float2negexp(qopt.value)
        ropts.suggest('NWCHEM', 'scf__thresh', conv, **kwgs)

    # qcdb/e_convergence --> nwchem/ccsd(t)__thresh or dft__convergence__energy
    qopt = ropts.scroll['QCDB']['E_CONVERGENCE']
    if qopt.disputed():
        conv = conv_float2negexp(qopt.value)
        ropts.suggest('NWCHEM', 'ccsd__thresh', conv, **kwgs)
        ropts.suggest('NWCHEM', 'dft__convergence__energy', conv, **kwgs)

    # qcdb/freeze_core --> nwchem/[mp2|ccsd|tce]__freeze
    fcae = ropts.scroll["QCDB"]["FREEZE_CORE"].value
    if fcae is True:
        ropts.suggest("NWCHEM", "mp2__freeze__atomic", True, **kwgs)
        ropts.suggest("NWCHEM", "ccsd__freeze__atomic", True, **kwgs)
        ropts.suggest("NWCHEM", "tce__freeze__atomic", True, **kwgs)
    elif fcae is False:
        ropts.suggest("NWCHEM", "mp2__freeze", 0, **kwgs)
        ropts.suggest("NWCHEM", "ccsd__freeze", 0, **kwgs)
        ropts.suggest("NWCHEM", "tce__freeze", 0, **kwgs)
