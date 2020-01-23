from typing import Dict, List


def nwchem_list() -> List:
    """Return an array of NWChem methods with energies. Appended
    to procedures[qcdb.energy()].

    """
    val = []
    val.append('nwchem')
    val.append('nwc-scf')
    val.append('nwc-hf')
    val.append('nwc-mp2')
    val.append('nwc-mp3')
    val.append('nwc-mp4')
    val.append('nwc-direct_mp2') #untested
    val.append('nwc-rimp2') #untested
    val.append('nwc-dft')
    val.append('nwc-cc2')
    val.append('nwc-ccsd')
    val.append('nwc-ccsdt')
    val.append('nwc-ccsdtq')
    val.append('nwc-ccsd[t]')
    val.append('nwc-ccsd(t)')
    val.append('nwc-eaccsd')
    val.append('nwc-qcisd')
    val.append('nwc-cisd')
    val.append('nwc-cisdt')
    val.append('nwc-cisdtq')
    val.append('nwc-lccd')
    val.append('nwc-lccsd')
    val.append('nwc-ccd')
    val.append('nwc-lr-ccsd')
    val.append('nwc-eom-ccsd')  #untested
    val.append('nwc-eom-ccsdt')  #untested
    val.append('nwc-eom-ccsdtq')  #untested
    val.extend(_dft_functionals_list())
    val.append('nwc-tddft')
    val.append('nwc-mcscf')

    return val

def nwchem_properties_list() -> List:
    """Return an array of NWChem methods with properties. Appended to
    procedures ['properties'].
    """

    val = []
    val.append('nwc-scf')
    val.append('nwc-hf')
    val.append('nwc-dft')

def nwchem_gradient_list() -> List:
    """Return an array of NWChem methods with energies. Appended
    to procedures['gradient'].

    """
    val = []
    #    val.append('nwchem')
    val.append('nwc-scf')
    val.append('nwc-hf')
    val.append('nwc-mp2')
    val.append('nwc-direct_mp2') #untested
    val.append('nwc-rimp2') #untested
    val.append('nwc-dft')
    val.append('nwc-ccsd')
    val.append('nwc-ccsdt')
    val.append('nwc-ccsdtq')
    val.append('nwc-ccsd(t)')
    val.append('nwc-lccd')
    val.append('nwc-eom-ccsd')
    val.append('nwc-eom-ccsdt')
    val.append('nwc-eom-ccsdtq')
    val.extend(_dft_functionals_list())
    val.append('nwc-tddft')

    return val


def nwchem_hessian_list() -> List:
    """Return an array of NWChem methods with energies. Appended to procedures
    ['hessian']"""

    val = []
    val.append('nwc-scf')
    val.append('nwc-hf')
    val.append('nwc-dft')

    return []


def _dft_functionals_list() -> List:
    """Return an array of available method of dft functionals through nwchem interface. 

    """
    val = []
    val.append('nwc-b3lyp')
    val.append('nwc-pbe0')
    val.append('nwc-m05')
    val.append('nwc-m05-2x')
    val.append('nwc-ft97')
    val.append('nwc-hcth')
    val.append('nwc-hcth120')
    val.append('nwc-dldf')
    val.append('nwc-b2plyp')
    val.append('nwc-blyp')
    val.append('nwc-bp86')
    val.append('nwc-pbe')
    val.append('nwc-b97-0')
    val.append('nwc-b97-1')
    val.append('nwc-b97-2')
    val.append('nwc-hcthp14')
    val.append('nwc-hcth407p')
    val.append('nwc-b1b95')
    val.append('nwc-m06')
    val.append('nwc-m06-2x')
    val.append('nwc-m06-hf')
    val.append('nwc-m08-hx')
    val.append('nwc-m08-so')
    val.append('nwc-m11')
    val.append('nwc-m11-l')
    val.append('nwc-mpw1b95')
    val.append('nwc-mpw1k')
    val.append('nwc-mpwb1k')
    val.append('nwc-pw6b95')
    val.append('nwc-pwb6k')
    val.append('nwc-b97-gga1')
    val.append('nwc-b97')
    val.append('nwc-tpssh')
    val.append('nwc-bhandh')
    val.append('nwc-cam-b3lyp')
    val.append('nwc-bop')
    val.append('nwc-hf')
    val.append('nwc-pw91')
    val.append('nwc-b1lyp')
    val.append('nwc-b1pw91')
    val.append('nwc-b3lyp5')
    val.append('nwc-b86bpbe')
    val.append('nwc-b97-1p')
    val.append('nwc-b97-3')
    val.append('nwc-bhandhlyp')
    val.append('nwc-bhhlyp')
    val.append('nwc-lrc-wpbe')
    val.append('nwc-lrc-wpbeh')
    val.append('nwc-mpw1pw')
    val.append('nwc-mpwlyp1m')
    val.append('nwc-mpwlyp1w')
    val.append('nwc-pbe0-13')
    val.append('nwc-pbeh')
    val.append('nwc-pbelyp1w')
    val.append('nwc-pw86pbe')
    val.append('nwc-tpsslyp1w')
    val.append('nwc-xlyp')
    val.append('nwc-wpbe')
    val.append('nwc-bhlyp')
    val.append('nwc-hcth407')
    val.append('nwc-pbeop')
    val.append('nwc-b97-d')
    val.append('nwc-cft97')
    val.append('nwc-op')
    val.append('nwc-optx')
    val.append('nwc-b98')
    val.append('nwc-acm')
    val.append('nwc-xtpss03')
    val.append('nwc-bb1k')
    val.append('nwc-vs98')
    val.append('nwc-m06-l')
    val.append('nwc-hcth147')
    
    return val


def nwchem_qcvar_list() -> Dict:
    """Return a dict with keys of most NWChem methods and values of dicts
    with the PSI Variables returned by those methods. Used by cbs()
    wrapper to avoid unnecessary computations in compound methods.
    Result is appended to ``VARH``.
        
    """
    VARH = {}
    VARH['nwc-scf'] = {'nwc-scf': 'HF TOTAL ENERGY'}
    VARH['nwc-hf'] = {'nwc-hf': 'HF TOTAL ENERGY'}
    VARH['nwc-mp2'] = {'nwc-hf': 'HF TOTAL ENERGY', 'nwc-mp2': 'MP2 TOTAL ENERGY'}
    VARH['nwc-mp3'] = {'nwc-hf': 'HF TOTAL ENERGY', 'nwc-mp2': 'MP2 TOTAL ENERGY', 'nwc-mp3': 'MP3 TOTAL ENERGY'}
    VARH['nwc-mp4'] = {
        'nwc-hf': 'HF TOTAL ENERGY',
        'nwc-mp2': 'MP2 TOTAL ENERGY',
        'nwc-mp3': 'MP3 TOTAL ENERGY',
        'nwc-mp4': 'MP4 TOTAL ENERGY'
    }
    VARH['nwc-lccd'] = {'nwc-hf': 'HF TOTAL ENERGY', 'nwc-lccd': 'LCCD TOTAL ENERGY'}
    VARH['nwc-cisd'] = {'nwc-hf': 'HF TOTAL ENERGY', 'nwc-mp2': 'MP2 TOTAL ENERGY', 'nwc-cisd': 'CISD TOTAL ENERGY'}
    VARH['nwc-cisdt'] = {
        'nwc-hf': 'HF TOTAL ENERGY',
        'nwc-mp2': 'MP2 TOTAL ENERGY',
        'nwc-cisd': 'CISD TOTAL ENERGY',
        'nwc-cisdt': 'CISDT TOTAL ENERGY'
    }
    VARH['nwc-cisdtq'] = {
        'nwc-hf': 'HF TOTAL ENERGY',
        'nwc-mp2': 'MP2 TOTAL ENERGY',
        'nwc-cisd': 'CISD TOTAL ENERGY',
        'nwc-cisdt': 'CISDT TOTAL ENERGY',
        'nwc-cisdtq': 'CISDTQ TOTAL ENERGY'
    }
    VARH['nwc-ccsd'] = {'nwc-hf': 'HF TOTAL ENERGY', 'nwc-mp2': 'MP2 TOTAL ENERGY', 'nwc-ccsd': 'CCSD TOTAL ENERGY'}

    VARH['nwc-ccsdt'] = {'nwc-hf': 'HF TOTAL ENERGY', 'nwc-ccsdt': 'CCSDT TOTAL ENERGY'}

    VARH['nwc-ccsdtq'] = {'nwc-hf': 'HF TOTAL ENERGY', 'nwc-ccsdtq': 'CCSDTQ TOTAL ENERGY'}

    VARH['nwc-ccsd(t)'] = {
        'nwc-hf': 'HF TOTAL ENERGY',
        'nwc-mp2': 'MP2 TOTAL ENERGY',
        'nwc-ccsd': 'CCSD TOTAL ENERGY',
        'nwc-ccsd(t)': 'CCSD(T) TOTAL ENERGY'
    }

    # VARH['nwc-dft'] = {
    #                        'nwc-dfttot' : 'DFT TOTAL ENERGY'}

    VARH['nwc-eom-ccsd'] = {
        'nwc-hf': 'HF TOTAL ENERGY',
        # 'nwc-dfttot' : 'DFT TOTAL ENERGY',
        'nwc-ccsd': 'CCSD TOTAL ENERGY'
    }

    return VARH
