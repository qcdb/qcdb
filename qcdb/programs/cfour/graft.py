from typing import Dict, List


def cfour_list() -> List:
    """Return an array of Cfour methods with energies. Appended
    to procedures['energy'].

    """
    val = []
    val.append('c4-cfour')
    val.append('c4-scf')
    val.append('c4-hf')
    val.append('c4-mp2')
    val.append('c4-mp3')
    val.append('c4-mp4(sdq)')
    val.append('c4-mp4')
    val.append("c4-lccd")
    val.append("c4-lccsd")
    val.append('c4-cc2')
    val.append("c4-ccd")
    val.append('c4-ccsd')
    val.append('c4-ccsd-dboc')
    val.append('c4-cc3')
    val.append("c4-ccsd+t(ccsd)")
    val.append('c4-ccsd(t)')
    val.append('c4-ccsdt')
    val.append('c4-ccsdt(q)')
    val.append('c4-ccsdtq')
    return val


def cfour_gradient_list() -> List:
    """Return an array of Cfour methods with analytical gradients.
    Appended to procedures['gradient'].

    """
    val = []
    val.append('cfour')
    val.append('c4-scf')
    val.append('c4-hf')
    val.append('c4-mp2')
    val.append('c4-mp3')
    val.append('c4-mp4(sdq)')
    val.append('c4-mp4')
    val.append('c4-cc2')
    val.append("c4-ccd")
    val.append('c4-ccsd')
    val.append('c4-cc3')
    val.append('c4-ccsd(t)')
    val.append('c4-ccsdt')
    val.append('c4-ccsdt(q)')
    val.append('c4-ccsdtq')
    return val


def cfour_hessian_list() -> List:
    """Return an array of Cfour methods with analytical Hessians.
    Appended to procedures['hessian'].

    """
    val = []
    val.append('cfour')
    val.append('c4-scf')
    val.append('c4-hf')
    val.append('c4-mp2')
    val.append('c4-ccd')
    val.append('c4-ccsd')
    val.append('c4-ccsd(t)')
    return val


def cfour_qcvar_list() -> Dict:
    """Return a dict with keys of most Cfour methods and values of dicts
    with the PSI Variables returned by those methods. Used by cbs()
    wrapper to avoid unnecessary computations in compound methods.
    Result is appended to ``VARH``.

    """
    VARH = {}
    VARH['c4-scf'] = {'c4-scf': 'SCF TOTAL ENERGY'}
    VARH['c4-hf'] = {'c4-hf': 'HF TOTAL ENERGY'}
    VARH['c4-mp2'] = {'c4-hf': 'HF TOTAL ENERGY', 'c4-mp2': 'MP2 TOTAL ENERGY'}
    VARH['c4-mp3'] = {
        'c4-hf': 'HF TOTAL ENERGY',
        'c4-mp2': 'MP2 TOTAL ENERGY',
        'c4-mp2.5': 'MP2.5 TOTAL ENERGY',
        'c4-mp3': 'MP3 TOTAL ENERGY'
    }
    VARH['c4-mp4(sdq)'] = {
        'c4-hf': 'HF TOTAL ENERGY',
        'c4-mp2': 'MP2 TOTAL ENERGY',
        'c4-mp2.5': 'MP2.5 TOTAL ENERGY',
        'c4-mp3': 'MP3 TOTAL ENERGY',
        'c4-mp4(sdq)': 'MP4(SDQ) TOTAL ENERGY'
    }
    VARH['c4-mp4'] = {
        'c4-hf': 'HF TOTAL ENERGY',
        'c4-mp2': 'MP2 TOTAL ENERGY',
        'c4-mp2.5': 'MP2.5 TOTAL ENERGY',
        'c4-mp3': 'MP3 TOTAL ENERGY',
        'c4-mp4(sdq)': 'MP4(SDQ) TOTAL ENERGY',
        'c4-mp4': 'MP4(SDTQ) TOTAL ENERGY'
    }
    VARH['c4-cc2'] = {'c4-hf': 'HF TOTAL ENERGY', 'c4-mp2': 'MP2 TOTAL ENERGY', 'c4-cc2': 'CC2 TOTAL ENERGY'}
    VARH['c4-ccsd'] = {'c4-hf': 'HF TOTAL ENERGY', 'c4-mp2': 'MP2 TOTAL ENERGY', 'c4-ccsd': 'CCSD TOTAL ENERGY'}
    VARH['c4-cc3'] = {'c4-hf': 'HF TOTAL ENERGY', 'c4-mp2': 'MP2 TOTAL ENERGY', 'c4-cc3': 'CC3 TOTAL ENERGY'}
    VARH['c4-ccsd(t)'] = {
        'c4-hf': 'HF TOTAL ENERGY',
        'c4-mp2': 'MP2 TOTAL ENERGY',
        'c4-ccsd': 'CCSD TOTAL ENERGY',
        'c4-ccsd(t)': 'CCSD(T) TOTAL ENERGY'
    }
    VARH['c4-ccsdt'] = {
        'c4-hf': 'HF TOTAL ENERGY',
        'c4-mp2': 'MP2 TOTAL ENERGY',
        'c4-ccsd': 'CCSD TOTAL ENERGY',
        'c4-ccsdt': 'CCSDT TOTAL ENERGY'
    }

    return VARH
