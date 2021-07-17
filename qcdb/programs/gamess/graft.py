from typing import Dict, List


def gamess_list() -> List:
    """Return an array of GAMESS methods with energies. Appended
    to procedures['energy'].
        """
    val = []
    val.append('gamess')
    val.append('gms-scf')
    val.append('gms-hf')
    val.append('gms-mp2')
    val.append("gms-cisd")
    val.append("gms-lccd")
    val.append("gms-ccd")
    val.append('gms-ccsd')
    val.append('gms-ccsd+t(ccsd)')
    val.append('gms-ccsd(t)')
    val.append('gms-cr-ccl')
    val.append('gms-ccsd(tq)')
    val.append('gms-fci')
    #        val.append('gms-dft')
    #        val.append('gms-efp')

    val.append("gms-pbe")
    val.append("gms-b3lyp")
    val.append("gms-b3lyp5")

    return val


def gamess_gradient_list() -> List:
    """Return an array of GAMESS methods with gradients. Appended
    to procedures['gradient'].
        """
    val = []
    val.append('gms-hf')
    val.append('gms-scf')
    val.append('gms-mp2')

    val.append("gms-pbe")
    val.append("gms-b3lyp")
    val.append("gms-b3lyp5")

    return val


def gamess_hessian_list() -> List:
    """Return an array of GAMESS methods with Hessians. Appended
    to procedures['hessian'].
        """
    val = []
    val.append('gms-hf')

    return val


def gamess_qcvar_list() -> Dict:
    """Return a dict with keys of most GAMESS methods and values of dicts
    with the PSI Variables returned by those methods. Used by cbs()
    wrapper to avoid unnecessary computations in compound methods.
    Result is appended to ``VARH``.

    """
    VARH = {}
    VARH['gms-scf'] = {'gms-scf': 'SCF TOTAL ENERGY'}
    VARH['gms-hf'] = {'gms-hf': 'HF TOTAL ENERGY'}
    VARH['gms-mp2'] = {'gms-hf': 'HF TOTAL ENERGY', 'gms-mp2': 'MP2 TOTAL ENERGY'}
    VARH['gms-ccsd'] = {'gms-hf': 'HF TOTAL ENERGY', 'gms-mp2': 'MP2 TOTAL ENERGY', 'gms-ccsd': 'CCSD TOTAL ENERGY'}
    VARH['gms-ccsd(t)'] = {
        'gms-hf': 'HF TOTAL ENERGY',
        'gms-mp2': 'MP2 TOTAL ENERGY',
        'gms-ccsd': 'CCSD TOTAL ENERGY',
        'gms-ccsd(t)': 'CCSD(T) TOTAL ENERGY'
    }
    return VARH
