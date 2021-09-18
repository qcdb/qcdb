"""
from https://github.com/psi4/psi4/blob/master/tests/psithon1/input.dat
Spectroscopic constants of H2, and the full ci cc-pVTZ level of theory

"""
import pprint

import pytest

import qcdb

from .utils import *


@pytest.mark.parametrize(
    "qcp",
    [
        pytest.param("p4-", marks=using("psi4")),
        pytest.param("gms-", marks=using("gamess")),
    ],
)
def test_diatomic(qcp):

    h2 = qcdb.set_molecule(
        """
    H
    H 1 R
    """
    )

    energies = []
    rvals = [0.65, 0.7, 0.75, 0.8, 0.85]
    qcdb.set_keywords(
        {
            "basis": "cc-pvtz",
            "d_convergence": 12,
            "e_convergence": 12,
            # "r_convergence": 12,
        }
    )

    for r in rvals:
        h2.R = r
        energies.append(qcdb.energy(qcp + "fci"))

    # Since h2 is the active molecule it will be used by default in diatomic.anharmonicity
    # However, if you need to provide the routine a molecule pass it as the third parameter:
    #   phys_const = diatomic.anharmonicity(rvals, energies, h2)

    # phys_consts = diatomic.anharmonicity(rvals, energies)
    phys_consts = qcdb.diatomic(rvals, energies, molecule=h2)
    ref_we = 4412.731844941288
    ref_ae = 3.280703358397913
    ref_wexe = 144.6025772908216
    ref_Be = 60.66938330300022
    ref_r0 = 0.752814273047763
    ref_De = 0.045872631987045
    ref_re = 0.742567407914979
    ref_B0 = 59.02903162380126
    ref_nu = 4123.526690359645
    assert compare_values(ref_re, phys_consts["re"], 5, "Equilibrium bond length")
    assert compare_values(ref_r0, phys_consts["r0"], 5, "Zero-point corrected bond length")
    assert compare_values(ref_Be, phys_consts["Be"], 5, "Equilibrium rotational constant")
    assert compare_values(ref_B0, phys_consts["B0"], 5, "Zero-point corrected rotational constant")
    assert compare_values(ref_we, phys_consts["we"], 4, "Harmonic vibrational frequency")
    assert compare_values(ref_wexe, phys_consts["wexe"], 4, "Anharmonicity")
    assert compare_values(ref_nu, phys_consts["nu"], 4, "Anharmonic vibrational frequency")
    assert compare_values(ref_ae, phys_consts["ae"], 5, "Vibration-rotation interaction constant")
    assert compare_values(ref_De, phys_consts["De"], 5, "Quartic centrifugal distortion constant")
