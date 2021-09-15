qcvardefs = {}
#.. include:: autodoc_abbr_options_c.rst
#
#.. _`apdx:qcvariables_alpha`:
#
#QC Variables by Alpha
#=====================
#
#.. note:: Lowercase letters in QCVariable names represent portions of
#   the variable name that vary by root number, calculation order, etc.
#   See text for fuller description.

_dip_components = ['X', 'Y', 'Z']
_quad_components = ['XX', 'XY', 'XZ', 'YY', 'YZ', 'ZZ']

qcvardefs['(T) CORRECTION ENERGY'] = {
    'units': 'Eh',
    'glossary': """
   The coupled-cluster perturbative triples correction.
"""
}

qcvardefs['A-(T) CORRECTION ENERGY'] = {
    'units': 'Eh',
    'glossary': """
   The coupled-cluster asymmetric perturbative triples correction.
   Identical to the "(AT)" and the "Lambda-CCSD(T)" correction.
"""
}

qcvardefs['T(CCSD) CORRECTION ENERGY'] = {
    'units': 'Eh',
    'glossary': """
   The coupled-cluster triples correction evaluated with CCSD amplitudes.
   Identical to the "[T]" bracket T correction.
"""
}

qcvardefs["QCISD(T) CORRECTION ENERGY"] = {
    "units": "Eh",
    'glossary': """
   The quadratice configuration interaction singles and doubles triples correction.
"""
}

qcvardefs['(Q) CORRECTION ENERGY'] = {
    'units': 'Eh',
    'glossary': """
   The coupled-cluster perturbative quadruples correction.
"""
}

qcvardefs['[Q] CORRECTION ENERGY'] = {
    'units': 'Eh',
    'glossary': """
   The coupled-cluster bracket perturbative quadruples correction.
"""
}

#.. qcvar:: AAA (T) CORRECTION ENERGY
#   AAB (T) CORRECTION ENERGY
#   ABB (T) CORRECTION ENERGY
#   BBB (T) CORRECTION ENERGY
#
#   Spin components of the UHF-based coupled-cluster perturbative triples correction [H].

#.. qcvar:: ACPF DIPOLE X
#   ACPF DIPOLE Y
#   ACPF DIPOLE Z
#
#   The three components of the dipole [Debye] for the
#   averaged coupled-pair functional level of theory.
#
#.. qcvar:: ACPF QUADRUPOLE XX
#   ACPF QUADRUPOLE XY
#   ACPF QUADRUPOLE XZ
#   ACPF QUADRUPOLE YY
#   ACPF QUADRUPOLE YZ
#   ACPF QUADRUPOLE ZZ
#
#   The six components of the quadrupole [Debye Ang] for the
#   averaged coupled-pair functional level of theory.
#
#.. qcvar:: ACPF TOTAL ENERGY
#   ACPF CORRELATION ENERGY
#
#   The total electronic energy [H] and correlation energy component [H]
#   for the averaged coupled-pair functional level of theory.
#
#.. qcvar:: AQCC DIPOLE X
#   AQCC DIPOLE Y
#   AQCC DIPOLE Z
#
#   The three components of the dipole [Debye] for the
#   averaged quadratic coupled-cluster level of theory.
#
#.. qcvar:: AQCC QUADRUPOLE XX
#   AQCC QUADRUPOLE XY
#   AQCC QUADRUPOLE XZ
#   AQCC QUADRUPOLE YY
#   AQCC QUADRUPOLE YZ
#   AQCC QUADRUPOLE ZZ

#   The six components of the quadrupole [Debye Ang] for the
#   averaged quadratic coupled-cluster level of theory.
#
#.. qcvar:: AQCC TOTAL ENERGY
#   AQCC CORRELATION ENERGY
#
#   The total electronic energy [H] and correlation energy component [H]
#   for the averaged quadratic coupled-cluster level of theory.
#.. qcvar:: BRUECKNER CONVERGED
#
#   Value 1 (0) when the Brueckner orbitals have (have not) converged.
#
#.. qcvar:: CBS TOTAL ENERGY
#   CBS CORRELATION ENERGY
#   CBS REFERENCE ENERGY
#
#   The total electronic energy [H] and its breakdown into reference total
#   energy [H] and correlation correction components [H] for the compound
#   method requested through cbs().
#
#.. qcvar:: CC ROOT n DIPOLE X
#   CC ROOT n DIPOLE Y
#   CC ROOT n DIPOLE Z
#
#   The three components of the dipole [Debye] for the requested
#   coupled cluster level of theory and root *n* (number starts at GS = 0).
#
#.. qcvar:: CC ROOT n QUADRUPOLE XX
#   CC ROOT n QUADRUPOLE XY
#   CC ROOT n QUADRUPOLE XZ
#   CC ROOT n QUADRUPOLE YY
#   CC ROOT n QUADRUPOLE YZ
#   CC ROOT n QUADRUPOLE ZZ
#
#   The six components of the quadrupole [Debye Ang] for the requested
#   coupled cluster level of theory and root *n* (numbering starts at GS = 0).
#
#.. qcvar:: CC ROOT n TOTAL ENERGY
#   CC ROOT n CORRELATION ENERGY
#
#   The total electronic energy [H] and correlation energy component [H]
#   for the requested coupled cluster level of theory and root
#   *n* (numbering starts at GS = 0).
#
#.. qcvar:: CC TOTAL ENERGY
#   CC CORRELATION ENERGY
#
for qcv in [
    "CC T1 DIAGNOSTIC",
    "CC D1 DIAGNOSTIC",
    "CC NEW D1 DIAGNOSTIC",
    "CC D2 DIAGNOSTIC",
]:
    qcvardefs[qcv] = {
        "units": "",
        "glossary":
        """
       Diagnostic of multireference character."""
}

#
#.. qcvar:: CC2 TOTAL ENERGY
#   CC2 CORRELATION ENERGY
#   CC3 TOTAL ENERGY
#   CC3 CORRELATION ENERGY
#   CC4 TOTAL ENERGY
#   CC4 CORRELATION ENERGY
#   CCnn TOTAL ENERGY
#   CCnn CORRELATION ENERGY
#
#   The total electronic energy [H] and correlation energy component [H]
#   for the requested approximate coupled-cluster (CC2, CC3, up to CC\ *nn*)
#   level of theory.
#
#.. qcvar:: CC DIPOLE X
#   CC DIPOLE Y
#   CC DIPOLE Z
#
#   The three components of the dipole [Debye] for the requested
#   coupled cluster level of theory and root.

#.. qcvar:: CC2 DIPOLE POLARIZABILITY @ xNM
#
#   The dipole polarizability [au] calculated at the CC2 level
#   for a given (x) wavelength, (x) rounded to nearest integer.
#
#.. qcvar:: CC2 SPECIFIC ROTATION (LEN) @ xNM
#
#   The specific rotation [deg/(dm (g/cm^3))] calculated at the CC2 level in the
#   length gauge for a given (x) wavelength, (x) rounded to nearest integer.
#
#.. qcvar:: CC2 SPECIFIC ROTATION (VEL) @ xNM
#
#   The specific rotation [deg/(dm (g/cm^3))] calculated at the CC2 level in the
#   velocity gauge for a given (x) wavelength, (x) rounded to nearest integer.
#
#.. qcvar:: CC2 SPECIFIC ROTATION (MVG) @ xNM
#
#   The specific rotation [deg/(dm (g/cm^3))] calculated at the CC2 level in the
#   modified velocity gauge for a given (x) wavelength, (x) rounded to nearest integer.
#
#.. qcvar:: CC QUADRUPOLE XX
#   CC QUADRUPOLE XY
#   CC QUADRUPOLE XZ
#   CC QUADRUPOLE YY
#   CC QUADRUPOLE YZ
#   CC QUADRUPOLE ZZ
#
#   The six components of the quadrupole [Debye Ang] for the requested
#   coupled cluster level of theory and root.
#
#.. qcvar:: CCSD TOTAL ENERGY
#   CCSD CORRELATION ENERGY
#   CCSDT TOTAL ENERGY
#   CCSDT CORRELATION ENERGY
#   CCSDTQ TOTAL ENERGY
#   CCSDTQ CORRELATION ENERGY
#   CCn TOTAL ENERGY
#   CCn CORRELATION ENERGY
#
#   The total electronic energy [H] and correlation energy component [H]
#   for the requested full coupled-cluster (CCSD, CCSDT, up to CC\ *n*)
#   level of theory.
#
#.. qcvar:: CCSD(T) TOTAL ENERGY
#   CCSD(T) CORRELATION ENERGY
#   A-CCSD(T) TOTAL ENERGY
#   A-CCSD(T) CORRELATION ENERGY
#   CCSDT(Q) TOTAL ENERGY
#   CCSDT(Q) CORRELATION ENERGY
#   CC(n-1)(n) TOTAL ENERGY
#   CC(n-1)(n) CORRELATION ENERGY
#
#   The total electronic energy [H] and correlation energy component [H]
#   for the perturbatively corrected coupled-cluster (CCSD(T), a-CCSD(T), CCSDT(Q),
#   up to CC(\ *n*\ -1)(\ *n*\ ) level of theory.
#
#.. qcvar:: CCSDT-1a TOTAL ENERGY
#   CCSDT-1a CORRELATION ENERGY
#   CCSDTQ-1a TOTAL ENERGY
#   CCSDTQ-1a CORRELATION ENERGY
#   CCn-1a TOTAL ENERGY
#   CCn-1a CORRELATION ENERGY
#
#   The total electronic energy [H] and correlation energy component [H]
#   for the approximate coupled-cluster (CCSD(T)-1a, CCSDT(Q)-1a,
#   up to CC\ *n*\ -1a) level of theory.
#
#.. qcvar:: CCSDT-1b TOTAL ENERGY
#   CCSDT-1b CORRELATION ENERGY
#   CCSDTQ-1b TOTAL ENERGY
#   CCSDTQ-1b CORRELATION ENERGY
#   CCn-1b TOTAL ENERGY
#   CCn-1b CORRELATION ENERGY
#
#   The total electronic energy [H] and correlation energy component [H]
#   for the approximate coupled-cluster (CCSD(T)-1b, CCSDT(Q)-1b,
#   up to CC\ *n*\ -1b) level of theory.
#
#.. qcvar:: CCSDT-3 TOTAL ENERGY
#   CCSDT-3 CORRELATION ENERGY
#   CCSDTQ-3 TOTAL ENERGY
#   CCSDTQ-3 CORRELATION ENERGY
#   CCn-3 TOTAL ENERGY
#   CCn-3 CORRELATION ENERGY
#
#   The total electronic energy [H] and correlation energy component [H]
#   for the approximate coupled-cluster (CCSD(T)-3, CCSDT(Q)-3,
#   up to CC\ *n*\ -3) level of theory.
#
#.. qcvar:: CCSD(T)_L TOTAL ENERGY
#   CCSD(T)_L CORRELATION ENERGY
#   CCSDT(Q)_L TOTAL ENERGY
#   CCSDT(Q)_L CORRELATION ENERGY
#   CC(n-1)(n)_L TOTAL ENERGY
#   CC(n-1)(n)_L CORRELATION ENERGY
#
#   The total electronic energy [H] and correlation energy component [H]
#   for the approximate coupled-cluster (CCSD(T)_L, CCSDT(Q)_L,
#   up to CC(\ *n*\ -1)(\ *n*\ )L level of theory.
#
#.. qcvar:: CCSD DIPOLE POLARIZABILITY @ xNM
#
#   The dipole polarizability [au] calculated at the CCSD level
#   for a given (x) wavelength, (x) rounded to nearest integer.

#.. qcvar:: CCSD SPECIFIC ROTATION (LEN) @ xNM
#
#   The specific rotation [deg/(dm (g/cm^3))] calculated at the CCSD level in the
#   length gauge for a given (x) wavelength, (x) rounded to nearest integer.
#
#.. qcvar:: CCSD SPECIFIC ROTATION (VEL) @ xNM
#
#   The specific rotation [deg/(dm (g/cm^3))] calculated at the CCSD level in the
#   velocity gauge for a given (x) wavelength, (x) rounded to nearest integer.
#
#.. qcvar:: CCSD SPECIFIC ROTATION (MVG) @ xNM
#
#   The specific rotation [deg/(dm (g/cm^3))] calculated at the CCSD level in the
#   modified velocity gauge for a given (x) wavelength, (x) rounded to nearest integer.
#
#.. qcvar:: CEPA(0) DIPOLE X
#   CEPA(0) DIPOLE Y
#   CEPA(0) DIPOLE Z
#
#   The three components of the dipole [Debye] for the
#   coupled electron pair approximation variant 0 level of theory.
#
#.. qcvar:: CEPA(0) QUADRUPOLE XX
#   CEPA(0) QUADRUPOLE XY
#   CEPA(0) QUADRUPOLE XZ
#   CEPA(0) QUADRUPOLE YY
#   CEPA(0) QUADRUPOLE YZ
#   CEPA(0) QUADRUPOLE ZZ
#
#   The six components of the quadrupole [Debye Ang] for the
#   coupled electron pair approximation variant 0 level of theory.

#.. qcvar:: CEPA(0) TOTAL ENERGY
#   CEPA(0) CORRELATION ENERGY
#   CEPA(1) TOTAL ENERGY
#   CEPA(1) CORRELATION ENERGY
#   CEPA(2) TOTAL ENERGY
#   CEPA(2) CORRELATION ENERGY
#   CEPA(3) TOTAL ENERGY
#   CEPA(3) CORRELATION ENERGY
#
#   The total electronic energy [H] and correlation energy component [H]
#   for the requested variant of coupled electron pair approximation level of theory.
#
#.. qcvar:: CFOUR ERROR CODE
#
#   The non-zero return value from a Cfour execution.
#
#.. qcvar:: CI DIPOLE X
#   CI DIPOLE Y
#   CI DIPOLE Z
#
#   The three components of the dipole [Debye] for the requested
#   configuration interaction level of theory and root.
#
#.. qcvar:: CI QUADRUPOLE XX
#   CI QUADRUPOLE XY
#   CI QUADRUPOLE XZ
#   CI QUADRUPOLE YY
#   CI QUADRUPOLE YZ
#   CI QUADRUPOLE ZZ

#   The six components of the quadrupole [Debye Ang] for the requested
#   configuration interaction level of theory and root.
#
#.. qcvar:: CI ROOT n -> ROOT m DIPOLE X
#   CI ROOT n -> ROOT m DIPOLE Y
#   CI ROOT n -> ROOT m DIPOLE Z
#
#   The three components of the transition dipole [Debye] between roots *n*
#   and *m* for the requested configuration interaction level of theory.
#
#.. qcvar:: CI ROOT n -> ROOT m QUADRUPOLE XX
#   CI ROOT n -> ROOT m QUADRUPOLE XY
#   CI ROOT n -> ROOT m QUADRUPOLE XZ
#   CI ROOT n -> ROOT m QUADRUPOLE YY
#   CI ROOT n -> ROOT m QUADRUPOLE YZ
#   CI ROOT n -> ROOT m QUADRUPOLE ZZ
#
#   The three components of the transition quadrupole [Debye Ang] between
#   roots *n* and *m* for the requested configuration interaction level of
#   theory.
#
#.. qcvar:: CI ROOT n DIPOLE X
#   CI ROOT n DIPOLE Y
#   CI ROOT n DIPOLE Z
#
#   The three components of the dipole [Debye] for the requested
#   configuration interaction level of theory and root *n*.
#
#.. qcvar:: CI ROOT n QUADRUPOLE XX
#   CI ROOT n QUADRUPOLE XY
#   CI ROOT n QUADRUPOLE XZ
#   CI ROOT n QUADRUPOLE YY
#   CI ROOT n QUADRUPOLE YZ
#   CI ROOT n QUADRUPOLE ZZ
#
#   The six components of the quadrupole [Debye Ang] for the requested
#   configuration interaction level of theory and root *n*.

for root in range(4):
    qcvardefs[f'CI ROOT {root} TOTAL ENERGY'] = {
        'units':
        'Eh',
        'glossary':
        """
       The total electronic energy [H] and correlation energy component [H]
       for the requested configuration interaction level of theory and root
       *n* (numbering starts at 0).
    """
    }

    qcvardefs[f'CI ROOT {root} CORRELATION ENERGY'] = {
        'units':
        'Eh',
        'glossary':
        """
       The total electronic energy [H] and correlation energy component [H]
       for the requested configuration interaction level of theory and root
       *n* (numbering starts at 0).
    """
    }

#.. qcvar:: CI STATE-AVERAGED TOTAL ENERGY
#   CI STATE-AVERAGED CORRELATION ENERGY
#
#   The total electronic energy [H] and correlation energy component [H]
#   for state-averaged CI/CASSCF levels of theory.

qcvardefs['CI TOTAL ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    """
   The total electronic energy [H]
   for the requested configuration interaction level of theory and root.
"""
}

qcvardefs['CI CORRELATION ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    """
   The correlation energy component [H]
   for the requested configuration interaction level of theory and root.
"""
}

#.. qcvar:: CISD DIPOLE X
#   CISD DIPOLE Y
#   CISD DIPOLE Z
#
#   The three components of the dipole [Debye] for the
#   configuration interaction singles and doubles level of theory and root.
#
#.. qcvar:: CISD QUADRUPOLE XX
#   CISD QUADRUPOLE XY
#   CISD QUADRUPOLE XZ
#   CISD QUADRUPOLE YY
#   CISD QUADRUPOLE YZ
#   CISD QUADRUPOLE ZZ
#
#   The six components of the quadrupole [Debye Ang] for the
#   configuration interaction singles and doubles level of theory and root.

#.. qcvar:: CISD TOTAL ENERGY
#   CISD CORRELATION ENERGY
#   CISDT TOTAL ENERGY
#   CISDT CORRELATION ENERGY
#   CISDTQ CORRELATION ENERGY
#   CISDTQ TOTAL ENERGY
#   CIn CORRELATION ENERGY
#   CIn TOTAL ENERGY
#
#   The total electronic energy [H] and correlation energy component [H]
#   for the labeled configuration interaction level of theory and root.
#   *n* is CI order for *n* > 4.
#
#.. qcvar:: CP-CORRECTED 2-BODY INTERACTION ENERGY
#
#   The interaction energy [H] considering only two-body interactions,
#   computed with counterpoise correction.
#   Related variable `UNCP-CORRECTED 2-BODY INTERACTION ENERGY <UNCP-CORRECTED2-BODYINTERACTIONENERGY>`.
#
#   .. math:: E_{\text{IE}} = E_{dimer} - \sum_{monomer}^{n}{E_{monomer}^{\text{CP}}}

#qcvardefs['CISD TOTAL ENERGY'] = {
#    'units':
#    'Eh',
#    'glossary':
#    """
#        The total energy for a configuration interaction with singles and doubles calculation.
#        """
#}
#
#qcvardefs['CISD CORRELATION ENERGY'] = {
#    'units':
#    'Eh',
#    'glossary':
#    """
#        The correlation energy corresponding to configuration interaction with single and doubles calculation.
#        """
#}

qcvardefs['CISDT TOTAL ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    """
        The total energy for a configuration interaction with singles, doubles, and triples calculation.
        """
}

qcvardefs['CISDT CORRELATION ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    """
        The correlation energy for a configuration interaction with singles, doubles, and triples calculation.
        """
}

qcvardefs['CISDTQ TOTAL ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    """
        The total energy for a configuration interaction with singles, doubles, triples, and quadruples calculation.
        """
}

qcvardefs['CISDTQ CORRELATION ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    """
        The correlation energy for a configuration interaction with singles, doubles, triples, and quadruples calculation.
        """
}

#qcvardefs['QCISD TOTAL ENERGY'] = {
#    'units':
#    'Eh',
#    'glossary':
#    """
#        The total energy for a quadratic configuration interaction with singles and doubles calculation.
#        """
#}
#
#qcvardefs['QCISD CORRELATION ENERGY'] = {
#    'units':
#    'Eh',
#    'glossary':
#    """
#        The correlation energy for a quadratic configuration interaction with singles and doubles calculation.
#        """
#}

qcvardefs['CURRENT CORRELATION ENERGY'] = {
    'units': 'Eh',
    'glossary': """
   The correlation energy corresponding to the :qcvar:`CURRENT ENERGY` variable.
"""
}

qcvardefs['CURRENT ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    """
   The total electronic energy of the most recent stage of a
   calculation (frequently overwritten). This is the quantity tracked by
   the geometry optimizer.
"""
}

qcvardefs['CURRENT DIPOLE'] = {
    'units': 'e a0',
    'dimension': '(3,)',
    'glossary': """
   The total dipole of the most recent stage of a calculation (frequently overwritten).
"""
}

qcvardefs['CURRENT GRADIENT'] = {
    'units':
    'Eh/a0',
    'dimension':
    '({nat}, 3)',
    'glossary':
    """
   The total electronic gradient of the most recent stage of a
   calculation (frequently overwritten). This is the quantity tracked by
   the geometry optimizer.
"""
}

qcvardefs['CURRENT DIPOLE GRADIENT'] = {
    'units':
    'Eh a0/u',  # = [(e a0/a0)^2/u]
    'dimension':
    '(3 * {nat}, 3)',  # T?
    'glossary':
    """
   The derivative of the dipole with respect to nuclear perturbations
   as a degree-of-freedom by dipole component array.
"""
}

qcvardefs['CURRENT HESSIAN'] = {
    'units': 'Eh/a0/a0',
    'dimension': '(3 * {nat}, 3 * {nat})',
    'glossary': """
   The total electronic Hessian of the most recent stage of a
   calculation.
"""
}

qcvardefs['CURRENT REFERENCE ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    """
   The total electronic energy of the reference stage corresponding to
   the :qcvar:`CURRENT ENERGY` variable.
"""
}

def define_prop_qcvars(mtd, extra=''):
    mtd = mtd.upper()

    qcvardefs[f'{mtd} DIPOLE'] = {
            'units':
            'e a0',
            'dimension': '(3,)',
            'glossary':
            """
           The total dipole for the {} level of theory
           and root *n* (number starts at GS = 0). {}
    """.format(mtd, extra)
    }

    qcvardefs[f'{mtd} QUADRUPOLE'] = {
            'units':
            'e a0^2',
            'dimension': ('(3,3)'),
            'glossary':
            """
           The total quadrupole for the {} level of theory
           and root *n* (number starts at GS = 0). {}
    """.format(mtd, extra)
    }


def define_transition_prop_qcvars(mtd, i, j, extra=''):
    mtd = mtd.upper()

    qcvardefs[f'{mtd} ROOT {i} -> ROOT {j} DIPOLE'] = {
            'units':
            'e a0',
            'dimension': '(3,)',
            'glossary':
            f"""
           The transition dipole array between roots {i} and {j} for the {mtd} level of theory (number starts at GS = 0). {extra}"""
    }

    qcvardefs[f'{mtd} ROOT {i} -> ROOT {j} QUADRUPOLE'] = {
            'units':
            'e a0^2',
            'dimension': ('(3,3)'),
            'glossary':
            """
           The redundant transition quadrupole between roots {i} and {j} for the {mtd} level of theory (number starts at GS = 0). {extra}"""
    }

    qcvardefs[f"{mtd} ROOT {i} -> ROOT {j} OSCILLATOR STRENGTH (LEN)"] = {
            "units": "",
            "glossary": 
   """The oscillator strength in length or velocity gauge of named method
   from ground state to root m in h symmetry (if available). DFT
   functional labeled if canonical."""
    }



def define_scf_qcvars(mtd, is_dft=True, extra='', is_dh=False, do_grad=False):
    global qcvardefs

    if mtd + ' TOTAL ENERGY' not in qcvardefs:
        qcvardefs['{} TOTAL ENERGY'.format(mtd)] = {
            'units':
            'Eh',
            'glossary':
            """
           The total electronic energy for the {} level of theory. {}
        """.format(mtd, extra)
        }

    if do_grad and mtd + ' TOTAL GRADIENT' not in qcvardefs:
        qcvardefs[f"{mtd} TOTAL GRADIENT"] = {
            "units": "Eh/a0",
            "dimension": "({nat}, 3)",
            "glossary": f"""The total electronic gradient for the {mtd} DFT level of theory. {extra}""",
        }

    if is_dft:
        qcvardefs['{} FUNCTIONAL TOTAL ENERGY'.format(mtd)] = {
            'units':
            'Eh',
            'glossary':
            r"""
           The total electronic energy for the underlying functional of the
           requested DFT method {}, without any dispersion correction. {}
        """.format(mtd, extra)
        }

    if is_dft and is_dh:
        qcvardefs['{} DOUBLE-HYBRID CORRECTION ENERGY'.format(mtd)] = {
            'units':
            'Eh',
            'glossary':
            r"""
           The scaled MP2 correlation energy correction appended to an
           underlying functional {}. {}
        """.format(mtd, extra)
        }

    qcvardefs[f'{mtd} DIPOLE'] = {
            'units':
            'e a0',
            'dimension': '(3,)',
            'glossary':
            """
           The total dipole for the {} level of theory. {}
    """.format(mtd, extra)
    }

    qcvardefs[f'{mtd} QUADRUPOLE'] = {
            'units':
            'e a0^2',
            'dimension': ('(3,3)'),
            'glossary':
            """
           The total quadrupole for the {} level of theory. {}
    """.format(mtd, extra)
    }


def define_spin_qcvars(mtd, description, do_spin=False, do_grad=False, extra=""):
    global qcvardefs

    if mtd + ' TOTAL ENERGY' not in qcvardefs:
        qcvardefs['{} TOTAL ENERGY'.format(mtd)] = {
            'units':
            'Eh',
            'glossary': """The total electronic energy for the {} level of theory. {}""".format(description, extra)
        }

    if mtd + ' CORRELATION ENERGY' not in qcvardefs:
        qcvardefs['{} CORRELATION ENERGY'.format(mtd)] = {
            'units':
            'Eh',
            'glossary': """The correlation energy for the {} level of theory. {}""".format(description, extra)
        }

    if do_grad and mtd + " TOTAL GRADIENT" not in qcvardefs:
        qcvardefs[f"{mtd} TOTAL GRADIENT"] = {
            "units": "Eh/a0",
            "dimension": "({nat}, 3)",
            "glossary": f"""The total electronic gradient for the {description} level of theory. {extra}""",
        }

    if do_grad and mtd + " TOTAL HESSIAN" not in qcvardefs:
        qcvardefs[f"{mtd} TOTAL HESSIAN"] = {
            "units": "Eh/a0/a0",
            "dimension": "(3 * {nat}, 3 * {nat})",
            "glossary": f"""The total electronic Hessian for the {description} level of theory. {extra}""",
        }

    if do_spin and f"{mtd} SAME-SPIN CORRELATION ENERGY" not in qcvardefs:
        qcvardefs[f"{mtd} SAME-SPIN CORRELATION ENERGY"] = {
            'units': 'Eh',
            'glossary':
            rf"""The unscaled portion of the {mtd} correlation energy from same-spin or triplet doubles correlations."""
        }

    if do_spin and f"{mtd} OPPOSITE-SPIN CORRELATION ENERGY" not in qcvardefs:
        qcvardefs[f"{mtd} OPPOSITE-SPIN CORRELATION ENERGY"] = {
            'units': 'Eh',
            'glossary':
            rf"""The unscaled portion of the {mtd} correlation energy from opposite-spin or singlet doubles correlations."""
        }

    if do_spin and f"{mtd} SINGLES ENERGY" not in qcvardefs:
        qcvardefs[f"{mtd} SINGLES ENERGY"] = {
            'units': 'Eh',
            'glossary':
            rf"""The singles portion of the {mtd} correlation energy. Zero except in ROHF."""
        }

    if do_spin and f"{mtd} DOUBLES ENERGY" not in qcvardefs:
        qcvardefs[f"{mtd} DOUBLES ENERGY"] = {
            'units': 'Eh',
            'glossary':
            rf"""The doubles portion of the {mtd} correlation energy including same-spin and opposite-spin correlations."""
        }


def define_dashd_qcvars(fctl, dashes, extra=''):
    for dd in dashes:
        qcvardefs['{}-{} DISPERSION CORRECTION ENERGY'.format(fctl.upper(), dd.upper())] = {
            'units':
            'Eh',
            'glossary':
            """
   The dispersion correction defined for appending to underlying functional
   {} when a DFT-D method is requested. {}
    """.format(fctl, extra)
        }

        qcvardefs['{}-{} DISPERSION CORRECTION GRADIENT'.format(fctl.upper(), dd.upper())] = {
            'units':
            'Eh',
            'dimension':
            '({nat}, 3)',
            'glossary':
            """
   The gradient to the dispersion correction defined for appending to underlying functional
   {} when a DFT-D method is requested. {}
    """.format(fctl, extra)
        }

        qcvardefs['{}-{} TOTAL ENERGY'.format(fctl.upper(), dd.upper())] = {
            'units':
            'Eh',
            'glossary':
            r"""
           The total electronic energy for the underlying functional of the
           requested DFT method {}, with dispersion correction. {}
        """.format(fctl, extra)
        }


def define_ex_transition_qcvars(max_root, mtd, sym):
    mtd = mtd.upper()
    sym = []

    for i in range(max_root):
        qcvardefs['EOM-{} ROOT 0 -> ROOT {} EXCITATION ENERGY - {} SYMMETRY'.format(mtd, i + 1, sym)] = {
            'units':
            'Eh',
            'glossary':
            """
                The excitation energy of EOM-{} in {} symmetry from 0 to root {}""".format(mtd, sym, i)
        }

    for i in range(max_root):
        qcvardefs['EOM-{} ROOT 0 -> ROOT {} TOTAL ENERGY - {} SYMMETRY'.format(mtd, i + 1, sym)] = {
            'units':
            'Eh',
            'glossary':
            """
                The total energy of EOM-{} in {} symmetry from 0 to root {}""".format(mtd, sym, i)
        }

def define_tddft_roots_qcvars(max_root, syms):
    for i in range(max_root):
        for sym in syms:
            qcvardefs[f"TDDFT ROOT {i + 1} EXCITATION ENERGY - {sym} SYMMETRY"] = {
                'units' : 'Eh',
                'glossary' : """ The excitation energy of time-dependent DFT in {} symmetry from 0 to root {}""".format(sym, i+1)
                }

            qcvardefs[f"TDDFT ROOT {i + 1} EXCITED STATE ENERGY - {sym} SYMMETRY"] = {
                'units' : 'Eh',
                'glossary' : """The excited state energy of time dependent DFT from root 0 to root {} in {} symmetry""".format(i+1, sym)
                }
        #define_prop_qcvars(f"TDDFT ROOT {i + 1}")
        define_transition_prop_qcvars("TDDFT", 0, i+1)


#TRACELESS QUADRUPOLE POLARIZABILITY XZYY
#.. qcvar:: db_name DATABASE MEAN ABSOLUTE DEVIATION
#
#   The mean absolute deviation [\ |kcalpermol|\ ] of the requested method
#   *name* from the stored reference values for the requested reactions in
#   database *db_name*. If no reference is available, this will be a large
#   and nonsensical value.
#
#   .. math:: \frac{1}{n}\sum_{rxn}^{n}{| \textsf{\textsl{name}}_{rxn}-\text{REF}_{rxn} | }
#
#.. qcvar:: db_name DATABASE MEAN SIGNED DEVIATION
#
#   The mean deviation [\ |kcalpermol|\ ] of the requested method *name*
#   from the stored reference values for the requested reactions in
#   database *db_name*. If no reference is available, this will be a large
#   and nonsensical value.
#
#   .. math:: \frac{1}{n}\sum_{rxn}^{n}{\textsf{\textsl{name}}_{rxn}-\text{REF}_{rxn}}
#
#.. qcvar:: db_name DATABASE ROOT-MEAN-SQUARE SIGNED DEVIATION
#
#   The rms deviation [\ |kcalpermol|\ ] of the requested method *name*
#   from the stored reference values for the requested reactions in
#   database *db_name*. If no reference is available, this will be a large
#   and nonsensical value.
#
#   .. math:: \sqrt{\frac{1}{n}\sum_{rxn}^{n}{(\textsf{\textsl{name}}_{rxn}-\text{REF}_{rxn})^2}}

qcvardefs['DFT FUNCTIONAL TOTAL ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The total electronic energy for the underlying functional of the
   requested DFT method, without any dispersion correction; the first four
   terms in Eq. :eq:`SCFterms` or :eq:`DFTterms`. Quantity
   :math:`E_{\text{FCTL}}` in Eqs.  :eq:`SCFterms` and :eq:`DFTterms`.
   Unless the method includes a dispersion correction, this quantity is
   equal to :qcvar:`SCF TOTAL ENERGY`.
"""
}

qcvardefs['DFT TOTAL ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The total electronic energy for the requested DFT method,
   :math:`E_{\text{DFT}}` in Eq. :eq:`DFTterms`.
???
   .. math::
   ???:nowrap:
   ???:label: DFTterms
???
   ???   \begin{align*}
   ???      E_{\text{DFT}} & = E_{NN} + E_{1e^-} + E_{2e^-} + E_{xc} + E_{\text{-D}} + E_{\text{DH}} \\
   ???                     & = E_{\text{FCTL}} + E_{\text{-D}} + E_{\text{DH}} \\
   ???                     & = E_{\text{SCF}} + E_{\text{DH}}
   ???   \end{align*}
???
   Unless the method is a DFT double-hybrid, this quantity is equal to
   :qcvar:`SCF TOTAL ENERGY`. If the method is neither a
   double-hybrid, nor dispersion corrected, this quantity is equal to
   :qcvar:`DFT FUNCTIONAL TOTAL ENERGY`.
"""
}

qcvardefs["DFT TOTAL GRADIENT"] = {
    "units": "Eh/a0",
    "glossary": r"""
    The total electronic gradient for the requested DFT method.
"""
}

qcvardefs['DFT VV10 ENERGY'] = {
    'units': 'Eh',
    'glossary': r"""
   The functional energy contribution to the total SCF energy (DFT only).
"""
}

qcvardefs['DFT XC ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The functional energy contribution [H] to the total SCF energy (DFT only).
   Quantity :math:`E_{xc}` in Eqs. :eq:`SCFterms` and :eq:`DFTterms`.
"""
}

qcvardefs['DISPERSION CORRECTION ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The dispersion correction appended to an underlying functional
   when a DFT-D method is requested. Quantity :math:`E_{\text{-D}}`
   in Eqs. :eq:`SCFterms` and :eq:`DFTterms`.
"""
}

qcvardefs['DISPERSION CORRECTION GRADIENT'] = {
    'units':
    'Eh/a0',
    'dimension':
    '({nat}, 3)',
    'glossary':
    r"""
   The gradient to the dispersion correction appended to an underlying functional
   when a DFT-D method is requested. Quantity :math:`E_{\text{-D}}`
"""
}

qcvardefs['DOUBLE-HYBRID CORRECTION ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The scaled MP2 correlation energy correction [H] appended to an
   underlying functional when a DH-DFT method is requested.
   Quantity :math:`E_{\text{DH}}` in Eq. :eq:`DFTterms`.
"""
}

qcvardefs['FCI TOTAL ENERGY'] = {
    'units': 'Eh',
    'glossary': r"""
   The total electronic energy
   for the full configuration interaction level of theory.
"""
}

qcvardefs['FCI CORRELATION ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The electronic correlation energy component [H]
   for the full configuration interaction level of theory.
"""
}

qcvardefs['HF TOTAL ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The total electronic energy for the Hartree--Fock method, without
   any dispersion correction; the first three (or four, since
   :math:`E_{xc} = 0`) terms in Eq. :eq:`SCFterms`. Quantity :math:`E_{\text{HF}}`
   in Eq. :eq:`SCFterms`.
"""
}

qcvardefs['HF TOTAL GRADIENT'] = {
    'units': 'Eh/a0',
    'dimension': '({nat}, 3)',
    'glossary': """
   The total electronic gradient of the Hartree--Fock method.
"""
}

qcvardefs['HF DIPOLE GRADIENT'] = {
    'units':
    'Eh a0/u',  # = [(e a0/a0)^2/u]
    'dimension':
    '(3 * {nat}, 3)',  # T?
    'glossary':
    """
   The derivative of the Hartree--Fock method dipole with respect to nuclear perturbations
   as a degree-of-freedom by dipole component array.
"""
}

qcvardefs['HF TOTAL HESSIAN'] = {
    'units': 'Eh/a0/a0',
    'dimension': '(3 * {nat}, 3 * {nat})',
    'glossary': """
   The total electronic energy for the Hartree-Fock method.
"""
}

qcvardefs['DMRG-SCF TOTAL ENERGY'] = {
    'units': 'Eh',
    'glossary': """
   The total DMRG total electonic energy. Not unique b/c oribital spaces
"""
}

qcvardefs['DMRG-CASPT2 TOTAL ENERGY'] = {
    'units': 'Eh',
    'glossary': """
   The total DMRG plus CASPT2 total electonic energy. Not unique b/c orbital spaces.
"""
}

#.. qcvar:: LCC2 (+LMP2) TOTAL ENERGY
#
#   The total electronic energy [H] for the local CC2 level of theory.
#
#.. qcvar:: LCCSD (+LMP2) TOTAL ENERGY
#
#   The total electronic energy [H] for the local CCSD level of theory.

qcvardefs['MP2 TOTAL ENERGY'] = {
    'units': 'Eh',
    'glossary': r"""
   The total electronic energy
   for the MP2 level of theory.
"""
}

qcvardefs['MP2 TOTAL GRADIENT'] = {
    'units': 'Eh/a0',
    'dimension': '({nat}, 3)',
    'glossary': """
   The total electronic gradient of the MP2 method.
"""
}

qcvardefs['MP2 TOTAL HESSIAN'] = {
    'units': 'Eh/a0/a0',
    'dimension': '(3 * {nat}, 3 * {nat})',
    'glossary': """
   The total electronic Hessian of the MP2 method.
"""
}

qcvardefs['MP2 CORRELATION ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The correlation energy component
   for the MP2 level of theory.

.. qcvar:: MP2 CORRELATION ENERGY

   The MP2 correlation energy for the requested DFT method,
   :math:`E_{\text{MP2corl}}` in Eq. :eq:`MP2corl`.
???
   .. math::
   ???:nowrap:
   ???:label: MP2corl
???
   ???   \begin{align*}
   ???      E_{\text{MP2corl}} & = E_{\text{S}} + E_{\text{SS}} + E_{\text{OS}} \\
   ???                         & = E_{\text{S}} + E_{\text{D}}
   ???   \end{align*}
???
"""
#                               & = E_{NN} + E_{1e^-} + E_{2e^-} + E_{xc} + E_{\text{-D}} + E_{\text{DH}} \\
#                               & = E_{\text{FCTL}} + E_{\text{-D}} + E_{\text{DH}} \\
#                               & = E_{\text{SCF}} + E_{\text{DH}}

#   Unless the method is a DFT double-hybrid, this quantity is equal to
#   :qcvar:`SCF TOTAL ENERGY`. If the method is neither a
#   double-hybrid, nor dispersion corrected, this quantity is equal to
#   :qcvar:`DFT FUNCTIONAL TOTAL ENERGY`.
}

qcvardefs['MP2 SAME-SPIN CORRELATION ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The unscaled portion of the MP2 correlation energy
   from same-spin or triplet doubles correlations.

   canonical_corl(os_scale=1, ss_scale=1) = singles + os_scale * (tot_corl - ss_corl) + ss_scale * ss_corl
   :math:`E_{\text{SS}}` in Eq. :eq:`MP2corl`.
"""
}

qcvardefs['MP2 OPPOSITE-SPIN CORRELATION ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The unscaled portion of the MP2 correlation energy
   from opposite-spin or singlet doubles correlations.
   :math:`E_{\text{OS}}` in Eq. :eq:`MP2corl`.
"""
}

qcvardefs['MP2 SINGLES ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The singles portion of the MP2 correlation energy.
   Zero except in ROHF.
   :math:`E_{\text{S}}` in Eq. :eq:`MP2corl`.
"""
}

qcvardefs['MP2 DOUBLES ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The doubles portion of the MP2 correlation energy
   including same-spin and opposite-spin correlations.
   :math:`E_{\text{D}}` in Eq. :eq:`MP2corl`.
"""
}

qcvardefs['SCS-MP2 CORRELATION ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The MP2-like correlation energy by reweighting
   MP2 DOUBLES ENERGY by 6/5 opposite-spin
   and 1/3 same-spin contributions, with any singles
   carried along.
"""
}

qcvardefs['SCS-MP2 TOTAL ENERGY'] = {
    'units': 'Eh',
    'glossary': r"""
   The total energy built from SCS-MP2 CORRELATION ENERGY
   and reference.
"""
}

qcvardefs['CUSTOM SCS-MP2 CORRELATION ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   Changeable quantity. The MP2-like correlation
   energy by any reweighting of SAME-SPIN or
   OPPOSITE-SPIN components. Depending on weights,
   this may equal any of MP2, SCS-MP2, SCS(N)-MP2,
   etc. quantities.
"""
}

qcvardefs['CUSTOM SCS-MP2 TOTAL ENERGY'] = {
    'units': 'Eh',
    'glossary': r"""
   The total energy built from
   CUSTOM SCS-MP2 CORRELATION ENERGY and reference.
"""
}

qcvardefs['SCS(N)-MP2 CORRELATION ENERGY'] = {
    'units':
    'Eh',
    'doi':
    '10.1021/ct6002737',
    'glossary':
    r"""
   The MP2-like correlation energy by reweighting
   MP2 DOUBLES ENERGY by 0 opposite-spin
   and 1.76 same-spin contributions, with any singles
   carried along.
"""
}

qcvardefs['SCS(N)-MP2 TOTAL ENERGY'] = {
    'units': 'Eh',
    'doi': '10.1021/ct6002737',
    'glossary': r"""
   The total energy built from SCS(N)-MP2 CORRELATION ENERGY
   and reference.
"""
}

qcvardefs['SCS-MP2-VDW CORRELATION ENERGY'] = {
    'units':
    'Eh',
    'doi':
    '10.1080/00268970802641242',
    'glossary':
    r"""
   The MP2-like correlation energy by reweighting
   MP2 DOUBLES ENERGY by 1.28 opposite-spin
   and 0.50 same-spin contributions, with any singles
   carried along.
"""
}

qcvardefs['SCS-MP2-VDW TOTAL ENERGY'] = {
    'units': 'Eh',
    'doi':
    '10.1080/00268970802641242',
    'glossary': r"""
   The total energy built from
   SCS-MP2-VDW CORRELATION ENERGY and reference.
"""
}

qcvardefs['CUSTOM D2 DISPERSION CORRECTION ENERGY'] = {
    'units': 'Eh',
    'glossary': r"""
   Label for D2-formula dispersion correction when
   parameters match no functional.
"""
}

qcvardefs['CUSTOM D2 DISPERSION CORRECTION GRADIENT'] = {
    'units': 'Eh',
    'dimension': '({nat}, 3)',
    'glossary': r"""
   Label for D2-formula dispersion correction gradient when
   parameters match no functional.
"""
}

qcvardefs['MP2.5 TOTAL ENERGY'] = {
    'units': 'Eh',
    'glossary': r"""
   The total electronic energy for the MP2.5 level of theory.
"""
}

qcvardefs['MP2.5 CORRELATION ENERGY'] = {
    'units': 'Eh',
    'glossary': r"""
   The correlation energy component for the MP2.5 level of theory.
???
   .. math::
   ???:nowrap:
   ???:label: MP2p5corl
???
   ???   \begin{align*}
   ???      E_{\text{MP2.5corl}} & = E_{\text{S}} + E_{\text{SS}} + E_{\text{OS}} \\
   ???                         & = E_{\text{S}} + E_{\text{D}}
   ???   \end{align*}
???
"""
}

qcvardefs['MP2.5 SINGLES ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The singles portion of the MP2.5 correlation energy.
   Zero except in ROHF.
   :math:`E_{\text{S}}` in Eq. :eq:`MP2p5corl`.
"""
}

qcvardefs['MP2.5 DOUBLES ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The doubles portion of the MP2.5 correlation energy
   including same-spin and opposite-spin correlations.
   :math:`E_{\text{D}}` in Eq. :eq:`MP2p5corl`.
"""
}

qcvardefs['MP2.5 SAME-SPIN CORRELATION ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The unscaled portion of the MP2.5 correlation energy
   from same-spin or triplet doubles correlations.

   canonical_corl(os_scale=1, ss_scale=1) = singles + os_scale * (tot_corl - ss_corl) + ss_scale * ss_corl
   :math:`E_{\text{SS}}` in Eq. :eq:`MP2p5corl`.
"""
}

qcvardefs['MP2.5 OPPOSITE-SPIN CORRELATION ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The unscaled portion of the MP2.5 correlation energy
   from opposite-spin or singlet doubles correlations.
   :math:`E_{\text{OS}}` in Eq. :eq:`MP2p5corl`.
"""
}

qcvardefs['MP3 CORRECTION ENERGY'] = {
    'units': 'Eh',
    'glossary': r"""
   The correlation energy difference between 2nd and 3-order Perturbation theory.
"""
}

qcvardefs['MP4 TOTAL ENERGY'] = {
    'units': 'Eh',
    'glossary': r"""
   The total electronic energy for 4-order Perturbation theory.
"""
}

qcvardefs['MP4 CORRELATION ENERGY'] = {
    'units': 'Eh',
    'glossary': r"""
   The correlation energy component for 4-order Perturbation theory.
"""
}

qcvardefs['MP4 CORRECTION ENERGY'] = {
    'units': 'Eh',
    'glossary': r"""
   The correlation energy difference between 3rd and 4th-order Perturbation theory.
"""
}

qcvardefs["MP4(T) CORRECTION ENERGY"] = {
    "units": "Eh",
    "glossary": r"""
   The MP4 triples component. Difference between MP4 and MP4(SDQ).
"""
}

qcvardefs['MP5 TOTAL ENERGY'] = {
    'units': 'Eh',
    'glossary': r"""
   The total electronic energy for 5-order Perturbation theory.
"""
}

qcvardefs['MP5 CORRELATION ENERGY'] = {
    'units': 'Eh',
    'glossary': r"""
   The correlation energy component for 5-order Perturbation theory.
"""
}

qcvardefs['MP6 TOTAL ENERGY'] = {
    'units': 'Eh',
    'glossary': r"""
   The total electronic energy for 6-order Perturbation theory.
"""
}

qcvardefs['MP6 CORRELATION ENERGY'] = {
    'units': 'Eh',
    'glossary': r"""
   The correlation energy component for 6-order Perturbation theory.
"""
}

qcvardefs['CCSD TOTAL ENERGY'] = {
    'units': 'Eh',
    'glossary': r"""
   The total electronic energy
   for the coupled cluster singles and doubles level of theory.
"""
}

qcvardefs['CCSD TOTAL GRADIENT'] = {
    'units': 'Eh/a0',
    'dimension': '({nat}, 3)',
    'glossary': """
   The total electronic gradient
    for the coupled cluster singles and doubles level of theory.
"""
}

qcvardefs['CCSD TOTAL HESSIAN'] = {
    'units': "Eh/a0/a0",
    "dimension": "(3 * {nat}, 3 * {nat})",
    'glossary': """
   The total electronic Hessian
    for the coupled cluster singles and doubles level of theory.
"""
}

qcvardefs['CCSD CORRELATION ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The correlation energy component
   for the coupled cluster singles and doubles level of theory.

.. qcvar:: CCSD CORRELATION ENERGY

   The CCSD correlation energy for the requested DFT method,
   :math:`E_{\text{CCSDcorl}}` in Eq. :eq:`CCSDcorl`.
???
   .. math::
   ???:nowrap:
   ???:label: CCSDcorl
???
   ???   \begin{align*}
   ???      E_{\text{CCSDcorl}} & = E_{\text{S}} + E_{\text{SS}} + E_{\text{OS}} \\
   ???                         & = E_{\text{S}} + E_{\text{D}}
   ???   \end{align*}
???
"""
}

qcvardefs['CCSD DBOC ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   A correction to the Born-Oppenheimer Approximation, calculated
   at the coupled cluster singles and doubles level of theory.
"""
}

qcvardefs['CUSTOM SCS-CCSD CORRELATION ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   Changeable quantity. The CCSD-like correlation
   energy by any reweighting of SAME-SPIN or
   OPPOSITE-SPIN components. Depending on weights,
   this may equal any of CCSD, SCS-CCSD,
   etc. quantities.
"""
}

qcvardefs['CUSTOM SCS-CCSD TOTAL ENERGY'] = {
    'units': 'Eh',
    'glossary': r"""
   The total energy built from
   CUSTOM SCS-CCSD CORRELATION ENERGY and reference.
"""
}

qcvardefs['CCSDT (PBE) TOTAL ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The total electronic energy
   for the coupled cluster singles, doubles, and triples level of theory.
"""
}

qcvardefs['CCSDT (PBE) CORRELATION ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The correlation energy component
   for the coupled cluster singles, doubles, and triples level of theory.

.. qcvar:: CCSDT CORRELATION ENERGY

   The CCSDT correlation energy for the requested DFT method,
   :math:`E_{\text{CCSDTcorl}}` in Eq. :eq:`CCSDTcorl`.
???
   .. math::
   ???:nowrap:
   ???:label: CCSDTcorl
???
   ???   \begin{align*}
   ???      E_{\text{CCSDTcorl}} & = E_{\text{S}} + E_{\text{SS}} + E_{\text{OS}} \\
   ???                         & = E_{\text{S}} + E_{\text{D}}
   ???   \end{align*}
???
"""
}

qcvardefs['CCSDTQ TOTAL ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The total electronic energy
   for the coupled cluster singles, doubles, triples, and quadruples level of theory.

   .. qcvar:: CCSDTQ TOTAL ENERGY
"""
}

qcvardefs['CCSDTQ TOTAL GRADIENT'] = {
    'units': 'Eh/a0',
    'dimension': '({nat}, 3)',
    'glossary': """
   The total electronic gradient
   for the coupled cluster singles, doubles, triples, and quadruples level of theory.
"""
}

qcvardefs['CCSDTQ CORRELATION ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The correlation energy component
   for the coupled cluster singles, doubles, triples, and quadruples level of theory.

.. qcvar:: CCSDTQ CORRELATION ENERGY

   The CCSDTQ correlation energy for the requested DFT method,
   :math:`E_{\text{CCSDTQcorl}}` in Eq. :eq:`CCSDTQcorl`.
???
   .. math::
   ???:nowrap:
   ???:label: CCSDTQcorl
???
   ???   \begin{align*}
   ???      E_{\text{CCSDTQcorl}} & = E_{\text{S}} + E_{\text{SS}} + E_{\text{OS}} \\
   ???                         & = E_{\text{S}} + E_{\text{D}}
   ???   \end{align*}
???
"""
}

qcvardefs['CCSDTQ SAME-SPIN CORRELATION ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The unscaled portion of the CCSDTQ correlation energy
   from same-spin or triplet doubles correlations.
"""
}

qcvardefs['CCSDTQ OPPOSITE-SPIN CORRELATION ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The unscaled portion of the CCSDTQ correlation energy
   from opposite-spin or singlet doubles correlations.
"""
}

qcvardefs['CCSDTQ SINGLES ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The singles portion of the CCSDTQ correlation energy.
   Zero except in ROHF.
   :math:`E_{\text{S}}` in Eq. :eq:`CCSDTQcorl`.
"""
}

qcvardefs['CCSDTQ DOUBLES ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The doubles portion of the CCSDTQ correlation energy
   including same-spin and opposite-spin correlations.
   :math:`E_{\text{D}}` in Eq. :eq:`CCSDTQcorl`.
"""
}


qcvardefs['CCSD SAME-SPIN CORRELATION ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The unscaled portion of the CCSD correlation energy
#   from same-spin or triplet doubles correlations.

#   canonical_corl(os_scale=1, ss_scale=1) = singles + os_scale * (tot_corl - ss_corl) + ss_scale * ss_corl
#   :math:`E_{\text{SS}}` in Eq. :eq:`CCSDcorl`.
"""
}

qcvardefs['CCSD OPPOSITE-SPIN CORRELATION ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
#   The unscaled portion of the CCSD correlation energy
#   from opposite-spin or singlet doubles correlations.
#   :math:`E_{\text{OS}}` in Eq. :eq:`CCSDcorl`.
"""
}

qcvardefs['CCSD SINGLES ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The singles portion of the CCSD correlation energy.
   Zero except in ROHF.
   :math:`E_{\text{S}}` in Eq. :eq:`CCSDcorl`.
"""
}

qcvardefs['CCSD DOUBLES ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The doubles portion of the CCSD correlation energy
   including same-spin and opposite-spin correlations.
   :math:`E_{\text{D}}` in Eq. :eq:`CCSDcorl`.
"""
}

qcvardefs['CCD TOTAL ENERGY'] = {
    'units': 'Eh',
    'glossary': r"""
   The total electronic energy
   for the coupled cluster doubles level of theory.
"""
}

qcvardefs['CCD CORRELATION ENERGY'] = {
    'units': 'Eh',
    'glossary': r"""
   The correlation energy component for the coupled cluster doubles level of theory.
"""
}

qcvardefs['CR-CCSD(T) TOTAL ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The total electronic energy
   for the completely renomalized coupled cluster singles and doubles plus parentheses perturbative connected triples level of theory.
"""
}

qcvardefs['CR-CCSD(T) CORRELATION ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The correlation energy component
   for the completely renomalized coupled cluster singles and doubles plus parentheses perturbative connected triples level of theory.
"""
}

#qcvardefs['CCSD[T] TOTAL ENERGY'] = {
#    'units':
#    'Eh',
#    'glossary':
#    r"""
#   The total electronic energy
#   for the coupled cluster singles and doubles plus bracket perturbative triples level of theory.
#"""
#}
#
#qcvardefs['CCSD[T] CORRELATION ENERGY'] = {
#    'units':
#    'Eh',
#    'glossary':
#    r"""
#   The correlation energy component
#   for the coupled cluster singles and doubles plus bracket perturbative triples level of theory.
#"""
#}

qcvardefs['CR-CCSD[T] TOTAL ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The total electronic energy
   for the completely renomalized coupled cluster singles and doubles plus bracket perturbative triples level of theory.
"""
}

qcvardefs['CR-CCSD[T] CORRELATION ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The correlation energy component
   for the completely renomalized coupled cluster singles and doubles plus bracket perturbative triples level of theory.
"""
}

qcvardefs['CCSDT[Q] TOTAL ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The total electronic energy
   for the coupled cluster singles doubles triples plus bracket perturbative quadruples level of theory.
"""
}

qcvardefs['CCSDT[Q] CORRELATION ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The correlation energy component
   for the coupled cluster singles and doubles and triples plus bracket perturbative quadruples level of theory.
"""
}

qcvardefs['CCSDT(Q) TOTAL ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The total electronic energy
   for the coupled cluster singles doubles triples plus perturbative quadruples level of theory.
"""
}

qcvardefs['CCSDT(Q) TOTAL GRADIENT'] = {
    'units': 'Eh/a0',
    'dimension': '({nat}, 3)',
    'glossary': """
   The total electronic gradient
   for the coupled cluster singles, doubles, and triples, plus perturbative quadruples level of theory.
"""
}

qcvardefs['CCSDT(Q) CORRELATION ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The correlation energy component
   for the coupled cluster singles and doubles and triples plus perturbative quadruples level of theory.
"""
}

qcvardefs['LCCD TOTAL ENERGY'] = {
    'units': 'Eh',
    'glossary': r"""
   The total energy for linearized coupled cluster doubles level of theory.
"""
}

qcvardefs['LCCD CORRELATION ENERGY'] = {
    'units': 'Eh',
    'glossary': r"""
   The correlation energy for linearized coupled cluster doubles level of theory.
"""
}

qcvardefs['LCCSD TOTAL ENERGY'] = {
    'units': 'Eh',
    'glossary': r"""
   The total energy for linearized coupled cluster singles and doubles level of theory.
"""
}

qcvardefs['LCCSD CORRELATION ENERGY'] = {
    'units': 'Eh',
    'glossary': r"""
   The correlation energy for linearized coupled cluster singles and doubles level of theory.
"""
}

qcvardefs['CR-CC(2,3),A TOTAL ENERGY'] = {'units': 'Eh', 'glossary': r"""
   The total energy for.
"""}

qcvardefs['CR-CC(2,3),A CORRELATION ENERGY'] = {'units': 'Eh', 'glossary': r"""
   The correlation energy for.
"""}

qcvardefs['CR-CC(2,3) TOTAL ENERGY'] = {'units': 'Eh', 'glossary': r"""
   The total energy for.
"""}

qcvardefs['CR-CC(2,3) CORRELATION ENERGY'] = {'units': 'Eh', 'glossary': r"""
   The correlation energy for.
"""}

qcvardefs['NUCLEAR REPULSION ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The nuclear repulsion energy contribution [H] to the total SCF energy.
   Quantity :math:`E_{NN}` in Eq. :eq:`SCFterms`.
???
   .. math:: E_{NN} = \sum_{i, j<i}^{N_{atom}}\frac{Z_i Z_j}{\lvert\mathbf{R}_i - \mathbf{R}_j\rvert}
   ???:label: ENN
???
"""
}  # TODO EFP?

#.. qcvar:: OCEPA(0) TOTAL ENERGY
#   OCEPA(0) CORRELATION ENERGY
#
#   The total electronic energy [H] and correlation energy component [H]
#   for the orbital-optimized CEPA(0) level of theory.
#
#.. qcvar:: OMP2 TOTAL ENERGY
#   OMP2 CORRELATION ENERGY
#
#   The total electronic energy [H] and correlation energy component [H]
#   for the orbital-optimized MP2 level of theory.
#
#.. qcvar:: OMP3 TOTAL ENERGY
#   OMP3 CORRELATION ENERGY
#
#   The total electronic energy [H] and correlation energy component [H]
#   for the orbital-optimized MP3 level of theory.

qcvardefs['ONE-ELECTRON ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The one-electron energy contribution [H] to the total SCF energy.
   Quantity :math:`E_{1e^-}` in Eq. :eq:`SCFterms`.
"""
}

qcvardefs['PCM POLARIZATION ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
    The Mutual polarization between the quantum chemical region and the classical polarizable continuum.
"""
}

qcvardefs['PE ENERGY'] = {'units': 'Eh', 'glossary': r"""
    The polarizable embedding energy.
"""}

#.. qcvar:: QCISD TOTAL ENERGY
#   QCISD CORRELATION ENERGY
#
#   The total electronic energy [H] and correlation energy component [H]
#   for the quadratic configuration interaction singles and doubles level
#   of theory.
#
#.. qcvar:: QCISD(T) TOTAL ENERGY
#   QCISD(T) CORRELATION ENERGY
#
#   The total electronic energy [H] and correlation energy component [H]
#   for the quadratic configuration interaction singles and doubles with
#   perturbative triples correction level of theory.
#
#.. qcvar:: SAPT DISP ENERGY
#   SAPT ELST ENERGY
#   SAPT EXCH ENERGY
#   SAPT IND ENERGY
#
#   Respectively, the dispersion, electrostatics, exchange, and induction
#   components of the total electronic interaction energy [H] for the the
#   requested SAPT level of theory. The sum of these four components yields
#   :qcvar:`SAPT TOTAL ENERGY`.
#
#.. qcvar:: SAPT TOTAL ENERGY
#
#   The total electronic interaction energy [H] for the requested SAPT
#   level of theory.
#
#.. qcvar:: SAPT0 TOTAL ENERGY
#   SSAPT0 TOTAL ENERGY
#   SAPT2 TOTAL ENERGY
#   SAPT2+ TOTAL ENERGY
#   SAPT2+(3) TOTAL ENERGY
#   SAPT2+3 TOTAL ENERGY
#
#   The total electronic interaction energy [H] for the labeled SAPT level
#   of theory.
#
#.. qcvar:: SAPT2+(CCD) TOTAL ENERGY
#   SAPT2+(3)(CCD) TOTAL ENERGY
#   SAPT2+3(CCD) TOTAL ENERGY
#
#   The total electronic interaction energy [H] for the labeled SAPT level
#   of theory that incorporates coupled-cluster dispersion.
#
#.. qcvar:: SAPT2+DMP2 TOTAL ENERGY
#   SAPT2+(3)DMP2 TOTAL ENERGY
#   SAPT2+3DMP2 TOTAL ENERGY
#   SAPT2+(CCD)DMP2 TOTAL ENERGY
#   SAPT2+(3)(CCD)DMP2 TOTAL ENERGY
#   SAPT2+3(CCD)DMP2 TOTAL ENERGY
#
#   The total electronic interaction energy [H] for the labeled SAPT level
#   of theory that incorporates MP2 induction correction.
#
#.. qcvar:: SCF DIPOLE X
#   SCF DIPOLE Y
#   SCF DIPOLE Z
#
#   The three components of the SCF dipole [Debye].

qcvardefs['SCF ITERATIONS'] = {'units': '', 'glossary': r"""
   The number of iterations in final? SCF set.
"""}

qcvardefs['SCF ITERATION ENERGY'] = {
    'units': 'Eh',
    'glossary': r"""
    The total SCF energy for the last completed iteration.
"""
}

qcvardefs['CCSD ITERATIONS'] = {'units': '', 'glossary': r"""
   The number of iterations in the CCSD set.
"""}

#.. qcvar:: SCF QUADRUPOLE XX
#   SCF QUADRUPOLE XY
#   SCF QUADRUPOLE XZ
#   SCF QUADRUPOLE YY
#   SCF QUADRUPOLE YZ
#   SCF QUADRUPOLE ZZ
#
#   The six components of the SCF quadrupole [Debye Ang].

qcvardefs['SCF TOTAL ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The total electronic energy of the SCF stage of the calculation.
   The :samp:`{method} CORRELATION ENERGY` variables from subsequent stages of a
   calculation are often the corresponding :samp:`{method} TOTAL ENERGY`
   variables less this quantity. Constructed from Eq. :eq:`SCFterms`,
   where this quantity is :math:`E_{\text{SCF}}`.
???
   .. math::
   ???:nowrap:
   ???:label: SCFterms
???
   ???   \begin{align*}
   ???      E_{\text{SCF}} & = E_{NN} + E_{1e^-} + E_{2e^-} + E_{xc} + E_{\text{-D}} \\
   ???                     & = E_{\text{FCTL/HF}} + E_{\text{-D}}
   ???   \end{align*}
???
   Unless the method includes a dispersion correction, this quantity is equal to :qcvar:`HF TOTAL ENERGY` (for HF) or
   :qcvar:`DFT FUNCTIONAL TOTAL ENERGY` (for
   DFT). Unless the method is a DFT double-hybrid, this quantity is equal
   to :qcvar:`DFT TOTAL ENERGY`.
"""
}

qcvardefs['SCF TOTAL GRADIENT'] = {
    'units': 'Eh/a0',
    'dimension': '({nat}, 3)',
    'glossary': """
   The total electronic gradient of the SCF stage of a calculation.
   May be HF or DFT.
"""
}

qcvardefs['SCF DIPOLE GRADIENT'] = {
    'units':
    'Eh a0/u',  # = [(e a0/a0)^2/u]
    'dimension':
    '(3 * {nat}, 3)',  # T?
    'glossary':
    """
   The derivative of the SCF dipole with respect to nuclear perturbations
   as a degree-of-freedom by dipole component array.
"""
}

qcvardefs['SCF TOTAL HESSIAN'] = {
    'units': 'Eh/a0/a0',
    'dimension': '(3 * {nat}, 3 * {nat})',
    'glossary': """
   The total electronic Hessian of the SCF stage of a calculation.
   May be HF or DFT.
"""
}

qcvardefs['TWO-ELECTRON ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The two-electron energy contribution [H] to the total SCF energy.
   Quantity :math:`E_{2e^-}` in Eq. :eq:`SCFterms`.
"""
}

qcvardefs['N ALPHA ELECTRONS'] = {
    'units': "",
    'glossary':
    r"""
The number of alpha electrons.
"""
}

qcvardefs['N BETA ELECTRONS'] = {
    'units': "",
    'glossary':
    r"""
  The number of beta electrons.
"""
}

qcvardefs["N BASIS FUNCTIONS"] = {
    "units": "",
    "glossary": r"""The number of basis functions"""}

qcvardefs["N MOLECULAR ORBITALS"] = {
    "units": "",
    "glossary": r"""The number of molecular orbitals"""}

qcvardefs['HOMO'] = {
    'units':
    'Eh a0/u', 
    'dimension':
    '(1 , 1)', 
    'glossary':
    """
    Highest occupied molecular orbitals
    """
}

qcvardefs['LUMO'] = {
    'units':
    'Eh a0/u',  
    'dimension':
    '(1 , 1)',  
    'glossary':
    """
    Lowest unoccupied molecular orbitals
    """
}

qcvardefs['N ATOMS'] = {
    'units':
    'Eh',
    'glossary':
    r"""
The number of atoms
"""
}

qcvardefs['N MO'] = {
    'units':
    'Eh',
    'glossary':
    r"""
The number of molecular orbitals
"""
}

qcvardefs['N BASIS'] = {
    'units':
    'Eh',
    'glossary':
    r"""
The number of molecular orbitals
"""
}

qcvardefs['MCSCF TOTAL ENERGY'] = {
        'units': 'Eh',
        'glossary': """ 
        The total energy for a MultiConfiguration Self-Consistent Field energy calculation.
        """
}

qcvardefs['DLDF-DAS2009 DISPERSION CORRECTION ENERGY'] = {
    'units': 'Eh',
    'glossary': r"""
   disp correction attaching to DLDF+D09 ORPHAN
"""
}

qcvardefs['WB97-CHG DISPERSION CORRECTION ENERGY'] = {
    'units': 'Eh',
    'glossary': r"""
   disp correction attaching to DLDF+D09 ORPHAN
"""
}

qcvardefs['B97-D TOTAL ENERGY'] = {
    'units': 'Eh',
    'glossary': r"""
   functional energy for B97-D w/ disp correction ORPHAN
"""
}

qcvardefs['B97-D FUNCTIONAL TOTAL ENERGY'] = {
    'units': 'Eh',
    'glossary': r"""
   functional energy for B97-D w/o disp correction ORPHAN
"""
}

# why not normal for this?
qcvardefs['B97-0 TOTAL ENERGY'] = {
    'units': 'Eh',
    'glossary': r"""
   functional energy for original hybrid B97-0 w/o disp correction ORPHAN
"""
}

qcvardefs['B97-0 FUNCTIONAL TOTAL ENERGY'] = {
    'units': 'Eh',
    'glossary': r"""
   functional energy for original hybrid B97-0 w/o disp correction ORPHAN
"""
}

qcvardefs['GROUND-STATE SYMMETRY'] = {
    'units': None,
    'glossary': r"""
        Ground state symmetry value of an excited state calculation.
        """
}

qcvardefs["CBS TOTAL ENERGY"] = {
    'units': 'Eh',
    "glossary": r"""
   The total electronic energy for the compound method requested through cbs().
"""
}

qcvardefs["CBS CORRELATION ENERGY"] = {
    'units': 'Eh',
    "glossary": r"""
   The correlation correction energy for the compound method requested through cbs().
"""
}

qcvardefs["CBS REFERENCE ENERGY"] = {
    'units': 'Eh',
    "glossary": r"""
   The reference total energy for the compound method requested through cbs().
"""
}

qcvardefs["GRID ELECTRONS ALPHA"] = {
    "units": "",
    "glossary": r"""
   The number of alpha electrons integrated by the xc quadrature grid.
"""
}

qcvardefs["GRID ELECTRONS BETA"] = {
    "units": "",
    "glossary": r"""
   The number of beta electrons integrated by the xc quadrature grid.
"""
}

qcvardefs["GRID ELECTRONS TOTAL"] = {
    "units": "",
    "glossary": r"""
   The number of total electrons integrated by the xc quadrature grid.
"""
}

## TODO kill off
#qcvardefs["-D ENERGY"] = {"units": "Eh", "glossary": ""}
#qcvardefs["C4-MP2 TOTAL GRADIENT"] = {"units": "Eh", "glossary": ""}
#qcvardefs["C4-MP3 TOTAL GRADIENT"] = {"units": "Eh", "glossary": ""}
#qcvardefs["C4-CCD TOTAL GRADIENT"] = {"units": "Eh", "glossary": ""}
#qcvardefs["C4-LCCD TOTAL GRADIENT"] = {"units": "Eh", "glossary": ""}
#qcvardefs["C4-LCCSD TOTAL GRADIENT"] = {"units": "Eh", "glossary": ""}
#qcvardefs["C4-CCSD TOTAL GRADIENT"] = {"units": "Eh", "glossary": ""}
#qcvardefs["C4-CCSD(T) TOTAL GRADIENT"] = {"units": "Eh", "glossary": ""}
#qcvardefs["C4-CCSDT-1A TOTAL GRADIENT"] = {"units": "Eh", "glossary": ""}
#qcvardefs["C4-CCSDT-1B TOTAL GRADIENT"] = {"units": "Eh", "glossary": ""}
#qcvardefs["C4-CCSDT-2 TOTAL GRADIENT"] = {"units": "Eh", "glossary": ""}
#qcvardefs["C4-CCSDT-3 TOTAL GRADIENT"] = {"units": "Eh", "glossary": ""}
#qcvardefs["C4-CCSDT TOTAL GRADIENT"] = {"units": "Eh", "glossary": ""}
#qcvardefs["MRCC TOTAL GRADIENT"] = {"units": "Eh", "glossary": ""}
#
#qcvardefs["C4-MP2 TOTAL HESSIAN"] = {"units": "Eh", "glossary": ""}
#qcvardefs["C4-CCD TOTAL HESSIAN"] = {"units": "Eh", "glossary": ""}
#qcvardefs["C4-CCSD TOTAL HESSIAN"] = {"units": "Eh", "glossary": ""}
#qcvardefs["C4-CCSD(T) TOTAL HESSIAN"] = {"units": "Eh", "glossary": ""}
#qcvardefs["C4-CCSDT-1A TOTAL HESSIAN"] = {"units": "Eh", "glossary": ""}
#qcvardefs["C4-CCSDT-1B TOTAL HESSIAN"] = {"units": "Eh", "glossary": ""}
#qcvardefs["C4-CCSDT-2 TOTAL HESSIAN"] = {"units": "Eh", "glossary": ""}
#qcvardefs["C4-CCSDT-3 TOTAL HESSIAN"] = {"units": "Eh", "glossary": ""}
#qcvardefs["C4-CCSDT TOTAL HESSIAN"] = {"units": "Eh", "glossary": ""}
#qcvardefs["MRCC TOTAL HESSIAN"] = {"units": "Eh", "glossary": ""}
## end TODO kill off

#.. qcvar:: UNCP-CORRECTED 2-BODY INTERACTION ENERGY
#
#   The interaction energy [H] considering only two-body interactions,
#   computed without counterpoise correction.
#   Related variable :qcvar:`CP-CORRECTED 2-BODY INTERACTION ENERGY <CP-CORRECTED2-BODYINTERACTIONENERGY>`.
#
#   .. math:: E_{\text{IE}} = E_{dimer} - \sum_{monomer}^{n}{E_{monomer}^{\text{unCP}}}
#
#.. qcvar:: ZAPTn TOTAL ENERGY
#   ZAPTn CORRELATION ENERGY
#
#   The total electronic energy [H] and correlation energy component [H]
#   for the labeled Z-averaged perturbation theory level.
#   *n* is ZAPT perturbation order.
define_scf_qcvars('HF', is_dft=False)
define_scf_qcvars('B3LYP', is_dh=False, do_grad=True)
define_scf_qcvars('B3LYP5', is_dh=False, do_grad=True)
define_scf_qcvars('SCF',
                  is_dft=False,
                  extra=' This is a generic HF/DFT quantity and not necessarily aligned across different calcs.')
define_scf_qcvars('B2PLYP', is_dh=True)
define_scf_qcvars('DSD-PBEP86', is_dh=True)
define_scf_qcvars('WPBE', is_dh=False)
define_scf_qcvars('WB97', is_dh=False)
define_scf_qcvars('PBE', is_dh=False, do_grad=True)
define_scf_qcvars('CAM-B3LYP', is_dh=False)
define_scf_qcvars('DLDF+D09', is_dh=False)
define_scf_qcvars('WB97X', is_dh=False)
define_scf_qcvars('PBE0-2', is_dh=True)
define_scf_qcvars('PBE0', is_dh=False)
define_scf_qcvars('SVWN', is_dh=False)
define_scf_qcvars('WB97X-D', is_dh=False)
define_scf_qcvars('PW91', is_dh=False)
define_scf_qcvars('BLYP', is_dh=False)
define_scf_qcvars('PW86PBE', is_dh=False)
define_scf_qcvars('FT97', is_dh=False)
define_scf_qcvars('BOP', is_dh=False)
define_scf_qcvars('MPWPW', is_dh=False)
define_scf_qcvars('SOGGA11', is_dh=False)
define_scf_qcvars('BP86', is_dh=False)
define_scf_qcvars('B86BPBE', is_dh=False)
define_scf_qcvars('PW6B95', is_dh=False)
define_scf_qcvars('MN15', is_dh=False)
define_scf_qcvars('SOGGA11-X', is_dh=False)
define_scf_qcvars('B2GPPLYP', is_dh=True)
define_scf_qcvars('PTPSS', is_dh=True)
define_scf_qcvars('PWPB95', is_dh=True)
define_scf_qcvars('DSD-BLYP', is_dh=True),
define_scf_qcvars('PBE0-DH', is_dh=True)
#define_scf_qcvars('B97-D', is_dh=False)
define_dashd_qcvars('bp86', dashes=['d2', 'd3', 'd3(bj)', 'd3m', 'd3m(bj)'])
define_dashd_qcvars('b3lyp', dashes=['d2', 'd3', 'd3(bj)', 'd3m', 'd3m(bj)'])
define_dashd_qcvars('b3lyp5', dashes=['d2', 'd3', 'd3(bj)', 'd3m', 'd3m(bj)'])
define_dashd_qcvars('b2plyp', dashes=['d2', 'd3', 'd3(bj)', 'd3m', 'd3m(bj)'])
define_dashd_qcvars('pbe', dashes=['d2', 'd3', 'd3(bj)', 'd3m', 'd3m(bj)'])
define_dashd_qcvars('b97', dashes=['d2', 'd3', 'd3(bj)', 'd3m', 'd3m(bj)'])
define_dashd_qcvars('blyp', dashes=['d2', 'd3', 'd3(bj)', 'd3m', 'd3m(bj)'])
define_dashd_qcvars('pbe0', dashes=['d2', 'd3', 'd3(bj)', 'd3m', 'd3m(bj)'])
define_dashd_qcvars('wpbe', dashes=['d2', 'd3', 'd3(bj)', 'd3m', 'd3m(bj)'])
#define_dashd_qcvars('wb97x', dashes=['d'])
define_ex_transition_qcvars(4, 'ccsd', ['B1', 'A1', 'A2', 'B2'])
define_tddft_roots_qcvars(10, ['B1U' , 'AG','B2U', 'B3U', 'AU'])
define_prop_qcvars('ccsd')
define_prop_qcvars('cc')  # TODO reconsider
define_spin_qcvars("LCCD", description="linearized coupled cluster doubles", do_spin=True, do_grad=True)
define_spin_qcvars("LCCSD", description="linearized coupled cluster singles and doubles", do_spin=True)
define_spin_qcvars("CEPA(0)", description="coupled electron pair approximation, variant 0", do_spin=True)
define_spin_qcvars("CCD", description="coupled cluster doubles", do_spin=True, do_grad=True)
define_spin_qcvars("MP3", description="3rd-order Moller--Plesset perturbation theory", do_spin=True, do_grad=True)
define_spin_qcvars("MP4(SDQ)", description="4rd-order Moller--Plesset perturbation theory without triples excitations", do_spin=False, do_grad=False)
define_spin_qcvars("MP4(SDTQ)", description="4rd-order Moller--Plesset perturbation theory", do_spin=False, do_grad=False)
define_spin_qcvars("CCSD+T(CCSD)", description="coupled cluster singles and doubles with triples evaluated at converged CCSD amplitudes", do_spin=False, do_grad=False)
define_spin_qcvars("CCSDT", description="coupled cluster singles, doubles, and triples excitations", do_spin=True, do_grad=True)
define_spin_qcvars("CCSDT-1A", description="coupled cluster singles, doubles, and triples excitations at approximation 1a", do_spin=False, do_grad=True)
define_spin_qcvars("CCSDT-1B", description="coupled cluster singles, doubles, and triples excitations at approximation 1b", do_spin=False, do_grad=True)
define_spin_qcvars("CCSDT-2", description="coupled cluster singles, doubles, and triples excitations at approximation 2", do_spin=False, do_grad=True)
define_spin_qcvars("CCSDT-3", description="coupled cluster singles, doubles, and triples excitations at approximation 3", do_spin=False, do_grad=True)
define_spin_qcvars("A-CCSD(T)", description="coupled cluster singles, doubles, and asymmetric perturbative triples excitations. Also known as Lambda-CCSD(T) or CCSD(AT).", do_spin=False, do_grad=True)
define_spin_qcvars("CCSD(T)", description="coupled cluster singles, doubles, and perturbative triples excitations.", do_spin=False, do_grad=True)
define_spin_qcvars("CISD", description="configuration interaction with singles and doubles", do_spin=True, do_grad=False)
define_spin_qcvars("QCISD", description="quadratic configuration interaction singles and doubles", do_spin=True, do_grad=False)
define_spin_qcvars("QCISD(T)", description="quadratic configuration interaction singles and doubles with perturbative triples", do_spin=False, do_grad=False)
