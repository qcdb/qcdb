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

.. qcvar:: (T) CORRECTION ENERGY
"""
}

qcvardefs['[T] CORRECTION ENERGY'] = {
    'units': 'Eh',
    'glossary': """
   The coupled-cluster bracket perturbative triples correction.

.. qcvar:: [T] CORRECTION ENERGY
"""
}

#.. qcvar:: A-(T) CORRECTION ENERGY
#
#   The coupled-cluster asymmetric perturbative triples correction [H].
#
#.. qcvar:: MP4(T) CORRECTION ENERGY
#
#   The MP4 triples component [H]. Quantity is second right-hand term in
#   Eq. :eq:`MP4terms`.
#
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
#.. qcvar:: CC T1 DIAGNOSTIC
#   CC D1 DIAGNOSTIC
#   CC NEW D1 DIAGNOSTIC
#   CC D2 DIAGNOSTIC
#
#   Diagnostic of multireference character.
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

qcvardefs['CISD TOTAL ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    """
        The total energy for a configuration interaction with singles and doubles calculation.
        """
}

qcvardefs['CISD CORRELATION ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    """
        The correlation energy corresponding to configuration interaction with single and doubles calculation.
        """
}

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

qcvardefs['QCISD TOTAL ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    """
        The total energy for a quadratic configuration interaction with singles and doubles calculation.
        """
}

qcvardefs['QCISD CORRELATION ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    """
        The correlation energy for a quadratic configuration interaction with singles and doubles calculation.
        """
}

qcvardefs['CURRENT CORRELATION ENERGY'] = {
    'units': 'Eh',
    'glossary': """
   The correlation energy corresponding to the :qcvar:`CURRENT ENERGY <CURRENTENERGY>` variable.
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
   the :qcvar:`CURRENT ENERGY <CURRENTENERGY>` variable.
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


def define_scf_qcvars(mtd, is_dft=True, extra='', is_dh=False):
    global qcvardefs

    if mtd + 'TOTAL ENERGY' not in qcvardefs:
        qcvardefs['{} TOTAL ENERGY'.format(mtd)] = {
            'units':
            'Eh',
            'glossary':
            """
           The total electronic energy for the {} level of theory. {}
        """.format(mtd, extra)
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

def define_tddft_roots_qcvars(max_root, sym):
    sym = []

    for i in range(max_root):
        qcvardefs['TDDFT ROOT {} EXCITATION ENERGY - {} SYMMETRY'.format(i+1, sym)] = {
                'units' : 'Eh',
                'glossary' : """ The excitation energy of time-dependent DFT in {} symmetry from 0 to root {}""".format(sym, i+1)
                }

        qcvardefs['TDDFT ROOT {} EXCITED STATE ENERGY - {} SYMMETRY'.format(i+1, sym)] = {
                'units' : 'Eh',
                'glossary' : """The excited state energy of time dependent DFT from root 0 to root {} in {} symmetry""".format(i+1, sym)
                }

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
   equal to :qcvar:`SCF TOTAL ENERGY <SCFTOTALENERGY>`.
"""
}

qcvardefs['DFT TOTAL ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The total electronic energy for the requested DFT method,
   :math:`E_{\text{DFT}}` in Eq. :eq:`DFTterms`.

   .. math::
      :nowrap:
      :label: DFTterms

         \begin{align*}
            E_{\text{DFT}} & = E_{NN} + E_{1e^-} + E_{2e^-} + E_{xc} + E_{\text{-D}} + E_{\text{DH}} \\
                           & = E_{\text{FCTL}} + E_{\text{-D}} + E_{\text{DH}} \\
                           & = E_{\text{SCF}} + E_{\text{DH}}
         \end{align*}

   Unless the method is a DFT double-hybrid, this quantity is equal to
   :qcvar:`SCF TOTAL ENERGY <SCFTOTALENERGY>`. If the method is neither a
   double-hybrid, nor dispersion corrected, this quantity is equal to
   :qcvar:`DFT FUNCTIONAL TOTAL ENERGY <DFTFUNCTIONALTOTALENERGY>`.
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
    'units': 'Eh/a0/a0',
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

   .. math::
      :nowrap:
      :label: MP2corl

         \begin{align*}
            E_{\text{MP2corl}} & = E_{\text{S}} + E_{\text{SS}} + E_{\text{OS}} \\
                               & = E_{\text{S}} + E_{\text{D}} \\
#                               & = E_{NN} + E_{1e^-} + E_{2e^-} + E_{xc} + E_{\text{-D}} + E_{\text{DH}} \\
#                               & = E_{\text{FCTL}} + E_{\text{-D}} + E_{\text{DH}} \\
#                               & = E_{\text{SCF}} + E_{\text{DH}}
         \end{align*}

#   Unless the method is a DFT double-hybrid, this quantity is equal to
#   :qcvar:`SCF TOTAL ENERGY <SCFTOTALENERGY>`. If the method is neither a
#   double-hybrid, nor dispersion corrected, this quantity is equal to
#   :qcvar:`DFT FUNCTIONAL TOTAL ENERGY <DFTFUNCTIONALTOTALENERGY>`.

"""
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

qcvardefs[''] = {'units': 'Eh', 'glossary': r"""
"""}

qcvardefs[''] = {'units': 'Eh', 'glossary': r"""
"""}

qcvardefs[''] = {'units': 'Eh', 'glossary': r"""
"""}

qcvardefs[''] = {'units': 'Eh', 'glossary': r"""
"""}

qcvardefs[''] = {'units': 'Eh', 'glossary': r"""
"""}

qcvardefs[''] = {'units': 'Eh', 'glossary': r"""
"""}

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

#.. qcvar:: MP3 TOTAL ENERGY
#   MP3 CORRELATION ENERGY
#
#   The total electronic energy [H] and correlation energy component [H]
#   for the MP3 level of theory.
#
#.. qcvar:: MP4(SDQ) TOTAL ENERGY
#   MP4(SDQ) CORRELATION ENERGY
#
#   The total electronic energy [H] and correlation energy component [H]
#   for the MP4 singles, doubles, quadruples level of theory.  Quantity
#   :qcvar:`MP4(SDQ) CORRELATION ENERGY <MP4(SDQ)CORRELATIONENERGY>` is
#   first right-hand term in Eq. :eq:`MP4terms`.
#
#.. qcvar:: MP4 TOTAL ENERGY
#   MP4 CORRELATION ENERGY
#   MP4(SDTQ) TOTAL ENERGY
#   MP4(SDTQ) CORRELATION ENERGY
#
#   The total electronic energy [H] and correlation energy component [H]
#   for the full MP4 level of theory. Quantity :qcvar:`MP4 CORRELATION
#   ENERGY <MP4CORRELATIONENERGY>` / :qcvar:`MP4(SDTQ) CORRELATION ENERGY
#   <MP4(SDTQ)CORRELATIONENERGY>` is left-hand term in Eq. :eq:`MP4terms`.
#
#   .. math:: E_{\text{MP4}} = E_{\text{MP4(SDQ)}} + E_{\text{MP4(T)}}
#      :label: MP4terms
#
#.. qcvar:: MPn TOTAL ENERGY
#   MPn CORRELATION ENERGY
#
#   The total electronic energy [H] and correlation energy component [H]
#   for the labeled |MollerPlesset| perturbation theory level.
#   *n* is MP perturbation order.

qcvardefs['MP3 TOTAL ENERGY'] = {
    'units': 'Eh',
    'glossary': r"""
   The total electronic energy for 3-order Perturbation theory.
"""
}

qcvardefs['MP3 CORRELATION ENERGY'] = {
    'units': 'Eh',
    'glossary': r"""
   The correlation energy component for 3-order Perturbation theory.
"""
}

qcvardefs['MP3 CORRECTION ENERGY'] = {
    'units': 'Eh',
    'glossary': r"""
   The correlation energy difference between 2nd and 3-order Perturbation theory.
"""
}

qcvardefs['MP3 SINGLES ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The singles portion of the MP3 correlation energy.
   Zero except in ROHF.
   :math:`E_{\text{S}}` in Eq. :eq:`MP3corl`.
"""
}

qcvardefs['MP3 DOUBLES ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The doubles portion of the MP3 correlation energy
   including same-spin and opposite-spin correlations.
   :math:`E_{\text{D}}` in Eq. :eq:`MP3corl`.
"""
}

qcvardefs['MP3 SAME-SPIN CORRELATION ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The unscaled portion of the MP3 correlation energy
   from same-spin or triplet doubles correlations.

   canonical_corl(os_scale=1, ss_scale=1) = singles + os_scale * (tot_corl - ss_corl) + ss_scale * ss_corl
   :math:`E_{\text{SS}}` in Eq. :eq:`MP3corl`.
"""
}

qcvardefs['MP3 OPPOSITE-SPIN CORRELATION ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The unscaled portion of the MP3 correlation energy
   from opposite-spin or singlet doubles correlations.
   :math:`E_{\text{OS}}` in Eq. :eq:`MP3corl`.
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
    'units': 'Eh/a0/a0',
    'dimension': '({nat}, 3)',
    'glossary': """
   The total electronic gradient
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

   .. math::
      :nowrap:
      :label: CCSDcorl

         \begin{align*}
#            E_{\text{CCSDcorl}} & = E_{\text{S}} + E_{\text{SS}} + E_{\text{OS}} \\
#                               & = E_{\text{S}} + E_{\text{D}} \\
         \end{align*}

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

   .. math:: 
      :nowrap:
      :label: CCSDTcorl

         \begin{align*}
#            E_{\text{CCSDTcorl}} & = E_{\text{S}} + E_{\text{SS}} + E_{\text{OS}} \\
#                               & = E_{\text{S}} + E_{\text{D}} \\
         \end{align*}

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

   .. math:: 
      :nowrap:
      :label: CCSDTQcorl

         \begin{align*}
#            E_{\text{CCSDTQcorl}} & = E_{\text{S}} + E_{\text{SS}} + E_{\text{OS}} \\
#                               & = E_{\text{S}} + E_{\text{D}} \\
         \end{align*}

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

qcvardefs['CCSD(T) TOTAL ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The total electronic energy
   for the coupled cluster singles and doubles plus perturbative triples level of theory.
"""
}

qcvardefs['CCSD(T) CORRELATION ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The correlation energy component
   for the coupled cluster singles and doubles plus perturbative triples level of theory.
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

qcvardefs['CCSDT TOTAL ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The total electronic energy 
   for the coupled cluster singles and doubles plus perturbative triples level of theory.
"""
}

qcvardefs['CCSDT CORRELATION ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The correlation energy component
   for the coupled cluster singles and doubles plus perturbative triples level of theory.
"""
}

qcvardefs['CCSD[T] TOTAL ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The total electronic energy
   for the coupled cluster singles and doubles plus bracket perturbative triples level of theory.
"""
}

qcvardefs['CCSD[T] CORRELATION ENERGY'] = {
    'units':
    'Eh',
    'glossary':
    r"""
   The correlation energy component
   for the coupled cluster singles and doubles plus bracket perturbative triples level of theory.
"""
}

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

   .. math:: E_{NN} = \sum_{i, j<i}^{N_{atom}}\frac{Z_i Z_j}{|\mathbf{R}_i - \mathbf{R}_j|}
      :label: ENN
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
#   :qcvar:`SAPT TOTAL ENERGY <SAPTTOTALENERGY>`.
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
   The :qcvar:`CORRELATION ENERGY` variables from subsequent stages of a
   calculation are often the corresponding :qcvar:`TOTAL ENERGY`
   variables less this quantity. Constructed from Eq. :eq:`SCFterms`,
   where this quantity is :math:`E_{\text{SCF}}`.

   .. math::
      :nowrap:
      :label: SCFterms

         \begin{align*}
            E_{\text{SCF}} & = E_{NN} + E_{1e^-} + E_{2e^-} + E_{xc} + E_{\text{-D}} \\
                           & = E_{\text{FCTL/HF}} + E_{\text{-D}}
         \end{align*}

   Unless the method includes a dispersion correction, this quantity is equal to :qcvar:`HF TOTAL ENERGY <HFTOTALENERGY>` (for HF) or
   :qcvar:`DFT FUNCTIONAL TOTAL ENERGY <DFTFUNCTIONALTOTALENERGY>` (for
   DFT). Unless the method is a DFT double-hybrid, this quantity is equal
   to :qcvar:`DFT TOTAL ENERGY <DFTTOTALENERGY>`.
"""
}

qcvardefs['SCF TOTAL GRADIENT'] = {
    'units': 'Eh/a0/a0',
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
    'units':
    'Eh',
    'glossary':
    r"""
The number of alpha electrons
"""
}

qcvardefs['N BETA ELECTRONS'] = {
    'units':
    'Eh',
    'glossary':
    r"""
  The number of beta electrons 
"""
}

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
define_scf_qcvars('B3LYP')
define_scf_qcvars('B3LYP5')
define_scf_qcvars('SCF',
                  is_dft=False,
                  extra=' This is a generic HF/DFT quantity and not necessarily aligned across different calcs.')
define_scf_qcvars('B2PLYP', is_dh=True)
define_scf_qcvars('DSD-PBEP86', is_dh=True)
define_scf_qcvars('WPBE', is_dh=False)
define_scf_qcvars('WB97', is_dh=False)
define_scf_qcvars('PBE', is_dh=False)
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
