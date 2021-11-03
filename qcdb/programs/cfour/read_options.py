from ...keywords import AliasKeyword, Keyword, Keywords, parsers


def load_cfour_keywords(options: Keywords) -> None:
    options.add(
        "cfour",
        Keyword(
            keyword="TRANSLATE_PSI4",
            default=True,
            validator=parsers.boolean,
            glossary="Will translate Psi4 options to CFOUR counterparts",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="ABCDTYPE",
            default="STANDARD",
            validator=parsers.enum("STANDARD AOBASIS"),
            glossary="Specifies the way the ⟨ab||cd⟩ molecular orbital integrals are handled in post-MP2 calculations.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="ACTIVE_ORBI",
            default="",
            validator= lambda x: str(x),
            glossary="Specifies the active orbitals used in a TCSCF calculation and has to be used in combination with the keyword CORE_ORBITALS.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="ANHARMONIC",
            default="OFF",
            validator=parsers.enum("CUBIC VPT2 FULLQUARTIC VIBROT OFF"),
            glossary="Specifies treatment of anharmonc effects by calculating cubic and/or quartic force fields.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="ANH_ALGORITHM",
            default="STANDARD",
            validator=parsers.enum("STANDARD PARALLEL"),
            glossary="Specifies which algorithm is used for ANHARM=VIBROT, ANHARM=VPT2, and ANHARM=FULLQUARTIC calculations.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="ANH_POINTS",
            default=1,
            validator=parsers.positive_integer,
            glossary="Specifies the number of points used for the numerical differentiation in the computation of cubic and quartic force constants based on analytically evaluated Hessians.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="ANH_DERIVATIVES",
            default="SECOND",
            validator=parsers.enum("FIRST SECOND"),
            glossary="Specifies whether the anharmonic force field is calculated using analytic gradients or analyic Hessians.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="ANH_STEPSIZE",
            default=50000,
            validator=parsers.positive_integer,
            glossary="Controls the stepsize used in anharmonic force field calculations. ",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="ANH_SYMMETRY",
            default="ABELIAN",
            validator=parsers.enum("ABELIAN NONABELIAN"),
            glossary="Specifies whether nonabelian symmetry is to be exploited in determining displacements for ANHARM=VIBROT or VPT2 calculations. ",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="AO_LADDERS",
            default="SINGLEPASS",
            validator=parsers.enum("MULTIPASS SINGLEPASS"),
            glossary="Can be used to control the algorithm used by CFOUR when terms involving ⟨ab||cd⟩ molecular orbital integrals are calculated in the atomic orbital basis.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="AV_SCF",
            default=False,
            validator=parsers.boolean,
            glossary="Experimental Use! Requests and averaged SCF over two states.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="BASIS",
            default="SPECIAL",
            validator=parsers.enum("STO-3G 3-21G 4-31G 6-31G 6-31G* 6-31G** 6-311G 6-311G* 6-311G** DZ DZP TZ TZP TZ2P PVDZ PVTZ \
                        PVQZ PV5Z PV6Z PCVDZ PCVTZ PCVQZ PCV5Z PCV6Z AUG-PVDZ AUG-PVTZ AUG-PVTZ AUG-PVQZ AUG-PV5Z \
                        AUG-PV6Z D-AUG-PVDZ D-AUG-PVTZ D-AUG-PVQZ D-AUG-PV5Z D-AUG-PV6Z cc-pVDZ cc-pVTZ cc-pVQZ \
                        cc-pV5Z cc-pV6Z cc-pCVDZ cc-pCVTZ cc-pCVQZ cc-pCV5Z cc-pCV6Z PWCVDZ PWCVTZ PWCVQZ PWCV5Z \
                        PWCV6Z PwCVDZ PwCVTZ PwCVQZ PwCV5Z PwCV6Z svp dzp tzp tzp2p qz2p pz3d2f 13s9p4d3f WMR ANO0 \
                        ANO1 ANO2 EVEN_TEMPERED SPECIAL"),
            glossary="Specifies the AO basis used in the calculation.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="BRUCK_CONV",
            default=4,
            validator=parsers.positive_integer,
            glossary="Specifies the convergence criterion in Brueckner based CC calculations. ",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="BRUECKNER",
            default=False,
            validator=parsers.boolean,
            glossary="Specifies whether Brueckner orbitals are to be determined for the specified CC method.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="CALC_LEVEL",
            default="SCF",
            validator=lambda x: str(x).upper(), #parsers.enum("SCF HF MBPT(2) MP2 MBPT(3) MP3 SDQ-MBPT(4) SDQ-MP4 MBPT(4) MP4 CCD CCSD CCSD(T) CCSDT-1 \
                        #CCSDT-1b CCSDT-2 CCSDT-3 CCSDT-4 CCSDT CC2 CC3 QCISD QCISD(T) CID CISD UCC(4) B-CCD CCSD[T] CCSD+T(CCSD)"),
            glossary="Defines the level of calculation to be performed. ",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="CACHE_RECS",
            default=10,
            validator=parsers.positive_integer,
            glossary="The number of records held in the i/o cache used by the post-SCF programs.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="CC_CONV",
            default=7,
            validator=parsers.positive_integer,
            glossary="Specifies the convergence criterion for the CC amplitude equations.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="CC_EXPORDER",
            default=5,
            validator=parsers.positive_integer,
            glossary="Specifies the maximum number of expansion vectors used in the iterative subspace to enhance convergence in the solution of the CC equations.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="CC_EXTRAPOLATION",
            default="DIIS",
            validator=parsers.enum("RLE DIIS NOJACOBI OFF"),
            glossary="Specifies the type of convergence acceleration used to solve the CC equations.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="CC_MAXCYC",
            default=50,
            validator=parsers.positive_integer,
            glossary="Specifies the maximum number of iterations in solving the CC amplitude equations.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="CC_PROGRAM",
            default="VCC",
            validator=parsers.enum("VCC ECC NCC MRCC SACC EXTERNAL"),
            glossary="Specifies which CC program is used.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="CHARGE",
            default=0,
            validator=parsers.integer,
            glossary="Specifies the molecular charge.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="CHOLESKY",
            default=False,
            validator=parsers.boolean,
            glossary="Specifies whether Cholesky decomposition of two-electron integrals is used.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="CIS_CONV",
            default=5,
            validator=parsers.positive_integer,
            glossary="Specifies the convergence threshold (as 10^−N) for CIS calculations.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="CONTINUUM",
            default="NONE",
            validator=parsers.enum("NONE VIRTUAL DVIRTUAL OCCUPIED DOCCUPIED"),
            glossary="Signifies that one or more 'continuum' orbitals should be added to the calculation.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="CONTRACTION",
            default="GENERAL",
            validator=parsers.enum("SEGMENTED GENERAL UNCONTRACTED"),
            glossary="specifies the contraction scheme used by the integral and integral derivative program.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="CONVERGENCE",
            default=4,
            validator=parsers.positive_integer,
            glossary="Specifies convergence criterion for geometry optimization.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="COORDINATES",
            default="INTERNAL",
            validator=parsers.enum("INTERNAL CARTESIAN XYZINT"),
            glossary="Specifies the type of coordinates used in the input file ZMAT.",
        ),
    )

    options.add(
            "cfour",
            Keyword(
                keyword="CORE_ORBITALS",
                default="",
                validator= lambda x: str(x),
                glossary="Specifies the core orbitals used in a TCSCF calculation and has to be used in combination with the keyword ACTIVE_ORBI.",
            ),
        )

    options.add(
        "cfour",
        Keyword(
            keyword="CPHF_CONVER",
            default=12,
            validator=parsers.positive_integer,
            glossary="Specifies the convergence criterion for the iterative solution of the CPHF and Z-vector equations.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="CPHF_MAXCYC",
            default=64,
            validator=parsers.positive_integer,
            glossary="Specifies the maximum number of cycles allowed for the solution of the CPHF- and/or Z-vector equations.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="CURVILINEAR",
            default=True,
            validator=parsers.boolean,
            glossary="Specifies whether or not Hessian matrix is transformed (nonlinearly) to curvilinear internal coordinates.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="DBOC",
            default=False,
            validator=parsers.boolean,
            glossary="Specifies whether the diagonal Born-Oppenheimer correction (DBOC) to the energy is evaluated.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="DCT",
            default=False,
            validator=parsers.boolean,
            glossary="Specifies whether the Dipole Coupling Tensor (DCT) is calculated.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="DERIV_LEVEL",
            default="ZERO",
            validator=parsers.enum("ZERO FIRST SECOND"),
            glossary="Specifies whether or not energy derivatives are to be calculated and if so whether first or second derivatives are computed.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="DIFF_TYPE",
            default="RELAXED",
            validator=parsers.enum("RELAXED UNRELAXED"),
            glossary="Specifies whether orbital-relaxed or orbital-unrelaxed derivatives are computed in the CC calculation.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="DROPMO",
            default=[],
            validator=lambda x: x,  # should be int list
            glossary="Specifies which molecular orbitals will be dropped from the post-SCF calculation. ",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="ECP",
            default=False,
            validator=parsers.boolean,
            glossary="Specifies whether effective core potentials (pseudopotentials) are used.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="EIGENVECTOR",
            default=1,
            validator=parsers.positive_integer,
            glossary="Specifies which eigenvector of the totally symmetric part of the block-factored Hessian is to be followed uphill in a transition state search.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="EL_ANHARM",
            default=False,
            validator=parsers.boolean,
            glossary="Experimental use! Requests the evaluation of electrical anharmonicities.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="ESTATE_CONV",
            default=5,
            validator=parsers.positive_integer,
            glossary="Specifies the threshold used in converging CC-LR/EOM-CC calculations. ",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="EOM_NONIT",
            default=False,
            validator=parsers.boolean,
            glossary="Controls whether non-iterative triples corrections are applied after various types of EOM-CCSD calculation.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="EOM_NSTATES",
            default="DAVIDSON",
            validator=parsers.enum("DAVIDSON MULTIROOT"),
            glossary="For experimental use only. Selects the iterative diagonalization algorithm for the EOMEE calculations. ",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="ESTATE_MAXCYC",
            default=20,
            validator=parsers.positive_integer,
            glossary="The maximum number of expansion vectors used in the solution of EOMCC equations.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="ESTATE_PROP",
            default="OFF",
            validator=parsers.enum("OFF EXPECTATION UNRELAXED RESPONSE"),
            glossary="This keyword applies only to EOM-CC calculations and specifies whether any excited or ionized state one-electron properties are to be calculated.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="ESTATE_SYM",
            default="",
            validator=lambda x: str(x),
            glossary="Specifies the number of excited states which are to be determined in each irreducible representation of the computational subgroup.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="ESTATE_TRANS",
            default="OFF",
            validator=parsers.enum("OFF EXPECTATION"),
            glossary="Specifies whether just the excitation energies or in addition transition moments are calculated. ",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="EVAL_HESS",
            default=0,
            validator=parsers.nonnegative_integer,
            glossary="Tells the program, in the course of a geometry optimization, to calculate the Hessian explicitly every N cycles.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="EXCITATION",
            default=0,
            validator=parsers.nonnegative_integer,
            glossary="Specifies in CC calculations using mrcc the excitation level if the calculation level has been chosen as CC(n), CI(n), or CCn(n).''",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="EXCITE",
            default="NONE",
            validator=parsers.enum("NONE EOMEE EOMIP EOMEA"),
            glossary="Specifies the type of EOM-CC/LR-CC treatment to be performed. ",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="EXTERN_POT",
            default=False,
            validator=parsers.boolean,
            glossary="Allows the use of external pointcharges.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="FC_FIELD",
            default=0,
            validator=parsers.integer,
            glossary="Specifies the strength of a Fermi-Contact pertubation as required for finite-field calculations of spin densities and the FC contributions to indirect spin-spin coupling constants.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="FD_CALCTYPE",
            default="GRADONLY",
            validator=parsers.enum("GRADONLY ENERONLY MIXED"),
            glossary="Specifies the algorithm used to compute the harmonic force constants in finite-difference calculations.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="FD_IRREPS",
            default="",
            validator=lambda x: str(x),
            glossary="Requests that only vibrational frequencies of certain symmetry types are evaluated in a VIBRATION=FINDIF calculation.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="FD_PROJECT",
            default=True,
            validator=parsers.boolean,
            glossary="Specifies whether or not rotational degrees of freedoms are projected out from the symmetry-adapted coordinates in a finite difference calculations.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="FD_STEPSIZE",
            default=5,
            validator=parsers.positive_integer,
            glossary="Specifies the step length in mass-weighted coordinates (in 10**(-4) amu**0.5 bohr) used in generating the force constant matrix by finite difference of Cartesian gradients.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="FD_USEGROUP",
            default="FULL",
            validator=parsers.enum("FULL COMP"),
            glossary="In finite difference calculations using the FINDIF option, this keyword specifies the point group to be used in generating the symmetry-adapted vibrational coordinates.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="FILE_RECSIZ",
            default=2048,
            validator=parsers.positive_integer,
            glossary="This specifies the physical length (in integer words) of the records used in the word-addressable direct access files used by CFOUR.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="FILE_STRIPE",
            default="0/0/0/0/0",
            validator=lambda x: str(x),
            glossary="This option allows the splitting of files.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="FINITE_PERTURBATION",
            default=0,
            validator=parsers.integer,
            glossary="Specifies the field strength for a perturbation",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="FOCK",
            default="AO",
            validator=parsers.enum("PK AO DIRECT"),
            glossary="This option is used to control the algorithm used for construction of the Fock matrix in SCF calculations.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="FREQ_ALGORITHM",
            default="STANDARD",
            validator=parsers.enum("STANDARD PARALLEL"),
            glossary="Experimental use.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="FROZEN_CORE",
            default=False,
            validator=parsers.boolean,
            glossary="Specifies whether in the correlation treatment all electron or only the valence electrons are considered. ",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="FROZEN_VIRT",
            default=False,
            validator=parsers.boolean,
            glossary="Specifies whether in the correlation treatment all virtual orbitals or only a subset of virtual orbitals are used.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="GAMMA_ABCD",
            default="DISK",
            validator=parsers.enum("DISK DIRECT"),
            glossary="Used to control the handling and storage of two-particle density matrix elements with four virtual indices.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="GENBAS_1",
            default="",
            validator=lambda x: str(x),
            glossary="This keyword applies only to Hydrogen and Helium atoms and specifies the number of contracted Gaussian functions per shell.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="GENBAS_2",
            default="",
            validator=lambda x: str(x),
            glossary="This keyword performs the same function as GENBAS_1, but applies to second-row atoms.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="GENBAS_3",
            default="",
            validator=lambda x: str(x),
            glossary="This keyword performs the same function as GENBAS_1, but applies to third-row atoms.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="GENBAS_4",
            default="",
            validator=lambda x: str(x),
            glossary="This keyword performs the same function as GENBAS_1, but applies to fourth-row atoms.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="GEO_CONV",
            default=5,
            validator=parsers.positive_integer,
            glossary="Specifies the convergence criterion for geometry optimization.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="GEO_MAXSTEP",
            default=300,
            validator=parsers.positive_integer,
            glossary="Specifies largest step (in millbohr) which is allowed in geometry optimizations.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="GEO_METHOD",
            default="SINGLE_POINT",
            validator=parsers.enum("NR RFA TS MANR SINGLE_POINT ENERONLY"),
            glossary="Specifies the used geometry optimization methods.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="GEO_MAXCYC",
            default=50,
            validator=parsers.positive_integer,
            glossary="Specifies convergence criterion for geometry optimization.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="GIAO",
            default=False,
            validator=parsers.boolean,
            glossary="Specifies whether gauge-including atomic orbitals (GIAOs, London atomic orbitals) are used.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="GRID",
            default="OFF",
            validator=parsers.enum("OFF CARTESIAN INTERNAL QUADRATURE"),
            glossary="Keyword used to control type of grid calculation.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="GUESS",
            default="MOREAD",
            validator=parsers.enum("MOREAD CORE"),
            glossary="Where the initial SCF eigenvectors are read from.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="HBAR",
            default=False,
            validator=parsers.boolean,
            glossary="This keyword determines which action is taken by the linear response program.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="HFSTABILITY",
            default="OFF",
            validator=parsers.enum("OFF ON FOLLOW"),
            glossary="Control analysis of the stability of RHF, ROHF and UHF wavefunctions, as well as a possible search for a lower SCF solution.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="INCORE",
            default="OFF",
            validator=parsers.enum("OFF ALL PARTIAL"),
            glossary="This keyword can be used to significantly reduce disk i/o, and should be implemented very soon. ",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="INPUT_MRCC",
            default=True,
            validator=parsers.boolean,
            glossary="Specifies whether an input for mrcc is written.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="INTEGRALS",
            default="VMOL",
            validator=parsers.enum("VMOL ARGOS"),
            glossary="This keyword defines what type of integral input will be written by JODA.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="JODA_PRINT",
            default=0,
            validator=parsers.nonnegative_integer,
            glossary="Controls amount of debug printing performed by Joda.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="LINEQ_CONV",
            default=7,
            validator=parsers.positive_integer,
            glossary="Convergence threshold for linear equations controlled by LINEQ_TYPE.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="LINEQ_TYPE",
            default="DIIS",
            validator=parsers.enum("POPLE DIIS"),
            glossary="Determines the algorithm used to solve linear equations.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="LINEQ_MAXCY",
            default=50,
            validator=parsers.positive_integer,
            glossary="The maximum number of iterations in all linear CC equations.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="LOCK_ORBOCC",
            default=False,
            validator=parsers.boolean,
            glossary="This keyword is used by the SCF program to determine if the orbital occupancy (by symmetry block) is allowed to change in the course of the calculation.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="MAXSTEP",
            default=300,
            validator=parsers.positive_integer,
            glossary="Specifies largest step (in millibohr) which is allowed in geometry optimizations.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="MEMORY_SIZE",
            default=100000000,
            validator=parsers.positive_integer,  #parsers.parse_memory,
            glossary="Specifies the amount of core memory used.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="MEM_UNIT",
            default="INTEGERWORDS",
            validator=parsers.enum("INTEGERWORDS KB MB GB TB"),
            glossary="Specifies the units in which the amount of requested core memory is given.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="METHOD",
            default="SINGLE_POINT",
            validator=parsers.enum("NR RFA TS MANR SINGLE_POINT"),
            glossary="Specifies the geometry optimization strategy.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="MRCC",
            default="OFF",
            validator=parsers.enum("OFF MK"),
            glossary="Specifies the type of MRCC calculation.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="MULTIPLICITY",
            default=1,
            validator=parsers.positive_integer,
            glossary="Specifies the spin multiplicity.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="NACOUPLING",
            default="OFF",
            validator=parsers.enum("OFF ON NACV LVC"),
            glossary="Calculation of non-adiabatic coupling.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="NEGEVAL",
            default="ABORT",
            validator=parsers.enum("ABORT SWITCH RFA"),
            glossary="Specifies what to do if negative eigenvalues are encountered in the totally symmetric Hessian during an NR or MANR geometry-optimization search.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="NEWNORM",
            default=False,
            validator=parsers.boolean,
            glossary="All components of spherical AO’s are normalized to 1.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="NONHF",
            default=False,
            validator=parsers.boolean,
            glossary="Specifies whether the reference function used in the correlation energy calculation satisfies the (spin-orbital) HF equations or not.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="NTOP_TAMP",
            default=15,
            validator=parsers.positive_integer,
            glossary="Specifies how many t amplitudes will be printed for each spin case and excitation level.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="OCCUPATION",
            default="",
            validator=lambda x: x,  # should be list(s)
            glossary="Specifies the orbital occupancy of the reference function in terms of the occupation numbers of the orbitals and their irreducible representations.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="OPEN-SHELL",
            default="SPIN-ORBITAL",
            validator=parsers.enum("SPIN-ORBITAL SR-CC PSA-CC_FULL SR-CC_FULL TD-CC"),
            glossary="Specifies which kind of open-shell CC treatment is employed.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="OPT_MAXCYC",
            default=50,
            validator=parsers.positive_integer,
            glossary="Specifies the maximum allowed number of geometry optimization cycles.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="ORBITALS",
            default="STANDARD",
            validator=parsers.mixedenum("STANDARD SEMICANONICAL 0 1"),
            glossary="Specifies the type of molecular orbitals used in post-HF calculations.",
        ),
    )

    #PARALLEL
    #PARA_PRINT
    #PARA_INT

    options.add(
        "cfour",
        Keyword(
            keyword="PERT_ORB",
            default="STANDARD",
            validator=parsers.mixedenum("STANDARD CANONICAL IJ_CANONICAL 0 1"),
            glossary="Specifies the type of perturbed orbitals used in energy derivative calculations.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="PNO_USE",
            default=False,
            validator=parsers.boolean,
            glossary="Specifies whether pair natural orbitals (PNOs) are used in a local MP2 treament.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="PNO_THRESHOLD",
            default=7,
            validator=parsers.positive_integer,
            glossary="Specifies the cutoff for the selection of PNOs are used in a local MP2 treament.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="POINTS",
            default="DOUBLE",
            validator=parsers.enum("SINGLE DOUBLE"),
            glossary="Specifies either single or double sided numerical differentiation in the finite difference evaluation of the Hessian.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="PRINT",
            default=0,
            validator=parsers.nonnegative_integer,
            glossary="Controls the amount of printing in the energy and energy derivative calculation programs.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="PROPS",
            default="OFF",
            validator=parsers.enum("OFF FIRST_ORDER SECOND_ORDER NMR HYPERPOL DYN_HYP SHG OPT_REC VERDET"),
            glossary="Specifies whether and which molecular property is calculated.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="PROP_INTEGRAL",
            default="INTERNAL",
            validator=parsers.enum("INTERNAL EXTERNAL"),
            glossary="Allows storage of property integrals computed in xvdint on internal files.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="QC_ALG",
            default="FLM",
            validator=parsers.enum("FLM AH"),
            glossary="Specify which algorithm to use for QCSCF.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="QC_MAXCYC",
            default=100,
            validator=parsers.positive_integer,
            glossary="Specify the maximum number of cycles in a QCSCF calculations.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="QC_MAXSCFCYC",
            default=20,
            validator=parsers.positive_integer,
            glossary="Specify the maximum number of standard SCF cycles before the QCSCF optimizazion.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="QC_NOISE",
            default=False,
            validator=parsers.boolean,
            glossary="Specify whether to add a random perturbation to the MO rotation gradient at the beginning of the second-order optimization.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="QC_RTRUST",
            default=0,
            validator=parsers.nonnegative_integer,
            glossary="Trust radius for the second-order SCF optimization.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="QC_START",
            default=-1,
            validator=parsers.integer,
            glossary="Specify the convergence threshold for the preliminary SCF iterations before switching to quadratically convergent ones.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="QRHF_GENERAL",
            default="",
            validator=lambda x: str(x),
            glossary="The presence of this keyword specifies that a QRHF based CC calculation, or alternatively, an SCF calculation that uses the QRHFGUES option, is to be performed.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="QRHFGUES",
            default=False,
            validator=parsers.boolean,
            glossary="If this keyword is set to ON, then the QRHF orbitals specified by the QRHF_GENERAL, QRHF_ORBITAL and QRHF_SPIN keywords are used as a starting guess for a restarted SCF procedure.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="QRHF_ORBITAL",
            default=1,
            validator=lambda x: str(x),
            glossary="The value of this keyword gives the offset with respect to the default orbital for the orbital which will be depopulated (or populated) in QRHF-CC calculations.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="QRHF_SPIN",
            default=1,
            validator=parsers.intenum("1 2"),  #parsers.integer,
            glossary="Specifies the spin of the electrons modified by the QRHF_GENERAL and QRHF_ORBITAL keywords, where a value of 1 means alpha spin, while 2 corresponds to a beta electron.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="RAMAN_INT",
            default="OFF",
            validator=parsers.enum("ON DYN OFF"),
            glossary="ON (=1) requests a calculation of Raman intensities based on the geometrical derivatives of the static polarizability tensor, while DYN (=2) requests a calculation of Raman intensities based on the derivatives of the dynamical polarizability tensor.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="RAMAN_ORB",
            default="UNRELAXED",
            validator=parsers.enum("RELAXED UNRELAXED"),
            glossary="specifies whether Raman intensities are calculated with orbital relaxation with respect to the electric field perturbation.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="RDO",
            default=True,
            validator=parsers.boolean,
            glossary="Specifies whether or not relaxed density natural orbitals are to be computed.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="REFERENCE",
            default="RHF",
            validator=parsers.enum("RHF UHF ROHF TCSCF CASSCF"),
            glossary="Specifies the type of SCF calculation to be performed.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="RELATIVISTIC",
            default="OFF",
            validator=parsers.enum("OFF MVD1 MVd2 DPT2 SF-DPT4 DPT4 SF-DPT6 SFREE X2C1E DPT"),
            glossary="Specifies the treatment of relativistic effects.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="RELAX_DENS",
            default=False,
            validator=parsers.boolean,
            glossary="Specifies whether the relaxed density matrix is computed for correlated wave functions.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="RES_RAMAN",
            default=False,
            validator=parsers.boolean,
            glossary="This option can be used to convert an analytically calculated gradient vector to a particular normal coordinate representation.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="RESTART_CC",
            default=False,
            validator=parsers.boolean,
            glossary="Offers the possibilty to restart a CC calculation which stopped for various reasons, e.g. time limit, in the correlation part.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="ROT_EVEC",
            default=0,
            validator=parsers.nonnegative_integer,
            glossary="Specifies which eigenvector of the orbital rotation Hessian is to be used to rotate the original SCF orbitals.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="SACC_ANSATZ",
            default="COS-CC",
            validator=parsers.enum("COS-CC NORMAL"),
            glossary="Specifies the type of the spin-adapted open-shell CC ansatz for the doublet case.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="SACC_CALC",
            default="DENSITY",
            validator=parsers.enum("DENSITY PROPERTIES"),
            glossary="Defines the type or the stage of the SACC calculation.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="SACC_ORBS",
            default="STANDARD",
            validator=parsers.enum("STANDARD SACC_ORBS"),
            glossary="Specifies the type of orbitals to be used for SACC calculations.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="SACC_PROP",
            default="RELAXED",
            validator=parsers.enum("RELAXED UNRELAXED"),
            glossary="Specifies whether orbital-relaxed or unrelaxed properties will be computed using the SACC program.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="SAVE_INTS",
            default=False,
            validator=parsers.boolean,
            glossary="Tells CFOUR whether to delete large files (AO integrals and MOINTS file for now) when they are no longer needed.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="SCALE_ON",
            default="MAG(S)",
            validator=parsers.enum("MAG(S) MAX(S)"),
            glossary="Controls whether step scaling is based on the absolute step length or the largest individual step in the internal coordinate space",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="SCF_CONV",
            default=7,
            validator=parsers.positive_integer,
            glossary="Specifies the convergence criterion for the HF-SCF equations.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="SCF_DAMPING",
            default=1000,
            validator=parsers.nonnegative_integer,
            glossary="Controls the damping (in the first iterations (specified by SCF_EXPSTART via D(new) = D(old) + X/1000 * [D(new) - D(old)] with X as the value specified by the keyword.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="SCF_EXPORDER",
            default=6,
            validator=parsers.positive_integer,
            glossary="Specifies the number (N) of density matrices to be used in the DIIS convergence acceleration procedure.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="SCF_EXPSTART",
            default=8,
            validator=parsers.positive_integer,
            glossary="Specifies the first iteration in which the DIIS convergence acceleration procedure is applied.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="SCF_EXTRAPOLATION",
            default=True,
            validator=parsers.boolean,
            glossary="Specifies whether or not the DIIS extrapolation is used to accelerate convergence of the SCF procedure.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="SCF_MAXCYC",
            default=150,
            validator=parsers.positive_integer,
            glossary="Specifies the maximum number of SCF iterations.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="SCF_PROG",
            default="SCF",
            validator=parsers.enum("SCF QCSCF DQCSCF"),
            glossary="Specifies the maximum number of SCF iterations.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="SD_FIELD",
            default=0,
            validator=parsers.integer,
            glossary="Specifies the strength of a spin-dipole pertubation as required for finite-field calculations of the SD contributions to indirect spin-spin coupling constants.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="SOPERT",
            default="OFF",
            validator=parsers.enum("OFF MKMRCC EMRCCSO"),
            glossary="Experimental use. Perturbative treatment of spin-orbit splittings in dublett-pi states via multireference coupled-cluster theory.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="SPHERICAL",
            default=True,
            validator=parsers.boolean,
            glossary="Specifies whether spherical harmonic (5d, 7f, 9g, etc.) or Cartesian (6d, 10f, 15g, etc.) basis functions are to be used.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="SPIN_FLIP",
            default=False,
            validator=parsers.boolean,
            glossary="Controls whether excitation energy calculations allow for a spin flip which changes the $M_s$ quantum number.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="SPIN_ORBIT",
            default="OFF",
            validator=parsers.enum("OFF ON SA-SOMF SOMF"),
            glossary="Experimental use. Requests calculation of one-electron spin-orbit integrals.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="SPIN_SCAL",
            default=False,
            validator=parsers.boolean,
            glossary="Requests the spin-component scaled variant of the MP2 approach.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="SPINROTATION",
            default=False,
            validator=parsers.boolean,
            glossary="Specifies whether nuclear spin-rotation tensors are computed within a NMR chemical shift calculation or not.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="SUBGROUP",
            default="DEFAULT",
            validator=parsers.enum("DEFAULT C1 C2 CS CI C2V C2H D2 D2H OFF"),
            glossary="Specifies an Abelian subgroup to be used in a calculation.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="SUBGRPAXIS",
            default="DEFAULT",
            validator=parsers.enum("DEFAULT C1 C2 CS CI C2V C2H D2 D2H OFF"),
            glossary="Specifies an Abelian subgroup to be used in a calculation.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="SYM_CHECK",
            default=True,
            validator=parsers.boolean,
            glossary="In principle can be used to force the SCF to converge a solution for which the density matrix transforms as the totally symmetric representation of the point group.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="SYMMETRY",
            default=True,
            validator=parsers.boolean,
            glossary="Specifies what subgroup of the full point group is to be used in the energy and/or gradient calculation",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="TAMP_SUM",
            default=5,
            validator=parsers.nonnegative_integer,
            glossary="Specifies how often the largest t amplitudes are to be printed.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="THERMOCHEMISTRY",
            default="ON",
            validator=parsers.enum("OFF ON VERBOSE"),
            glossary="Specifies whether to calculate finite-temperature thermodynamic corrections after a frequency calculation.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="TOL_CHOLESKY",
            default=4,
            validator=parsers.integer,
            glossary="Specifies the threshold to be used in Cholesky decomposition of the two-electron integrals.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="TRANS_INV",
            default="USE",
            validator=parsers.enum("USE IGNORE"),
            glossary="Specifies whether or not translational invariance is exploited in geometrical derivative calculations.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="TREAT_PERT",
            default="SIMULTANEOUS",
            validator=parsers.enum("SIMULTANEOUS SEQUENTIAL"),
            glossary="Specifies whether in a correlated NMR chemical shift calculations all perturbations are treated at once or sequentially.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="T3_EXTRAPOL",
            default=False,
            validator=parsers.boolean,
            glossary="Specifies whether the T3 amplitudes are included.",
        ),
    )
    options.add(
        "cfour",
        Keyword(
            keyword="UIJ_THRESHOLD",
            default=25,
            validator=parsers.positive_integer,
            glossary="Specifies the threshold value (given as an integer) for the treatment of CPHF coefficients in second derivative calculations using perturbed canonical orbitals.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="UNITS",
            default="ANGSTROM",
            validator=parsers.enum("ANGSTROM BOHR"),
            glossary="Specifies the units used for molecular geometry input.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="UPDATE_HESSIAN",
            default=True,
            validator=parsers.boolean,
            glossary="Specifies whether or not the Hessian update is carried out.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="VIBRATION",
            default="NO",
            validator=parsers.enum("NO ANALYTIC FINDIF EXACT"),
            glossary="Specifies whether (harmonic) vibrational frequencies are calculated or not.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="VNATORB",
            default=False,
            validator=parsers.boolean,
            glossary="This keyword specifies whether virtual natural orbitals are to be used or not.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="VTRAN",
            default="FULL/PARTIAL",
            validator=parsers.enum("FULL/PARTIAL FULL PARTIAL"),
            glossary="This keyword defines what type of integral transformation is to be performed in the program VTRAN.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="XFIELD",
            default=0,
            validator=parsers.nonnegative_integer,  #float?
            glossary="Specifies the X-component of an external electric field. ",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="XFORM_TOL",
            default=11,
            validator=parsers.positive_integer,
            glossary="The tolerance for storing transformed integrals.",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="YFIELD",
            default=0,
            validator=parsers.nonnegative_integer,
            glossary="Specifies the Y-component of an external electric field. ",
        ),
    )

    options.add(
        "cfour",
        Keyword(
            keyword="ZFIELD",
            default=0,
            validator=parsers.nonnegative_integer,
            glossary="specifies the Z-component of an external electric field.",
        ),
    )

