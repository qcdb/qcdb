from ...keywords import Keyword, Keywords, parsers


def load_gamess_keywords(options: Keywords) -> None:

    # $SYSTEM

    options.add(
        'gamess',
        Keyword(
            keyword='system__mwords',
            default=400,
            validator=parsers.positive_integer,
            glossary=
            """Maximum replicated memory which your job can use, on every core. Units of 1000^2 words (not 1024^2).""")
    )

    options.add(
        'gamess',
        Keyword(
            keyword='system__memddi',
            default=0,
            validator=parsers.nonnegative_integer,
            glossary=
            """Grand total memory needed for the distributed data interface (DDI) storage, given in units of 1000^2 words.""")
    )

    # $CONTRL

    options.add(
        'gamess',
        Keyword(keyword='contrl__scftyp',
                default='rhf',
                validator=parsers.enum("RHF UHF ROHF GVB MCSCF NONE"),
                glossary=""""""))

    options.add(
        'gamess',
        Keyword(keyword='contrl__dfttyp', default='none', validator=parsers.enum("NONE B3LYP PBE0"), glossary=""""""))

    options.add('gamess',
                Keyword(keyword='contrl__mplevl', default=2, validator=parsers.intenum("0 2"), glossary=""""""))

    options.add(
        'gamess',
        Keyword(
            keyword='contrl__cctyp',
            default='none',
            validator=parsers.enum(
                "NONE LCCD CCD CCSD CCSD(T) R-CC CR-CC CR-CCL CCSD(TQ) CR-CC(Q) EOM-CCSD CR-EOM CR-EOML IP-EOM2 IP-EOM3A EA-EOM2 EA-EOM3A"
            ),
            glossary=""""""))

    options.add(
        'gamess',
        Keyword(keyword='contrl__runtyp',
                default='energy',
                validator=parsers.enum("ENERGY GRADIENT HESSIAN"),
                glossary=""""""))

    options.add(
        'gamess',
        Keyword(keyword='contrl__cityp',
                default='none',
                validator=parsers.enum("NONE CIS SFCIS ALDET ORMAS FSOCI GENCI GUGA"),
                glossary=""""""))

    options.add(
        'gamess',
        Keyword(
            keyword='contrl__maxit',
            default=30,
            validator=parsers.positive_integer,  # psi4 maxiter
            glossary=""""""))

    options.add('gamess', Keyword(keyword='contrl__icharg', default=0, validator=parsers.integer, glossary=""""""))

    options.add('gamess',
                Keyword(keyword='contrl__mult', default=1, validator=parsers.positive_integer, glossary=""""""))

    options.add(
        'gamess',
        Keyword(keyword='contrl__coord',
                default='unique',
                validator=parsers.enum("UNIQUE HINT PRINAXIS ZMT ZMTMPC FRAGONLY CART"),
                glossary=""""""))

    options.add(
        'gamess',
        Keyword(keyword='contrl__units',
                default='angs',
                validator=parsers.enum("ANGS BOHR"),
                glossary="""Distance units, any angles must be in degrees."""))

    options.add(
        'gamess',
        Keyword(keyword='contrl__ispher',
                default=-1,
                validator=parsers.intenum("-1 0 1"),
                glossary="""Spherical Harmonics option."""))

    options.add(
        'gamess',
        Keyword(keyword='basis__gbasis',
                default='sto',
                validator=parsers.enum("STO N21 N31 N311 G3L G3LX"),
                glossary=""""""))

    options.add(
        "gamess",
        Keyword(keyword="contrl__numgrd",
                default=False,
                validator=parsers.boolean,
                glossary="""Flag to allow numerical differentiation of the energy."""))

    # $SCF

    options.add(
        'gamess',
        Keyword(
            keyword='scf__conv',
            default=1.e-6,
            validator=parsers.parse_convergence,  # psi4 d_convergence
            glossary="""SCF density convergence criteria"""))

    options.add(
        'gamess',
        Keyword(keyword='scf__ethrsh',
                default=0.5,
                validator=lambda x: float(x),
                glossary="""Energy error threshol for initiating DIIS."""))

    options.add(
        'gamess',
        Keyword(keyword='scf__dirscf',
                default=False,
                validator=parsers.boolean,
                glossary="""Do activate a direct SCF calculation?"""))

    # $MP2

    options.add(
        'gamess',
        Keyword(keyword='mp2__nacore',
                default=0,
                validator=parsers.nonnegative_integer,
                glossary="""Omits the first n occupied orbitals from the calculation. Default is chemical core."""))

    options.add(
        'gamess',
        Keyword(keyword='mp2__nbcore',
                default=0,
                validator=parsers.nonnegative_integer,
                glossary="""same as |gamess__mp2__nacore| for beta orbitals of UHF. Generally equals nacore."""))

    options.add(
        'gamess',
        Keyword(
            keyword='mp2__ospt',
            default="ZAPT",
            validator=parsers.enum("ZAPT RMP"),
            glossary="""Selects open shell spin-restricted perturbation. Only applies when gamess_mp2__scftype=ROHF."""
        ))

    # $CCINP

    options.add(
        'gamess',
        Keyword(keyword='ccinp__ncore',
                default=0,
                validator=parsers.nonnegative_integer,
                glossary="""Omits the first n occupied orbitals from the calculation. Default is chemical core."""))

    options.add(
        'gamess',
        Keyword(keyword='ccinp__iconv',
                default=7,
                validator=parsers.positive_integer,
                glossary="""Convergence criterion for the cluster amplitudes."""))

    options.add(
        'gamess',
        Keyword(
            keyword='ccinp__kmicro',
            default=6,
            validator=parsers.nonnegative_integer,
            glossary=
            """Perfor DIIS extrapolation of the open shell CCSD every n iterations. Enter 0 to avoid using the DIIS converger."""
        ))

    options.add(
        'gamess',
        Keyword(
            keyword='ccinp__maxcc',
            default=30,
            validator=parsers.positive_integer,
            glossary="""Max iterations for CCSD iterations, also ROHF left CC vector solver."""
        ))

    # $DFT

    options.add(
        'gamess',
        Keyword(
            keyword='dft__nrad',  # psi4 dft_angular_points
            default=96,
            validator=parsers.positive_integer,
            glossary="""Number of radial points in the Euler-MacLaurin quadrature."""))

    options.add(
        'gamess',
        Keyword(
            keyword='dft__nleb',  # psi4 dft_spherical_points
            default=302,
            validator=parsers.intenum("86 110 146 170 194 302 350 434 590 770 974 1202 1454 1730 2030"),
            glossary="""Number of angular points in the Lebedev grids."""))

    #options.add_int("GAMESS_EOMINP_NSTATE", 1);

    # $CIDET

    options.add(
        'gamess',
        Keyword(keyword='cidet__ncore',
                default=0,
                validator=parsers.nonnegative_integer,
                glossary="""Total number of orbitals doubly occupied in all determinants."""))

    options.add(
        'gamess',
        Keyword(
            keyword='cidet__nact',
            default=1,  #None,
            validator=parsers.positive_integer,
            glossary="""Total number of active orbitals."""))

    options.add(
        'gamess',
        Keyword(
            keyword='cidet__nels',
            default=2,  #None,
            validator=parsers.positive_integer,
            glossary="""Total number of active electrons."""))


    # <<< FORCE

    options.add(
        "gamess",
        Keyword(
            keyword="force__method",
            default="analytic",
            validator=parsers.enum("ANALYTIC SEMINUM FULLNUM"),
            glossary="""Chooses between fully analytic, numerical differention of analytic first derivatives, or double differentiation of energies."""))

