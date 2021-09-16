from . import parsers
from .keywords import Keyword, Keywords


def load_qcdb_keywords(options: Keywords) -> None:

    options.add(
        "qcdb",
        Keyword(
            keyword="translate_qcdb",
            default=True,
            validator=parsers.boolean,
            glossary="assert always translated. should qcdb or program defaults prevail?",
        ),
    )

    options.add(
        "qcdb",
        Keyword(  # true global
            keyword="memory",
            default="700 mb",
            validator=parsers.parse_memory,
            glossary="Total memory allocation in bytes.",
        ),
    )

    options.add(
        "qcdb",
        Keyword(  # true global
            keyword="basis", default="", validator=lambda x: x.upper(), glossary="Primary orbital basis set."
        ),
    )

    options.add(
        "qcdb",
        Keyword(
            keyword="e_convergence",
            default=1.0e-6,
            validator=parsers.parse_convergence,
            glossary="Convergence criterion for energy.",
        ),
    )

    options.add(
        "qcdb",
        Keyword(  # derived shorthand global
            keyword="scf__e_convergence",
            default=1.0e-6,
            validator=parsers.parse_convergence,
            glossary="Convergence criterion for SCF energy.",
        ),
    )

    options.add(
        "qcdb",
        Keyword(  # derived shorthand global
            keyword="scf__d_convergence",
            default=1.0e-6,
            validator=parsers.parse_convergence,
            glossary="Convergence criterion for SCF density (psi: which is defined as the RMS value of the orbital gradient.",
        ),
    )

    options.add(
        "qcdb",
        Keyword(  # true global
            keyword="puream",
            default=True,
            validator=parsers.sphcart,
            glossary="""Do use pure angular momentum basis functions?
  If not explicitly set, the default comes from the basis set.
  **Cfour Interface:** Keyword translates into |cfour__cfour_spherical|.""",
        ),
    )

    options.add(
        "qcdb",
        Keyword(
            keyword="reference",  # TODO don't want higher and local
            default="",
            validator=lambda x: x.upper(),  # TODO
            glossary="""Reference wavefunction type.
    **Cfour Interface:** Keyword translates into |cfour__cfour_reference|.""",
        ),
    )
    # options.add_str("REFERENCE", "RHF", "RHF ROHF UHF CUHF RKS UKS")

    options.add(
        "qcdb",
        Keyword(
            keyword="scf__reference",
            default="",
            validator=lambda x: x.upper(),  # TODO
            glossary="""Reference wavefunction type.
    **Cfour Interface:** Keyword translates into |cfour__cfour_reference|.""",
        ),
    )
    # options.add_str("REFERENCE", "RHF", "RHF ROHF UHF CUHF RKS UKS")

    options.add(
        "qcdb",
        Keyword(
            keyword="scf_type",  # TODO ditto, 2-leveled
            default="",
            validator=lambda x: x.upper(),  # TODO
            glossary="""What algorithm to use for the SCF computation. See Table :ref:`SCF
    Convergence & Algorithm <table:conv_scf>` for default algorithm for
    different calculation types.""",
        ),
    )
    # options.add_str("SCF_TYPE", "PK", "DIRECT DF PK OUT_OF_CORE CD GTFOCK");

    options.add(
        "qcdb",
        Keyword(
            keyword="scf__scf_type",
            default="",
            validator=lambda x: x.upper(),  # TODO
            glossary="""What algorithm to use for the SCF computation. See Table :ref:`SCF
    Convergence & Algorithm <table:conv_scf>` for default algorithm for
    different calculation types.""",
        ),
    )
    # options.add_str("SCF_TYPE", "PK", "DIRECT DF PK OUT_OF_CORE CD GTFOCK");

    options.add(
        "qcdb",
        Keyword(
            keyword="scf__maxiter",
            default=100,
            validator=parsers.positive_integer,
            glossary="""Maximum number of iterations.
    **Cfour Interface:** Keyword translates into |cfour__cfour_scf_maxcyc|.""",
        ),
    )

    options.add(
        "qcdb",
        Keyword(
            keyword="scf__damping_percentage",
            default=0.0,
            validator=parsers.percentage,
            glossary=r"""The amount (percentage) of damping to apply to the early density updates.
        0 will result in a full update, 100 will completely stall the update. A
        value around 20 (which corresponds to 20\% of the previous iteration's
        density being mixed into the current density)
        could help to solve problems with oscillatory convergence.""",
        ),
    )

    options.add(
        "qcdb",
        Keyword(
            keyword="mp2__mp2_type",
            default="",
            validator=lambda x: x.upper(),  # TODO
            glossary="""What algorithm to use for MP2 computation. See 
        :ref:`Cross-module Redundancies <table:managedmethods>` for details.""",
        ),
    )
    # options.add_str("MP2_TYPE", "DF", "DF CONV CD");

    options.add(
        "qcdb",
        Keyword(
            keyword="ss_scale",
            default=1.0,
            validator=lambda x: float(x),
            glossary="""Custom same-spin scaling value.""",
        ),
    )

    options.add(
        "qcdb",
        Keyword(
            keyword="os_scale",
            default=1.0,
            validator=lambda x: float(x),
            glossary="""Custom opposite-spin scaling value.""",
        ),
    )

    options.add(
        "qcdb",
        Keyword(
            keyword="writer_file_label",
            default="",
            validator=lambda x: x,
            glossary="""Base filename for text files written by QCDB, such as the
  MOLDEN output file, the Hessian file, the internal coordinate file, etc.
""",
        ),
    )

    options.add(
        "qcdb",
        Keyword(
            keyword="geom_maxiter",
            default=20,
            validator=parsers.positive_integer,
            glossary="""Maximum number of geometry optimization steps.""",
        ),
    )

    options.add(
        "qcdb",
        Keyword(
            keyword="g_convergence",
            default="qchem",
            validator=parsers.enum(
                "QCHEM MOLPRO GAU GAU_LOOSE GAU_TIGHT INTERFRAG_TIGHT GAU_VERYTIGHT TURBOMOLE CFOUR NWCHEM_LOOSE"
            ),
            glossary="""Set of optimization criteria. Specification of any MAX_*_G_CONVERGENCE
      or RMS_*_G_CONVERGENCE options will append to overwrite the criteria set here
      unless |optking__flexible_g_convergence| is also on.      See Table :ref:`Geometry Convergence <table:optkingconv>` for details.""",
        ),
    )

    options.add(
        "qcdb",
        Keyword(
            keyword="max_force_g_convergence",
            default=3.0e-4,
            validator=parsers.parse_convergence,
            glossary="""Convergence criterion for geometry optmization: maximum force
      (internal coordinates, atomic units).""",
        ),
    )

    options.add(
        "qcdb",
        Keyword(
            keyword="qc_module",
            default="",
            validator=lambda x: x.lower(),
            glossary="""Select alternate implementations.""",
        ),
    )
    # /*- Convergence criterion for geometry optmization: rms force
    # (internal coordinates, atomic units). -*/
    # options.add_double("RMS_FORCE_G_CONVERGENCE", 3.0e-4);
    # /*- Convergence criterion for geometry optmization: maximum energy change. -*/
    # options.add_double("MAX_ENERGY_G_CONVERGENCE", 1.0e-6);
    # /*- Convergence criterion for geometry optmization: maximum displacement
    # (internal coordinates, atomic units). -*/
    # options.add_double("MAX_DISP_G_CONVERGENCE", 1.2e-3);
    # /*- Convergence criterion for geometry optmization: rms displacement
    # (internal coordinates, atomic units). -*/
    # options.add_double("RMS_DISP_G_CONVERGENCE", 1.2e-3);
    # /*- Even if a user-defined threshold is set, allow for normal, flexible convergence criteria -*/
    # options.add_bool("FLEXIBLE_G_CONVERGENCE", false);

    options.add(
        "qcdb",
        Keyword(
            keyword="freeze_core",
            default=False,
            validator=parsers.boolean,
            glossary="""Handle core orbitals in correlated calculation.""",
        ),
    )


#    options.add('qcdb', Keyword(
#            keyword='',
#            default=,
#            validator=,
#            glossary="""."""))

#    options.add('qcdb', Keyword(
#            keyword='',
#            default=,
#            validator=,
#            glossary="""."""))

#    options.add('qcdb', Keyword(
#            keyword='',
#            default=,
#            validator=,
#            glossary="""."""))

#    options.add('qcdb', Keyword(
#            keyword='',
#            default=,
#            validator=,
#            glossary="""."""))
