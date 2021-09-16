import collections

from ...keywords import Keyword, Keywords


def load_psi4_keywords(options: Keywords) -> None:

    opts = _query_options_defaults_from_psi()

    def p4_validator(val):
        try:
            nuval = val.upper()
        except AttributeError:
            nuval = val
        return nuval

    for m in opts:
        for o, v in opts[m].items():
            if m == "GLOBALS":
                keyword = o
            else:
                keyword = m + "__" + o

            options.add("psi4", Keyword(keyword=keyword, default=v["value"], validator=p4_validator))


def _basic_validator(val):
    try:
        nuval = val.upper()
    except AttributeError:
        nuval = val
    return nuval


def load_cfour_keywords_from_psi4(options: Keywords) -> None:

    opts = _query_options_defaults_from_psi()
    opts = opts["CFOUR"]

    for o, v in opts.items():
        if o.startswith("CFOUR_"):
            options.add("cfour", Keyword(keyword=o[6:], default=v["value"], validator=_basic_validator))


def _query_options_defaults_from_psi(changedOnly=False):
    """Function to return a string of commands to replicate the
    current state of user-modified options. Used to capture C++
    options information for distributed (sow/reap) input files.

    .. caution:: Some features are not yet implemented. Buy a developer a coffee.

       - Need some option to get either all or changed

       - Need some option to either get dict or set string or psimod command list

       - command return doesn't revoke has_changed setting for unchanged with changedOnly=False

    """
    import psi4

    modules = [
        # PSI4 Modules
        "ADC",
        "CCENERGY",
        "CCEOM",
        "CCDENSITY",
        "CCLAMBDA",
        "CCHBAR",
        "CCRESPONSE",
        "CCSORT",
        "CCTRIPLES",
        "CLAG",
        "CPHF",
        "CIS",
        "DCT",
        "DETCI",
        "DFMP2",
        "DFTSAPT",
        "FINDIF",
        "FNOCC",
        "LMP2",
        "MCSCF",
        "MINTS",
        "MRCC",
        "OCC",
        "OPTKING",
        "PSIMRCC",
        "RESPONSE",
        "SAPT",
        "SCF",
        "STABILITY",
        "THERMO",
        "TRANSQT",
        "TRANSQT2",
        # External Modules
        "CFOUR",
    ]

    options = collections.defaultdict(dict)

    for opt in psi4.core.get_global_option_list():
        hoc = psi4.core.has_global_option_changed(opt)
        if hoc or not changedOnly:
            if opt in ["DFT_CUSTOM_FUNCTIONAL", "EXTERN"]:  # Feb 2017 hack
                continue
            val = psi4.core.get_global_option(opt)
            options["GLOBALS"][opt] = {"value": val, "has_changed": hoc}
        for module in modules:
            if psi4.core.option_exists_in_module(module, opt):
                hoc = psi4.core.has_option_changed(module, opt)
                if hoc or not changedOnly:
                    val = psi4.core.get_option(module, opt)
                    options[module][opt] = {"value": val, "has_changed": hoc}

    return options
