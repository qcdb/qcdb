from decimal import Decimal
from typing import Dict, Union

import numpy as np
from qcelemental import Datum

from ..exceptions import ValidationError
from .glossary import qcvardefs
from .identities import wfn_qcvars


def certify_and_datumize(
    dicary: Dict[str, Union[float, Decimal, np.ndarray]], *, plump: bool = False, nat: int = None
) -> Dict[str, Datum]:
    """Convert a QCVariable: data dictionary to QCVariable: Datum dictionary, filling in details from the glossary.

    Parameters
    ----------
    dicary : dict
        Dictionary where keys are strings (proper case – mostly UC) and values are float-like

    """
    calcinfo = []
    for pv, var in dicary.items():
        if pv in qcvardefs.keys():
            doi = qcvardefs[pv].get("doi", None)
            if plump and isinstance(var, np.ndarray) and var.ndim == 1:
                var = var.reshape(eval(qcvardefs[pv]["dimension"].format(nat=nat)))
            calcinfo.append(Datum(pv, qcvardefs[pv]["units"], var, doi=doi, glossary=qcvardefs[pv]["glossary"]))

        else:
            defined = sorted(qcvardefs.keys())
            raise ValidationError("Undefined QCvar!: {}\n{}".format(pv, "\n\t".join(defined)))

    return {info.label: info for info in calcinfo}


def build_out(rawvars: Dict[str, Datum], verbose: int = 1) -> None:
    """Apply standard QC identities to QCVariables `rawvars` to build more (e.g., correlation from total and HF energies).

    Dictionary `wfn_qcvars` had keys with names of QCVariables to be created and values with dictionary of two keys:
    `args`, the QCVariables that contribute tot the key and `func`, a functional (or lambda) to combine them. This
    function builds that key QCVariables if all the contributors are available in `rawvars`. Updates internally so
    multiple passes not needed.

    Parameters
    ----------
    verbose
        Controls print level. Per-var printing with >=2.

    Returns
    -------
    None
        But input dictionary `rawvars` is updated.

    """
    for action in wfn_qcvars():
        pvar = action["form"]
        buildline = """building {} {}""".format(pvar, "." * (50 - len(pvar)))

        data_rich_args = []

        for pv in action["args"]:
            if isinstance(pv, str):
                if pv in rawvars:
                    data_rich_args.append(rawvars[pv])
                    if verbose >= 3:
                        print(f"{pv=} {rawvars[pv]}")
                else:
                    if verbose >= 2:
                        print("""{}EMPTY, missing {}""".format(buildline, pv))
                    break
            else:
                data_rich_args.append(pv)
        else:
            result = action["func"](data_rich_args)
            if verbose >= 3:
                print(f"{result=}")
            # rawvars[pvar] = result
            # with data coming from file --> variable, looks more precise than it is. hack
            rawvars.__setitem__(pvar, result, 6)
            if verbose >= 1:
                print("""{}SUCCESS""".format(buildline))

            if pvar == "CURRENT CORRELATION ENERGY" and abs(float(rawvars[pvar])) < 1.0e-16:
                rawvars.pop(pvar)
