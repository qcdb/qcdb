import numpy as np

from qcelemental import Datum

from ..exceptions import *
from ..pdict import PreservingDict
from .glossary import qcvardefs

from .psivardefs import wfn_psivars

try:
    basestring
except NameError:
    basestring = str


def certify(dicary, plump=False, nat=None):
    """

    Parameters
    ----------
    dicary : dict
        Dictionary where keys are strings (proper case – mostly UC) and values are float-like

    """
    calcinfo = []
    for pv, var in dicary.items():
        if pv in qcvardefs.keys():
            doi = qcvardefs[pv].get('doi', None)
            if plump and isinstance(var, np.ndarray) and var.ndim == 1:
                var = var.reshape(eval(qcvardefs[pv]['dimension'].format(nat=nat)))
            calcinfo.append(Datum(pv, qcvardefs[pv]['units'], var, doi=doi, glossary=qcvardefs[pv]['glossary']))

        else:
            defined = sorted(qcvardefs.keys())
            raise ValidationError('Undefined QCvar!: {}\n{}'.format(pv, '\n\t'.join(defined)))

    return {info.label: info for info in calcinfo}


def build_out(rawvars, verbose=1):
    """

    Parameters
    ----------
    verbose : int, optional
        Controls print level. Per-var printing with >=2.

    Returns
    -------
    None
        But input dictionary `rawvars` is updated.

    """
    for action in wfn_psivars():
        pvar = action['form']
        buildline = """building {} {}""".format(pvar, '.' * (50 - len(pvar)))

        data_rich_args = []

        for pv in action['args']:
            if isinstance(pv, basestring):
                if pv in rawvars:
                    data_rich_args.append(rawvars[pv])
                else:
                    if verbose >= 2:
                        print("""{}EMPTY, missing {}""".format(buildline, pv))
                    break
            else:
                data_rich_args.append(pv)
        else:
            result = action['func'](data_rich_args)
            #rawvars[pvar] = result
            # with data coming from file --> variable, looks more precise than it is. hack
            rawvars.__setitem__(pvar, result, 6)
            if verbose >= 1:
                print("""{}SUCCESS""".format(buildline))


def expand_qcvars(qcvars, qvdefs, verbose=1):
    """Dictionary *qvdefs* has keys with names of PsiVariables to be
    created and values with dictionary of two keys: 'args', the
    PsiVariables that contribute to the key and 'func', a function (or
    lambda) to combine them. This function builds those PsiVariables if
    all the contributors are available. Helpful printing is available when
    PRINT > 2.

    """
    for pvar, action in qvdefs.items():
        if verbose >= 2:
            print("""building %s %s""" % (pvar, '.' * (50 - len(pvar))), end='')

        data_rich_args = []

        for pv in action['args']:
            if isinstance(pv, basestring):
                if pv in qcvars:
                    data_rich_args.append(qcvars[pv])
                else:
                    if verbose >= 2:
                        print("""EMPTY, missing {}""".format(pv))
                    break
            else:
                data_rich_args.append(pv)
        else:
            result = action['func'](data_rich_args)
            core.set_variable(pvar, result)
            if verbose >= 2:
                print("""SUCCESS""")


# probably move to driver_helpers
def get_variable_details(qcvar, qcvars=None):
    from ..driver import pe

    if qcvars is None:
        qcvars = pe.active_qcvars

    capslink = {k.lower(): k for k in qcvardefs.keys()}

    print(qcvars[capslink[qcvar.lower()]])
    print(qcvardefs[capslink[qcvar.lower()]])

