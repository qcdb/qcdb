from ..datastructures import *
from ..exceptions import *
from ..pdict import PreservingDict
from .glossary import qcvardefs

from .psivardefs import wfn_psivars

try:
    basestring
except NameError:
    basestring = str


def certify_qcvars(dicary):
    """

    Parameters
    ----------
    dicary : dict
        Dictionary where keys are strings (proper case – mostly UC) and values are float-like

    """
    calcinfo = []
    for pv, var in dicary.items():
        print('PV', pv, var)
        if pv in qcvardefs.keys():
            calcinfo.append(QCAspect(pv, qcvardefs[pv]['units'], var, ''))
        else:
            raise ValidationError('Undefined QCvar!: {}'.format(pv))

    return {info.lbl: info for info in calcinfo}
    #return calcinfo


def fill_in(rawvars, verbose=2):
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
        if verbose >= 2:
            print("""building %s %s""" % (pvar, '.' * (50 - len(pvar))), end='')

        data_rich_args = []

        for pv in action['args']:
            if isinstance(pv, basestring):
                if pv in rawvars:
                    data_rich_args.append(rawvars[pv])
                else:
                    if verbose >= 2:
                        print("""EMPTY, missing {}""".format(pv))
                    break
            else:
                data_rich_args.append(pv)
        else:
            result = action['func'](data_rich_args)
            rawvars[pvar] = result
            if verbose >= 2:
                print("""SUCCESS""")


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

