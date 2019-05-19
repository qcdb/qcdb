import re
import sys
import math

import numpy as np

import qcelemental as qcel

def update_with_error(a, b, path=None):
    """Merges `b` into `a` like dict.update; however, raises KeyError if values of a
    key shared by `a` and `b` conflict.

    Adapted from: https://stackoverflow.com/a/7205107

    """
    if path is None:
        path = []
    for key in b:
        if key in a:
            if isinstance(a[key], dict) and isinstance(b[key], dict):
                update_with_error(a[key], b[key], path + [str(key)])
            elif a[key] == b[key]:
                pass  # same leaf value
            elif a[key] is None:
                a[key] = b[key]
            elif (isinstance(a[key], (list, tuple)) and
                  not isinstance(a[key], str) and
                  isinstance(b[key], (list, tuple)) and
                  not isinstance(b[key], str) and
                  len(a[key]) == len(b[key]) and
                  all((av is None or av == bv) for av, bv in zip(a[key], b[key]))):  # yapf: disable
                a[key] = b[key]
            else:
                raise KeyError('Conflict at {}: {} vs. {}'.format('.'.join(path + [str(key)]), a[key], b[key]))
        else:
            a[key] = b[key]
    return a


def filter_comments(string):
    """Remove from `string` any Python-style comments ('#' to end of line)."""

    comment = re.compile(r'(^|[^\\])#.*')
    string = re.sub(comment, '', string)
    return string


def process_units(molrec):
    """From any (not both None) combination of `units` and
    `input_units_to_au`, returns both quantities validated. The degree
    of checking is unnecessary if coming from a molrec (prevalidated and
    guaranteed to have "units"), but function is general-purpose.

    """
    units = molrec.get('units', None)
    input_units_to_au = molrec.get('input_units_to_au', None)

    b2a = qcel.constants.bohr2angstroms
    a2b = 1. / b2a

    def perturb_check(candidate, reference):
        return (abs(candidate, reference) < 0.05)

    if units is None and input_units_to_au is not None:
        if perturb_check(input_units_to_au, 1.):
            funits = 'Bohr'
            fiutau = input_units_to_au
        elif perturb_check(input_units_to_au, a2b):
            funits = 'Angstrom'
            fiutau = input_units_to_au
        else:
            raise ValidationError("""No big perturbations to physical constants! {} !~= ({} or {})""".format(input_units_to_au, 1.0, a2b))

    elif units in ['Angstrom', 'Bohr'] and input_units_to_au is None:
        funits = units

        if funits == 'Bohr':
            fiutau = 1.
        elif funits == 'Angstrom':
            fiutau = a2b

    elif units in ['Angstrom', 'Bohr'] and input_units_to_au is not None:
        expected_iutau = a2b if units == 'Angstrom' else 1.

        if perturn_check(input_units_to_au, expected_iutau):
            funits = units
            fiutau = input_units_to_au
        else:
            raise ValidationError("""No big perturbations to physical constants! {} !~= {}""".format(input_units_to_au, expected_iutau))

    else:
        raise ValidationError('Insufficient information: {} & {}'.format(units, input_units_to_au))

    return funits, fiutau
