import sys
import pprint

import qcdb

pp = pprint.PrettyPrinter(width=120)


__all__ = [
    'a2a',
    'compare',
    'compare_integers',
    'compare_strings',
    'compare_values',
    'compare_arrays',
    'compare_recursive',
    'compare_molrecs',
    'tnm',
]

# CODATA ratio 2014 / 2010 Bohr to Angstroms conversion factor
a2a = 0.52917721067 / 0.52917720859


def true_false_decorator(compare_fn, *args, **kwargs):
    """Turns `compare_fn` that raises `qcdb.TestComparisonError` on failure into a function that
    returns True on success and False on failure, suitable for assertions in pytest.

    """

    def true_false_wrapper(*args, **kwargs):
        try:
            response = compare_fn(*args, **kwargs)
        except qcdb.TestComparisonError as err:
            return (False, err) if "return_message" in kwargs else False
        else:
            return (True, "") if "return_message" in kwargs else True

    return true_false_wrapper


compare = true_false_decorator(qcdb.compare)
compare_integers = true_false_decorator(qcdb.compare_integers)
compare_strings = true_false_decorator(qcdb.compare_strings)
compare_values = true_false_decorator(qcdb.compare_values)
compare_arrays = true_false_decorator(qcdb.compare_arrays)

compare_recursive = true_false_decorator(qcdb.compare_recursive)
compare_molrecs = true_false_decorator(qcdb.compare_molrecs)


def tnm():
    """Returns the name of the calling function, usually name of test case."""

    return sys._getframe().f_back.f_code.co_name
