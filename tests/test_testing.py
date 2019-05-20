import pytest

from functools import partial

import numpy as np

import qcdb

pytestmark = pytest.mark.quick

_arrs = {
    'a1234_14': np.arange(4),
    'blip14': np.arange(4) + [0., 0.02, 0.005, 0.02],
    'a1234_22': np.arange(4).reshape((2, 2)),
    'blip22': (np.arange(4) + [0., 0.02, 0.005, 0.02]).reshape((2, 2)),
    'iblip14': np.arange(4) + [0, 1, 0, 1],
    'iblip22': (np.arange(4) + [0, 1, 0, 1]).reshape((2, 2)),
}

_dcts = {
    'ell': {'a': _arrs['a1234_14'], 'b': {'ba': _arrs['a1234_14'], 'bb': _arrs['a1234_22']}},
    'elll': {'a': _arrs['a1234_14'], 'b': {'ba': _arrs['a1234_14'], 'bb': _arrs['a1234_22'], 'bc': 4}},
    'ellnone': {'a': _arrs['a1234_14'], 'b': {'ba': _arrs['a1234_14'], 'bb': _arrs['a1234_22'], 'bc': None}},
    'ellshort': {'a': np.arange(3), 'b': {'ba': _arrs['a1234_14'], 'bb': _arrs['a1234_22']}},
    'blipell': {'a': _arrs['blip14'], 'b': {'ba': _arrs['a1234_14'], 'bb': _arrs['blip22']}},
}  # yapf: disable


@pytest.mark.parametrize(
    "fn,args,kwargs",
    [
        # scalar int
        (qcdb.compare_integers, [1, 1, 'labeled'], {}),
        (qcdb.compare_integers, [1, 1], {'verbose': 2}),
        (qcdb.compare_strings, ['a', 'a', 'labeled'], {}),
        (qcdb.compare_strings, ['a', 'a'], {'verbose': 2}),
        (qcdb.compare, [1, 1, 'labeled'], {}),
        (qcdb.compare, [1, 1], {'quiet': True}),

        # array int
        (qcdb.compare, [_arrs['a1234_14'], _arrs['a1234_14'], 'labeled'], {}),
        (qcdb.compare, [_arrs['a1234_14'], _arrs['a1234_14']], {'quiet': True}),

        # scalar float
        (qcdb.compare_values, [4.0, 4.001, 2, 'qcdb api'], {}),
        (qcdb.compare_values, [4.0, 4.0000001, 'qcel api'], {}),
        (qcdb.compare_values, [4.0, 4.001, 2], {}),
        (qcdb.compare_values, [4.0, 4.001], {'atol': 1.e-2}),
        (qcdb.compare_values, [4.0, 4.0000001], {}),

        # array float
        (qcdb.compare_values, [_arrs['a1234_22'], _arrs['a1234_22'], 'labeled'], {}),
        (qcdb.compare_values, [_arrs['a1234_22'], _arrs['a1234_22']], {'quiet': True}),
        (qcdb.compare_values, [_arrs['a1234_22'], _arrs['blip22'], 0.1, 'labeled'], {}),
        (qcdb.compare_values, [_arrs['a1234_22'], _arrs['blip22'], 0.1], {'quiet': True}),
        (qcdb.compare_values, [_arrs['a1234_22'], _arrs['blip22'], 'labeled'], {'atol': 0.1}),
        (qcdb.compare_values, [_arrs['a1234_22'], _arrs['blip22']], {'atol': 0.1, 'quiet': True}),
        (qcdb.compare_arrays, [[-1.2, 1.2], [-1.2, 1.20000000002]], {}),
        (qcdb.compare_arrays, [[-1.2, 1.2], [-1.2, 1.20000000002], 6], {}),

        # dicts
        (qcdb.compare_recursive, [_dcts['ell'], _dcts['ell'], 'labeled'], {}),
        (qcdb.compare_recursive, [_dcts['ell'], _dcts['ell']], {}),

    ])  # yapf: disable
def test_qcdb_compare_true(fn, args, kwargs):
    assert fn(*args, **kwargs)


@pytest.mark.parametrize(
    "fn,args,kwargs",
    [
        # scalar int
        (qcdb.compare_integers, [1, 2, 'labeled'], {}),
        (qcdb.compare_integers, [1, 2], {'verbose': 2}),
        (qcdb.compare_strings, ['a', 'b', 'labeled'], {}),
        (qcdb.compare_strings, ['a', 'b'], {'verbose': 2}),
        (qcdb.compare, [1, 2, 'labeled'], {}),
        (qcdb.compare, [1, 2], {'quiet': True}),

        # array int
        (qcdb.compare, [_arrs['a1234_14'], _arrs['iblip14'], 'labeled'], {}),
        (qcdb.compare, [_arrs['a1234_14'], _arrs['iblip14']], {'quiet': True}),

        # scalar float
        (qcdb.compare_values, [4.0, 4.1, 2, 'qcdb api'], {}),
        (qcdb.compare_values, [4.0, 4.0001, 'qcel api'], {}),
        (qcdb.compare_values, [4.0, 4.1, 2], {}),
        (qcdb.compare_values, [4.0, 4.1], {'atol': 1.e-2}),
        (qcdb.compare_values, [4.0, 4.0001], {}),
        (qcdb.compare_values, [4.0, 4.001, 4], {'atol': 1.e-1}),  # arg trumps kwarg

        # array float
        (qcdb.compare_values, [_arrs['a1234_22'], _arrs['blip22'], 'labeled'], {}),
        (qcdb.compare_values, [_arrs['a1234_22'], _arrs['blip22']], {'quiet': True}),
        (qcdb.compare_values, [_arrs['a1234_22'], _arrs['blip22'], 7, 'labeled'], {}),
        (qcdb.compare_values, [_arrs['a1234_22'], _arrs['blip22'], 7], {'quiet': True}),
        (qcdb.compare_values, [_arrs['a1234_22'], _arrs['blip22'], 'labeled'], {'atol': 1e-7}),
        (qcdb.compare_values, [_arrs['a1234_22'], _arrs['blip22']], {'atol': 1e-7, 'quiet': True}),
        (qcdb.compare_arrays, [[-1.2, 1.2], [-1.2, 1.2002]], {}),
        (qcdb.compare_arrays, [[-1.2, 1.2], [-1.2, 1.20000000002], 12], {}),
        (qcdb.compare_arrays, [[-1.2, 1.2], [-1.2, 1.2002, 2.4]], {}),

        # dicts
        (qcdb.compare_recursive, [_dcts['elll'], _dcts['ell']], {}),

    ])  # yapf: disable
def test_qcdb_compare_raise(fn, args, kwargs):
    with pytest.raises(qcdb.TestComparisonError):
        fn(*args, **kwargs)


@pytest.mark.parametrize(
    "fn,args,kwargs",
    [
        (qcdb.compare_recursive, [_dcts['ell'], _dcts['ell'], 4], {}),
        (qcdb.compare_matrices, [None, None], {}),
        (qcdb.compare_dicts, [None, None], {}),
    ])  # yapf: disable
def test_compare_upgrade(fn, args, kwargs):
    with pytest.raises(qcdb.UpgradeHelper):
        fn(*args, **kwargs)


def _true_false_handler(passfail, label, message, return_message=False, quiet=False):
    print(f"""    {label:.<66}{'PASSED' if passfail else 'FAILED'}""")
    return passfail


_tf_compare_integers = partial(qcdb.testing._psi4_compare_integers, return_handler=_true_false_handler)


def test_alt_handler_compare_true():
    assert _tf_compare_integers(1, 1) is True


def test_alt_handler_compare_false():
    assert _tf_compare_integers(1, 2) is False
