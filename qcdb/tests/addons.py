import pytest

from qcelemental.util import parse_version, which, which_import
from qcengine.testing import _programs as _programs_qcng


# Figure out what is imported
_programs = {
    # non-QC
    # "networkx": which_import("networkx", return_bool=True),
    # QC
    # "adcc": which_import("adcc", return_bool=True),
}


def has_program(name):
    if name in _programs:
        return _programs[name]
    elif name in _programs_qcng:
        return _programs_qcng[name]
    else:
        raise KeyError(f"Program {name} not registered with QCDB testing.")


_using_cache = {}


def using(program):

    if program == "geometric":
        return pytest.mark.skipif(True, reason="QCDB written for a local pre-QCSchema geometric interface")

    if program not in _using_cache:
        import_message = f"Not detecting module {program}. Install package if necessary to enable tests."
        skip = pytest.mark.skipif(has_program(program) is False, reason=import_message)
        _using_cache[program] = skip

    return _using_cache[program]
