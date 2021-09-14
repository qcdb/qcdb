import re
import pprint

from qcelemental.molparse import regex
import qcdb

import pytest
from .utils import *

NUMBER = r"(?x:" + regex.NUMBER + ")"


@pytest.mark.parametrize(
    "program,keywords",
    [
        pytest.param("cfour", {}, marks=using("cfour")),
        pytest.param("gamess", {}, marks=using("gamess")),
        pytest.param("nwchem", {}, marks=using("nwchem")),
        pytest.param("psi4", {"psi4_scf_type": "pk"}, marks=using("psi4")),
    ],
)
def test_scf_conv(program, keywords):
    prefix = qcdb.util.program_prefix(program)

    qcdb.set_keywords(keywords)

    hcn2 = qcdb.Molecule(
        """
H 0.000000 0.000000 3.409664
C 0.000000 0.000000 2.340950
N 0.000000 0.000000 1.185959
--
H 0.000000 0.000000 -0.148222
C 0.000000 0.000000 -1.222546
N 0.000000 0.000000 -2.379698
"""
    )

    ene, wfn = qcdb.energy(prefix + "hf/6-311G*", molecule=hcn2, return_wfn=True)
    pprint.pprint(wfn)

    assert wfn["success"] is True
    assert compare_values(-185.74334367736, ene, 6, "ene")

    ans = _harvest_scf_convergence(wfn["stdout"])
    print(ans)

    # assert 0, ans


# FAILED qcdb/tests/test_keyword_mode.py::test_scf_conv[cfour] - AssertionError: 1e-07
# FAILED qcdb/tests/test_keyword_mode.py::test_scf_conv[gamess] - AssertionError: 1e-05
# FAILED qcdb/tests/test_keyword_mode.py::test_scf_conv[nwchem] - AssertionError: 0.0001
# FAILED qcdb/tests/test_keyword_mode.py::test_scf_conv[psi4] - AssertionError: (1e-06, 1e-06)


def _harvest_scf_convergence(txt: str) -> float:
    # Psi4
    mobj = re.search(
        # fmt: off
        r"^\s+" + r"Energy threshold" + r"\s+=\s+" + NUMBER + r"\s*" +
        r"^\s+" + r"Density threshold" + r"\s+=\s+" + NUMBER + r"\s*$",
        # fmt: on
        txt,
        re.MULTILINE,
    )
    if mobj:
        e_conv = float(mobj.group(1))
        d_conv = float(mobj.group(2))
        return d_conv

    # CFOUR
    mobj = re.search(
        # fmt: off
        r"^\s+" + "SCF convergence tolerance:" + r"\s+" + r"10\*\*\(-\s*" + r"(?P<conv>\d+)" + r"\s*\)\s*$",
        # fmt: on
        txt,
        re.MULTILINE,
    )
    if mobj:
        scf_conv = 10 ** (-1 * int(mobj.group("conv")))
        return scf_conv

    # GAMESS
    mobj = re.search(
        # fmt: off
        r"^\s+" + r"DENSITY MATRIX CONV=\s+" + NUMBER + r"\s*$",
        # fmt: on
        txt,
        re.MULTILINE,
    )
    if mobj:
        d_conv = float(mobj.group(1))
        return d_conv

    # NWChem
    mobj = re.search(
        # fmt: off
        r"^\s+" + r"Convergence threshold" + r"\s*:\s*" + NUMBER + r"\s*" +
        r"^\s+" + r"Maximum no. of iterations" + r"\s*:\s*" + NUMBER + r"\s*" +
        r"^\s+" + r"Final Fock-matrix accuracy" + r"\s*:\s*" + NUMBER + r"\s*$",
        # fmt: on
        txt,
        re.MULTILINE,
    )
    if mobj:
        d_conv = float(mobj.group(1))
        return d_conv
