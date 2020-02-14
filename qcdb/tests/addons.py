import pytest

from qcengine.testing import using

__all__ = [
    'using_cfour',
    'using_dftd3',
    'using_gamess',
    'using_geometric',
    'using_nwchem',
    'using_psi4',
]

using_geometric = pytest.mark.skipif(True, reason='QCDB written for a local pre-QCSchema interface')
using_cfour = using("cfour")
using_dftd3 = using("dftd3")
using_gamess = using("gamess")
using_nwchem = using("nwchem")
using_psi4 = using("psi4")
