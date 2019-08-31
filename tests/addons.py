import pytest

from qcengine.testing import using_cfour, using_dftd3, using_gamess, using_nwchem, using_psi4

__all__ = [
    'using_cfour',
    'using_dftd3',
    'using_gamess',
    'using_geometric',
    'using_nwchem',
    'using_psi4',
]

using_geometric = pytest.mark.skipif(True, reason='QCDB written for a local pre-QCSchema interface')
