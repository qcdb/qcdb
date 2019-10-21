from qcelemental.util import which, which_import

from .. import data_dir
from ..keywords import Keywords
from ..keywords.read_options import load_qcdb_keywords
from ..molecule import Molecule
from ..programs.gamess.read_options import load_gamess_keywords
from ..programs.nwchem.read_options import load_nwchem_keywords
from ..programs.psi4.read_options import load_cfour_keywords_from_psi4, load_psi4_keywords


def clean_options() -> None:
    """Initialize the empty global keywords object."""

    global nu_options
    nu_options = Keywords()


def load_options() -> None:
    """Initialize the global keywords object with program defaults."""
    global nu_options

    load_program_options(nu_options)


def load_program_options(options: Keywords) -> None:
    """Initialize program keywords with defaults onto `options`."""

    load_qcdb_keywords(options)
    if which('xcfour') and which_import('psi4'):
        load_cfour_keywords_from_psi4(options)
    if which('nwchem'):
        load_nwchem_keywords(options)
    if which('rungms'):
        load_gamess_keywords(options)
    if which('psi4') and which_import('psi4'):
        load_psi4_keywords(options)
    if which_import('resp_qcdb'):
        import resp_qcdb
        resp.load_keywords(nu_options)


# here liveth the options _during_ function calls
nu_options = None
clean_options()

# here liveth the molecule between function calls
active_molecule = Molecule("""H 0 0 0\nH 0.74 0 0""")

# here liveth the QCVariables when not attached to jobrec
active_qcvars = {}
