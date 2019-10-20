import collections

from qcelemental.util import which, which_import

from .. import data_dir
from ..keywords import Keywords
from ..molecule import Molecule
from ..keywords.read_options import load_qcdb_defaults
from ..programs.psi4.options import load_cfour_defaults_from_psi4, load_psi4_defaults
from ..programs.nwchem.options import load_nwchem_defaults
from ..programs.gamess.options import load_gamess_defaults


def clean_options():
    global nu_options
    nu_options = Keywords()


def load_options():
    global nu_options

    load_program_options(nu_options)


def load_program_options(options):
    """Initialize program options defaults onto `options`."""

    load_qcdb_defaults(options)
    if which('xcfour') and which_import('psi4'):
        load_cfour_defaults_from_psi4(options)
    if which('nwchem'):
        load_nwchem_defaults(options)
    if which('rungms'):
        load_gamess_defaults(options)
    if which('psi4') and which_import('psi4'):
        load_psi4_defaults(options)
    if which_import('resp_qcdb'):
        import resp_qcdb
        resp.load_defaults(nu_options)


# here liveth the options _during_ function calls
nu_options = None
clean_options()

# here liveth the molecule between function calls
active_molecule = Molecule("""H 0 0 0\nH 0.74 0 0""")

# here liveth the QCVariables when not attached to jobrec
active_qcvars = {}
