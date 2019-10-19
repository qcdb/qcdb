import collections

from qcelemental.util import which, which_import

from .. import data_dir
from ..molecule import Molecule
from ..keyword.read_options import load_qcdb_defaults
from ..keyword import Keywords
from ..programs.psi4.options import load_cfour_defaults_from_psi4, load_psi4_defaults
from ..programs.nwchem.options import load_nwchem_defaults
from ..programs.gamess.options import load_gamess_defaults


def clean_nu_options():
    global nu_options
    nu_options = Keywords()


def load_nu_options():
    global nu_options

    load_options(nu_options)
    #print('OPTIONS LOADED')
    #print(nu_options)


def load_options(options):
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


nu_options = None
clean_nu_options()

def clean_options():
    global active_options
    active_options = collections.defaultdict(lambda: collections.defaultdict(dict))
    active_options['GLOBALS']

# here liveth the options _during_ function calls
active_options = None
clean_options()

# here liveth the molecule between function calls
active_molecule = Molecule("""H 0 0 0\nH 0.74 0 0""")

# here liveth the QCVariables when not attached to jobrec
active_qcvars = {}

