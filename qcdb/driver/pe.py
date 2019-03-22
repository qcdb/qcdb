import collections

from .. import data_dir
from ..molecule import Molecule
from ..moptions.read_options2 import RottenOptions, load_qcdb_defaults
from ..iface_psi4.options import load_cfour_defaults_from_psi4, load_psi4_defaults, load_nwchem_defaults_from_psi4
from ..intf_gamess.options import load_gamess_defaults


def clean_nu_options():
    global nu_options
    nu_options = RottenOptions()


def load_nu_options():
    global nu_options

    load_qcdb_defaults(nu_options)
    load_cfour_defaults_from_psi4(nu_options)
    load_nwchem_defaults_from_psi4(nu_options)
    load_gamess_defaults(nu_options)
    load_psi4_defaults(nu_options)
    try:
        import resp
    except ImportError:
        print('import resp_qcdb FAIL')
    else:
        print('resp at', resp.__file__)
        resp.load_defaults(nu_options)
    #print('OPTIONS LOADED')
    #print(nu_options)


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

