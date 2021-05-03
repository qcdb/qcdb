"""Module to facilitate quantum chemical computations on chemical
databases. Contains Molecule class and physical constants from psi4 suite.

isort:skip_file
"""
__author__ = 'Lori A. Burns'

# Figure out psidatadir: envvar trumps staged/installed
import os
qcdb_module_loc = os.path.dirname(os.path.abspath(__file__))
pymod = os.path.normpath(os.sep.join(['@PYMOD_INSTALL_LIBDIR@', '@CMAKE_INSTALL_LIBDIR@', 'qcdb']))
if pymod.startswith(os.sep + os.sep):
    pymod = pymod[1:]
pymod_dir_step = os.sep.join(['..'] * pymod.count(os.sep))
data_dir = os.sep.join([qcdb_module_loc, pymod_dir_step, '@CMAKE_INSTALL_DATADIR@', 'qcdb'])

if 'PSIDATADIR' in os.environ.keys():
    data_dir = os.path.expanduser(os.environ['PSIDATADIR'])
elif 'CMAKE_INSTALL_DATADIR' in data_dir:
    data_dir = os.sep.join([os.path.abspath(os.path.dirname(__file__)), '..', 'share', 'qcdb'])

data_dir = os.path.abspath(data_dir)
data_dir = qcdb_module_loc
if not os.path.isdir(data_dir):
    raise KeyError('Unable to read the data folder - check the PSIDATADIR environment variable'
                   '      Current value of PSIDATADIR is {}'.format(data_dir))

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

from .driver import energy, properties, hessian, gradient, frequency
from .driver import optking, geometric
from .driver import vpt2
from .driver.cbs_driver import cbs
from .driver.cbs_helpers import *
from .driver.driver_helpers import get_variable, has_variable, print_variables, variable
from .driver.driver_helpers import set_keywords, set_options, get_active_options
from .driver.driver_helpers import set_molecule, activate
from .driver.yaml import yaml_run
from .exceptions import *

from .basisset import BasisSet, basishorde
from .keywords import AliasKeyword, Keyword, Keywords
from .molecule import Molecule
from . import vib

## Load items that are useful to access from an input file
from .vib import compare_vibinfos
from .testing import *

