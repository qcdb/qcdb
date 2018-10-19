import os
import subprocess

import pytest


def _which(command):
    # environment is $PSIPATH:$PATH, less any None values
    lenv = {'PATH': ':'.join([os.path.abspath(x) for x in os.environ.get('PSIPATH', '').split(':') if x != '']) +
                    ':' + os.environ.get('PATH')}
    lenv = {k: v for k, v in lenv.items() if v is not None}

    # thanks, http://stackoverflow.com/a/11270665
    try:
        from subprocess import DEVNULL  # py33
    except ImportError:
        DEVNULL = open(os.devnull, 'wb')

    try:
        dashout = subprocess.Popen(command, stdout=DEVNULL, stderr=subprocess.STDOUT, env=lenv)
    except OSError as e:
        return False
    else:
        return True


def _plugin_import(plug):
    import sys
    if sys.version_info >= (3, 4):
        from importlib import util
        plug_spec = util.find_spec(plug)
    else:
        import pkgutil
        plug_spec = pkgutil.find_loader(plug)
    if plug_spec is None:
        return False
    else:
        return True


def is_psi4_new_enough(version_feature_introduced):
    if not _plugin_import('psi4'):
        return False
    import psi4
    from pkg_resources import parse_version
    return parse_version(psi4.__version__) >= parse_version(version_feature_introduced)


#def is_numpy_new_enough(version_feature_introduced):
#    if not _plugin_import('numpy'):
#        return False
#    import numpy
#    from pkg_resources import parse_version
#    return parse_version(numpy.version.version) >= parse_version(version_feature_introduced)


#using_scipy = pytest.mark.skipif(_plugin_import('scipy') is False,
#                                reason='Not detecting module scipy. Install package if necessary and add to envvar PYTHONPATH')

using_psi4 = pytest.mark.skipif(_plugin_import('psi4') is False,
                                 reason='Not detecting module psi4. Install package and add to envvar PYTHONPATH')

#using_psi4_libxc = pytest.mark.skipif(is_psi4_new_enough("1.2a1.dev100") is False,
#                                reason="Psi4 does not include DFT rewrite to use Libxc. Update to development head")
#
#using_psi4_efpmints = pytest.mark.skipif(is_psi4_new_enough("1.2a1.dev507") is False,
#                                reason="Psi4 does not include EFP integrals in mints. Update to development head")
#
#using_psi4_python_integral_deriv = pytest.mark.skipif(is_psi4_new_enough("1000") is False,
#                                reason="Psi4 does not include derivatives of integrals exported to python. Update to development head")

using_psi4_molrec = pytest.mark.skipif(is_psi4_new_enough("1.2a1.dev999") is False,
                                reason="Psi4 does not use the new Molecule parsing. Update to development head")

#using_numpy_113 = pytest.mark.skipif(is_numpy_new_enough("1.13.0") is False,
#                                reason='NumPy does not include 1.13 features. Update package and add to envvar PYTHONPATH')
#
#using_matplotlib = pytest.mark.skipif(_plugin_import('matplotlib') is False,
#                                reason='Note detecting module matplotlib. Install package if necessary and add to envvar PYTHONPATH')

using_pylibefp = pytest.mark.skipif(_plugin_import('pylibefp') is False,
                                reason='Not detecting module pylibefp. Install package if necessary and add to envvar PYTHONPATH')

using_cfour = pytest.mark.skipif(_which('xcfour') is False,
                                reason='Not detecting executable xcfour. Install program if necessary and add to envvar PATH')

using_nwchem = pytest.mark.skipif(_which('nwchem') is False,
                                reason='Not detecting executable nwchem. Install program if necessary and add to envvar PATH')

using_geometric = pytest.mark.skipif(_plugin_import('geometric') is False,
                                reason='Not detecting module geomeTRIC. Install package if necessary and add to envvar PYTHONPATH')

using_dftd3 = pytest.mark.skipif(_which('dftd3') is False,
                                reason='Not detecting executable dftd3. Install package if necessary and add to envvar PATH')
