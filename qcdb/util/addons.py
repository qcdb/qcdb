

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

    return (plug_spec is not None)


def _is_psi4_new_enough(version_feature_introduced):
    if not _plugin_import('psi4'):
        return False
    import psi4
    from pkg_resources import parse_version
    return parse_version(psi4.__version__) >= parse_version(version_feature_introduced)


def _is_numpy_new_enough(version_feature_introduced):
    if not _plugin_import('numpy'):
        return False
    import numpy
    from pkg_resources import parse_version
    return parse_version(numpy.version.version) >= parse_version(version_feature_introduced)

_addons_ = {
    'psi4': _plugin_import('psi4'),
    #"ambit": _CMake_to_Py_boolean("@ENABLE_ambit@"),
    #"chemps2": _CMake_to_Py_boolean("@ENABLE_CheMPS2@"),
    #"dkh": _CMake_to_Py_boolean("@ENABLE_dkh@"),
    #"libefp": _CMake_to_Py_boolean("@ENABLE_libefp@"),
    #"erd": _CMake_to_Py_boolean("@ENABLE_erd@"),
    #"gdma": _CMake_to_Py_boolean("@ENABLE_gdma@"),
    #"pcmsolver": _CMake_to_Py_boolean("@ENABLE_PCMSolver@"),
    #"simint": _CMake_to_Py_boolean("@ENABLE_simint@"),
    #"dftd3": _which("dftd3"),
    "cfour": _which("xcfour"),
    #"mrcc": _which("dmrcc"),
    #"gcp": _which("gcp"),
    #"v2rdm_casscf": _plugin_import("v2rdm_casscf"),
    #"forte": _plugin_import("forte"),
    #"snsmp2": _plugin_import("snsmp2"),
}

def addons(request=None):
    """Returns boolean of whether Add-On *request* is available to Psi4,
    either compiled in or searchable in $PSIPATH:$PATH, as relevant. If
    *request* not passed, returns list of available Add-Ons.

    """
    if request is None:
        return sorted([k for k, v in _addons_.items() if v])
    return _addons_[request.lower()]


