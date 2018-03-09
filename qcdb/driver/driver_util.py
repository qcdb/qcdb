from ..exceptions import *
from ..util import yes, no, der0th, der1st, der2nd, find_approximate_string_matches
from .proc_table import procedures


def kwargs_lower(kwargs):
    """Sanitize user's `kwargs`.
    * Rebuilds and returns `kwargs` dictionary with all keys made lowercase.
    * Should be called by every function that could be called directly by the user.
    * Turns boolean-like values into actual booleans.
    * Turns dertype values into ints.
    * Turns values lowercase if sensible.

    """
    caseless_kwargs = {}
    for key, value in kwargs.items():
        lkey = key.lower()
        if lkey in ['subset', 'banner']:  # only kw for which case matters
            lvalue = value
        else:
            try:
                lvalue = value.lower()
            except (AttributeError, KeyError):
                lvalue = value

        if lkey in ['irrep', 'check_bsse', 'linkage', 'bsse_type']:
            caseless_kwargs[lkey] = lvalue

        elif 'dertype' in lkey:
            if der0th.match(str(lvalue)):
                caseless_kwargs[lkey] = 0
            elif der1st.match(str(lvalue)):
                caseless_kwargs[lkey] = 1
            elif der2nd.match(str(lvalue)):
                caseless_kwargs[lkey] = 2
            else:
                raise KeyError('Derivative type key ({}) was not recognized'.format(key))

        elif yes.match(str(lvalue)):
            caseless_kwargs[lkey] = True

        elif no.match(str(lvalue)):
            caseless_kwargs[lkey] = False

        else:
            caseless_kwargs[lkey] = lvalue

    return caseless_kwargs


def get_package(lowername, kwargs, default_package='psi4'):
    pkgprefix = {'p4-': 'psi4',
                 'c4-': 'cfour',
                 'd3-': 'dftd3',
                }
    for k, v in pkgprefix.items():
        if lowername.startswith(k):
            package = v
            break
    else:
        package = kwargs.get('package', default_package)

    return package


def get_package2(lowername, user_package=None, default_package='psi4'):
    pkgprefix = {'p4-': 'psi4',
                 'c4-': 'cfour',
                 'd3-': 'dftd3',
                }
    for k, v in pkgprefix.items():
        if lowername.startswith(k):
            package = v
            break
    else:
        package = user_package if user_package else default_package

    return package


def find_derivative_type(ptype, method_name, user_dertype, user_package):
    r"""
    Figures out the derivative type (0, 1, 2) for a given method_name. Will
    first use user default and then the highest available derivative type for
    a given method.
    """
    if ptype not in ['gradient', 'hessian']:
        raise ValidationError("find_derivative_type: ptype must either be gradient or hessian: {}".format(ptype))

    dertype = "(auto)"
    package = get_package2(method_name, user_package=user_package)
#   MOVED   #    jobrec = procedures['energy'][package][lowername](lowername, molecule=molecule, options=pe.nu_options, ptype='energy', **kwargs)

    # If user type is None, try to find the highest derivative
    if user_dertype is None:
        if (ptype == 'hessian') and (method_name in procedures['hessian'][package]):
            dertype = 2
            # Will need special logic if we ever have managed Hessians
        elif method_name in procedures['gradient'][package]:
            dertype = 1
# LAB oh, deah!
            #if procedures['gradient'][package][method_name].__name__.startswith('select_'):
            #    try:
            #        procedures['gradient'][method_name](method_name, probe=True)
            #    except ManagedMethodError:
            #        dertype = 0
        elif method_name in procedures['energy'][package]:
            dertype = 0
    else:
        # Quick sanity check. Only *should* be able to be None or int, but hey, kids today...
        if not isinstance(user_dertype, int):
            raise ValidationError("find_derivative_type: user_dertype ({}) should only be None or int!".format(user_dertype))
        dertype = user_dertype

    # TODO consider generalization for btwn-pkg

    # Summary validation
    if (dertype == 2) and (method_name in procedures['hessian'][package]):
        pass
    elif (dertype == 1) and (method_name in procedures['gradient'][package]):
        pass
    elif (dertype == 0) and (method_name in procedures['energy'][package]):
        pass
    else:
        alternatives = ''
        alt_method_name = find_approximate_string_matches(method_name, procedures['energy'][package].keys(), 2)
        if len(alt_method_name) > 0:
            alternatives = """ Did you mean? {}""".format(' '.join(alt_method_name))

        raise ValidationError("""Derivative method 'name' ({}) and derivative level 'dertype' ({}) are not available.{}""".
            format(method_name, dertype, alternatives))

    return dertype

