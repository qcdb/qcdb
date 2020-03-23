from ..exceptions import ValidationError
from ..util import der0th, der1st, der2nd, find_approximate_string_matches, no, yes
from .proc_table import procedures

pkgprefix = {'p4-': 'psi4',
             'c4-': 'cfour',
             'd3-': 'dftd3',
             'nwc-': 'nwchem',
             'gms-': 'gamess',
            }

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
    for k, v in pkgprefix.items():
        if lowername.startswith(k):
            package = v
            break
    else:
        package = kwargs.get('package', default_package)

    return package


def get_package2(lowername, user_package=None, default_package='psi4'):
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
        elif method_name in procedures['properties'][package]:
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
    elif (dertype == 0) and (method_name in procedure['properties'][package]):
        pass
    else:
        alternatives = ''
        alt_method_name = find_approximate_string_matches(method_name, procedures['energy'][package].keys(), 2)
        if len(alt_method_name) > 0:
            alternatives = """ Did you mean? {}""".format(' '.join(alt_method_name))

        raise ValidationError("""Derivative method 'name' ({}) and derivative level 'dertype' ({}) are not available.{}""".
            format(method_name, dertype, alternatives))

    return dertype


def _process_displacement(derivfunc, method, molecule, displacement, n, ndisp, **kwargs):
    """A helper function to perform all processing for an individual finite
       difference computation.

       Parameters
       ----------
       derivfunc : func
           The function computing the target derivative.
       method : str
          A string specifying the method to be used for the computation.
       molecule: psi4.core.molecule or qcdb.molecule
          The molecule for the computation. All processing is handled internally.
          molecule must not be modified!
       displacement : dict
          A dictionary containing the necessary information for the displacement.
          See driver_findif/_geom_generator.py docstring for details.
       n : int
          The number of the displacement being computed, for print purposes.
       ndisp : int
           The total number of geometries, for print purposes.

       Returns
       -------
       wfn: :py:class:`~psi4.core.Wavefunction`
           The wavefunction computed.
    """
    import sys
    import numpy as np
    import psi4

    # print progress to file and screen
    psi4.core.print_out('\n')
    #p4util.banner('Loading displacement %d of %d' % (n, ndisp))
    print(""" %d""" % (n), end=('\n' if (n == ndisp) else ''))
    sys.stdout.flush()

    parent_group = molecule.point_group()
    clone = molecule.clone()
    clone.reinterpret_coordentry(False)
    clone.fix_orientation(True)

    # Load in displacement (flat list) into the active molecule
    geom_array = np.reshape(displacement["geometry"], (-1, 3))
    clone.set_geometry(psi4.core.Matrix.from_array(geom_array))

    # If the user insists on symmetry, weaken it if some is lost when displacing.
    if molecule.symmetry_from_input():
        disp_group = clone.find_highest_point_group()
        new_bits = parent_group.bits() & disp_group.bits()
        new_symm_string = qcdb.PointGroup.bits_to_full_name(new_bits)
        clone.reset_point_group(new_symm_string)

    # clean possibly necessary for n=1 if its irrep (unsorted in displacement list) different from initial G0 for freq
    psi4.core.clean()

    # Perform the derivative calculation
    derivative, wfn = derivfunc(method, return_wfn=True, molecule=clone, **kwargs)
    displacement["energy"] = wfn['qcvars']['CURRENT ENERGY'].data

    # If we computed a first or higher order derivative, set it.
    if derivfunc.__name__ == 'gradient':
        displacement["gradient"] = wfn['qcvars']['CURRENT GRADIENT'].data

    # clean may be necessary when changing irreps of displacements
    psi4.core.clean()

    return wfn
