import re

from . import pe
from . import driver_util
from .. import moptions
#from . import pe
#from .driver import options

def _cbs_wrapper_methods(**kwargs):
    cbs_method_kwargs = ['scf_wfn', 'corl_wfn', 'delta_wfn']
    cbs_method_kwargs += ['delta%d_wfn' % x for x in range(2, 6)]

    cbs_methods = []
    for method in cbs_method_kwargs:
        if method in kwargs:
            cbs_methods.append(kwargs[method])
    return cbs_methods


def _parse_cbs_gufunc_string(method_name):
    """Process a 'mtd/bas + D:bas2'-like `method_name` into list of methods and bases."""

    method_name_list = re.split( """\+(?![^\[\]]*\]|[^\(\)]*\))""", method_name)
    if len(method_name_list) > 2:
        raise ValidationError("CBS gufunc: Text parsing is only valid for a single delta, please use the CBS wrapper directly")

    method_list = []
    basis_list = []
    for num, method_str in enumerate(method_name_list):
        if (method_str.count("[") > 1) or (method_str.count("]") > 1):
            raise ValidationError("""CBS gufunc: Too many brackets given! %s """ % method_str)

        if method_str.count('/') != 1:
            raise ValidationError("""CBS gufunc: All methods must specify a basis with '/'. %s""" % method_str)

        if num > 0:
            method_str = method_str.strip()
            if method_str[:2].lower() != 'd:':
                raise ValidationError("""CBS gufunc: Delta method must start with 'D:'.""")
            else:
                method_str = method_str[2:]
        method, basis = method_str.split('/')
        method_list.append(method)
        basis_list.append(basis)
    return method_list, basis_list


def _cbs_gufunc(func, total_method_name, molecule, **kwargs):
    """Text-based wrapper of the CBS function."""

    # Catch kwarg issues
    print('\nINTO _cbs_gufunc', 'KW', kwargs)
    kwargs = driver_util.kwargs_lower(kwargs)
    return_wfn = kwargs.pop('return_wfn', False)
#    core.clean_variables()
    user_dertype = kwargs.pop('dertype', None)
    cbs_verbose = kwargs.pop('cbs_verbose', False)
    ptype = kwargs.pop('ptype', None)

#    # Make sure the molecule the user provided is the active one
#    molecule = kwargs.pop('molecule') #, core.get_active_molecule())
    molecule.update_geometry()

    # Sanitize total_method_name
    label = total_method_name
    total_method_name = total_method_name.lower().replace(' ', '')

    # Split into components
    method_list, basis_list = _parse_cbs_gufunc_string(total_method_name)

    # Single energy call?
    single_call = len(method_list) == 1
    single_call &= '[' not in basis_list[0]
    single_call &= ']' not in basis_list[0]

    if single_call:
        method_name = method_list[0]
        basis = basis_list[0]

        print('gufunc pre opt', pe.nu_options)
        #print('gufunc pre opt', pe.active_options)
        # Save some global variables so we can reset them later
        #optstash = options.OptionsState(['BASIS'])
        #core.set_global_option('BASIS', basis)
        #pe.active_options['GLOBALS']['BASIS']['value'] = basis
        #pe.nu_options.scroll['QCDB']['BASIS'].value = basis
        pe.nu_options.require('QCDB', 'BASIS', basis, accession=1234)

        #print('\n gufunc_calling', func, method_name, return_wfn, 'OPT', pe.active_options, 'KW', kwargs)
        print('\n gufunc_calling', func, method_name, return_wfn, 'OPT', pe.nu_options, 'KW', kwargs)
        ptype_value, wfn = func(method_name, return_wfn=True, molecule=molecule, **kwargs)
#        core.clean()

#        optstash.restore()

        if return_wfn:
            return (ptype_value, wfn)
        else:
            return ptype_value

    # If we are not a single call, let CBS wrapper handle it!
    cbs_kwargs = {}
    cbs_kwargs['ptype'] = ptype
    cbs_kwargs['return_wfn'] = True
    cbs_kwargs['molecule'] = molecule
    cbs_kwargs['verbose'] = cbs_verbose

    # Find method and basis
    if method_list[0] in ['scf', 'hf', 'c4-scf', 'c4-hf']:
        cbs_kwargs['scf_wfn'] = method_list[0]
        cbs_kwargs['scf_basis'] = basis_list[0]
        if 'scf_scheme' in kwargs:
            cbs_kwargs['scf_scheme'] = kwargs['scf_scheme']
    else:
        cbs_kwargs['corl_wfn'] = method_list[0]
        cbs_kwargs['corl_basis'] = basis_list[0]
        if 'corl_scheme' in kwargs:
            cbs_kwargs['corl_scheme'] = kwargs['corl_scheme']

    if len(method_list) > 1:
        cbs_kwargs['delta_wfn'] = method_list[1]
        cbs_kwargs['delta_basis'] = basis_list[1]
        if 'delta_scheme' in kwargs:
            cbs_kwargs['delta_scheme'] = kwargs['delta_scheme']

    ptype_value, wfn = cbs(func, label, **cbs_kwargs)

    if return_wfn:
        return (ptype_value, wfn)
    else:
        return ptype_value

