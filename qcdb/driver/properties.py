"""Module with a *procedures* dictionary specifying available quantum
chemical methods and functions driving the main quantum chemical
functionality, namely single-point energies, geometry optimizations,
properties, and vibrational frequency calculations.

"""
import copy
import pprint

from ..keywords import register_kwds
from . import pe  # keep this at top of imports
from . import cbs_driver, driver_helpers, driver_util
from .proc_table import procedures

pp = pprint.PrettyPrinter(width=120)



#   def properties(*args, **kwargs):
#       r"""Function to compute various properties.
#
#       :aliases: prop()
#
#       :returns: none.
#
#       .. caution:: Some features are not yet implemented. Buy a developer a coffee.
#
#          - This function at present has a limited functionality.
#            Consult the keywords sections of other modules for further properties capabilities.
#
#       +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+
#       | Name               | Calls Method                                  | Reference      | Supported Properties                                          |
#       +====================+===============================================+================+===============================================================+
#       | scf                | Self-consistent field method(s)               | RHF/ROHF/UHF   | Listed :ref:`here <sec:oeprop>`                               |
#       +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+
#       | hf                 | HF Self-consistent field method(s)            | RHF/ROHF/UHF   | Listed :ref:`here <sec:oeprop>`                               |
#       +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+
#       | mp2                | MP2 with density fitting only (mp2_type df)   | RHF            | Listed :ref:`here <sec:oeprop>`                               |
#       +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+
#       | cc2                | 2nd-order approximate CCSD                    | RHF            | dipole, quadrupole, polarizability, rotation, roa_tensor      |
#       +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+
#       | ccsd               | Coupled cluster singles and doubles (CCSD)    | RHF            | dipole, quadrupole, polarizability, rotation, roa_tensor      |
#       +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+
#       | eom-cc2            | 2nd-order approximate EOM-CCSD                | RHF            | oscillator_strength, rotational_strength                      |
#       +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+
#       | eom-ccsd           | Equation-of-motion CCSD (EOM-CCSD)            | RHF            | oscillator_strength, rotational_strength                      |
#       +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+
#       | cisd, cisdt,       | Configuration interaction                     | RHF/ROHF       | Listed :ref:`here <sec:oeprop>`, transition_dipole,           |
#       | cisdt, cisdtq,     |                                               |                | transition_quadrupole                                         |
#       | ci5, ..., fci      |                                               |                |                                                               |
#       +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+
#       | casscf, rasscf     | Multi-configurational SCF                     | RHF/ROHF       | Listed :ref:`here <sec:oeprop>`, transition_dipole,           |
#       |                    |                                               |                | transition_quadrupole                                         |
#       +--------------------+-----------------------------------------------+----------------+---------------------------------------------------------------+
#
#       :type name: string
#       :param name: ``'ccsd'`` || etc.
#
#           First argument, usually unlabeled. Indicates the computational method
#           to be applied to the system.
#
#       :type properties: array of strings
#       :param properties: |dl| ``[]`` |dr| || ``['rotation', 'polarizability', 'oscillator_strength', 'roa']`` || etc.
#
#           Indicates which properties should be computed. Defaults to dipole and quadrupole.
#
#       :type molecule: :ref:`molecule <op_py_molecule>`
#       :param molecule: ``h2o`` || etc.
#
#           The target molecule, if not the last molecule defined.
#
#       :examples:
#
#       >>> # [1] Optical rotation calculation
#       >>> properties('cc2', properties=['rotation'])
#
#       """
@register_kwds(pe.nu_options)
def properties(name, **kwargs):
    r"""Function to compute the single-point electronic properties."""

    from . import load_proc_table
    kwargs = driver_util.kwargs_lower(kwargs)

    if 'options' in kwargs:
        driver_helpers.set_options(kwargs.pop('options'))

    # Bounce if name is function
    if hasattr(name, '__call__'):
        return name(properties, kwargs.pop('label', 'custom function'), ptype='energy', **kwargs)

    # Allow specification of methods to arbitrary order
    lowername = name.lower()
    lowername, level = driver_helpers._parse_arbitrary_order(lowername)
    if level:
        kwargs['level'] = level

    # Make sure the molecule the user provided is the active one
    molecule = kwargs.pop('molecule', driver_helpers.get_active_molecule())
    molecule.update_geometry()

    if len(pe.nu_options.scroll) == 0:
        #print('EMPTY OPT')
        pe.load_options()


#    # Bounce to CP if bsse kwarg
#    if kwargs.get('bsse_type', None) is not None:
#        return driver_nbody.nbody_gufunc(properties, name, ptype='energy', **kwargs)

    # Bounce to CBS if "method/basis" name
    if '/' in lowername:
        return cbs_driver._cbs_gufunc(properties, name, ptype='properties', molecule=molecule, **kwargs)

    # Commit to procedures['properties'] call hereafter
    return_wfn = kwargs.pop('return_wfn', False)
    pe.active_qcvars = {}

#    #for precallback in hooks['properties']['pre']:
#    #    precallback(lowername, **kwargs)
#
#    optstash = driver_util._set_convergence_criterion('properties', lowername, 6, 8, 6, 8, 6)
#
#    # Before invoking the procedure, we rename any file that should be read.
#    # This is a workaround to do restarts with the current PSI4 capabilities
#    # before actual, clean restarts are put in there
#    # Restartfile is always converted to a single-element list if
#    # it contains a single string
#    if 'restart_file' in kwargs:
#        restartfile = kwargs['restart_file']  # Option still available for procedure-specific action
#        if restartfile != list(restartfile):
#            restartfile = [restartfile]
#        # Rename the files to be read to be consistent with psi4's file system
#        for item in restartfile:
#            name_split = re.split(r'\.', item)
#            filenum = name_split[len(name_split) - 1]
#            try:
#                filenum = int(filenum)
#            except ValueError:
#                filenum = 32  # Default file number is the checkpoint one
#            psioh = core.IOManager.shared_object()
#            psio = core.IO.shared_object()
#            filepath = psioh.get_file_path(filenum)
#            namespace = psio.get_default_namespace()
#            pid = str(os.getpid())
#            prefix = 'psi'
#            targetfile = filepath + prefix + '.' + pid + '.' + namespace + '.' + str(filenum)
#            shutil.copy(item, targetfile)

#PR    print('QWER', pe.nu_options.print_changed())
    package = driver_util.get_package(lowername, kwargs)
    #for k, v in pkgprefix.items():
    #    if lowername.startswith(k):
    #        package = v
    #        break
    #else:
    #    package = kwargs.get('package', 'psi4')
    #print('\nENE calling', 'procedures', package, lowername, 'with', lowername, molecule, pe.nu_options, kwargs)
    #jobrec = procedures['properties'][package][lowername](lowername, molecule=molecule, options=pe.active_options, **kwargs)
    jobrec = procedures['properties'][package][lowername](lowername, molecule=molecule, options=pe.nu_options, ptype='properties', **kwargs)

#    for postcallback in hooks['properties']['post']:
#        postcallback(lowername, wfn=wfn, **kwargs)
#
#    optstash.restore()
    #jobrec.pop('raw_output')  # just to moderate printint to screen
#PR    pp.pprint(jobrec)
    pe.active_qcvars = copy.deepcopy(jobrec['qcvars'])

   # if return_wfn:  # TODO current properties safer than wfn.energy() for now, but should be revisited

#        # TODO place this with the associated call, very awkward to call this in other areas at the moment
#        if lowername in ['efp', 'mrcc', 'dmrg', 'psimrcc']:
#            core.print_out("\n\nWarning! %s does not have an associated derived wavefunction." % name)
#            core.print_out("The returned wavefunction is the incoming reference wavefunction.\n\n")
#        elif 'sapt' in lowername:
#            core.print_out("\n\nWarning! %s does not have an associated derived wavefunction." % name)
#            core.print_out("The returned wavefunction is the dimer SCF wavefunction.\n\n")

    #    return (float(jobrec['qcvars']['CURRENT ENERGY'].data), jobrec)
    #else:
    #    return float(jobrec['qcvars']['CURRENT ENERGY'].data)
        # float() is for decimal.Decimal
