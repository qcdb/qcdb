"""Module with a *procedures* dictionary specifying available quantum
chemical methods and functions driving the main quantum chemical
functionality, namely single-point energies, geometry optimizations,
properties, and vibrational frequency calculations.

"""
import copy
import pprint

from ..molecule import Molecule
from . import driver_helpers, driver_util, pe
from .gradient import gradient

pp = pprint.PrettyPrinter(width=120)

#   import numpy as np


# def optimize(name, **kwargs):
#    r"""Function to perform a geometry optimization.
#
#    :aliases: opt()
#
#    :returns: *float* |w--w| Total electronic energy of optimized structure in Hartrees.
#
#    :returns: (*float*, :py:class:`~psi4.core.Wavefunction`) |w--w| energy and wavefunction when **return_wfn** specified.
#
#    :raises: psi4.OptimizationConvergenceError if |optking__geom_maxiter| exceeded without reaching geometry convergence.
#
#    :PSI variables:
#
#    .. hlist::
#       :columns: 1
#
#       * :qcvar:`CURRENT ENERGY <CURRENTENERGY>`
#
#    :type name: string
#    :param name: ``'scf'`` || ``'mp2'`` || ``'ci5'`` || etc.
#
#        First argument, usually unlabeled. Indicates the computational method
#        to be applied to the database. May be any valid argument to
#        :py:func:`~driver.energy`.
#
#    :type molecule: :ref:`molecule <op_py_molecule>`
#    :param molecule: ``h2o`` || etc.
#
#        The target molecule, if not the last molecule defined.
#
#    :type return_wfn: :ref:`boolean <op_py_boolean>`
#    :param return_wfn: ``'on'`` || |dl| ``'off'`` |dr|
#
#        Indicate to additionally return the :py:class:`~psi4.core.Wavefunction`
#        calculation result as the second element (after *float* energy) of a tuple.
#
#    :type return_history: :ref:`boolean <op_py_boolean>`
#    :param return_history: ``'on'`` || |dl| ``'off'`` |dr|
#
#        Indicate to additionally return dictionary of lists of geometries,
#        energies, and gradients at each step in the optimization.
#
#    :type func: :ref:`function <op_py_function>`
#    :param func: |dl| ``gradient`` |dr| || ``energy`` || ``cbs``
#
#        Indicates the type of calculation to be performed on the molecule.
#        The default dertype accesses ``'gradient'`` or ``'energy'``, while
#        ``'cbs'`` performs a multistage finite difference calculation.
#        If a nested series of python functions is intended (see :ref:`sec:intercalls`),
#        use keyword ``opt_func`` instead of ``func``.
#
#    :type mode: string
#    :param mode: |dl| ``'continuous'`` |dr| || ``'sow'`` || ``'reap'``
#
#        For a finite difference of energies optimization, indicates whether
#        the calculations required to complete the
#        optimization are to be run in one file (``'continuous'``) or are to be
#        farmed out in an embarrassingly parallel fashion
#        (``'sow'``/``'reap'``). For the latter, run an initial job with
#        ``'sow'`` and follow instructions in its output file. For maximum
#        flexibility, ``return_wfn`` is always on in ``'reap'`` mode.
#
#    :type dertype: :ref:`dertype <op_py_dertype>`
#    :param dertype: ``'gradient'`` || ``'energy'``
#
#        Indicates whether analytic (if available) or finite difference
#        optimization is to be performed.
#
#    :type hessian_with: string
#    :param hessian_with: ``'scf'`` || ``'mp2'`` || etc.
#
#        Indicates the computational method with which to perform a hessian
#        analysis to guide the geometry optimization.
#
#    .. warning:: Optimizations where the molecule is specified in Z-matrix format
#       with dummy atoms will result in the geometry being converted to a Cartesian representation.
#
#    .. note:: Analytic gradients area available for all methods in the table
#        below. Optimizations with other methods in the energy table proceed
#        by finite differences.
#
#    .. _`table:grad_gen`:
#
#    +-------------------------+---------------------------------------------------------------------------------------------------------------+
#    | name                    | calls method                                                                                                  |
#    +=========================+===============================================================================================================+
#    | efp                     | efp-only optimizations                                                                                        |
#    +-------------------------+---------------------------------------------------------------------------------------------------------------+
#    | scf                     | Hartree--Fock (HF) or density functional theory (DFT) :ref:`[manual] <sec:scf>`                               |
#    +-------------------------+---------------------------------------------------------------------------------------------------------------+
#    | hf                      | HF self consistent field (SCF) :ref:`[manual] <sec:scf>`                                                      |
#    +-------------------------+---------------------------------------------------------------------------------------------------------------+
#    | dct                     | density cumulant functional theory :ref:`[manual] <sec:dct>`                                                  |
#    +-------------------------+---------------------------------------------------------------------------------------------------------------+
#    | mp2                     | 2nd-order |MollerPlesset| perturbation theory (MP2) :ref:`[manual] <sec:dfmp2>` :ref:`[details] <tlmp2>`      |
#    +-------------------------+---------------------------------------------------------------------------------------------------------------+
#    | mp3                     | 3rd-order |MollerPlesset| perturbation theory (MP3) :ref:`[manual] <sec:occ_nonoo>` :ref:`[details] <tlmp3>`  |
#    +-------------------------+---------------------------------------------------------------------------------------------------------------+
#    | mp2.5                   | average of MP2 and MP3 :ref:`[manual] <sec:occ_nonoo>` :ref:`[details] <tlmp25>`                              |
#    +-------------------------+---------------------------------------------------------------------------------------------------------------+
#    | omp2                    | orbital-optimized second-order MP perturbation theory :ref:`[manual] <sec:occ_oo>`                            |
#    +-------------------------+---------------------------------------------------------------------------------------------------------------+
#    | omp3                    | orbital-optimized third-order MP perturbation theory :ref:`[manual] <sec:occ_oo>`                             |
#    +-------------------------+---------------------------------------------------------------------------------------------------------------+
#    | omp2.5                  | orbital-optimized MP2.5 :ref:`[manual] <sec:occ_oo>`                                                          |
#    +-------------------------+---------------------------------------------------------------------------------------------------------------+
#    | lccd                    | Linear CCD :ref:`[manual] <sec:occ_nonoo>` :ref:`[details] <tllccd>`                                          |
#    +-------------------------+---------------------------------------------------------------------------------------------------------------+
#    | olccd                   | orbital optimized LCCD :ref:`[manual] <sec:occ_oo>`                                                           |
#    +-------------------------+---------------------------------------------------------------------------------------------------------------+
#    | ccd                     | coupled cluster doubles  (CCD) :ref:`[manual] <sec:occ_nonoo>`                                                |
#    +-------------------------+---------------------------------------------------------------------------------------------------------------+
#    | ccsd                    | coupled cluster singles and doubles (CCSD) :ref:`[manual] <sec:cc>` :ref:`[details] <tlccsd>`                 |
#    +-------------------------+---------------------------------------------------------------------------------------------------------------+
#    | ccsd(t)                 | CCSD with perturbative triples (CCSD(T)) :ref:`[manual] <sec:cc>` :ref:`[details] <tlccsdt>`                  |
#    +-------------------------+---------------------------------------------------------------------------------------------------------------+
#    | eom-ccsd                | equation of motion (EOM) CCSD :ref:`[manual] <sec:eomcc>`                                                     |
#    +-------------------------+---------------------------------------------------------------------------------------------------------------+
#
#    .. _`table:grad_scf`:
#
#
#    .. include:: ../autodoc_dft_opt.rst
#
#    .. include:: ../cfour_table_grad.rst
#
#
#    :examples:
#
#    >>> # [1] Analytic hf optimization
#    >>> optimize('hf')
#
#    >>> # [2] Finite difference mp5 optimization with gradient
#    >>> #     printed to output file
#    >>> e, wfn = opt('mp5', return_wfn='yes')
#    >>> wfn.gradient().print_out()
#
#    >>> # [3] Forced finite difference hf optimization run in
#    >>> #     embarrassingly parallel fashion
#    >>> optimize('hf', dertype='energy', mode='sow')
#
#    >>> # [4] Can automatically perform complete basis set extrapolations
#    >>> optimize('MP2/cc-pV([D,T]+d)Z')
#
#    >>> # [5] Can automatically perform delta corrections that include extrapolations
#    >>> # even with a user-defined extrapolation formula. See sample inputs named
#    >>> # cbs-xtpl* for more examples of this input style
#    >>> optimize("MP2/aug-cc-pv([d,t]+d)z + d:ccsd(t)/cc-pvdz", corl_scheme=myxtplfn_2)
#
#    >>> # [6] Get info like geometry, gradient, energy back after an
#    >>> #     optimization fails. Note that the energy and gradient
#    >>> #     correspond to the last optimization cycle, whereas the
#    >>> #     geometry (by default) is the anticipated *next* optimization step.
#    >>> try:
#    >>>     optimize('hf/cc-pvtz')
#    >>> except psi4.OptimizationConvergenceError as ex:
#    >>>     next_geom_coords_as_numpy_array = np.asarray(ex.wfn.molecule().geometry())
#
#    """
def optking(name, **kwargs):
    from . import load_proc_table

    kwargs = driver_util.kwargs_lower(kwargs)

    if hasattr(name, "__call__"):
        lowername = name
    else:
        lowername = name.lower()

    return_wfn = kwargs.pop("return_wfn", False)

    import psi4

    psi4.core.clean()
    psi4.core.clean_options()
    # psi4.set_output_file("pytest_output.dat", True)

    #    return_history = kwargs.pop('return_history', False)
    #    if return_history:
    #        # Add wfn once the deep copy issues are worked out
    #        step_energies      = []
    #        step_gradients     = []
    #        step_coordinates   = []
    #
    #    # For CBS wrapper, need to set retention on INTCO file
    #    if custom_gradient or ('/' in lowername):
    #        psi4.core.IOManager.shared_object().set_specific_retention(1, True)

    #    if kwargs.get('bsse_type', None) is not None:
    #        raise ValidationError("Optimize: Does not currently support 'bsse_type' arguements")
    #
    #    full_hess_every = core.get_option('OPTKING', 'FULL_HESS_EVERY')
    #    steps_since_last_hessian = 0
    #
    #    if custom_gradient and core.has_option_changed('OPTKING', 'FULL_HESS_EVERY'):
    #        raise ValidationError("Optimize: Does not support custom Hessian's yet.")
    #    else:
    #        hessian_with_method = kwargs.get('hessian_with', lowername)
    #
    #    # are we in sow/reap mode?
    #    opt_mode = kwargs.get('mode', 'continuous').lower()
    #    if opt_mode not in ['continuous', 'sow', 'reap']:
    #        raise ValidationError("""Optimize execution mode '%s' not valid.""" % (opt_mode))
    #
    #    optstash = p4util.OptionsState(
    #        ['OPTKING', 'INTRAFRAG_STEP_LIMIT'],
    #        ['FINDIF', 'HESSIAN_WRITE'],
    #        ['OPTKING', 'CART_HESS_READ'],
    #        ['SCF', 'GUESS_PERSIST'],  # handle on behalf of cbs()
    #        ['SCF', 'GUESS'])

    iopt = kwargs.get("opt_iter", 1)

    # Make sure the molecule the user provided is the active one
    molecule = kwargs.pop("molecule", driver_helpers.get_active_molecule())
    molecule.update_geometry()

    if len(pe.nu_options.scroll) == 0:
        print("EMPTY OPT")
        pe.load_options()

    #    # If we are freezing cartesian, do not orient or COM
    #    if core.get_local_option("OPTKING", "FROZEN_CARTESIAN"):
    #        molecule.fix_orientation(True)
    #        molecule.fix_com(True)
    #    molecule.update_geometry()
    #
    #    # Shifting the geometry so need to copy the active molecule
    moleculeclone = molecule.clone()
    #
    #    initial_sym = moleculeclone.schoenflies_symbol()
    #    while iopt <= core.get_option('OPTKING', 'GEOM_MAXITER'):

    while iopt <= pe.nu_options.scroll["QCDB"]["GEOM_MAXITER"].value:
        #        current_sym = moleculeclone.schoenflies_symbol()
        #        if initial_sym != current_sym:
        #            raise ValidationError("""Point group changed! (%s <-- %s) You should restart """
        #                                  """using the last geometry in the output, after """
        #                                  """carefully making sure all symmetry-dependent """
        #                                  """input, such as DOCC, is correct.""" %
        #                                  (current_sym, initial_sym))
        #        kwargs['opt_iter'] = iopt
        #
        #        # Use orbitals from previous iteration as a guess
        #        #   set within loop so that can be influenced by fns to optimize (e.g., cbs)
        #        if (iopt > 1) and (opt_mode == 'continuous') and (not core.get_option('SCF', 'GUESS_PERSIST')):
        #            core.set_local_option('SCF', 'GUESS', 'READ')
        #
        #        # Before computing gradient, save previous molecule and wavefunction if this is an IRC optimization
        #        if (iopt > 1) and (core.get_option('OPTKING', 'OPT_TYPE') == 'IRC'):
        #            old_thisenergy = core.get_variable('CURRENT ENERGY')

        # Compute the gradient - preserve opt data despite core.clean calls in gradient
        #        psi4.core.IOManager.shared_object().set_specific_retention(1, True)
        G, jobrec = gradient(lowername, return_wfn=True, molecule=moleculeclone, **kwargs)
        thisenergy = float(jobrec["qcvars"]["CURRENT ENERGY"].data)

        #        # above, used to be getting energy as last of energy list from gradient()
        #        # thisenergy below should ultimately be testing on wfn.energy()
        #
        #        # Record optimization steps
        #        # Add wavefunctions later
        #        if return_history:
        #            step_energies.append(thisenergy)
        #            step_coordinates.append(moleculeclone.geometry())
        #            step_gradients.append(G.clone())
        #
        #        # S/R: Quit after getting new displacements or if forming gradient fails
        #        if opt_mode == 'sow':
        #            return (0.0, None)
        #        elif opt_mode == 'reap' and thisenergy == 0.0:
        #            return (0.0, None)

        psi4.core.set_variable("CURRENT ENERGY", thisenergy)
        psi4.core.set_gradient(psi4.core.Matrix.from_array(G))

        #        # S/R: Move opt data file from last pass into namespace for this pass
        #        if opt_mode == 'reap' and iopt != 0:
        #            core.IOManager.shared_object().set_specific_retention(1, True)
        #            core.IOManager.shared_object().set_specific_path(1, './')
        #            if 'opt_datafile' in kwargs:
        #                restartfile = kwargs.pop('opt_datafile')
        #                shutil.copy(restartfile, p4util.get_psifile(1))
        #
        #        # opt_func = kwargs.get('opt_func', kwargs.get('func', energy))
        #        # if opt_func.__name__ == 'complete_basis_set':
        #        #     core.IOManager.shared_object().set_specific_retention(1, True)
        #
        #        if full_hess_every > -1:
        #            core.set_global_option('HESSIAN_WRITE', True)
        #
        #        # compute Hessian as requested; frequency wipes out gradient so stash it
        #        if ((full_hess_every > -1) and (iopt == 1)) or (steps_since_last_hessian + 1 == full_hess_every):
        #            G = core.get_gradient()  # TODO
        #            core.IOManager.shared_object().set_specific_retention(1, True)
        #            core.IOManager.shared_object().set_specific_path(1, './')
        #            frequencies(hessian_with_method, **kwargs)
        #            steps_since_last_hessian = 0
        #            core.set_gradient(G)
        #            core.set_global_option('CART_HESS_READ', True)
        #        elif (full_hess_every == -1) and core.get_global_option('CART_HESS_READ') and (iopt == 1):
        #            pass
        #            # Do nothing; user said to read existing hessian once
        #        else:
        #            core.set_global_option('CART_HESS_READ', False)
        #            steps_since_last_hessian += 1

        popts = {}
        for k, v in pe.nu_options.scroll["QCDB"].items():
            if v.disputed() and k.endswith("G_CONVERGENCE"):
                popts[k] = v.value

        for k, v in pe.nu_options.scroll["PSI4"].items():
            if v.disputed() and k.endswith("G_CONVERGENCE"):
                popts[k] = v.value
        psi4.driver.p4util.python_helpers.set_options(popts)

        # Take step. communicate to/from/within optking through legacy_molecule
        psi4.core.set_legacy_molecule(psi4.core.Molecule.from_dict(moleculeclone.to_dict()))
        optking_rval = psi4.core.optking()
        moleculeclone = Molecule.from_dict(psi4.core.get_legacy_molecule().to_dict())
        moleculeclone.update_geometry()
        if optking_rval == psi4.core.PsiReturnType.EndLoop:
            #            # if this is the end of an IRC run, set wfn, energy, and molecule to that
            #            # of the last optimized IRC point
            #            if core.get_option('OPTKING', 'OPT_TYPE') == 'IRC':
            #                thisenergy = old_thisenergy
            print("Optimizer: Optimization complete!")
            #            core.print_out('\n    Final optimized geometry and variables:\n')
            #            moleculeclone.print_in_input_format()
            #            # Check if user wants to see the intcos; if so, don't delete them.
            #            if core.get_option('OPTKING', 'INTCOS_GENERATE_EXIT') == False:
            #                if core.get_option('OPTKING', 'KEEP_INTCOS') == False:
            psi4.core.opt_clean()
            #                    core.opt_clean()
            # Changing environment to optimized geometry as expected by user
            molecule.set_geometry(moleculeclone.geometry())
            #            for postcallback in hooks['optimize']['post']:
            #                postcallback(lowername, wfn=wfn, **kwargs)
            #            core.clean()
            psi4.core.clean()

            #            # S/R: Clean up opt input file
            #            if opt_mode == 'reap':
            #                with open('OPT-master.in', 'wb') as fmaster:
            #                    fmaster.write('# This is a psi4 input file auto-generated from the gradient() wrapper.\n\n'.encode('utf-8'))
            #                    fmaster.write('# Optimization complete!\n\n'.encode('utf-8'))
            #
            #            # Cleanup binary file 1
            #            if custom_gradient or ('/' in lowername):
            #                psi4.core.IOManager.shared_object().set_specific_retention(1, False)
            #
            #            optstash.restore()
            #
            #            if return_history:
            #                history = { 'energy'        : step_energies ,
            #                            'gradient'      : step_gradients ,
            #                            'coordinates'   : step_coordinates,
            #                          }

            #            if return_wfn and return_history:
            #                return (thisenergy, wfn, history)
            #            elif return_wfn and not return_history:
            #                return (thisenergy, wfn)
            #            elif return_history and not return_wfn:
            #                return (thisenergy, history)
            if return_wfn:
                return (thisenergy, jobrec)
            else:
                return thisenergy

        elif optking_rval == psi4.core.PsiReturnType.Failure:
            print("Optimizer: Optimization failed!")
            #            if (core.get_option('OPTKING', 'KEEP_INTCOS') == False):
            #                core.opt_clean()
            psi4.core.opt_clean()
            molecule.set_geometry(moleculeclone.geometry())
            psi4.core.clean()
            #            optstash.restore()
            #            raise OptimizationConvergenceError("""geometry optimization""", iopt - 1, wfn)
            return thisenergy

        print("\n    Structure for next step:\n")
        print(moleculeclone)
        #        moleculeclone.print_in_input_format()

        #        # S/R: Preserve opt data file for next pass and switch modes to get new displacements
        #        if opt_mode == 'reap':
        #            kwargs['opt_datafile'] = p4util.get_psifile(1)
        #            kwargs['mode'] = 'sow'

        iopt += 1

    #    if core.get_option('OPTKING', 'INTCOS_GENERATE_EXIT') == False:
    #        if core.get_option('OPTKING', 'KEEP_INTCOS') == False:
    #            core.opt_clean()
    psi4.core.IOManager.shared_object().set_specific_retention(1, False)
    psi4.core.clean()
    psi4.core.opt_clean()


#    optstash.restore()
#    raise OptimizationConvergenceError("""geometry optimization""", iopt - 1, wfn)

optimize = optking


def geometric(name, **kwargs):
    from . import load_proc_table

    kwargs = driver_util.kwargs_lower(kwargs)

    if hasattr(name, "__call__"):
        pass
    else:
        name.lower()

    return_wfn = kwargs.pop("return_wfn", False)

    # Make sure the molecule the user provided is the active one
    molecule = kwargs.pop("molecule", driver_helpers.get_active_molecule())
    molecule.update_geometry()

    if len(pe.nu_options.scroll) == 0:
        print("EMPTY OPT")
        pe.load_options()

    import yaml

    tricin = {}
    tricin["molecule"] = molecule.to_dict(np_out=False, force_units="Angstrom")
    cereal = {}
    cereal["driver"] = gradient
    cereal["method"] = name
    cereal["options"] = {}  # outrageous
    qwer = {"template": yaml.dump(cereal), "molecule": molecule.to_dict(np_out=False, force_units="Angstrom")}
    tricin["template"] = yaml.dump(qwer)
    tricin["argparse"] = ["--qcdb", "tmp_tricin.yaml"]  # prepare_options_for_geometric({'qcdb': True}).split()

    import geometric

    otricrec = geometric.json_run(tricin)
    # pp.pprint(otricrec)

    ngeom = otricrec["coords"].reshape((-1, 3))
    molecule.set_geometry(ngeom)

    # going to do an extra grad to get E, G, wfn until better tric comm
    G, jobrec = gradient(name, return_wfn=True, **kwargs)
    E = jobrec["qcvars"]["CURRENT ENERGY"].data

    # jobrec['qcvars'] = {info.lbl: info for info in calcinfo}
    pp.pprint(jobrec)
    pe.active_qcvars = copy.deepcopy(jobrec["qcvars"])

    if return_wfn:
        return (E, jobrec)
    else:
        return E


def prepare_options_for_geometric(dicary):

    sysargv = " ".join(["--" + k + " " + str(v) for k, v in dicary.items()])
    sysargv += " tmp_tricin.yaml"
    print("sysargv", sysargv)

    return sysargv


#    jrec = {}
#    jrec['molecule'] = mymol.to_dict(np_out=False)
#    jrec['template'] = tricin
#
#    pp.pprint(jrec)
#
#    with open('trictest2.json', 'w') as handle:
#        json.dump(jrec, handle)
#
#    #jrec.pop('template')
#
#    import yaml
#    with open('trictest3.yaml', 'w') as handle:
#        yaml.dump(jrec, handle)


# molecule:
#  elbl: ['', '', '', '', '', '', '', '', '']
#  elea: [16, 1, 1, 16, 1, 1, 16, 1, 1]
#  elem: [O, H, H, O, H, H, O, H, H]
#  elez: [8, 1, 1, 8, 1, 1, 8, 1, 1]
#  fix_com: false
#  fix_orientation: false
#  fragment_charges: [0.0]
#  fragment_multiplicities: [1]
#  fragment_separators: []
#  geom: [0.2647072205488637, -1.064967384509428, -0.932720644153007, 1.0840372205488638,
#    -1.5303473845094284, -0.764380644153007, 0.5264072205488637, -0.15103738450942847,
#    -1.0444306441530071, -1.6898327794511359, -0.5409773845094284, 0.9230093558469931,
#    -0.9943027794511362, -0.8286873845094284, 0.33165935584699274, -1.6763127794511363,
#    0.4142226154905714, 0.8626393558469929, 1.2732772205488636, 1.5375726154905711,
#    0.08292935584699303, 1.419267220548864, 2.2003126154905708, -0.5921206441530071,
#    2.0508472205488633, 0.9806526154905715, 0.044609355846992965]
#  mass: [15.99491461956, 1.00782503207, 1.00782503207, 15.99491461956, 1.00782503207,
#    1.00782503207, 15.99491461956, 1.00782503207, 1.00782503207]
#  molecular_charge: 0.0
#  molecular_multiplicity: 1
#  real: [true, true, true, true, true, true, true, true, true]
#  units: Angstrom
# template: "\ndriver: !!python/name:qcdb.gradient\nmethod: hf\noptions:\n  basis: cc-pvdz\n"


# PYTHONPATH=/home/psilocaluser/miniconda3/envs/idp35p4/lib/python3.5/site-packages:/home/psilocaluser/gits/geomeTRIC:/home/psilocaluser/miniconda3/envs/tric35/lib/python3.5/site-packages/:$PYTHONPATH /home/psilocaluser/gits/geomeTRIC/geometric/optimize.py --qcdb trictest3.yaml
