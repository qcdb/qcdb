#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2017 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

"""Module with a *procedures* dictionary specifying available quantum
chemical methods and functions driving the main quantum chemical
functionality, namely single-point energies, geometry optimizations,
properties, and vibrational frequency calculations.

"""
from __future__ import print_function
from __future__ import absolute_import
#   import os
#   import re
#   import sys
#   import shutil

import copy
import pprint
pp = pprint.PrettyPrinter(width=120)

#   import numpy as np

from . import pe
from . import driver_util
from . import driver_helpers
from . import cbs_driver
##   from psi4.driver import driver_nbody
##   from psi4.driver import p4util
from .proc_table import procedures
##   from psi4.driver.procrouting import *
##   from psi4.driver.p4util.exceptions import *
##   # never import wrappers or aliases into this file


def gradient(name, **kwargs):
#       r"""Function complementary to :py:func:~driver.optimize(). Carries out one gradient pass,
#       deciding analytic or finite difference.
#   
#       :returns: :py:class:`~psi4.core.Matrix` |w--w| Total electronic gradient in Hartrees/Bohr.
#   
#       :returns: (:py:class:`~psi4.core.Matrix`, :py:class:`~psi4.core.Wavefunction`) |w--w| gradient and wavefunction when **return_wfn** specified.
#   
#       :examples:
#   
#       >>> # [1] Single-point dft gradient getting the gradient
#       >>> #     in file, core.Matrix, and np.array forms
#       >>> set gradient_write on
#       >>> G, wfn = gradient('b3lyp-d', return_wfn=True)
#       >>> wfn.gradient().print_out()
#       >>> np.array(G)
#   
#       """
    from . import endorsed_plugins
    kwargs = driver_util.kwargs_lower(kwargs)
    text = ''
   
#       # Bounce to CP if bsse kwarg (someday)
#       if kwargs.get('bsse_type', None) is not None:
#           raise ValidationError("Gradient: Cannot specify bsse_type for gradient yet.")
   
    # Figure out what kind of gradient this is
    if hasattr(name, '__call__'):
        if name.__name__ in ['cbs', 'complete_basis_set']:
            gradient_type = 'cbs_wrapper'
#       else:
#           # Bounce to name if name is non-CBS function
#           gradient_type = 'custom_function'
    elif '/' in name:
        gradient_type = 'cbs_gufunc'
    else:
        gradient_type = 'conventional'
   
##    lowername = name.lower()
##    package = driver_util.get_package(lowername, kwargs)

    # Figure out lowername, dertype, and func
    # If we have analytical gradients we want to pass to our wrappers, otherwise we want to run
    # finite-diference energy or cbs energies
#       # TODO MP5/cc-pv[DT]Z behavior unkown due to "levels"
    user_dertype = kwargs.pop('dertype', None)
    if gradient_type == 'custom_function':
        pass
#        if user_dertype is None:
#            dertype = 0
#            core.print_out("\nGradient: Custom function passed in without a defined dertype, assuming fd-energy based gradient.\n")
#        else:
#            core.print_out("\nGradient: Custom function passed in with a dertype of %d\n" % user_dertype)
#            dertype = user_dertype
#
#        if dertype == 1:
#            return name(gradient, kwargs.pop('label', 'custom function'), ptype='gradient', **kwargs)
#        else:
#            optstash = driver_util._set_convergence_criterion('energy', 'scf', 8, 10, 8, 10, 8)
#            lowername = name

#    elif gradient_type == 'cbs_wrapper':
#        cbs_methods = driver_cbs._cbs_wrapper_methods(**kwargs)
#        dertype = min([_find_derivative_type('gradient', method, user_dertype) for method in cbs_methods])
#        if dertype == 1:
#            # Bounce to CBS (directly) in pure-gradient mode if name is CBS and all parts have analytic grad. avail.
#            return name(gradient, kwargs.pop('label', 'custom function'), ptype='gradient', **kwargs)
#        else:
#            optstash = driver_util._set_convergence_criterion('energy', cbs_methods[0], 8, 10, 8, 10, 8)
#            lowername = name
#            # Pass through to G by E

    elif gradient_type == 'cbs_gufunc':
        cbs_methods = cbs_driver._parse_cbs_gufunc_string(name.lower())[0]
        dertype = min(driver_util.find_derivative_type('gradient', method, user_dertype, kwargs.get('package', None)) for method in cbs_methods)
        lowername = name.lower()
        molecule = kwargs.pop('molecule', driver_helpers.get_active_molecule())
        if dertype == 1:
            # Bounce to CBS in pure-gradient mode if "method/basis" name and all parts have analytic grad. avail.
            return cbs_driver._cbs_gufunc(gradient, name, ptype='gradient', molecule=molecule, **kwargs)
#        else:
#            # Set method-dependent scf convergence criteria (test on procedures['energy'] since that's guaranteed)
#            optstash = driver_util._set_convergence_criterion('energy', cbs_methods[0], 8, 10, 8, 10, 8)

    else:
        # Allow specification of methods to arbitrary order
        lowername = name.lower()
#        lowername, level = driver_util._parse_arbitrary_order(lowername)
#        if level:
#            kwargs['level'] = level

#        # Prevent methods that do not have associated gradients
#        if lowername in energy_only_methods:
#            raise ValidationError("gradient('%s') does not have an associated gradient" % name)

        dertype = driver_util.find_derivative_type('gradient', lowername, user_dertype, kwargs.get('package', None))

#        # Set method-dependent scf convergence criteria (test on procedures['energy'] since that's guaranteed)
#        optstash = driver_util._set_convergence_criterion('energy', lowername, 8, 10, 8, 10, 8)

    # Commit to procedures[] call hereafter
  #  lowername = name.lower()
    return_wfn = kwargs.pop('return_wfn', False)
    package = driver_util.get_package2(lowername, kwargs.get('package', None))
#    core.clean_variables()
#
#    # no analytic derivatives for scf_type cd
#    if core.get_option('SCF', 'SCF_TYPE') == 'CD':
#        if (dertype == 1):
#            raise ValidationError("""No analytic derivatives for SCF_TYPE CD.""")

    # Make sure the molecule the user provided is the active one
    molecule = kwargs.pop('molecule', driver_helpers.get_active_molecule())
    molecule.update_geometry()

    if len(pe.nu_options.scroll) == 0:
        print('EMPTY OPT')
        pe.load_nu_options()

#    # S/R: Mode of operation- whether finite difference opt run in one job or files farmed out
#    opt_mode = kwargs.get('mode', 'continuous').lower()
#    if opt_mode == 'continuous':
#        pass
#    elif opt_mode == 'sow':
#        if dertype == 1:
#            raise ValidationError("""Optimize execution mode 'sow' not valid for analytic gradient calculation.""")
#    elif opt_mode == 'reap':
#        opt_linkage = kwargs.get('linkage', None)
#        if opt_linkage is None:
#            raise ValidationError("""Optimize execution mode 'reap' requires a linkage option.""")
#    else:
#        raise ValidationError("""Optimize execution mode '%s' not valid.""" % (opt_mode))

    # Does dertype indicate an analytic procedure both exists and is wanted?
    if dertype == 1:
        text += """qcdb.gradient() will perform analytic gradient computation.\n"""

        # Perform the gradient calculation
        jobrec = procedures['gradient'][package][lowername](lowername, molecule=molecule, options=pe.nu_options, ptype='gradient', **kwargs)

        pp.pprint(jobrec)
        pe.active_qcvars = copy.deepcopy(jobrec['qcvars'])

#        optstash.restore()
        if return_wfn:
            return (jobrec['qcvars']['CURRENT GRADIENT'].data, jobrec)
        else:
            return jobrec['qcvars']['CURRENT GRADIENT'].data

    else:
        raise FeatureNotImplemented("""gradient(dertype=0)""")
#        core.print_out("""gradient() will perform gradient computation by finite difference of analytic energies.\n""")
#
#        opt_iter = kwargs.get('opt_iter', 1)
#        if opt_iter is True:
#            opt_iter = 1
#
#        if opt_iter == 1:
#            print('Performing finite difference calculations')
#
#        # Shifting the geometry so need to copy the active molecule
#        moleculeclone = molecule.clone()
#
#        # Obtain list of displacements
#        # print("about to generate displacements")
#        displacements = core.fd_geoms_1_0(moleculeclone)
#        # print(displacements)
#        ndisp = len(displacements)
#        # print("generated displacments")
#
#        # This version is pretty dependent on the reference geometry being last (as it is now)
#        print(""" %d displacements needed ...""" % (ndisp), end='')
#        energies = []
#
#        # S/R: Write instructions for sow/reap procedure to output file and reap input file
#        if opt_mode == 'sow':
#            instructionsO = """\n    The optimization sow/reap procedure has been selected through mode='sow'. In addition\n"""
#            instructionsO += """    to this output file (which contains no quantum chemical calculations), this job\n"""
#            instructionsO += """    has produced a number of input files (OPT-%s-*.in) for individual components\n""" % (str(opt_iter))
#            instructionsO += """    and a single input file (OPT-master.in) with an optimize(mode='reap') command.\n"""
#            instructionsO += """    These files may look very peculiar since they contain processed and pickled python\n"""
#            instructionsO += """    rather than normal input. Follow the instructions in OPT-master.in to continue.\n\n"""
#            instructionsO += """    Alternatively, a single-job execution of the gradient may be accessed through\n"""
#            instructionsO += """    the optimization wrapper option mode='continuous'.\n\n"""
#            core.print_out(instructionsO)
#
#            instructionsM = """\n#    Follow the instructions below to carry out this optimization cycle.\n#\n"""
#            instructionsM += """#    (1)  Run all of the OPT-%s-*.in input files on any variety of computer architecture.\n""" % (str(opt_iter))
#            instructionsM += """#       The output file names must be as given below.\n#\n"""
#            for rgt in range(ndisp):
#                pre = 'OPT-' + str(opt_iter) + '-' + str(rgt + 1)
#                instructionsM += """#             psi4 -i %-27s -o %-27s\n""" % (pre + '.in', pre + '.out')
#            instructionsM += """#\n#    (2)  Gather all the resulting output files in a directory. Place input file\n"""
#            instructionsM += """#         OPT-master.in into that directory and run it. The job will be minimal in\n"""
#            instructionsM += """#         length and give summary results for the gradient step in its output file.\n#\n"""
#            if opt_iter == 1:
#                instructionsM += """#             psi4 -i %-27s -o %-27s\n#\n""" % ('OPT-master.in', 'OPT-master.out')
#            else:
#                instructionsM += """#             psi4 -a -i %-27s -o %-27s\n#\n""" % ('OPT-master.in', 'OPT-master.out')
#            instructionsM += """#    After each optimization iteration, the OPT-master.in file is overwritten so return here\n"""
#            instructionsM += """#    for new instructions. With the use of the psi4 -a flag, OPT-master.out is not\n"""
#            instructionsM += """#    overwritten and so maintains a history of the job. To use the (binary) optimizer\n"""
#            instructionsM += """#    data file to accelerate convergence, the OPT-master jobs must run on the same computer.\n\n"""
#
#            with open('OPT-master.in', 'wb') as fmaster:
#                fmaster.write('# This is a psi4 input file auto-generated from the gradient() wrapper.\n\n'.encode('utf-8'))
#                fmaster.write(p4util.format_molecule_for_input(moleculeclone).encode('utf-8'))
#                fmaster.write(p4util.format_options_for_input().encode('utf-8'))
#                p4util.format_kwargs_for_input(fmaster, lmode=2, return_wfn=True, dertype=dertype, **kwargs)
#                fmaster.write(("""retE, retwfn = optimize('%s', **kwargs)\n\n""" % (lowername)).encode('utf-8'))
#                fmaster.write(instructionsM.encode('utf-8'))
#
#        for n, displacement in enumerate(displacements):
#            rfile = 'OPT-%s-%s' % (opt_iter, n + 1)
#
#            # Build string of title banner
#            banners = ''
#            banners += """core.print_out('\\n')\n"""
#            banners += """p4util.banner(' Gradient %d Computation: Displacement %d ')\n""" % (opt_iter, n + 1)
#            banners += """core.print_out('\\n')\n\n"""
#
#            if opt_mode == 'continuous':
#
#                # print progress to file and screen
#                core.print_out('\n')
#                p4util.banner('Loading displacement %d of %d' % (n + 1, ndisp))
#                print(""" %d""" % (n + 1), end=('\n' if (n + 1 == ndisp) else ''))
#                sys.stdout.flush()
#
#                # Load in displacement into the active molecule
#                moleculeclone.set_geometry(displacement)
#
#                # Perform the energy calculation
#                E, wfn = energy(lowername, return_wfn=True, molecule=moleculeclone, **kwargs)
#                energies.append(core.get_variable('CURRENT ENERGY'))
#
#            # S/R: Write each displaced geometry to an input file
#            elif opt_mode == 'sow':
#                moleculeclone.set_geometry(displacement)
#
#                # S/R: Prepare molecule, options, and kwargs
#                with open('%s.in' % (rfile), 'wb') as freagent:
#                    freagent.write('# This is a psi4 input file auto-generated from the gradient() wrapper.\n\n'.encode('utf-8'))
#                    freagent.write(p4util.format_molecule_for_input(moleculeclone).encode('utf-8'))
#                    freagent.write(p4util.format_options_for_input().encode('utf-8'))
#                    p4util.format_kwargs_for_input(freagent, **kwargs)
#
#                    # S/R: Prepare function call and energy save
#                    freagent.write(("""electronic_energy = energy('%s', **kwargs)\n\n""" % (lowername)).encode('utf-8'))
#                    freagent.write(("""core.print_out('\\nGRADIENT RESULT: computation %d for item %d """ % (os.getpid(), n + 1)).encode('utf-8'))
#                    freagent.write("""yields electronic energy %20.12f\\n' % (electronic_energy))\n\n""".encode('utf-8'))
#
#            # S/R: Read energy from each displaced geometry output file and save in energies array
#            elif opt_mode == 'reap':
#                exec(banners)
#                core.set_variable('NUCLEAR REPULSION ENERGY', moleculeclone.nuclear_repulsion_energy())
#                energies.append(p4util.extract_sowreap_from_output(rfile, 'GRADIENT', n, opt_linkage, True))
#
#        # S/R: Quit sow after writing files. Initialize skeleton wfn to receive grad for reap
#        if opt_mode == 'sow':
#            optstash.restore()
#            if return_wfn:
#                return (None, None)  # any point to building a dummy wfn here?
#            else:
#                return None
#        elif opt_mode == 'reap':
#            core.set_variable('CURRENT ENERGY', energies[-1])
#            wfn = core.Wavefunction.build(molecule, core.get_global_option('BASIS'))
#
#        # Compute the gradient; last item in 'energies' is undisplaced
#        core.set_local_option('FINDIF', 'GRADIENT_WRITE', True)
#        G = core.fd_1_0(molecule, energies)
#        G.print_out()
#        wfn.set_gradient(G)
#
#        optstash.restore()
#
#        if return_wfn:
#            return (wfn.gradient(), wfn)
#        else:
#            return wfn.gradient()
   
   
#def optimize(name, **kwargs):
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
#       * :psivar:`CURRENT ENERGY <CURRENTENERGY>`
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
#    | dcft                    | density cumulant functional theory :ref:`[manual] <sec:dcft>`                                                 |
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
#    kwargs = p4util.kwargs_lower(kwargs)
#
#    if hasattr(name, '__call__'):
#        lowername = name
#        custom_gradient = True
#    else:
#        lowername = name.lower()
#        custom_gradient = False
#
#    return_wfn = kwargs.pop('return_wfn', False)
#
#    return_history = kwargs.pop('return_history', False)
#    if return_history:
#        # Add wfn once the deep copy issues are worked out
#        step_energies      = []
#        step_gradients     = []
#        step_coordinates   = []
#
#    # For CBS wrapper, need to set retention on INTCO file
#    if custom_gradient or ('/' in lowername):
#        core.IOManager.shared_object().set_specific_retention(1, True)
#
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
#
#    n = kwargs.get('opt_iter', 1)
#
#    # Make sure the molecule the user provided is the active one
#    molecule = kwargs.pop('molecule', core.get_active_molecule())
#
#    # If we are freezing cartesian, do not orient or COM
#    if core.get_local_option("OPTKING", "FROZEN_CARTESIAN"):
#        molecule.fix_orientation(True)
#        molecule.fix_com(True)
#    molecule.update_geometry()
#
#    # Shifting the geometry so need to copy the active molecule
#    moleculeclone = molecule.clone()
#
#    initial_sym = moleculeclone.schoenflies_symbol()
#    while n <= core.get_option('OPTKING', 'GEOM_MAXITER'):
#        current_sym = moleculeclone.schoenflies_symbol()
#        if initial_sym != current_sym:
#            raise ValidationError("""Point group changed! (%s <-- %s) You should restart """
#                                  """using the last geometry in the output, after """
#                                  """carefully making sure all symmetry-dependent """
#                                  """input, such as DOCC, is correct.""" %
#                                  (current_sym, initial_sym))
#        kwargs['opt_iter'] = n
#
#        # Use orbitals from previous iteration as a guess
#        #   set within loop so that can be influenced by fns to optimize (e.g., cbs)
#        if (n > 1) and (opt_mode == 'continuous') and (not core.get_option('SCF', 'GUESS_PERSIST')):
#            core.set_local_option('SCF', 'GUESS', 'READ')
#
#        # Before computing gradient, save previous molecule and wavefunction if this is an IRC optimization
#        if (n > 1) and (core.get_option('OPTKING', 'OPT_TYPE') == 'IRC'):
#            old_thisenergy = core.get_variable('CURRENT ENERGY')
#
#        # Compute the gradient
#        G, wfn = gradient(lowername, return_wfn=True, molecule=moleculeclone, **kwargs)
#        thisenergy = core.get_variable('CURRENT ENERGY')
#
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
#
#        core.set_gradient(G)
#
#        # S/R: Move opt data file from last pass into namespace for this pass
#        if opt_mode == 'reap' and n != 0:
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
#        if ((full_hess_every > -1) and (n == 1)) or (steps_since_last_hessian + 1 == full_hess_every):
#            G = core.get_gradient()  # TODO
#            core.IOManager.shared_object().set_specific_retention(1, True)
#            core.IOManager.shared_object().set_specific_path(1, './')
#            frequencies(hessian_with_method, **kwargs)
#            steps_since_last_hessian = 0
#            core.set_gradient(G)
#            core.set_global_option('CART_HESS_READ', True)
#        elif (full_hess_every == -1) and core.get_global_option('CART_HESS_READ') and (n == 1):
#            pass
#            # Do nothing; user said to read existing hessian once
#        else:
#            core.set_global_option('CART_HESS_READ', False)
#            steps_since_last_hessian += 1
#
#        # Take step. communicate to/from/within optking through legacy_molecule
#        core.set_legacy_molecule(moleculeclone)
#        optking_rval = core.optking()
#        moleculeclone = core.get_legacy_molecule()
#        moleculeclone.update_geometry()
#        if optking_rval == core.PsiReturnType.EndLoop:
#            # if this is the end of an IRC run, set wfn, energy, and molecule to that
#            # of the last optimized IRC point
#            if core.get_option('OPTKING', 'OPT_TYPE') == 'IRC':
#                thisenergy = old_thisenergy
#            print('Optimizer: Optimization complete!')
#            core.print_out('\n    Final optimized geometry and variables:\n')
#            moleculeclone.print_in_input_format()
#            # Check if user wants to see the intcos; if so, don't delete them.
#            if core.get_option('OPTKING', 'INTCOS_GENERATE_EXIT') == False:
#                if core.get_option('OPTKING', 'KEEP_INTCOS') == False:
#                    core.opt_clean()
#            # Changing environment to optimized geometry as expected by user
#            molecule.set_geometry(moleculeclone.geometry())
#            for postcallback in hooks['optimize']['post']:
#                postcallback(lowername, wfn=wfn, **kwargs)
#            core.clean()
#
#            # S/R: Clean up opt input file
#            if opt_mode == 'reap':
#                with open('OPT-master.in', 'wb') as fmaster:
#                    fmaster.write('# This is a psi4 input file auto-generated from the gradient() wrapper.\n\n'.encode('utf-8'))
#                    fmaster.write('# Optimization complete!\n\n'.encode('utf-8'))
#
#            # Cleanup binary file 1
#            if custom_gradient or ('/' in lowername):
#                core.IOManager.shared_object().set_specific_retention(1, False)
#
#            optstash.restore()
#
#            if return_history:
#                history = { 'energy'        : step_energies ,
#                            'gradient'      : step_gradients ,
#                            'coordinates'   : step_coordinates,
#                          }
#
#            if return_wfn and return_history:
#                return (thisenergy, wfn, history)
#            elif return_wfn and not return_history:
#                return (thisenergy, wfn)
#            elif return_history and not return_wfn:
#                return (thisenergy, history)
#            else:
#                return thisenergy
#
#        elif optking_rval == core.PsiReturnType.Failure:
#            print('Optimizer: Optimization failed!')
#            if (core.get_option('OPTKING', 'KEEP_INTCOS') == False):
#                core.opt_clean()
#            molecule.set_geometry(moleculeclone.geometry())
#            core.clean()
#            optstash.restore()
#            raise OptimizationConvergenceError("""geometry optimization""", n - 1, wfn)
#            return thisenergy
#
#        core.print_out('\n    Structure for next step:\n')
#        moleculeclone.print_in_input_format()
#
#        # S/R: Preserve opt data file for next pass and switch modes to get new displacements
#        if opt_mode == 'reap':
#            kwargs['opt_datafile'] = p4util.get_psifile(1)
#            kwargs['mode'] = 'sow'
#
#        n += 1
#
#    if core.get_option('OPTKING', 'INTCOS_GENERATE_EXIT') == False:
#        if core.get_option('OPTKING', 'KEEP_INTCOS') == False:
#            core.opt_clean()
#
#    optstash.restore()
#    raise OptimizationConvergenceError("""geometry optimization""", n - 1, wfn)
#
#
