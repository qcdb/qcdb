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

from .. import moptions
from . import pe
from . import driver_util
from . import driver_helpers
from . import cbs_driver
#from psi4.driver import driver_nbody
from .proc_table import procedures


@moptions.register_opts(pe.nu_options)
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
   
    if 'options' in kwargs:
        driver_helpers.set_options(kwargs.pop('options'))

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

    if len(pe.nu_options.scroll) == 0:
        print('EMPTY OPT')
        pe.load_nu_options()


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

        #print('GRADIENT() JOBREC (j@io) <<<')
        #pp.pprint(jobrec)
        #print('>>>')

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
   
