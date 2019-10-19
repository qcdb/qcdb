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
import copy
import pprint
pp = pprint.PrettyPrinter(width=120)

from . import pe  # keep this at top of imports
from ..keyword import register_kwds
from . import driver_util
from . import driver_helpers
from . import cbs_driver
#from psi4.driver import driver_nbody
from . proc_table import procedures

   
#def properties(name, molecule, **kwargs):
#    r"""Function to compute the single-point electronic properties."""

#       :returns: *float* |w--w| Total electronic properties in Hartrees. SAPT & EFP return interaction energy.
#   
#       :returns: (*float*, :py:class:`~psi4.core.Wavefunction`) |w--w| properties and wavefunction when **return_wfn** specified.
#   
#       :PSI variables:
#   
#       .. hlist::
#          :columns: 1
#   
#          * :psivar:`CURRENT ENERGY <CURRENTENERGY>`
#          * :psivar:`CURRENT REFERENCE ENERGY <CURRENTREFERENCEENERGY>`
#          * :psivar:`CURRENT CORRELATION ENERGY <CURRENTCORRELATIONENERGY>`
#   
#       :type name: string
#       :param name: ``'scf'`` || ``'mp2'`` || ``'ci5'`` || etc.
#   
#           First argument, usually unlabeled. Indicates the computational method
#           to be applied to the system.
#   
#       :type molecule: :ref:`molecule <op_py_molecule>`
#       :param molecule: ``h2o`` || etc.
#   
#           The target molecule, if not the last molecule defined.
#   
#       :type return_wfn: :ref:`boolean <op_py_boolean>`
#       :param return_wfn: ``'on'`` || |dl| ``'off'`` |dr|
#   
#           Indicate to additionally return the :py:class:`~psi4.core.Wavefunction`
#           calculation result as the second element (after *float* properties) of a tuple.
#   
#       :type restart_file: string
#       :param restart_file: ``['file.1, file.32]`` || ``./file`` || etc.
#   
#           Binary data files to be renamed for calculation restart.
#   
#       .. _`table:properties_gen`:
#   
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | name                    | calls method                                                                                                  |
#       +=========================+===============================================================================================================+
#       | efp                     | effective fragment potential (EFP) :ref:`[manual] <sec:libefp>`                                               |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | scf                     | Hartree--Fock (HF) or density functional theory (DFT) :ref:`[manual] <sec:scf>`                               |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | hf                      | HF self consistent field (SCF) :ref:`[manual] <sec:scf>`                                                      |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | hf3c                    | HF with dispersion, BSSE, and basis set corrections :ref:`[manual] <sec:gcp>`                                 |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | pbeh3c                  | PBEh with dispersion, BSSE, and basis set corrections :ref:`[manual] <sec:gcp>`                               |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | dcft                    | density cumulant functional theory :ref:`[manual] <sec:dcft>`                                                 |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | mp2                     | 2nd-order |MollerPlesset| perturbation theory (MP2) :ref:`[manual] <sec:dfmp2>` :ref:`[details] <tlmp2>`      |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | mp3                     | 3rd-order |MollerPlesset| perturbation theory (MP3) :ref:`[manual] <sec:occ_nonoo>` :ref:`[details] <tlmp3>`  |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | fno-mp3                 | MP3 with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                                  |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | mp2.5                   | average of MP2 and MP3 :ref:`[manual] <sec:occ_nonoo>` :ref:`[details] <tlmp25>`                              |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | mp4(sdq)                | 4th-order MP perturbation theory (MP4) less triples :ref:`[manual] <sec:fnompn>`                              |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | fno-mp4(sdq)            | MP4 (less triples) with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                   |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | mp4                     | full MP4 :ref:`[manual] <sec:fnompn>` :ref:`[details] <tlmp4>`                                                |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | fno-mp4                 | full MP4 with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                             |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | mp\ *n*                 | *n*\ th-order |MollerPlesset| (MP) perturbation theory :ref:`[manual] <sec:arbpt>`                            |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | zapt\ *n*               | *n*\ th-order z-averaged perturbation theory (ZAPT) :ref:`[manual] <sec:arbpt>`                               |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | omp2                    | orbital-optimized second-order MP perturbation theory :ref:`[manual] <sec:occ_oo>`                            |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | scs-omp2                | spin-component scaled OMP2 :ref:`[manual] <sec:occ_oo>`                                                       |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | scs(n)-omp2             | a special version of SCS-OMP2 for nucleobase interactions :ref:`[manual] <sec:occ_oo>`                        |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | scs-omp2-vdw            | a special version of SCS-OMP2 (from ethene dimers) :ref:`[manual] <sec:occ_oo>`                               |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | sos-omp2                | spin-opposite scaled OMP2 :ref:`[manual] <sec:occ_oo>`                                                        |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | sos-pi-omp2             | A special version of SOS-OMP2 for pi systems :ref:`[manual] <sec:occ_oo>`                                     |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | omp3                    | orbital-optimized third-order MP perturbation theory :ref:`[manual] <sec:occ_oo>`                             |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | scs-omp3                | spin-component scaled OMP3 :ref:`[manual] <sec:occ_oo>`                                                       |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | scs(n)-omp3             | a special version of SCS-OMP3 for nucleobase interactions :ref:`[manual] <sec:occ_oo>`                        |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | scs-omp3-vdw            | a special version of SCS-OMP3 (from ethene dimers) :ref:`[manual] <sec:occ_oo>`                               |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | sos-omp3                | spin-opposite scaled OMP3 :ref:`[manual] <sec:occ_oo>`                                                        |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | sos-pi-omp3             | A special version of SOS-OMP3 for pi systems :ref:`[manual] <sec:occ_oo>`                                     |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | omp2.5                  | orbital-optimized MP2.5 :ref:`[manual] <sec:occ_oo>`                                                          |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | lccsd, cepa(0)          | coupled electron pair approximation variant 0 :ref:`[manual] <sec:fnocepa>` :ref:`[details] <tllccsd>`        |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | fno-lccsd, fno-cepa(0)  | CEPA(0) with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                              |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | cepa(1)                 | coupled electron pair approximation variant 1 :ref:`[manual] <sec:fnocepa>`                                   |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | fno-cepa(1)             | CEPA(1) with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                              |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | cepa(3)                 | coupled electron pair approximation variant 3 :ref:`[manual] <sec:fnocepa>`                                   |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | fno-cepa(3)             | CEPA(3) with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                              |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | acpf                    | averaged coupled-pair functional :ref:`[manual] <sec:fnocepa>`                                                |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | fno-acpf                | ACPF with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                                 |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | aqcc                    | averaged quadratic coupled cluster :ref:`[manual] <sec:fnocepa>`                                              |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | fno-aqcc                | AQCC with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                                 |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | qcisd                   | quadratic CI singles doubles (QCISD) :ref:`[manual] <sec:fnocc>`                                              |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | fno-qcisd               | QCISD with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                                |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | lccd                    | Linear CCD :ref:`[manual] <sec:occ_nonoo>` :ref:`[details] <tllccd>`                                          |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | fno-lccd                | LCCD with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                                 |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | olccd                   | orbital optimized LCCD :ref:`[manual] <sec:occ_oo>`                                                           |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | cc2                     | approximate coupled cluster singles and doubles (CC2) :ref:`[manual] <sec:cc>`                                |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | ccd                     | coupled cluster doubles  (CCD) :ref:`[manual] <sec:occ_nonoo>`                                                |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | ccsd                    | coupled cluster singles and doubles (CCSD) :ref:`[manual] <sec:cc>` :ref:`[details] <tlccsd>`                 |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | bccd                    | Brueckner coupled cluster doubles (BCCD) :ref:`[manual] <sec:cc>`                                             |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | fno-ccsd                | CCSD with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                                 |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | qcisd(t)                | QCISD with perturbative triples :ref:`[manual] <sec:fnocc>`                                                   |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | fno-qcisd(t)            | QCISD(T) with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                             |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | ccsd(t)                 | CCSD with perturbative triples (CCSD(T)) :ref:`[manual] <sec:cc>` :ref:`[details] <tlccsdt>`                  |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | ccsd(at)                | CCSD with asymmetric perturbative triples (CCSD(AT)) :ref:`[manual] <sec:cc>` :ref:`[details] <tlccsdat>`     |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | bccd(t)                 | BCCD with perturbative triples :ref:`[manual] <sec:cc>`                                                       |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | fno-ccsd(t)             | CCSD(T) with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                              |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | cc3                     | approximate CC singles, doubles, and triples (CC3) :ref:`[manual] <sec:cc>`                                   |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | ccproperties                | **expert** full control over ccenergy module                                                                  |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | dfocc                   | **expert** full control over dfocc module                                                                     |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | cisd                    | configuration interaction (CI) singles and doubles (CISD) :ref:`[manual] <sec:ci>` :ref:`[details] <tlcisd>`  |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | fno-cisd                | CISD with frozen natural orbitals :ref:`[manual] <sec:fnocc>`                                                 |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | cisdt                   | CI singles, doubles, and triples (CISDT) :ref:`[manual] <sec:ci>`                                             |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | cisdtq                  | CI singles, doubles, triples, and quadruples (CISDTQ) :ref:`[manual] <sec:ci>`                                |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | ci\ *n*                 | *n*\ th-order CI :ref:`[manual] <sec:ci>`                                                                     |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | fci                     | full configuration interaction (FCI) :ref:`[manual] <sec:ci>`                                                 |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | detci                   | **expert** full control over detci module                                                                     |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | casscf                  | complete active space self consistent field (CASSCF)  :ref:`[manual] <sec:ci>`                                |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | rasscf                  | restricted active space self consistent field (RASSCF)  :ref:`[manual] <sec:ci>`                              |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | mcscf                   | multiconfigurational self consistent field (SCF) :ref:`[manual] <sec:psimrcc>`                                |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | psimrcc                 | Mukherjee multireference coupled cluster (Mk-MRCC) :ref:`[manual] <sec:psimrcc>`                              |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | dmrg-scf                | density matrix renormalization group SCF :ref:`[manual] <sec:chemps2>`                                        |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | dmrg-caspt2             | density matrix renormalization group CASPT2 :ref:`[manual] <sec:chemps2>`                                     |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | dmrg-ci                 | density matrix renormalization group CI :ref:`[manual] <sec:chemps2>`                                         |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | sapt0                   | 0th-order symmetry adapted perturbation theory (SAPT) :ref:`[manual] <sec:sapt>`                              |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | ssapt0                  | 0th-order SAPT with special exchange scaling :ref:`[manual] <sec:sapt>`                                       |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | fisapt0                 | 0th-order functional and/or intramolecular SAPT :ref:`[manual] <sec:fisapt>`                                  |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | sapt2                   | 2nd-order SAPT, traditional definition :ref:`[manual] <sec:sapt>`                                             |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | sapt2+                  | SAPT including all 2nd-order terms :ref:`[manual] <sec:sapt>`                                                 |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | sapt2+(3)               | SAPT including perturbative triples :ref:`[manual] <sec:sapt>`                                                |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | sapt2+3                 | SAPT including all 3rd-order terms :ref:`[manual] <sec:sapt>`                                                 |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | sapt2+(ccd)             | SAPT2+ with CC-based dispersion :ref:`[manual] <sec:sapt>`                                                    |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | sapt2+(3)(ccd)          | SAPT2+(3) with CC-based dispersion :ref:`[manual] <sec:sapt>`                                                 |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | sapt2+3(ccd)            | SAPT2+3 with CC-based dispersion :ref:`[manual] <sec:sapt>`                                                   |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | sapt2+dmp2              | SAPT including all 2nd-order terms and MP2 correction :ref:`[manual] <sec:sapt>`                              |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | sapt2+(3)dmp2           | SAPT including perturbative triples and MP2 correction :ref:`[manual] <sec:sapt>`                             |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | sapt2+3dmp2             | SAPT including all 3rd-order terms and MP2 correction :ref:`[manual] <sec:sapt>`                              |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | sapt2+(ccd)dmp2         | SAPT2+ with CC-based dispersion and MP2 correction :ref:`[manual] <sec:sapt>`                                 |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | sapt2+(3)(ccd)dmp2      | SAPT2+(3) with CC-based dispersion and MP2 correction :ref:`[manual] <sec:sapt>`                              |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | sapt2+3(ccd)dmp2        | SAPT2+3 with CC-based dispersion and MP2 correction :ref:`[manual] <sec:sapt>`                                |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | sapt0-ct                | 0th-order SAPT plus charge transfer (CT) calculation :ref:`[manual] <sec:saptct>`                             |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | sapt2-ct                | SAPT2 plus CT :ref:`[manual] <sec:saptct>`                                                                    |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | sapt2+-ct               | SAPT2+ plus CT :ref:`[manual] <sec:saptct>`                                                                   |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | sapt2+(3)-ct            | SAPT2+(3) plus CT :ref:`[manual] <sec:saptct>`                                                                |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | sapt2+3-ct              | SAPT2+3 plus CT :ref:`[manual] <sec:saptct>`                                                                  |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | sapt2+(ccd)-ct          | SAPT2+(CCD) plus CT :ref:`[manual] <sec:saptct>`                                                              |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | sapt2+(3)(ccd)-ct       | SAPT2+(3)(CCD) plus CT :ref:`[manual] <sec:saptct>`                                                           |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | sapt2+3(ccd)-ct         | SAPT2+3(CCD) plus CT :ref:`[manual] <sec:saptct>`                                                             |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | adc                     | 2nd-order algebraic diagrammatic construction (ADC) :ref:`[manual] <sec:adc>`                                 |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | eom-cc2                 | EOM-CC2 :ref:`[manual] <sec:eomcc>`                                                                           |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | eom-ccsd                | equation of motion (EOM) CCSD :ref:`[manual] <sec:eomcc>`                                                     |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#       | eom-cc3                 | EOM-CC3 :ref:`[manual] <sec:eomcc>`                                                                           |
#       +-------------------------+---------------------------------------------------------------------------------------------------------------+
#   
#       .. comment missing and why
#       .. comment a certain isapt --- marginally released
#       .. comment mrcc --- this is handled in its own table
#       .. comment psimrcc_scf --- convenience fn
#   
#       .. include:: ../autodoc_dft_properties.rst
#   
#       .. include:: ../mrcc_table_properties.rst
#   
#       .. include:: ../cfour_table_properties.rst
#   
#       :examples:
#   
#       >>> # [1] Coupled-cluster singles and doubles calculation with psi code
#       >>> properties('ccsd')
#   
#       >>> # [2] Charge-transfer SAPT calculation with scf projection from small into
#       >>> #     requested basis, with specified projection fitting basis
#       >>> set basis_guess true
#       >>> set df_basis_guess jun-cc-pVDZ-JKFIT
#       >>> properties('sapt0-ct')
#   
#       >>> # [3] Arbitrary-order MPn calculation
#       >>> properties('mp7')
#   
#       >>> # [4] Converge scf as singlet, then run detci as triplet upon singlet reference
#       >>> # Note that the integral transformation is not done automatically when detci is run in a separate step.
#       >>> molecule H2 {\n0 1\nH\nH 1 0.74\n}
#       >>> set basis cc-pVDZ
#       >>> set reference rohf
#       >>> scf_e, scf_wfn = properties('scf', return_wfn=True)
#       >>> H2.set_multiplicity(3)
#       >>> core.MintsHelper(scf_wfn.basisset()).integrals()
#       >>> properties('detci', ref_wfn=scf_wfn)
#   
#       >>> # [5] Run two CI calculations, keeping the integrals generated in the first one.
#       >>> molecule ne {\nNe\n}
#       >>> set basis cc-pVDZ
#       >>> cisd_e, cisd_wfn = properties('cisd', return_wfn=True)
#       >>> properties('fci', ref_wfn=cisd_wfn)
#   
#       >>> # [6] Can automatically perform complete basis set extrapolations
#       >>> properties("CCSD/cc-pV[DT]Z")
#   
#       >>> # [7] Can automatically perform delta corrections that include extrapolations
#       >>> # even with a user-defined extrapolation formula. See sample inputs named
#       >>> # cbs-xtpl* for more examples of this input style
#       >>> properties("MP2/aug-cc-pv([d,t]+d)z + d:ccsd(t)/cc-pvdz", corl_scheme=myxtplfn_2)
#   
#       """

@register_kwds(pe.nu_options)
def properties(name, **kwargs):
    r"""Function to compute the single-point electronic properties."""

    from . import endorsed_plugins
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
        pe.load_nu_options()


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
def properties(*args, **kwargs):
    r"""Function to compute various properties."""

    from . import endorsed_plugins
    kwargs = driver_util.kwargs_lower(kwargs)

    # Make sure the molecule the user provided is the active one
    molecule = kwargs.pop('molecule', driver_helpers.get_active_molecule())
    molecule.update_geometry()

    if len(pe.nu_options.scroll) == 0:
        #print('EMPTY OPT')
        pe.load_nu_options()

    # Allow specification of methods to arbitrary order
    lowername = args[0].lower()
    lowername, level = driver_helpers._parse_arbitrary_order(lowername)
    if level:
        kwargs['level'] = level

    # Commit to procedures['properties'] call hereafter
    return_wfn = kwargs.pop('return_wfn', False)
    pe.active_qcvars = {}

    properties = kwargs.get('properties', ['dipole', 'quadrupole', 'octupole']) #nwc offer octupole
    if len(args) > 1:
        properties += args[1:]
    kwargs['properties'] = list(set(properties))
#    kwargs['properties'] = p4util.drop_duplicates(properties)

    package = driver_util.get_package(lowername, kwargs)
#    optstash = driver_util._set_convergence_criterion('properties', lowername, 6, 10, 6, 10, 8)
#    wfn = procedures['properties'][lowername](lowername, **kwargs)
    jobrec = procedures['properties'][package][lowername](lowername, molecule=molecule, options=pe.nu_options, ptype='properties', **kwargs)

    pp.pprint(jobrec)
    pe.active_qcvars = copy.deepcopy(jobrec['qcvars'])

#    if return_wfn:
#        return (float(jobrec['qcvars']['CURRENT ENERGY'].data), jobrec)
#    else:
#        return float(jobrec['qcvars']['CURRENT ENERGY'].data)


#def gdma(wfn, datafile=""):
#    """Function to use wavefunction information in *wfn* and, if specified,
#    additional commands in *filename* to run GDMA analysis.
#
#    .. include:: ../autodoc_abbr_options_c.rst
#
#    .. versionadded:: 0.6
#
#    :returns: None
#
#    :type wfn: :py:class:`~psi4.core.Wavefunction`
#    :param wfn: set of molecule, basis, orbitals from which to generate DMA analysis
#
#    :type datafile: string
#    :param datafile: optional control file (see GDMA manual) to peform more complicated DMA
#                     analyses.  If this option is used, the File keyword must be set to read
#                     a filename.fchk, where filename is provided by |globals__writer_file_label| .
#
#    :examples:
#
#    >>> # [1] DMA analysis from MP2 wavefunction.  N.B. gradient must be requested to generate MP2 density.
#    >>> grad, wfn = gradient('mp2', return_wfn=True)
#    >>> gdma(wfn)
#
#    """
#    # Start by writing a G* checkpoint file, for the GDMA code to read in
#    fw = core.FCHKWriter(wfn)
#    molname = wfn.molecule().name()
#    prefix = core.get_writer_file_prefix(molname)
#    fchkfile = prefix + '.fchk'
#    fw.write(fchkfile)
#
#    if datafile:
#        commands = datafile
#    else:
#        densname = wfn.name()
#        if densname == "DFT":
#            densname = "SCF"
#        commands = 'psi4_dma_datafile.dma'
#        radii = core.get_option('GDMA', 'GDMA_RADIUS')
#        origin = core.get_option('GDMA', 'GDMA_ORIGIN')
#        with open(commands, 'w') as f:
#            f.write("File %s Density %s\n" % (fchkfile, densname))
#            f.write("Angstrom\n")
#            f.write("%s\n" % core.get_option('GDMA', 'GDMA_MULTIPOLE_UNITS'))
#            f.write("Multipoles\n")
#            if origin:
#                try:
#                    f.write("Origin %f %f %f\n" % (float(origin[0]), float(origin[1]), float(origin[2])))
#                except:
#                    raise ValidationError("The GDMA origin array should contain three entries: x, y, and z.")
#            f.write("Switch %f\n" % core.get_option('GDMA', 'GDMA_SWITCH'))
#            if radii:
#                f.write("Radius %s\n" % " ".join([str(r) for r in radii]))
#            f.write("Limit %d\n" % core.get_option('GDMA', 'GDMA_LIMIT'))
#            f.write("Start\n")
#            f.write("Finish\n")
#    core.run_gdma(wfn, commands)
#
#    os.remove(fchkfile)
#    # If we generated the DMA control file, we should clean up here
#    if not datafile:
#        os.remove(commands)
#
#
#def fchk(wfn, filename):
#    """Function to write wavefunction information in *wfn* to *filename* in
#    Gaussian FCHK format.
#
#    .. versionadded:: 0.6
#
#    :returns: None
#
#    :type filename: string
#    :param filename: destination file name for FCHK file
#
#    :type wfn: :py:class:`~psi4.core.Wavefunction`
#    :param wfn: set of molecule, basis, orbitals from which to generate fchk file
#
#    :examples:
#
#    >>> # [1] FCHK file for DFT calculation
#    >>> E, wfn = properties('b3lyp', return_wfn=True)
#    >>> fchk(wfn, 'mycalc.fchk')
#
#    """
#    fw = core.FCHKWriter(wfn)
#    fw.write(filename)
#
#
#def molden(wfn, filename=None, density_a=None, density_b=None, dovirtual=None):
#    """Function to write wavefunction information in *wfn* to *filename* in
#    molden format. Will write natural orbitals from *density* (MO basis) if supplied.
#    Warning! Most post-SCF Wavefunctions do not build the density as this is often
#    much more costly than the properties. In addition, the Wavefunction density attributes
#    (Da and Db) return the SO density and must be transformed to the MO basis
#    to use with this function.
#
#    .. versionadded:: 0.5
#       *wfn* parameter passed explicitly
#
#    :returns: None
#
#    :type wfn: :py:class:`~psi4.core.Wavefunction`
#    :param wfn: set of molecule, basis, orbitals from which to generate cube files
#
#    :type filename: string
#    :param filename: destination file name for MOLDEN file (optional)
#
#    :type density_a: :py:class:`~psi4.core.Matrix`
#    :param density_a: density in the MO basis to build alpha NO's from (optional)
#
#    :type density_b: :py:class:`~psi4.core.Matrix`
#    :param density_b: density in the MO basis to build beta NO's from, assumes restricted if not supplied (optional)
#
#    :type dovirtual: bool
#    :param dovirtual: do write all the MOs to the MOLDEN file (true) or discard the unoccupied MOs, not valid for NO's (false) (optional)
#
#    :examples:
#
#    >>> # [1] Molden file for DFT calculation
#    >>> E, wfn = properties('b3lyp', return_wfn=True)
#    >>> molden(wfn, 'mycalc.molden')
#
#    >>> # [2] Molden file for CI/MCSCF computation using NO roots
#    >>> E, wfn = properties('ci', return_wfn=True)
#    >>> molden(wfn, 'no_root1.molden', density_a=wfn.opdm(0, 0, "A", True))
#
#    >>> # [3] The following does NOT work, please see below
#    >>> E, wfn = properties('ccsd', return_wfn=True)
#    >>> molden(wfn, 'ccsd_no.molden', density_a=wfn.Da())
#
#    >>> # [4] This WILL work, note the transformation of Da (SO->MO)
#    >>> E, wfn = properties('ccsd', properties=['dipole'], return_wfn=True)
#    >>> Da_so = wfn.Da()
#    >>> Da_mo = Matrix.triplet(wfn.Ca(), Da_so, wfn.Ca(), True, False, False)
#    >>> molden(wfn, 'ccsd_no.molden', density_a=Da_mo)
#
#    """
#
#    if filename is None:
#        filename = core.get_writer_file_prefix(wfn.molecule().name()) + ".molden"
#
#    if dovirtual is None:
#        dovirt = bool(core.get_option("SCF", "MOLDEN_WITH_VIRTUAL"))
#
#    else:
#        dovirt = dovirtual
#
#    if density_a:
#        nmopi = wfn.nmopi()
#        nsopi = wfn.nsopi()
#
#        NO_Ra = core.Matrix("NO Alpha Rotation Matrix", nmopi, nmopi)
#        NO_occa = core.Vector(nmopi)
#        density_a.diagonalize(NO_Ra, NO_occa, core.DiagonalizeOrder.Descending)
#        NO_Ca = core.Matrix("Ca Natural Orbitals", nsopi, nmopi)
#        NO_Ca.gemm(False, False, 1.0, wfn.Ca(), NO_Ra, 0)
#
#        if density_b:
#            NO_Rb = core.Matrix("NO Beta Rotation Matrix", nmopi, nmopi)
#            NO_occb = core.Vector(nmopi)
#            density_b.diagonalize(NO_Rb, NO_occb, core.DiagonalizeOrder.Descending)
#            NO_Cb = core.Matrix("Cb Natural Orbitals", nsopi, nmopi)
#            NO_Cb.gemm(False, False, 1.0, wfn.Cb(), NO_Rb, 0)
#
#        else:
#            NO_occb = NO_occa
#            NO_Cb = NO_Ca
#
#        mw = core.MoldenWriter(wfn)
#        mw.write(filename, NO_Ca, NO_Cb, NO_occa, NO_occb, NO_occa, NO_occb, dovirt)
#
#    else:
#        try:
#            occa = wfn.occupation_a()
#            occb = wfn.occupation_b()
#        except AttributeError:
#            core.print_out("\n!Molden warning: This wavefunction does not have occupation numbers.\n"
#                           "Writing zero's for occupation numbers\n\n")
#            occa = core.Vector(wfn.nmopi())
#            occb = core.Vector(wfn.nmopi())
#
#        mw = core.MoldenWriter(wfn)
#        mw.write(filename, wfn.Ca(), wfn.Cb(), wfn.epsilon_a(), wfn.epsilon_b(), occa, occb, dovirt)
#
#
## Aliases
#opt = optimize
#freq = frequency
#frequencies = frequency
#prop = properties
