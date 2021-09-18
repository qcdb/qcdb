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
from .driver_nbody import nbody_gufunc
from .proc_table import procedures

pp = pprint.PrettyPrinter(width=120)


# def energy(name, molecule, **kwargs):
#    r"""Function to compute the single-point electronic energy."""

#       :returns: *float* |w--w| Total electronic energy in Hartrees. SAPT & EFP return interaction energy.
#
#       :returns: (*float*, :py:class:`~psi4.core.Wavefunction`) |w--w| energy and wavefunction when **return_wfn** specified.
#
#       :PSI variables:
#
#       .. hlist::
#          :columns: 1
#
#          * :qcvar:`CURRENT ENERGY <CURRENTENERGY>`
#          * :qcvar:`CURRENT REFERENCE ENERGY <CURRENTREFERENCEENERGY>`
#          * :qcvar:`CURRENT CORRELATION ENERGY <CURRENTCORRELATIONENERGY>`
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
#           calculation result as the second element (after *float* energy) of a tuple.
#
#       :type restart_file: string
#       :param restart_file: ``['file.1, file.32]`` || ``./file`` || etc.
#
#           Binary data files to be renamed for calculation restart.
#
#       .. _`table:energy_gen`:
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
#       | dct                     | density cumulant functional theory :ref:`[manual] <sec:dct>`                                                  |
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
#       | ccenergy                | **expert** full control over ccenergy module                                                                  |
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
#       .. include:: ../autodoc_dft_energy.rst
#
#       .. include:: ../mrcc_table_energy.rst
#
#       .. include:: ../cfour_table_energy.rst
#
#       :examples:
#
#       >>> # [1] Coupled-cluster singles and doubles calculation with psi code
#       >>> energy('ccsd')
#
#       >>> # [2] Charge-transfer SAPT calculation with scf projection from small into
#       >>> #     requested basis, with specified projection fitting basis
#       >>> set basis_guess true
#       >>> set df_basis_guess jun-cc-pVDZ-JKFIT
#       >>> energy('sapt0-ct')
#
#       >>> # [3] Arbitrary-order MPn calculation
#       >>> energy('mp7')
#
#       >>> # [4] Converge scf as singlet, then run detci as triplet upon singlet reference
#       >>> # Note that the integral transformation is not done automatically when detci is run in a separate step.
#       >>> molecule H2 {\n0 1\nH\nH 1 0.74\n}
#       >>> set basis cc-pVDZ
#       >>> set reference rohf
#       >>> scf_e, scf_wfn = energy('scf', return_wfn=True)
#       >>> H2.set_multiplicity(3)
#       >>> core.MintsHelper(scf_wfn.basisset()).integrals()
#       >>> energy('detci', ref_wfn=scf_wfn)
#
#       >>> # [5] Run two CI calculations, keeping the integrals generated in the first one.
#       >>> molecule ne {\nNe\n}
#       >>> set basis cc-pVDZ
#       >>> cisd_e, cisd_wfn = energy('cisd', return_wfn=True)
#       >>> energy('fci', ref_wfn=cisd_wfn)
#
#       >>> # [6] Can automatically perform complete basis set extrapolations
#       >>> energy("CCSD/cc-pV[DT]Z")
#
#       >>> # [7] Can automatically perform delta corrections that include extrapolations
#       >>> # even with a user-defined extrapolation formula. See sample inputs named
#       >>> # cbs-xtpl* for more examples of this input style
#       >>> energy("MP2/aug-cc-pv([d,t]+d)z + d:ccsd(t)/cc-pvdz", corl_scheme=myxtplfn_2)
#
#       """


@register_kwds(pe.nu_options)
def energy(name, **kwargs):
    r"""Function to compute the single-point electronic energy."""

    from . import load_proc_table

    kwargs = driver_util.kwargs_lower(kwargs)

    if "options" in kwargs:
        driver_helpers.set_options(kwargs.pop("options"))

    # Bounce if name is function
    if hasattr(name, "__call__"):
        return name(energy, kwargs.pop("label", "custom function"), ptype="energy", **kwargs)

    # Allow specification of methods to arbitrary order
    lowername = name.lower()
    lowername, level = driver_helpers._parse_arbitrary_order(lowername)
    if level:
        kwargs["level"] = level

    # Make sure the molecule the user provided is the active one
    molecule = kwargs.pop("molecule", driver_helpers.get_active_molecule())
    molecule.update_geometry()

    if len(pe.nu_options.scroll) == 0:
        # print('EMPTY OPT')
        pe.load_options()

    # Bounce to CP if bsse kwarg
    if kwargs.get("bsse_type", None) is not None:
        # return nbody_driver.nbody_gufunc(energy, name, ptype="energy", molecule=molecule, **kwargs)
        return nbody_gufunc(energy, name, ptype="energy", molecule=molecule, **kwargs)

    # Bounce to CBS if "method/basis" name
    if "/" in lowername:
        return cbs_driver._cbs_gufunc(energy, name, ptype="energy", molecule=molecule, **kwargs)

    # Commit to procedures['energy'] call hereafter
    return_wfn = kwargs.pop("return_wfn", False)
    pe.active_qcvars = {}

    #    #for precallback in hooks['energy']['pre']:
    #    #    precallback(lowername, **kwargs)
    #
    #    optstash = driver_util._set_convergence_criterion('energy', lowername, 6, 8, 6, 8, 6)
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

    # PR    print('QWER', pe.nu_options.print_changed())
    package = driver_util.get_package(lowername, kwargs)
    # for k, v in pkgprefix.items():
    #    if lowername.startswith(k):
    #        package = v
    #        break
    # else:
    #    package = kwargs.get('package', 'psi4')
    # print('\nENE calling', 'procedures', package, lowername, 'with', lowername, molecule, pe.nu_options, kwargs)
    # jobrec = procedures['energy'][package][lowername](lowername, molecule=molecule, options=pe.active_options, **kwargs)
    jobrec = procedures["energy"][package][lowername](
        lowername, molecule=molecule, options=pe.nu_options, ptype="energy", **kwargs
    )

    #    for postcallback in hooks['energy']['post']:
    #        postcallback(lowername, wfn=wfn, **kwargs)
    #
    #    optstash.restore()
    # jobrec.pop('raw_output')  # just to moderate printint to screen
    # PR    pp.pprint(jobrec)
    pe.active_qcvars = copy.deepcopy(jobrec["qcvars"])

    if return_wfn:  # TODO current energy safer than wfn.energy() for now, but should be revisited

        #        # TODO place this with the associated call, very awkward to call this in other areas at the moment
        #        if lowername in ['efp', 'mrcc', 'dmrg', 'psimrcc']:
        #            core.print_out("\n\nWarning! %s does not have an associated derived wavefunction." % name)
        #            core.print_out("The returned wavefunction is the incoming reference wavefunction.\n\n")
        #        elif 'sapt' in lowername:
        #            core.print_out("\n\nWarning! %s does not have an associated derived wavefunction." % name)
        #            core.print_out("The returned wavefunction is the dimer SCF wavefunction.\n\n")

        return (float(jobrec["qcvars"]["CURRENT ENERGY"].data), jobrec)
    else:
        return float(jobrec["qcvars"]["CURRENT ENERGY"].data)
