import copy
import pprint
import re
import sys

import numpy as np
from qcelemental import Datum

from ..exceptions import ValidationError
from ..keywords import register_kwds
from ..qcvars import VARH
from ..util import banner
from . import driver_helpers, driver_util, pe
from .cbs_helpers import *

pp = pprint.PrettyPrinter(width=120)


def _cbs_wrapper_methods(**kwargs):
    cbs_method_kwargs = ["scf_wfn", "corl_wfn", "delta_wfn"]
    cbs_method_kwargs += ["delta%d_wfn" % x for x in range(2, 6)]

    cbs_methods = []
    for method in cbs_method_kwargs:
        if method in kwargs:
            cbs_methods.append(kwargs[method])
    return cbs_methods


def _parse_cbs_gufunc_string(method_name):
    """Process a 'mtd/bas + D:bas2'-like `method_name` into list of methods and bases."""

    method_name_list = re.split(r"""\+(?![^\[\]]*\]|[^\(\)]*\))""", method_name)
    if len(method_name_list) > 3:
        raise ValidationError(
            "CBS gufunc: Text parsing is only valid for a single dor double delta, please use the CBS wrapper directly"
        )

    method_list = []
    basis_list = []
    for num, method_str in enumerate(method_name_list):
        if (method_str.count("[") > 1) or (method_str.count("]") > 1):
            raise ValidationError("""CBS gufunc: Too many brackets given! %s """ % method_str)

        if method_str.count("/") != 1:
            raise ValidationError("""CBS gufunc: All methods must specify a basis with '/'. %s""" % method_str)

        if num > 0:
            method_str = method_str.strip()
            if method_str[:2].lower() != "d:":
                raise ValidationError("""CBS gufunc: Delta method must start with 'D:'.""")
            else:
                method_str = method_str[2:]
        method, basis = method_str.split("/")
        method_list.append(method)
        basis_list.append(basis)
    return method_list, basis_list


def _cbs_gufunc(func, total_method_name, molecule, **kwargs):
    """Text-based wrapper of the CBS function."""

    # Catch kwarg issues
    # print('\nINTO _cbs_gufunc', 'KW', kwargs)
    kwargs = driver_util.kwargs_lower(kwargs)
    return_wfn = kwargs.pop("return_wfn", False)
    #    core.clean_variables()
    kwargs.pop("dertype", None)
    cbs_verbose = kwargs.pop("cbs_verbose", False)
    ptype = kwargs.pop("ptype", None)

    #    # Make sure the molecule the user provided is the active one
    #    molecule = kwargs.pop('molecule') #, core.get_active_molecule())
    molecule.update_geometry()

    # Sanitize total_method_name
    label = total_method_name
    total_method_name = total_method_name.lower().replace(" ", "")

    # Split into components
    method_list, basis_list = _parse_cbs_gufunc_string(total_method_name)

    # Single energy call?
    single_call = len(method_list) == 1
    single_call &= "[" not in basis_list[0]
    single_call &= "]" not in basis_list[0]

    if single_call:
        method_name = method_list[0]
        basis = basis_list[0]

        # print('gufunc pre opt')
        # print(pe.nu_options.print_changed())
        # print('gufunc pre opt', pe.active_options)
        # Save some global variables so we can reset them later
        # optstash = options.OptionsState(['BASIS'])
        # core.set_global_option('BASIS', basis)
        # pe.active_options['GLOBALS']['BASIS']['value'] = basis
        # pe.nu_options.scroll['QCDB']['BASIS'].value = basis
        pe.nu_options.require("QCDB", "BASIS", basis, accession=1234)

        # print('\n gufunc_calling', func, method_name, return_wfn, 'OPT', pe.active_options, 'KW', kwargs)
        # print('\n gufunc_calling', func, method_name, return_wfn, 'OPT', 'KW', kwargs)
        ptype_value, wfn = func(method_name, return_wfn=True, molecule=molecule, **kwargs)
        #        core.clean()

        #        optstash.restore()

        if return_wfn:
            return (ptype_value, wfn)
        else:
            return ptype_value

    # If we are not a single call, let CBS wrapper handle it!
    cbs_kwargs = {}
    cbs_kwargs["ptype"] = ptype
    cbs_kwargs["return_wfn"] = True
    cbs_kwargs["molecule"] = molecule
    cbs_kwargs["verbose"] = cbs_verbose

    # Find method and basis
    pkgmtd = method_list[0].split("-", 1)
    if method_list[0] in ["scf", "hf"] or (
        len(pkgmtd) > 1 and (pkgmtd[1] in ["scf", "hf"] and (pkgmtd[0] + "-" in driver_util.pkgprefix))
    ):
        # method_list[0][3:] in ['scf', 'hf'] and method_list[0][:3] in pkgprefix):
        # if method_list[0] in ['scf', 'hf', 'c4-scf', 'c4-hf', 'p4-scf', 'p4-hf']:
        cbs_kwargs["scf_wfn"] = method_list[0]
        cbs_kwargs["scf_basis"] = basis_list[0]
        if "scf_scheme" in kwargs:
            cbs_kwargs["scf_scheme"] = kwargs["scf_scheme"]
    else:
        cbs_kwargs["corl_wfn"] = method_list[0]
        cbs_kwargs["corl_basis"] = basis_list[0]
        if "corl_scheme" in kwargs:
            cbs_kwargs["corl_scheme"] = kwargs["corl_scheme"]

    if len(method_list) > 2:
        cbs_kwargs["delta2_wfn"] = method_list[2]
        cbs_kwargs["delta2_basis"] = basis_list[2]
        if "delta2_scheme" in kwargs:
            cbs_kwargs["delta2_scheme"] = kwargs["delta2_scheme"]
    if len(method_list) > 1:
        cbs_kwargs["delta_wfn"] = method_list[1]
        cbs_kwargs["delta_basis"] = basis_list[1]
        if "delta_scheme" in kwargs:
            cbs_kwargs["delta_scheme"] = kwargs["delta_scheme"]

    ptype_value, wfn = cbs(func, label, **cbs_kwargs)

    if return_wfn:
        return (ptype_value, wfn)
    else:
        return ptype_value


###################################
##  Start of Complete Basis Set  ##
###################################


@register_kwds(pe.nu_options)
def cbs(func, label, **kwargs):
    r"""Function to define a multistage energy method from combinations of
    basis set extrapolations and delta corrections and condense the
    components into a minimum number of calculations.

    :aliases: complete_basis_set()

    :returns: (*float*) -- Total electronic energy in Hartrees

    :PSI variables:

    .. hlist::
       :columns: 1

       * :qcvar:`CBS TOTAL ENERGY`
       * :qcvar:`CBS REFERENCE ENERGY`
       * :qcvar:`CBS CORRELATION ENERGY`
       * :qcvar:`CURRENT ENERGY`
       * :qcvar:`CURRENT REFERENCE ENERGY`
       * :qcvar:`CURRENT CORRELATION ENERGY`

    .. caution:: Some features are not yet implemented. Buy a developer a coffee.

       - No way to tell function to boost fitting basis size for all calculations.

       - No way to extrapolate def2 family basis sets

       - Need to add more extrapolation schemes

    As represented in the equation below, a CBS energy method is defined in several
    sequential stages (scf, corl, delta, delta2, delta3, delta4, delta5) covering treatment
    of the reference total energy, the correlation energy, a delta correction to the
    correlation energy, and a second delta correction, etc.. Each is activated by its
    stage_wfn keyword and is only allowed if all preceding stages are active.

    * Energy Methods
        The presence of a stage_wfn keyword is the indicator to incorporate
        (and check for stage_basis and stage_scheme keywords) and compute
        that stage in defining the CBS energy.

        The cbs() function requires, at a minimum, ``name='scf'`` and ``scf_basis``
        keywords to be specified for reference-step only jobs and ``name`` and
        ``corl_basis`` keywords for correlated jobs.

        The following energy methods have been set up for cbs().

        .. hlist::
           :columns: 5

           * scf
           * hf
           * mp2
           * mp2.5
           * mp3
           * mp4(sdq)
           * mp4
           * mp\ *n*
           * omp2
           * omp2.5
           * omp3
           * olccd
           * lccd
           * lccsd
           * cepa(0)
           * cepa(1)
           * cepa(3)
           * acpf
           * aqcc
           * qcisd
           * cc2
           * ccsd
           * fno-ccsd
           * bccd
           * cc3
           * qcisd(t)
           * ccsd(t)
           * fno-ccsd(t)
           * bccd(t)
           * cisd
           * cisdt
           * cisdtq
           * ci\ *n*
           * fci
           * mrccsd
           * mrccsd(t)
           * mrccsdt
           * mrccsdt(q)

    :type name: str
    :param name: ``'scf'`` || ``'ccsd'`` || etc.

        First argument, usually unlabeled. Indicates the computational method
        for the correlation energy, unless only reference step to be performed,
        in which case should be ``'scf'``. Overruled if stage_wfn keywords supplied.

    :type scf_wfn: str
    :param scf_wfn: |dl| ``'scf'`` |dr| || ``'c4-scf'`` || etc.

        Indicates the energy method for which the reference energy is to be
        obtained. Generally unnecessary, as 'scf' is *the* scf in |PSIfour| but
        can be used to direct lone scf components to run in |PSIfour| or Cfour
        in a mixed-program composite method.

    :type corl_wfn: str
    :param corl_wfn: ``'mp2'`` || ``'ccsd(t)'`` || etc.

        Indicates the energy method for which the correlation energy is to be
        obtained. Can also be specified with ``name`` or as the unlabeled
        first argument to the function.

    :type delta_wfn: str
    :param delta_wfn: ``'ccsd'`` || ``'ccsd(t)'`` || etc.

        Indicates the (superior) energy method for which a delta correction
        to the correlation energy is to be obtained.

    :type delta_wfn_lesser: str
    :param delta_wfn_lesser: |dl| ``corl_wfn`` |dr| || ``'mp2'`` || etc.

        Indicates the inferior energy method for which a delta correction
        to the correlation energy is to be obtained.

    :type delta2_wfn: str
    :param delta2_wfn: ``'ccsd'`` || ``'ccsd(t)'`` || etc.

        Indicates the (superior) energy method for which a second delta correction
        to the correlation energy is to be obtained.

    :type delta2_wfn_lesser: str
    :param delta2_wfn_lesser: |dl| ``delta_wfn`` |dr| || ``'ccsd(t)'`` || etc.

        Indicates the inferior energy method for which a second delta correction
        to the correlation energy is to be obtained.

    :type delta3_wfn: str
    :param delta3_wfn: ``'ccsd'`` || ``'ccsd(t)'`` || etc.

        Indicates the (superior) energy method for which a third delta correction
        to the correlation energy is to be obtained.

    :type delta3_wfn_lesser: str
    :param delta3_wfn_lesser: |dl| ``delta2_wfn`` |dr| || ``'ccsd(t)'`` || etc.

        Indicates the inferior energy method for which a third delta correction
        to the correlation energy is to be obtained.

    :type delta4_wfn: str
    :param delta4_wfn: ``'ccsd'`` || ``'ccsd(t)'`` || etc.

        Indicates the (superior) energy method for which a fourth delta correction
        to the correlation energy is to be obtained.

    :type delta4_wfn_lesser: str
    :param delta4_wfn_lesser: |dl| ``delta3_wfn`` |dr| || ``'ccsd(t)'`` || etc.

        Indicates the inferior energy method for which a fourth delta correction
        to the correlation energy is to be obtained.

    :type delta5_wfn: str
    :param delta5_wfn: ``'ccsd'`` || ``'ccsd(t)'`` || etc.

        Indicates the (superior) energy method for which a fifth delta correction
        to the correlation energy is to be obtained.

    :type delta5_wfn_lesser: str
    :param delta5_wfn_lesser: |dl| ``delta4_wfn`` |dr| || ``'ccsd(t)'`` || etc.

        Indicates the inferior energy method for which a fifth delta correction
        to the correlation energy is to be obtained.

    * Basis Sets
        Currently, the basis set set through ``set`` commands have no influence
        on a cbs calculation.

    :type scf_basis: :ref:`basis string <apdx:basisElement>`
    :param scf_basis: |dl| ``corl_basis`` |dr| || ``'cc-pV[TQ]Z'`` || ``'jun-cc-pv[tq5]z'`` || ``'6-31G*'`` || etc.

        Indicates the sequence of basis sets employed for the reference energy.
        If any correlation method is specified, ``scf_basis`` can default
        to ``corl_basis``.

    :type corl_basis: :ref:`basis string <apdx:basisElement>`
    :param corl_basis: ``'cc-pV[TQ]Z'`` || ``'jun-cc-pv[tq5]z'`` || ``'6-31G*'`` || etc.

        Indicates the sequence of basis sets employed for the correlation energy.

    :type delta_basis: :ref:`basis string <apdx:basisElement>`
    :param delta_basis: ``'cc-pV[TQ]Z'`` || ``'jun-cc-pv[tq5]z'`` || ``'6-31G*'`` || etc.

        Indicates the sequence of basis sets employed for the delta correction
        to the correlation energy.

    :type delta2_basis: :ref:`basis string <apdx:basisElement>`
    :param delta2_basis: ``'cc-pV[TQ]Z'`` || ``'jun-cc-pv[tq5]z'`` || ``'6-31G*'`` || etc.

        Indicates the sequence of basis sets employed for the second delta correction
        to the correlation energy.

    :type delta3_basis: :ref:`basis string <apdx:basisElement>`
    :param delta3_basis: ``'cc-pV[TQ]Z'`` || ``'jun-cc-pv[tq5]z'`` || ``'6-31G*'`` || etc.

        Indicates the sequence of basis sets employed for the third delta correction
        to the correlation energy.

    :type delta4_basis: :ref:`basis string <apdx:basisElement>`
    :param delta4_basis: ``'cc-pV[TQ]Z'`` || ``'jun-cc-pv[tq5]z'`` || ``'6-31G*'`` || etc.

        Indicates the sequence of basis sets employed for the fourth delta correction
        to the correlation energy.

    :type delta5_basis: :ref:`basis string <apdx:basisElement>`
    :param delta5_basis: ``'cc-pV[TQ]Z'`` || ``'jun-cc-pv[tq5]z'`` || ``'6-31G*'`` || etc.

        Indicates the sequence of basis sets employed for the fifth delta correction
        to the correlation energy.

    * Schemes
        Transformations of the energy through basis set extrapolation for each
        stage of the CBS definition. A complaint is generated if number of basis
        sets in stage_basis does not exactly satisfy requirements of stage_scheme.
        An exception is the default, ``'xtpl_highest_1'``, which uses the best basis
        set available. See :ref:`sec:cbs_xtpl` for all available schemes.

    :type scf_scheme: Callable
    :param scf_scheme: |dl| ``xtpl_highest_1`` |dr| || ``scf_xtpl_helgaker_3`` || etc.

        Indicates the basis set extrapolation scheme to be applied to the reference energy.
        Defaults to :py:func:`~scf_xtpl_helgaker_3` if three valid basis sets
        present in ``scf_basis``, :py:func:`~scf_xtpl_helgaker_2` if two valid basis
        sets present in ``scf_basis``, and :py:func:`~xtpl_highest_1` otherwise.

    :type corl_scheme: Callable
    :param corl_scheme: |dl| ``xtpl_highest_1`` |dr| || ``corl_xtpl_helgaker_2`` || etc.

        Indicates the basis set extrapolation scheme to be applied to the correlation energy.
        Defaults to :py:func:`~corl_xtpl_helgaker_2` if two valid basis sets
        present in ``corl_basis`` and :py:func:`~xtpl_highest_1` otherwise.

    :type delta_scheme: Callable
    :param delta_scheme: |dl| ``xtpl_highest_1`` |dr| || ``corl_xtpl_helgaker_2`` || etc.

        Indicates the basis set extrapolation scheme to be applied to the delta correction
        to the correlation energy.
        Defaults to :py:func:`~corl_xtpl_helgaker_2` if two valid basis sets
        present in ``delta_basis`` and :py:func:`~xtpl_highest_1` otherwise.

    :type delta2_scheme: Callable
    :param delta2_scheme: |dl| ``xtpl_highest_1`` |dr| || ``corl_xtpl_helgaker_2`` || etc.

        Indicates the basis set extrapolation scheme to be applied to the second delta correction
        to the correlation energy.
        Defaults to :py:func:`~corl_xtpl_helgaker_2` if two valid basis sets
        present in ``delta2_basis`` and :py:func:`~xtpl_highest_1` otherwise.

    :type delta3_scheme: Callable
    :param delta3_scheme: |dl| ``xtpl_highest_1`` |dr| || ``corl_xtpl_helgaker_2`` || etc.

        Indicates the basis set extrapolation scheme to be applied to the third delta correction
        to the correlation energy.
        Defaults to :py:func:`~corl_xtpl_helgaker_2` if two valid basis sets
        present in ``delta3_basis`` and :py:func:`~xtpl_highest_1` otherwise.

    :type delta4_scheme: Callable
    :param delta4_scheme: |dl| ``xtpl_highest_1`` |dr| || ``corl_xtpl_helgaker_2`` || etc.

        Indicates the basis set extrapolation scheme to be applied to the fourth delta correction
        to the correlation energy.
        Defaults to :py:func:`~corl_xtpl_helgaker_2` if two valid basis sets
        present in ``delta4_basis`` and :py:func:`~xtpl_highest_1` otherwise.

    :type delta5_scheme: Callable
    :param delta5_scheme: |dl| ``xtpl_highest_1`` |dr| || ``corl_xtpl_helgaker_2`` || etc.

        Indicates the basis set extrapolation scheme to be applied to the fifth delta correction
        to the correlation energy.
        Defaults to :py:func:`~corl_xtpl_helgaker_2` if two valid basis sets
        present in ``delta5_basis`` and :py:func:`~xtpl_highest_1` otherwise.

    :type molecule: :ref:`molecule <op_py_molecule>`
    :param molecule: ``h2o`` || etc.

        The target molecule, if not the last molecule defined.

    :examples:


    >>> # [1] replicates with cbs() the simple model chemistry scf/cc-pVDZ: set basis cc-pVDZ energy('scf')
    >>> cbs(name='scf', scf_basis='cc-pVDZ')

    >>> # [2] replicates with cbs() the simple model chemistry mp2/jun-cc-pVDZ: set basis jun-cc-pVDZ energy('mp2')
    >>> cbs(name='mp2', corl_basis='jun-cc-pVDZ')

    >>> # [3] DTQ-zeta extrapolated scf reference energy
    >>> cbs(name='scf', scf_basis='cc-pV[DTQ]Z', scf_scheme=scf_xtpl_helgaker_3)

    >>> # [4] DT-zeta extrapolated mp2 correlation energy atop a T-zeta reference
    >>> cbs(corl_wfn='mp2', corl_basis='cc-pv[dt]z', corl_scheme=corl_xtpl_helgaker_2)

    >>> # [5] a DT-zeta extrapolated coupled-cluster correction atop a TQ-zeta extrapolated mp2 correlation energy atop a Q-zeta reference (both equivalent)
    >>> cbs(corl_wfn='mp2', corl_basis='aug-cc-pv[tq]z', delta_wfn='ccsd(t)', delta_basis='aug-cc-pv[dt]z')
    >>> cbs(energy, wfn='mp2', corl_basis='aug-cc-pv[tq]z', corl_scheme=corl_xtpl_helgaker_2, delta_wfn='ccsd(t)', delta_basis='aug-cc-pv[dt]z', delta_scheme=corl_xtpl_helgaker_2)

    >>> # [6] a D-zeta ccsd(t) correction atop a DT-zeta extrapolated ccsd cluster correction atop a TQ-zeta extrapolated mp2 correlation energy atop a Q-zeta reference
    >>> cbs(name='mp2', corl_basis='aug-cc-pv[tq]z', corl_scheme=corl_xtpl_helgaker_2, delta_wfn='ccsd', delta_basis='aug-cc-pv[dt]z', delta_scheme=corl_xtpl_helgaker_2, delta2_wfn='ccsd(t)', delta2_wfn_lesser='ccsd', delta2_basis='aug-cc-pvdz')

    >>> # [7] cbs() coupled with database()
    >>> TODO database('mp2', 'BASIC', subset=['h2o','nh3'], symm='on', func=cbs, corl_basis='cc-pV[tq]z', corl_scheme=corl_xtpl_helgaker_2, delta_wfn='ccsd(t)', delta_basis='sto-3g')

    >>> # [8] cbs() coupled with optimize()
    >>> TODO optimize('mp2', corl_basis='cc-pV[DT]Z', corl_scheme=corl_xtpl_helgaker_2, func=cbs)

    """
    kwargs = driver_util.kwargs_lower(kwargs)
    return_wfn = kwargs.pop("return_wfn", False)
    verbose = kwargs.pop("verbose", 0)
    ptype = kwargs.pop("ptype")

    kwgs = {"accession": kwargs["accession"], "verbose": verbose}

    # Establish function to call (only energy makes sense for cbs)
    if ptype not in ["energy", "gradient", "hessian"]:
        raise ValidationError(
            """Wrapper complete_basis_set is unhappy to be calling function '%s' instead of 'energy'.""" % ptype
        )

    #    optstash = p4util.OptionsState(
    #        ['BASIS'],
    #        ['WFN'],
    #        ['WRITER_FILE_LABEL'])

    # Define some quantum chemical knowledge, namely what methods are subsumed in others

    do_scf = True
    do_corl = False
    do_delta = False
    do_delta2 = False
    do_delta3 = False
    do_delta4 = False
    do_delta5 = False

    user_writer_file_label = pe.nu_options.scroll["QCDB"][
        "WRITER_FILE_LABEL"
    ].value  # core.get_global_option('WRITER_FILE_LABEL')

    # Make sure the molecule the user provided is the active one
    molecule = kwargs.pop("molecule", driver_helpers.get_active_molecule())
    molecule.update_geometry()
    natom = molecule.natom()

    #    molstr = molecule.create_psi4_string_from_molecule()

    # Establish method for correlation energy
    cbs_corl_wfn = kwargs.pop("corl_wfn", "").lower()
    if cbs_corl_wfn:
        do_corl = True

    # Establish method for reference energy
    if do_corl and cbs_corl_wfn.startswith("c4-"):
        default_scf = "c4-hf"
    elif do_corl and cbs_corl_wfn.startswith("p4-"):
        default_scf = "p4-hf"
    elif do_corl and cbs_corl_wfn.startswith("nwc-"):
        default_scf = "nwc-hf"
    elif do_corl and cbs_corl_wfn.startswith("gms-"):
        default_scf = "gms-hf"
    else:
        default_scf = "hf"
    cbs_scf_wfn = kwargs.pop("scf_wfn", default_scf).lower()

    if do_scf:
        if cbs_scf_wfn not in VARH:
            raise ValidationError(
                """Requested SCF method '%s' is not recognized. Add it to VARH in wrapper.py to proceed."""
                % (cbs_scf_wfn)
            )

    # ... resume correlation logic
    if do_corl:
        if cbs_corl_wfn not in VARH:
            raise ValidationError(
                """Requested CORL method '%s' is not recognized. Add it to VARH in wrapper.py to proceed."""
                % (cbs_corl_wfn)
            )

        cbs_corl_wfn_lesser = kwargs.get("corl_wfn_lesser", cbs_scf_wfn).lower()
        if cbs_corl_wfn_lesser not in VARH:
            raise ValidationError(
                """Requested CORL method lesser '%s' is not recognized. Add it to VARH in wrapper.py to proceed."""
                % (cbs_delta_wfn_lesser)
            )

    # Establish method for delta correction energy
    if "delta_wfn" in kwargs:
        do_delta = True
        cbs_delta_wfn = kwargs["delta_wfn"].lower()
        if cbs_delta_wfn not in VARH:
            raise ValidationError(
                """Requested DELTA method '%s' is not recognized. Add it to VARH in wrapper.py to proceed."""
                % (cbs_delta_wfn)
            )

        # if 'delta_wfn_lesser' in kwargs:
        #    cbs_delta_wfn_lesser = kwargs['delta_wfn_lesser'].lower()
        #    if cbs_delta_wfn.split('-', 1)[0] != cbs_delta_wfn_lesser.split('-', 1)[0]:
        #        raise ValidationError("""Bad idea to mix programs.""")
        # else:
        #    prog = cbs_delta_wfn.split('-', 1)[0]
        #    mtd = cbs_corl_wfn.split('-', 1)[1]
        #    cbs_delta_wfn_lesser = prog + '-' + mtd

        cbs_delta_wfn_lesser = kwargs.get("delta_wfn_lesser", cbs_corl_wfn).lower()
        if cbs_delta_wfn_lesser not in VARH:
            raise ValidationError(
                """Requested DELTA method lesser '%s' is not recognized. Add it to VARH in wrapper.py to proceed."""
                % (cbs_delta_wfn_lesser)
            )

    # Establish method for second delta correction energy
    if "delta2_wfn" in kwargs:
        do_delta2 = True
        cbs_delta2_wfn = kwargs["delta2_wfn"].lower()
        if cbs_delta2_wfn not in VARH:
            raise ValidationError(
                """Requested DELTA2 method '%s' is not recognized. Add it to VARH in wrapper.py to proceed."""
                % (cbs_delta2_wfn)
            )

        # if 'delta2_wfn_lesser' in kwargs:
        #    cbs_delta2_wfn_lesser = kwargs['delta2_wfn_lesser'].lower()
        #    if cbs_delta2_wfn.split('-', 1)[0] != cbs_delta2_wfn_lesser.split('-', 1)[0]:
        #        raise ValidationError("""Bad idea to mix programs.""")
        # else:
        #    prog = cbs_delta2_wfn.split('-', 1)[0]
        #    mtd = cbs_delta_wfn.split('-', 1)[1]
        #    cbs_delta2_wfn_lesser = prog + '-' + mtd

        cbs_delta2_wfn_lesser = kwargs.get("delta2_wfn_lesser", cbs_delta_wfn).lower()
        if cbs_delta2_wfn_lesser not in VARH:
            raise ValidationError(
                """Requested DELTA2 method lesser '%s' is not recognized. Add it to VARH in wrapper.py to proceed."""
                % (cbs_delta2_wfn_lesser)
            )

    #    # Establish method for third delta correction energy
    #    if 'delta3_wfn' in kwargs:
    #        do_delta3 = True
    #        cbs_delta3_wfn = kwargs['delta3_wfn'].lower()
    #        if cbs_delta3_wfn not in VARH.keys():
    #            raise ValidationError("""Requested DELTA3 method '%s' is not recognized. Add it to VARH in wrapper.py to proceed.""" % (cbs_delta3_wfn))
    #
    #        cbs_delta3_wfn_lesser = kwargs.get('delta3_wfn_lesser', cbs_delta2_wfn).lower()
    #        if not (cbs_delta3_wfn_lesser in VARH.keys()):
    #            raise ValidationError("""Requested DELTA3 method lesser '%s' is not recognized. Add it to VARH in wrapper.py to proceed.""" % (cbs_delta3_wfn_lesser))
    #
    #    # Establish method for fourth delta correction energy
    #    if 'delta4_wfn' in kwargs:
    #        do_delta4 = True
    #        cbs_delta4_wfn = kwargs['delta4_wfn'].lower()
    #        if not (cbs_delta4_wfn in VARH.keys()):
    #            raise ValidationError('Requested DELTA4 method \'%s\' is not recognized. Add it to VARH in wrapper.py to proceed.' % (cbs_delta4_wfn))
    #
    #        if 'delta4_wfn_lesser' in kwargs:
    #            cbs_delta4_wfn_lesser = kwargs['delta4_wfn_lesser'].lower()
    #        else:
    #            cbs_delta4_wfn_lesser = cbs_delta3_wfn
    #        if not (cbs_delta4_wfn_lesser in VARH.keys()):
    #            raise ValidationError('Requested DELTA4 method lesser \'%s\' is not recognized. Add it to VARH in wrapper.py to proceed.' % (cbs_delta4_wfn_lesser))
    #
    #    # Establish method for fifth delta correction energy
    #    if 'delta5_wfn' in kwargs:
    #        do_delta5 = True
    #        cbs_delta5_wfn = kwargs['delta5_wfn'].lower()
    #        if not (cbs_delta5_wfn in VARH.keys()):
    #            raise ValidationError('Requested DELTA5 method \'%s\' is not recognized. Add it to VARH in wrapper.py to proceed.' % (cbs_delta5_wfn))
    #
    #        if 'delta5_wfn_lesser' in kwargs:
    #            cbs_delta5_wfn_lesser = kwargs['delta5_wfn_lesser'].lower()
    #        else:
    #            cbs_delta5_wfn_lesser = cbs_delta4_wfn
    #        if not (cbs_delta5_wfn_lesser in VARH.keys()):
    #            raise ValidationError('Requested DELTA5 method lesser \'%s\' is not recognized. Add it to VARH in wrapper.py to proceed.' % (cbs_delta5_wfn_lesser))

    # Check that user isn't skipping steps in scf + corl + delta + delta2 sequence
    if do_scf and not do_corl and not do_delta and not do_delta2 and not do_delta3 and not do_delta4 and not do_delta5:
        pass
    elif do_scf and do_corl and not do_delta and not do_delta2 and not do_delta3 and not do_delta4 and not do_delta5:
        pass
    elif do_scf and do_corl and do_delta and not do_delta2 and not do_delta3 and not do_delta4 and not do_delta5:
        pass
    elif do_scf and do_corl and do_delta and do_delta2 and not do_delta3 and not do_delta4 and not do_delta5:
        pass
    # elif do_scf and do_corl and do_delta and do_delta2 and do_delta3 and not do_delta4 and not do_delta5:
    #    pass
    # elif do_scf and do_corl and do_delta and do_delta2 and do_delta3 and do_delta4 and not do_delta5:
    #    pass
    # elif do_scf and do_corl and do_delta and do_delta2 and do_delta3 and do_delta4 and do_delta5:
    #    pass
    else:
        raise ValidationError(
            "Requested scf (%s) + corl (%s) + delta (%s) + delta2 (%s) + delta3 (%s) + delta4 (%s) + delta5 (%s) not valid. These steps are cummulative."
            % (do_scf, do_corl, do_delta, do_delta2, do_delta3, do_delta4, do_delta5)
        )

    # Establish list of valid basis sets for correlation energy
    if do_corl:
        if "corl_basis" in kwargs:
            BSTC, ZETC = _expand_bracketed_basis(kwargs["corl_basis"].lower(), molecule=molecule)
        else:
            raise ValidationError("""CORL basis sets through keyword '%s' are required.""" % ("corl_basis"))

    # Establish list of valid basis sets for scf energy
    if "scf_basis" in kwargs:
        BSTR, ZETR = _expand_bracketed_basis(kwargs["scf_basis"].lower(), molecule=molecule)
    elif do_corl:
        BSTR = BSTC[:]
        ZETR = ZETC[:]
    else:
        raise ValidationError(
            """SCF basis sets through keyword '%s' are required. Or perhaps you forgot the '%s'."""
            % ("scf_basis", "corl_wfn")
        )

    # Establish list of valid basis sets for delta correction energy
    if do_delta:
        if "delta_basis" in kwargs:
            BSTD, ZETD = _expand_bracketed_basis(kwargs["delta_basis"].lower(), molecule=molecule)
        else:
            raise ValidationError("""DELTA basis sets through keyword '%s' are required.""" % ("delta_basis"))

    # Establish list of valid basis sets for second delta correction energy
    if do_delta2:
        if "delta2_basis" in kwargs:
            BSTD2, ZETD2 = _expand_bracketed_basis(kwargs["delta2_basis"].lower(), molecule=molecule)
        else:
            raise ValidationError("""DELTA2 basis sets through keyword '%s' are required.""" % ("delta2_basis"))

    #    # Establish list of valid basis sets for third delta correction energy
    #    if do_delta3:
    #        if 'delta3_basis' in kwargs:
    #            BSTD3, ZETD3 = validate_bracketed_basis(kwargs['delta3_basis'].lower())
    #        else:
    #            raise ValidationError('DELTA3 basis sets through keyword \'%s\' are required.' % ('delta3_basis'))
    #
    #    # Establish list of valid basis sets for fourth delta correction energy
    #    if do_delta4:
    #        if 'delta4_basis' in kwargs:
    #            BSTD4, ZETD4 = validate_bracketed_basis(kwargs['delta4_basis'].lower())
    #        else:
    #            raise ValidationError('DELTA4 basis sets through keyword \'%s\' are required.' % ('delta4_basis'))
    #
    #    # Establish list of valid basis sets for fifth delta correction energy
    #    if do_delta5:
    #        if 'delta5_basis' in kwargs:
    #            BSTD5, ZETD5 = validate_bracketed_basis(kwargs['delta5_basis'].lower())
    #        else:
    #            raise ValidationError('DELTA5 basis sets through keyword \'%s\' are required.' % ('delta5_basis'))

    # Establish treatment for scf energy (validity check useless since python will catch it long before here)
    if (len(BSTR) == 3) and ("scf_basis" in kwargs):
        cbs_scf_scheme = scf_xtpl_helgaker_3
    elif (len(BSTR) == 2) and ("scf_basis" in kwargs):
        cbs_scf_scheme = scf_xtpl_helgaker_2
    elif (len(BSTR) == 1) and ("scf_basis" in kwargs):
        cbs_scf_scheme = xtpl_highest_1
    elif "scf_basis" in kwargs:
        raise ValidationError("""SCF basis sets of number %d cannot be handled.""" % (len(BSTR)))
    elif do_corl:
        cbs_scf_scheme = xtpl_highest_1
        BSTR = [BSTC[-1]]
        ZETR = [ZETC[-1]]
    if "scf_scheme" in kwargs:
        cbs_scf_scheme = kwargs["scf_scheme"]

    # Establish treatment for correlation energy
    if do_corl:
        if len(BSTC) == 2:
            cbs_corl_scheme = corl_xtpl_helgaker_2
        else:
            cbs_corl_scheme = xtpl_highest_1
        if "corl_scheme" in kwargs:
            cbs_corl_scheme = kwargs["corl_scheme"]

    # Establish treatment for delta correction energy
    if do_delta:
        if len(BSTD) == 2:
            cbs_delta_scheme = corl_xtpl_helgaker_2
        else:
            cbs_delta_scheme = xtpl_highest_1
        if "delta_scheme" in kwargs:
            cbs_delta_scheme = kwargs["delta_scheme"]

    # Establish treatment for delta2 correction energy
    if do_delta2:
        if len(BSTD2) == 2:
            cbs_delta2_scheme = corl_xtpl_helgaker_2
        else:
            cbs_delta2_scheme = xtpl_highest_1
        if "delta2_scheme" in kwargs:
            cbs_delta2_scheme = kwargs["delta2_scheme"]

    #    # Establish treatment for delta3 correction energy
    #    if do_delta3:
    #        if len(BSTD3) == 2:
    #            cbs_delta3_scheme = corl_xtpl_helgaker_2
    #        else:
    #            cbs_delta3_scheme = xtpl_highest_1
    #        if 'delta3_scheme' in kwargs:
    #            cbs_delta3_scheme = kwargs['delta3_scheme']
    #
    #    # Establish treatment for delta4 correction energy
    #    if do_delta4:
    #        if len(BSTD4) == 2:
    #            cbs_delta4_scheme = corl_xtpl_helgaker_2
    #        else:
    #            cbs_delta4_scheme = xtpl_highest_1
    #        if 'delta4_scheme' in kwargs:
    #            cbs_delta4_scheme = kwargs['delta4_scheme']
    #
    #    # Establish treatment for delta5 correction energy
    #    if do_delta5:
    #        if len(BSTD5) == 2:
    #            cbs_delta5_scheme = corl_xtpl_helgaker_2
    #        else:
    #            cbs_delta5_scheme = xtpl_highest_1
    #        if 'delta5_scheme' in kwargs:
    #            cbs_delta5_scheme = kwargs['delta5_scheme']

    # Build string of title banner
    cbsbanners = banner(" CBS Setup: {} ".format(label))
    #    cbsbanners = ''
    #    cbsbanners += """core.print_out('\\n')\n"""
    #    cbsbanners += """p4util.banner(' CBS Setup: %s ' % label)\n"""
    #    cbsbanners += """core.print_out('\\n')\n\n"""
    #    exec(cbsbanners)

    # Call schemes for each portion of total energy to 'place orders' for calculations needed
    d_fields = ["d_stage", "d_scheme", "d_basis", "d_wfn", "d_need", "d_coef", "d_energy", "d_gradient", "d_hessian"]
    f_fields = ["f_wfn", "f_basis", "f_zeta", "f_energy", "f_gradient", "f_hessian"]
    GRAND_NEED = []
    MODELCHEM = []
    if do_scf:
        NEED = _expand_scheme_orders(cbs_scf_scheme, BSTR, ZETR, cbs_scf_wfn, natom)
        GRAND_NEED.append(
            dict(
                zip(
                    d_fields,
                    ["scf", cbs_scf_scheme, _contract_bracketed_basis(BSTR), cbs_scf_wfn, NEED, +1, 0.0, None, None],
                )
            )
        )

    if do_corl:
        NEED = _expand_scheme_orders(cbs_corl_scheme, BSTC, ZETC, cbs_corl_wfn, natom)
        GRAND_NEED.append(
            dict(
                zip(
                    d_fields,
                    ["corl", cbs_corl_scheme, _contract_bracketed_basis(BSTC), cbs_corl_wfn, NEED, +1, 0.0, None, None],
                )
            )
        )

        NEED = _expand_scheme_orders(cbs_corl_scheme, BSTC, ZETC, cbs_corl_wfn_lesser, natom)
        GRAND_NEED.append(
            dict(
                zip(
                    d_fields,
                    [
                        "corl",
                        cbs_corl_scheme,
                        _contract_bracketed_basis(BSTC),
                        cbs_corl_wfn_lesser,
                        NEED,
                        -1,
                        0.0,
                        None,
                        None,
                    ],
                )
            )
        )

    if do_delta:
        NEED = _expand_scheme_orders(cbs_delta_scheme, BSTD, ZETD, cbs_delta_wfn, natom)
        GRAND_NEED.append(
            dict(
                zip(
                    d_fields,
                    [
                        "delta",
                        cbs_delta_scheme,
                        _contract_bracketed_basis(BSTD),
                        cbs_delta_wfn,
                        NEED,
                        +1,
                        0.0,
                        None,
                        None,
                    ],
                )
            )
        )

        NEED = _expand_scheme_orders(cbs_delta_scheme, BSTD, ZETD, cbs_delta_wfn_lesser, natom)
        GRAND_NEED.append(
            dict(
                zip(
                    d_fields,
                    [
                        "delta",
                        cbs_delta_scheme,
                        _contract_bracketed_basis(BSTD),
                        cbs_delta_wfn_lesser,
                        NEED,
                        -1,
                        0.0,
                        None,
                        None,
                    ],
                )
            )
        )

    if do_delta2:
        NEED = _expand_scheme_orders(cbs_delta2_scheme, BSTD2, ZETD2, cbs_delta2_wfn, natom)
        GRAND_NEED.append(
            dict(
                zip(
                    d_fields,
                    [
                        "delta2",
                        cbs_delta2_scheme,
                        _contract_bracketed_basis(BSTD2),
                        cbs_delta2_wfn,
                        NEED,
                        +1,
                        0.0,
                        None,
                        None,
                    ],
                )
            )
        )

        NEED = _expand_scheme_orders(cbs_delta2_scheme, BSTD2, ZETD2, cbs_delta2_wfn_lesser, natom)
        GRAND_NEED.append(
            dict(
                zip(
                    d_fields,
                    [
                        "delta2",
                        cbs_delta2_scheme,
                        _contract_bracketed_basis(BSTD2),
                        cbs_delta2_wfn_lesser,
                        NEED,
                        -1,
                        0.0,
                        None,
                        None,
                    ],
                )
            )
        )

    #    if do_delta3:
    #        NEED = call_function_in_1st_argument(cbs_delta3_scheme,
    #            mode='requisition', basisname=BSTD3, basiszeta=ZETD3, wfnname=cbs_delta3_wfn)
    #        GRAND_NEED.append(dict(zip(d_fields, ['delta3', cbs_delta3_scheme,
    #            reconstitute_bracketed_basis(NEED), cbs_delta3_wfn, NEED, +1, 0.0])))
    #
    #        NEED = call_function_in_1st_argument(cbs_delta3_scheme,
    #            mode='requisition', basisname=BSTD3, basiszeta=ZETD3, wfnname=cbs_delta3_wfn_lesser)
    #        GRAND_NEED.append(dict(zip(d_fields, ['delta3', cbs_delta3_scheme,
    #            reconstitute_bracketed_basis(NEED), cbs_delta3_wfn_lesser, NEED, -1, 0.0])))
    #
    #    if do_delta4:
    #        NEED = call_function_in_1st_argument(cbs_delta4_scheme,
    #            mode='requisition', basisname=BSTD4, basiszeta=ZETD4, wfnname=cbs_delta4_wfn)
    #        GRAND_NEED.append(dict(zip(d_fields, ['delta4', cbs_delta4_scheme,
    #            reconstitute_bracketed_basis(NEED), cbs_delta4_wfn, NEED, +1, 0.0])))
    #
    #        NEED = call_function_in_1st_argument(cbs_delta4_scheme,
    #            mode='requisition', basisname=BSTD4, basiszeta=ZETD4, wfnname=cbs_delta4_wfn_lesser)
    #        GRAND_NEED.append(dict(zip(d_fields, ['delta4', cbs_delta4_scheme,
    #            reconstitute_bracketed_basis(NEED), cbs_delta4_wfn_lesser, NEED, -1, 0.0])))
    #
    #    if do_delta5:
    #        NEED = call_function_in_1st_argument(cbs_delta5_scheme,
    #            mode='requisition', basisname=BSTD5, basiszeta=ZETD5, wfnname=cbs_delta5_wfn)
    #        GRAND_NEED.append(dict(zip(d_fields, ['delta5', cbs_delta5_scheme,
    #            reconstitute_bracketed_basis(NEED), cbs_delta5_wfn, NEED, +1, 0.0])))
    #
    #        NEED = call_function_in_1st_argument(cbs_delta5_scheme,
    #            mode='requisition', basisname=BSTD5, basiszeta=ZETD5, wfnname=cbs_delta5_wfn_lesser)
    #        GRAND_NEED.append(dict(zip(d_fields, ['delta5', cbs_delta5_scheme,
    #            reconstitute_bracketed_basis(NEED), cbs_delta5_wfn_lesser, NEED, -1, 0.0])))

    for stage in GRAND_NEED:
        for lvl in stage["d_need"].items():
            MODELCHEM.append(lvl[1])

    # Apply chemical reasoning to choose the minimum computations to run
    JOBS = MODELCHEM[:]

    addlremark = {"energy": "", "gradient": ", GRADIENT", "hessian": ", HESSIAN"}
    instructions = ""
    instructions += """    Naive listing of computations required.\n"""
    for mc in JOBS:
        instructions += """   %12s / %-24s for  %s%s\n""" % (
            mc["f_wfn"],
            mc["f_basis"],
            VARH[mc["f_wfn"]][mc["f_wfn"]],
            addlremark[ptype],
        )

    #     Remove duplicate modelchem portion listings
    for mc in MODELCHEM:
        dups = -1
        for indx_job, job in enumerate(JOBS):
            if (job["f_wfn"] == mc["f_wfn"]) and (job["f_basis"] == mc["f_basis"]):
                dups += 1
                if dups >= 1:
                    del JOBS[indx_job]

    #     Remove chemically subsumed modelchem portion listings
    if ptype == "energy":
        for mc in MODELCHEM:
            for wfn in VARH[mc["f_wfn"]]:
                for indx_job, job in enumerate(JOBS):
                    if (
                        (VARH[mc["f_wfn"]][wfn] == VARH[job["f_wfn"]][job["f_wfn"]])
                        and (mc["f_basis"] == job["f_basis"])
                        and not (mc["f_wfn"] == job["f_wfn"])
                    ):
                        del JOBS[indx_job]

    instructions += """\n    Enlightened listing of computations required.\n"""
    for mc in JOBS:
        instructions += """   %12s / %-24s for  %s%s\n""" % (
            mc["f_wfn"],
            mc["f_basis"],
            VARH[mc["f_wfn"]][mc["f_wfn"]],
            addlremark[ptype],
        )

    #     Expand listings to all that will be obtained
    JOBS_EXT = []
    for job in JOBS:
        for wfn in VARH[job["f_wfn"]]:
            JOBS_EXT.append(
                dict(
                    zip(
                        f_fields,
                        [
                            wfn,
                            job["f_basis"],
                            job["f_zeta"],
                            0.0,
                            np.zeros((natom, 3)),
                            np.zeros((3 * natom, 3 * natom)),
                        ],
                    )
                )
            )

    instructions += """\n    Full listing of computations to be obtained (required and bonus).\n"""
    for mc in JOBS_EXT:
        instructions += """   %12s / %-24s for  %s%s\n""" % (
            mc["f_wfn"],
            mc["f_basis"],
            VARH[mc["f_wfn"]][mc["f_wfn"]],
            addlremark[ptype],
        )
    print(instructions)

    #    psioh = core.IOManager.shared_object()
    #    psioh.set_specific_retention(constants.PSIF_SCF_MOS, True)
    # projection across point groups not allowed and cbs() usually a mix of symm-enabled and symm-tol calls
    #   needs to be communicated to optimize() so reset by that optstash
    #    core.set_local_option('SCF', 'GUESS_PERSIST', True)

    Njobs = 0
    # Run necessary computations
    for mc in JOBS:
        kwargs["name"] = mc["f_wfn"]

        # Build string of title banner
        cbsbanners = banner(
            " CBS Computation: {} / {}{} ".format(mc["f_wfn"].upper(), mc["f_basis"].upper(), addlremark[ptype])
        )
        print(cbsbanners)
        #        cbsbanners = ''
        #        cbsbanners += """core.print_out('\\n')\n"""
        #        cbsbanners += """p4util.banner(' CBS Computation: %s / %s%s ')\n""" % \
        #            (mc['f_wfn'].upper(), mc['f_basis'].upper(), addlremark[ptype])
        #        cbsbanners += """core.print_out('\\n')\n\n"""
        #        exec(cbsbanners)

        # Build string of molecule and commands that are dependent on the database
        pe.nu_options.require("QCDB", "BASIS", mc["f_basis"], **kwgs)
        pe.nu_options.require(
            "QCDB",
            "WRITER_FILE_LABEL",
            "-".join(
                filter(
                    None,
                    [
                        pe.nu_options.scroll["QCDB"]["WRITER_FILE_LABEL"].value,
                        mc["f_wfn"].lower(),
                        mc["f_basis"].lower(),
                    ],
                )
            )[:60],
            **kwgs,
        )

        #        commands = '\n'
        #        commands += """\ncore.set_global_option('BASIS', '%s')\n""" % (mc['f_basis'])
        #        commands += """core.set_global_option('WRITER_FILE_LABEL', '%s')\n""" % \
        #            (user_writer_file_label + ('' if user_writer_file_label == '' else '-') + mc['f_wfn'].lower() + '-' + mc['f_basis'].lower())
        #        exec(commands)

        # Make energy(), etc. call
        safe_kwargs = {k: v for k, v in kwargs.items() if k not in ["scf_scheme", "corl_scheme", "delta_scheme"]}
        response, jrec = func(molecule=molecule, return_wfn=True, **safe_kwargs)
        # response = func(molecule=molecule, **kwargs)
        if ptype == "energy":
            mc["f_energy"] = response

        elif ptype == "gradient":
            mc["f_gradient"] = response
            mc["f_energy"] = float(jrec["qcvars"]["CURRENT ENERGY"].data)
            if verbose > 1:
                print(np.array_str(mc["f_gradient"], max_line_width=120, precision=8, suppress_small=True))
                # mc['f_gradient'].print_out()

        elif ptype == "hessian":
            mc["f_hessian"] = response
            mc["f_energy"] = float(jrec["qcvars"]["CURRENT ENERGY"].data)
            mc["f_gradient"] = jrec["qcvars"]["CURRENT GRADIENT"].data  # risky if GRAD not computed
            if verbose > 1:
                print(np.array_str(mc["f_hessian"], max_line_width=120, precision=8, suppress_small=True))
        Njobs += 1
        if verbose > 1:
            print("\nCURRENT ENERGY: %14.16f\n" % mc["f_energy"])

        # Fill in energies for subsumed methods
        if ptype == "energy":
            for wfn in VARH[mc["f_wfn"]]:
                for job in JOBS_EXT:
                    if (wfn == job["f_wfn"]) and (mc["f_basis"] == job["f_basis"]):
                        # job['f_energy'] = core.get_variable(VARH[wfn][wfn])
                        job["f_energy"] = float(jrec["qcvars"][VARH[wfn][wfn]].data)

        if verbose > 1:
            jrec["qcvars"].print_variables()
        pe.active_qcvars = {}  # even want to use global qcvars?
        # core.clean_variables()
        # core.clean()

        # Copy data from 'run' to 'obtained' table
        for mce in JOBS_EXT:
            if (mc["f_wfn"] == mce["f_wfn"]) and (mc["f_basis"] == mce["f_basis"]):
                mce["f_energy"] = mc["f_energy"]
                mce["f_gradient"] = mc["f_gradient"]
                mce["f_hessian"] = mc["f_hessian"]

    #    psioh.set_specific_retention(constants.PSIF_SCF_MOS, False)

    # Build string of title banner
    cbsbanners = banner(" CBS Results: {}".format(label))
    print(cbsbanners)
    #    cbsbanners = ''
    #    cbsbanners += """core.print_out('\\n')\n"""
    #    cbsbanners += """p4util.banner(' CBS Results: %s ' % label)\n"""
    #    cbsbanners += """core.print_out('\\n')\n\n"""
    #    exec(cbsbanners)

    # Insert obtained energies into the array that stores the cbs stages
    for stage in GRAND_NEED:
        for lvl in stage["d_need"].items():
            MODELCHEM.append(lvl[1])

            for job in JOBS_EXT:
                pkgmtdlvl = lvl[1]["f_wfn"].split("-", 1)
                pkgmtdjob = job["f_wfn"].split("-", 1)

                # if lvl[1]['f_wfn'][:3] in pkgprefix:
                #    rawlvl = lvl[1]['f_wfn'][3:]
                if (pkgmtdlvl[0] + "-") in driver_util.pkgprefix:
                    rawlvl = pkgmtdlvl[1]
                else:
                    rawlvl = lvl[1]["f_wfn"]

                # if job['f_wfn'][:3] in pkgprefix:
                #    rawjob = job['f_wfn'][3:]
                if (pkgmtdjob[0] + "-") in driver_util.pkgprefix:
                    rawjob = pkgmtdjob[1]
                else:
                    rawjob = job["f_wfn"]

                # Dont ask
                # if (((lvl[1]['f_wfn'] == job['f_wfn']) or
                #     ((lvl[1]['f_wfn'][3:] == job['f_wfn']) and lvl[1]['f_wfn'].startswith('c4-')) or
                #     ((lvl[1]['f_wfn'] == job['f_wfn'][3:]) and job['f_wfn'].startswith('c4-')) or
                #     (('c4-' + lvl[1]['f_wfn']) == job['f_wfn']) or
                #     (lvl[1]['f_wfn'] == ('c4-' + job['f_wfn']))) and
                if rawlvl == rawjob and (lvl[1]["f_basis"] == job["f_basis"]):
                    lvl[1]["f_energy"] = job["f_energy"]
                    lvl[1]["f_gradient"] = job["f_gradient"]
                    lvl[1]["f_hessian"] = job["f_hessian"]

    # Make xtpl() call
    finalenergy = 0.0
    finalgradient = np.zeros((natom, 3))
    finalhessian = np.zeros((3 * natom, 3 * natom))
    for stage in GRAND_NEED:
        hiloargs = _contract_scheme_orders(stage["d_need"], "f_energy")
        stage["d_energy"] = stage["d_scheme"](**hiloargs)
        finalenergy += stage["d_energy"] * stage["d_coef"]

        if ptype == "gradient":
            hiloargs = _contract_scheme_orders(stage["d_need"], "f_gradient")
            stage["d_gradient"] = stage["d_scheme"](**hiloargs)
            finalgradient = stage["d_coef"] * np.copy(stage["d_gradient"]) + finalgradient

        elif ptype == "hessian":
            hiloargs = _contract_scheme_orders(stage["d_need"], "f_gradient")
            stage["d_gradient"] = stage["d_scheme"](**hiloargs)
            finalgradient = stage["d_coef"] * np.copy(stage["d_gradient"]) + finalgradient

            hiloargs = _contract_scheme_orders(stage["d_need"], "f_hessian")
            stage["d_hessian"] = stage["d_scheme"](**hiloargs)
            finalhessian = stage["d_coef"] * np.copy(stage["d_hessian"]) + finalhessian

    # Build string of results table
    table_delimit = "  " + "-" * 105 + "\n"
    tables = ""
    tables += """\n   ==> %s <==\n\n""" % ("Components")
    tables += table_delimit
    tables += """     %6s %20s %1s %-26s %3s %16s   %-s\n""" % (
        "",
        "Method",
        "/",
        "Basis",
        "Rqd",
        "Energy [Eh]",
        "Variable",
    )
    tables += table_delimit
    for job in JOBS_EXT:
        star = ""
        for mc in MODELCHEM:
            if (job["f_wfn"] == mc["f_wfn"]) and (job["f_basis"] == mc["f_basis"]):
                star = "*"
        tables += """     %6s %20s %1s %-27s %2s %16.8f   %-s\n""" % (
            "",
            job["f_wfn"],
            "/",
            job["f_basis"],
            star,
            job["f_energy"],
            VARH[job["f_wfn"]][job["f_wfn"]],
        )
    tables += table_delimit

    tables += """\n   ==> %s <==\n\n""" % ("Stages")
    tables += table_delimit
    tables += """     %6s %20s %1s %-27s %2s %16s   %-s\n""" % (
        "Stage",
        "Method",
        "/",
        "Basis",
        "Wt",
        "Energy [Eh]",
        "Scheme",
    )
    tables += table_delimit
    for stage in GRAND_NEED:
        tables += """     %6s %20s %1s %-27s %2d %16.8f   %-s\n""" % (
            stage["d_stage"],
            stage["d_wfn"],
            "/",
            stage["d_basis"],
            stage["d_coef"],
            stage["d_energy"],
            stage["d_scheme"].__name__,
        )
    tables += table_delimit

    tables += """\n   ==> %s <==\n\n""" % ("CBS")
    tables += table_delimit
    tables += """     %6s %20s %1s %-27s %2s %16s   %-s\n""" % (
        "Stage",
        "Method",
        "/",
        "Basis",
        "",
        "Energy [Eh]",
        "Scheme",
    )
    tables += table_delimit
    if do_scf:
        tables += """     %6s %20s %1s %-27s %2s %16.8f   %-s\n""" % (
            GRAND_NEED[0]["d_stage"],
            GRAND_NEED[0]["d_wfn"],
            "/",
            GRAND_NEED[0]["d_basis"],
            "",
            GRAND_NEED[0]["d_energy"],
            GRAND_NEED[0]["d_scheme"].__name__,
        )
    if do_corl:
        tables += """     %6s %20s %1s %-27s %2s %16.8f   %-s\n""" % (
            GRAND_NEED[1]["d_stage"],
            GRAND_NEED[1]["d_wfn"],
            "/",
            GRAND_NEED[1]["d_basis"],
            "",
            GRAND_NEED[1]["d_energy"] - GRAND_NEED[2]["d_energy"],
            GRAND_NEED[1]["d_scheme"].__name__,
        )
    if do_delta:
        tables += """     %6s %20s %1s %-27s %2s %16.8f   %-s\n""" % (
            GRAND_NEED[3]["d_stage"],
            GRAND_NEED[3]["d_wfn"] + " - " + GRAND_NEED[4]["d_wfn"],
            "/",
            GRAND_NEED[3]["d_basis"],
            "",
            GRAND_NEED[3]["d_energy"] - GRAND_NEED[4]["d_energy"],
            GRAND_NEED[3]["d_scheme"].__name__,
        )
    if do_delta2:
        tables += """     %6s %20s %1s %-27s %2s %16.8f   %-s\n""" % (
            GRAND_NEED[5]["d_stage"],
            GRAND_NEED[5]["d_wfn"] + " - " + GRAND_NEED[6]["d_wfn"],
            "/",
            GRAND_NEED[5]["d_basis"],
            "",
            GRAND_NEED[5]["d_energy"] - GRAND_NEED[6]["d_energy"],
            GRAND_NEED[5]["d_scheme"].__name__,
        )
    #    if do_delta3:
    #        tables += """     %6s %20s %1s %-27s %2s %16.8f   %-s\n""" % (GRAND_NEED[6]['d_stage'], GRAND_NEED[6]['d_wfn'] + ' - ' + GRAND_NEED[7]['d_wfn'],
    #                  '/', GRAND_NEED[6]['d_basis'], '', GRAND_NEED[6]['d_energy'] - GRAND_NEED[7]['d_energy'], GRAND_NEED[6]['d_scheme'].__name__)
    #    if do_delta4:
    #        tables += """     %6s %20s %1s %-27s %2s %16.8f   %-s\n""" % (GRAND_NEED[8]['d_stage'], GRAND_NEED[8]['d_wfn'] + ' - ' + GRAND_NEED[9]['d_wfn'],
    #                  '/', GRAND_NEED[8]['d_basis'], '', GRAND_NEED[8]['d_energy'] - GRAND_NEED[9]['d_energy'], GRAND_NEED[8]['d_scheme'].__name__)
    #    if do_delta5:
    #        tables += """     %6s %20s %1s %-27s %2s %16.8f   %-s\n""" % (GRAND_NEED[10]['d_stage'], GRAND_NEED[10]['d_wfn'] + ' - ' + GRAND_NEED[11]['d_wfn'],
    #                  '/', GRAND_NEED[10]['d_basis'], '', GRAND_NEED[10]['d_energy'] - GRAND_NEED[11]['d_energy'], GRAND_NEED[10]['d_scheme'].__name__)
    tables += """     %6s %20s %1s %-27s %2s %16.8f   %-s\n""" % ("total", "CBS", "", "", "", finalenergy, "")
    tables += table_delimit

    print(tables)

    calcinfo = []
    # core.set_variable('CBS REFERENCE ENERGY', GRAND_NEED[0]['d_energy'])
    # core.set_variable('CBS CORRELATION ENERGY', finalenergy - GRAND_NEED[0]['d_energy'])
    # core.set_variable('CBS TOTAL ENERGY', finalenergy)
    # core.set_variable('CURRENT REFERENCE ENERGY', GRAND_NEED[0]['d_energy'])
    # core.set_variable('CURRENT CORRELATION ENERGY', finalenergy - GRAND_NEED[0]['d_energy'])
    # core.set_variable('CURRENT ENERGY', finalenergy)
    calcinfo.append(Datum("CBS NUMBER", "", Njobs))

    # new skeleton wavefunction w/mol, highest-SCF basis (just to choose one), & not energy
    jobrec = {}
    # TODO hack
    jobrec["molecule"] = jrec["molecule"]
    #    basis = core.BasisSet.build(molecule, "ORBITAL", 'sto-3g')
    #    wfn = core.Wavefunction(molecule, basis)

    #    optstash.restore()

    if ptype == "energy":
        finalquantity = finalenergy
        calcinfo.append(Datum("CURRENT ENERGY", "Eh", finalquantity))

    elif ptype == "gradient":
        calcinfo.append(Datum("CURRENT ENERGY", "Eh", finalenergy))

        finalquantity = finalgradient
        calcinfo.append(Datum("CURRENT GRADIENT", "Eh/a0", finalquantity))
        if finalquantity.shape[0] < 20:
            print("CURRENT GRADIENT")
            print(finalquantity)

    elif ptype == "hessian":
        calcinfo.append(Datum("CURRENT ENERGY", "Eh", finalenergy))
        calcinfo.append(Datum("CURRENT GRADIENT", "Eh/a0", finalgradient))

        finalquantity = finalhessian
        calcinfo.append(Datum("CURRENT HESSIAN", "Eh/a0/a0", finalquantity))
        if finalquantity.shape[0] < 20:
            print("CURRENT GRADIENT")
            print(finalquantity)
            print("CURRENT HESSIAN")
            print(finalquantity)

    jobrec["qcvars"] = {info.label: info for info in calcinfo}
    pp.pprint(jobrec)
    pe.active_qcvars = copy.deepcopy(jobrec["qcvars"])

    if return_wfn:
        return (finalquantity, jobrec)
    else:
        return finalquantity


_lmh_labels = {
    1: ["HI"],
    2: ["LO", "HI"],
    3: ["LO", "MD", "HI"],
    4: ["LO", "MD", "M2", "HI"],
    5: ["LO", "MD", "M2", "M3", "HI"],
}


def _expand_scheme_orders(scheme, basisname, basiszeta, wfnname, natom):
    """Check that the length of *basiszeta* array matches the implied degree of
    extrapolation in *scheme* name. Return a dictionary of same length as
    basiszeta, with *basisname* and *basiszeta* distributed therein.

    """
    Nxtpl = len(basiszeta)

    if int(scheme.__name__.split("_")[-1]) != Nxtpl:
        raise ValidationError("""Call to '%s' not valid with '%s' basis sets.""" % (scheme.__name__, len(basiszeta)))

    f_fields = ["f_wfn", "f_basis", "f_zeta", "f_energy", "f_gradient", "f_hessian"]
    NEED = {}
    for idx in range(Nxtpl):
        NEED[_lmh_labels[Nxtpl][idx]] = dict(
            zip(
                f_fields,
                [wfnname, basisname[idx], basiszeta[idx], 0.0, np.zeros((natom, 3)), np.zeros((3 * natom, 3 * natom))],
            )
        )
    return NEED


def _contract_scheme_orders(needdict, datakey="f_energy"):
    """Prepared named arguments for extrapolation functions by
    extracting zetas and values (which one determined by *datakey*) out
    of *needdict* and returning a dictionary whose keys are contructed
    from _lmh_labels.

    """
    largs = {}
    largs["mtdname"] = needdict["HI"]["f_wfn"]
    Nxtpl = len(needdict)
    zlabels = _lmh_labels[Nxtpl]  # e.g., ['LO', 'HI']

    for zeta in range(Nxtpl):
        zlab = zlabels[zeta]  # e.g., LO
        largs["z" + zlab] = needdict[zlab]["f_zeta"]
        largs["value" + zlab] = needdict[zlab][datakey]

    return largs


_zeta_values = "dtq5678"
_zeta_val2sym = {k + 2: v for k, v in enumerate(_zeta_values)}
_zeta_sym2val = {v: k for k, v in _zeta_val2sym.items()}


#   MOVED   #zeta_values = ['d', 't', 'q', '5', '6', '7', '8']
#   MOVED   #zeta_val2sym = {k + 2: v for k, v in zip(range(7), zeta_values)}
#   MOVED   #zeta_sym2val = {v: k for k, v in zeta_val2sym.items()}


def _expand_bracketed_basis(basisstring, molecule=None):
    r"""Function to transform and validate basis series specification
    *basisstring* for cbs(). A basis set with no paired square brackets is
    passed through with zeta level 0 (e.g., '6-31+G(d,p)' is returned as
    [6-31+G(d,p)] and [0]). A basis set with square brackets is checked
    for sensible sequence and Dunning-ness and returned as separate basis
    sets (e.g., 'cc-pV[Q5]Z' is returned as [cc-pVQZ, cc-pV5Z] and [4,
    5]). This function checks that the basis is valid by trying to build
    the qcdb.BasisSet object for *molecule* or for H2 if None. Allows
    out-of-order zeta specification (e.g., [qtd]) and numeral for number
    (e.g., [23]) but not skipped zetas (e.g., [dq]) or zetas outside [2,
    8] or non-Dunning sets or non-findable .gbs sets.

    """
    import qcdb

    BSET = []
    ZSET = []
    legit_compound_basis = re.compile(
        r"^(?P<pre>.*cc-.*|def2-|.*pcs+eg-|.*)\[(?P<zeta>[dtq2345678,s1]*)\](?P<post>.*z.*|)$", re.IGNORECASE
    )
    pc_basis = re.compile(r".*pcs+eg-$", re.IGNORECASE)
    def2_basis = re.compile(r"def2-", re.IGNORECASE)
    zapa_basis = re.compile(r".*zapa.*", re.IGNORECASE)

    if legit_compound_basis.match(basisstring):
        basisname = legit_compound_basis.match(basisstring)
        # handle def2-svp* basis sets as double-zeta
        if def2_basis.match(basisname.group("pre")):
            bn_gz = basisname.group("zeta").replace("s", "d")
        # handle pc-n basis set polarisation -> zeta conversion
        elif pc_basis.match(basisname.group("pre")):
            bn_gz = basisname.group("zeta").replace("4", "5").replace("3", "4").replace("2", "3").replace("1", "2")
        else:
            bn_gz = basisname.group("zeta")
        # filter out commas and be forgiving of e.g., t5q or 3q
        zetas = [z for z in _zeta_values if (z in bn_gz or str(_zeta_values.index(z) + 2) in bn_gz)]
        for b in zetas:
            if ZSET and (int(ZSET[len(ZSET) - 1]) - _zeta_values.index(b)) != 1:
                raise ValidationError(
                    """Basis set '%s' has skipped zeta level '%s'."""
                    % (basisstring, _zeta_val2sym[_zeta_sym2val[b] - 1])
                )
            # reassemble def2-svp* properly instead of def2-dzvp*
            if def2_basis.match(basisname.group("pre")) and b == "d":
                BSET.append(basisname.group("pre") + "s" + basisname.group("post")[1:])
            # reassemble pc-n basis sets properly
            elif pc_basis.match(basisname.group("pre")):
                BSET.append(basisname.group("pre") + "{0:d}".format(_zeta_sym2val[b] - 1))
            # assemble nZaPa basis sets
            elif zapa_basis.match(basisname.group("post")):
                bzapa = b.replace("d", "2").replace("t", "3").replace("q", "4")
                BSET.append(basisname.group("pre") + bzapa + basisname.group("post"))
            else:
                BSET.append(basisname.group("pre") + b + basisname.group("post"))
            ZSET.append(_zeta_values.index(b) + 2)
    elif re.match(r".*\[.*\].*$", basisstring, flags=re.IGNORECASE):
        raise ValidationError(
            """Basis series '%s' invalid. Specify a basis series matching"""
            """ '*cc-*[dtq2345678,]*z*'. or 'def2-[sdtq]zvp*' or '*pcs[s]eg-[1234]' or '[1234567]ZaPa' """
            % (basisstring)
        )
    else:
        BSET.append(basisstring)
        ZSET.append(0)

    if molecule is None:
        molecule = Molecule("""\nH\nH 1 1.00\n""")

    for basis in BSET:
        try:
            qcdb.BasisSet.pyconstruct(molecule, "BASIS", basis)
        # qbs = BasisSet.pyconstruct(jobrec['molecule'], 'BASIS', jobrec['options']['BASIS'])
        except qcdb.BasisSetNotFound:
            sys.exc_info()[1]
            raise ValidationError("""Basis set '%s' not available for molecule.""" % (basis))

    return (BSET, ZSET)


def _contract_bracketed_basis(basisarray):
    r"""Function to reform a bracketed basis set string from a sequential series
    of basis sets *basisarray* (e.g, form 'cc-pv[q5]z' from array [cc-pvqz, cc-pv5z]).
    Used to print a nicely formatted basis set string in the results table.

    """
    if len(basisarray) == 1:
        return basisarray[0]

    else:
        zetaindx = [i for i in range(len(basisarray[0])) if basisarray[0][i] != basisarray[1][i]][0]
        ZSET = [bas[zetaindx] for bas in basisarray]

        pre = basisarray[0][:zetaindx]
        post = basisarray[0][zetaindx + 1 :]
        basisstring = pre + "[" + "".join(ZSET) + "]" + post
        return basisstring
