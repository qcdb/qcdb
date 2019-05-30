#import re
import math

import numpy as np

from ..exceptions import *
#from . import pe
#from . import driver_util
#from . import driver_helpers
#from .. import moptions
##from . import pe
##from .driver import options
#from .. import util
#from .. import qcvars


#_zeta_values = ['d', 't', 'q', '5', '6', '7', '8']
#_zeta_val2sym = {k + 2: v for k, v in zip(range(7), _zeta_values)}
#_zeta_sym2val = {v: k for k, v in _zeta_val2sym.items()}
_zeta_values = 'dtq5678'
_zeta_val2sym = {k + 2: v for k, v in enumerate(_zeta_values)}
_zeta_sym2val = {v: k for k, v in _zeta_val2sym.items()}


def xtpl_highest_1(mtdname, zHI, valueHI, verbose=1):
    r"""Scheme for total or correlation energies with a single basis or the highest
    zeta-level among an array of bases. Used by :py:func:`qcdb.cbs`.

    .. math:: E_{\textrm{total}}^X = E_{\textrm{total}}^X

    Parameters
    ----------
    mtdname : str
        Method name (e.g., 'mp2') used in summary printing.
    zHI : int
        Zeta number of the basis set.
    valueHI : float or numpy.ndarray
        Energy, gradient, or Hessian value at the basis set.
    verbose : int, optional
        Controls volume of printing.

    Returns
    -------
    float or numpy.ndarray
        Eponymous function applied to input zetas and values; type from `valueHI`.

    Examples
    --------
    >>> # [1] Fancy way to get HF/cc-pCVQZ
    >>> qcdb.energy(qcdb.cbs, scf_wfn='hf', scf_basis='cc-pcvqz', scf_scheme=qcdb.xtpl_highest_1)

    """
    if isinstance(valueHI, float):

        if verbose:
            # Output string with extrapolation parameters
            cbsscheme = ''
            cbsscheme += """\n   ==> {} <==\n\n""".format(mtdname.upper())
            cbsscheme += """   HI-zeta ({}) Energy:               {:16.12f}\n""".format(zHI, valueHI)

            print(cbsscheme)

        return valueHI

    elif isinstance(valueHI, np.ndarray):

        if verbose > 2:
            core.print_out("""   HI-zeta (%s) Total Energy:\n""" % (str(zHI)))
            valueHI.print_out()

        return valueHI


def scf_xtpl_helgaker_2(mtdname, zLO, valueLO, zHI, valueHI, alpha=1.63, verbose=1):
    r"""Extrapolation scheme for reference energies with two adjacent zeta-level bases.
    Used by :py:func:`qcdb.cbs`.
    `Halkier, Helgaker, Jorgensen, Klopper, & Olsen, Chem. Phys. Lett. 302 (1999) 437-446
    <https://doi.org/10.1016/S0009-2614(99)00179-7>`_

    .. math:: E_{\textrm{total}}^X = E_{\textrm{total}}^{\infty} + \beta e^{-\alpha X}, \alpha = 1.63

    Parameters
    ----------
    mtdname : str
        Method name (e.g., 'HF') used in summary printing.
    zLO : int
        Zeta number of the smaller basis set in 2-point extrapolation.
    valueLO : float or numpy.ndarray
        Energy, gradient, or Hessian value at the smaller basis set in 2-point
        extrapolation.
    zHI : int
        Zeta number of the larger basis set in 2-point extrapolation.
        Must be `zLO + 1`.
    valueHI : float or numpy.ndarray
        Energy, gradient, or Hessian value at the larger basis set in 2-point
        extrapolation.
    verbose : int, optional
        Controls volume of printing.
    alpha : float, optional
        Fitted 2-point parameter.

    Returns
    -------
    float or numpy.ndarray
        Eponymous function applied to input zetas and values; type from `valueLO`.

    Examples
    --------
    >>> # [1] Hartree-Fock extrapolation
    >>> qcdb.energy(qcdb.cbs, scf_wfn='hf', scf_basis='cc-pV[DT]Z', scf_scheme=qcdb.scf_xtpl_helgaker_2)

    """
    if type(valueLO) != type(valueHI):
        raise ValidationError("scf_xtpl_helgaker_2: Inputs must be of the same datatype! (%s, %s)"
                              % (type(valueLO), type(valueHI)))

    beta_division = 1 / (math.exp(-1 * alpha * zLO) * (math.exp(-1 * alpha) - 1))
    beta_mult = math.exp(-1 * alpha * zHI)

    if isinstance(valueLO, float):
        beta = (valueHI - valueLO) / (math.exp(-1 * alpha * zLO) * (math.exp(-1 * alpha) - 1))
        value = valueHI - beta * math.exp(-1 * alpha * zHI)

        if verbose:
            # Output string with extrapolation parameters
            cbsscheme = ''
            cbsscheme += """\n   ==> Helgaker 2-point SCF extrapolation for method: %s <==\n\n""" % (mtdname.upper())
            cbsscheme += """   LO-zeta (%s) Energy:               % 16.12f\n""" % (str(zLO), valueLO)
            cbsscheme += """   HI-zeta (%s) Energy:               % 16.12f\n""" % (str(zHI), valueHI)
            cbsscheme += """   Alpha (exponent) Value:           % 16.12f\n""" % (alpha)
            cbsscheme += """   Beta (coefficient) Value:         % 16.12f\n\n""" % (beta)

            name_str = "%s/(%s,%s)" % (mtdname.upper(), _zeta_val2sym[zLO].upper(), _zeta_val2sym[zHI].upper())
            cbsscheme += """   @Extrapolated """
            cbsscheme += name_str + ':'
            cbsscheme += " " * (18 - len(name_str))
            cbsscheme += """% 16.12f\n\n""" % value
            print(cbsscheme)

        return value

    elif isinstance(valueLO, np.ndarray):
        beta = (valueHI - valueLO) * beta_division
        value = valueHI - beta * beta_mult
        #beta.name = 'Helgaker SCF (%s, %s) beta' % (zLO, zHI)
        #value.name = 'Helgaker SCF (%s, %s) data' % (zLO, zHI)

        if verbose > 2:
            core.print_out("""\n   ==> Helgaker 2-point SCF extrapolation for method: %s <==\n\n""" % (mtdname.upper()))
            core.print_out("""   LO-zeta ({})""".format(zLO))
            core.print_out("""   LO-zeta Data""")
            print(valueLO)
            core.print_out("""   HI-zeta (%s)""" % str(zHI))
            core.print_out("""   HI-zeta Data""")
            valueHI.print_out()
            core.print_out("""   Extrapolated Data:\n""")
            value.print_out()
            core.print_out("""   Alpha (exponent) Value:          %16.8f\n""" % (alpha))
            core.print_out("""   Beta Data:\n""")
            beta.print_out()

        return value

    else:
        raise ValidationError("scf_xtpl_helgaker_2: datatype is not recognized '%s'." % type(valueLO))


def scf_xtpl_helgaker_3(mtdname, zLO, valueLO, zMD, valueMD, zHI, valueHI, verbose=1):
    r"""Extrapolation scheme for reference energies with three adjacent zeta-level bases.
    Used by :py:func:`qcdb.cbs`.
    `Halkier, Helgaker, Jorgensen, Klopper, & Olsen, Chem. Phys. Lett. 302 (1999) 437-446
    <https://doi.org/10.1016/S0009-2614(99)00179-7>`_

    .. math:: E_{\textrm{total}}^X = E_{\textrm{total}}^{\infty} + \beta e^{-\alpha X}

    Parameters
    ----------
    mtdname : str
        Method name (e.g., 'HF') used in summary printing.
    zLO : int
        Zeta number of the smaller basis set in 3-point extrapolation.
    valueLO : float or numpy.ndarray
        Energy, gradient, or Hessian value at the smaller basis set in 3-point
        extrapolation.
    zMD : int
        Zeta number of the medium basis set in 3-point extrapolation.
        Must be `zLO + 1`.
    valueMD : float or numpy.ndarray
        Energy, gradient, or Hessian value at the medium basis set in 3-point
        extrapolation.
    zHI : int
        Zeta number of the larger basis set in 3-point extrapolation.
        Must be `zLO + 2`.
    valueHI : float or numpy.ndarray
        Energy, gradient, or Hessian value at the larger basis set in 3-point
        extrapolation.
    verbose : int, optional
        Controls volume of printing.

    Returns
    -------
    float or numpy.ndarray
        Eponymous function applied to input zetas and values; type from `valueLO`.

    Examples
    --------
    >>> # [1] Hartree-Fock extrapolation
    >>> qcdb.energy(qcdb.cbs, scf_wfn='hf', scf_basis='cc-pV[DTQ]Z', scf_scheme=qcdb.scf_xtpl_helgaker_3)

    """
    if (type(valueLO) != type(valueMD)) or (type(valueMD) != type(valueHI)):
        raise ValidationError("scf_xtpl_helgaker_3: Inputs must be of the same datatype! (%s, %s, %s)"
                              % (type(valueLO), type(valueMD), type(valueHI)))

    if isinstance(valueLO, float):

        ratio = (valueHI - valueMD) / (valueMD - valueLO)
        alpha = -1 * math.log(ratio)
        beta = (valueHI - valueMD) / (math.exp(-1 * alpha * zMD) * (ratio - 1))
        value = valueHI - beta * math.exp(-1 * alpha * zHI)

        if verbose:
            # Output string with extrapolation parameters
            cbsscheme = ''
            cbsscheme += """\n   ==> Helgaker 3-point SCF extrapolation for method: %s <==\n\n""" % (mtdname.upper())
            cbsscheme += """   LO-zeta (%s) Energy:               % 16.12f\n""" % (str(zLO), valueLO)
            cbsscheme += """   MD-zeta (%s) Energy:               % 16.12f\n""" % (str(zMD), valueMD)
            cbsscheme += """   HI-zeta (%s) Energy:               % 16.12f\n""" % (str(zHI), valueHI)
            cbsscheme += """   Alpha (exponent) Value:           % 16.12f\n""" % (alpha)
            cbsscheme += """   Beta (coefficient) Value:         % 16.12f\n\n""" % (beta)

            name_str = "%s/(%s,%s,%s)" % (mtdname.upper(), _zeta_val2sym[zLO].upper(), _zeta_val2sym[zMD].upper(),
                                                             _zeta_val2sym[zHI].upper())
            cbsscheme += """   @Extrapolated """
            cbsscheme += name_str + ':'
            cbsscheme += " " * (18 - len(name_str))
            cbsscheme += """% 16.12f\n\n""" % value
            print(cbsscheme)

        return value

    elif isinstance(valueLO, np.ndarray):
        valueLO = np.array(valueLO)
        valueMD = np.array(valueMD)
        valueHI = np.array(valueHI)

        nonzero_mask = np.abs(valueHI) > 1.e-14
        top = (valueHI - valueMD)[nonzero_mask]
        bot = (valueMD - valueLO)[nonzero_mask]

        ratio = top / bot
        alpha = -1 * np.log(np.abs(ratio))
        beta = top / (np.exp(-1 * alpha * zMD) * (ratio - 1))
        np_value = valueHI.copy()
        np_value[nonzero_mask] -= beta * np.exp(-1 * alpha * zHI)
        np_value[~nonzero_mask] = 0.0
        return np_value

        ## Build and set from numpy routines
        #value = core.Matrix(*valueHI.shape)
        #value_view = np.asarray(value)
        #value_view[:] = np_value
        #return value

    else:
        raise ValidationError("scf_xtpl_helgaker_3: datatype is not recognized '%s'." % type(valueLO))


def corl_xtpl_helgaker_2(mtdname, zLO, valueLO, zHI, valueHI, verbose=1):
    r"""Extrapolation scheme for correlation energies with two adjacent zeta-level bases.
    Used by :py:func:`qcdb.cbs`.
    `Halkier, Helgaker, Jorgensen, Klopper, Koch, Olsen, & Wilson, Chem. Phys. Lett. 286 (1998) 243-252
    <https://doi.org/10.1016/S0009-2614(98)00111-0>`_

    .. math:: E_{\textrm{corl}}^X = E_{\textrm{corl}}^{\infty} + \beta X^{-3}

    Parameters
    ----------
    mtdname : str
        Method name (e.g., 'MP2') used in summary printing.
    zLO : int
        Zeta number of the smaller basis set in 2-point extrapolation.
    valueLO : float or numpy.ndarray
        Energy, gradient, or Hessian value at the smaller basis set in 2-point
        extrapolation.
    zHI : int
        Zeta number of the larger basis set in 2-point extrapolation.
        Must be `zLO + 1`.
    valueHI : float or numpy.ndarray
        Energy, gradient, or Hessian value at the larger basis set in 2-point
        extrapolation.
    verbose : int, optional
        Controls volume of printing.

    Returns
    -------
    float or numpy.ndarray
        Eponymous function applied to input zetas and values; type from `valueLO`.

    Examples
    --------
    >>> # [1] CISD extrapolation
    >>> qcdb.energy(qcdb.cbs, corl_wfn='cisd', corl_basis='cc-pV[DT]Z', corl_scheme=qcdb.corl_xtpl_helgaker_2)

    """
    if type(valueLO) != type(valueHI):
        raise ValidationError("corl_xtpl_helgaker_2: Inputs must be of the same datatype! (%s, %s)"
                              % (type(valueLO), type(valueHI)))

    if isinstance(valueLO, float):
        value = (valueHI * zHI ** 3 - valueLO * zLO ** 3) / (zHI ** 3 - zLO ** 3)
        beta = (valueHI - valueLO) / (zHI ** (-3) - zLO ** (-3))

#        final = valueSCF + value
        final = value
        if verbose:
            # Output string with extrapolation parameters
            cbsscheme = """\n\n   ==> Helgaker 2-point correlated extrapolation for method: %s <==\n\n""" % (mtdname.upper())
#            cbsscheme += """   HI-zeta (%1s) SCF Energy:           % 16.12f\n""" % (str(zHI), valueSCF)
            cbsscheme += """   LO-zeta (%s) Energy:               % 16.12f\n""" % (str(zLO), valueLO)
            cbsscheme += """   HI-zeta (%s) Energy:               % 16.12f\n""" % (str(zHI), valueHI)
#            cbsscheme += """   Beta (coefficient) Value:         % 16.12f\n""" % beta
            cbsscheme += """   Extrapolated Energy:              % 16.12f\n\n""" % value
            #cbsscheme += """   LO-zeta (%s) Correlation Energy:   % 16.12f\n""" % (str(zLO), valueLO)
            #cbsscheme += """   HI-zeta (%s) Correlation Energy:   % 16.12f\n""" % (str(zHI), valueHI)
            #cbsscheme += """   Beta (coefficient) Value:         % 16.12f\n""" % beta
            #cbsscheme += """   Extrapolated Correlation Energy:  % 16.12f\n\n""" % value

            name_str = "%s/(%s,%s)" % (mtdname.upper(), _zeta_val2sym[zLO].upper(), _zeta_val2sym[zHI].upper())
            cbsscheme += """   @Extrapolated """
            cbsscheme += name_str + ':'
            cbsscheme += " " * (19 - len(name_str))
            cbsscheme += """% 16.12f\n\n""" % final
            print(cbsscheme)

        return final

    elif isinstance(valueLO, np.ndarray):
        beta = (valueHI - valueLO) / (zHI ** (-3) - zLO ** (-3))
        #beta.name = 'Helgaker Corl (%s, %s) beta' % (zLO, zHI)

        value = (valueHI * zHI ** 3 - valueLO * zLO ** 3) / (zHI ** 3 - zLO ** 3)
        #value.name = 'Helgaker Corr (%s, %s) data' % (zLO, zHI)

        if verbose > 2:
            core.print_out("""\n   ==> Helgaker 2-point correlated extrapolation for """
                           """method: %s <==\n\n""" % (mtdname.upper()))
            core.print_out("""   LO-zeta (%s) Data\n""" % (str(zLO)))
            valueLO.print_out()
            core.print_out("""   HI-zeta (%s) Data\n""" % (str(zHI)))
            valueHI.print_out()
            core.print_out("""   Extrapolated Data:\n""")
            value.print_out()
            core.print_out("""   Beta Data:\n""")
            beta.print_out()

#        value.add(valueSCF)
        return value

    else:
        raise ValidationError("corl_xtpl_helgaker_2: datatype is not recognized '%s'." % type(valueLO))

