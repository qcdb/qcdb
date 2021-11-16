from typing import Dict, List

import numpy as np
#from psi4.driver.p4util.exceptions import *
from qcelemental import constants


def least_squares_fit_polynomial(xvals, fvals, localization_point, no_factorials=True, weighted=True, polynomial_order=4):
    """Performs and unweighted least squares fit of a polynomial, with specified order
       to an array of input function values (fvals) evaluated at given locations (xvals).
       See https://doi.org/10.1063/1.4862157, particularly eqn (7) for details. """
    xpts = np.array(xvals) - localization_point
    if weighted:
        R = 1.0
        p_nu = 1
        epsilon = 1e-3
        zvals = np.square(xpts/R)
        weights = np.exp(-zvals) / (zvals**p_nu + epsilon**p_nu)
    else:
        weights = None
    fit = np.polynomial.polynomial.polyfit(xpts, fvals, polynomial_order, w=weights)
    # Remove the 1/n! coefficients
    if no_factorials:
        scalefac = 1.0
        for n in range(2,polynomial_order+1):
            scalefac *= n
            fit[n] *= scalefac
    return fit


#def anharmonicity(rvals: List, energies: List, molecule, plot_fit: str = '') -> Dict:
def diatomic(rvals: List, energies: List, molecule, plot_fit: str = '') -> Dict:
    """Generates spectroscopic constants for a diatomic molecules.
       Fits a diatomic potential energy curve using a weighted least squares approach
       (c.f. https://doi.org/10.1063/1.4862157, particularly eqn. 7), locates the minimum
       energy point, and then applies second order vibrational perturbation theory to obtain spectroscopic
       constants.  Any number of points greater than 4 may be provided, and they should bracket the minimum.
       The data need not be evenly spaced, and can be provided in any order.  The data are weighted such that
       those closest to the minimum have highest impact.

       A dictionary with the following keys, which correspond to spectroscopic constants, is returned:

       :param rvals: The bond lengths (in Angstrom) for which energies are
           provided, of length at least 5 and equal to the length of the energies array

       :param energies: The energies (Eh) computed at the bond lengths in the rvals list

       :param plot_fit: A string describing where to save a plot of the harmonic and anharmonic fits, the
           inputted data points, re, r0 and the first few energy levels, if matplotlib
           is available.  Set to 'screen' to generate an interactive plot on the screen instead. If a filename is
           provided, the image type is determined by the extension; see matplotlib for supported file types.

       :returns: (*dict*) Keys: "re", "r0", "we", "wexe", "nu", "ZPVE(harmonic)", "ZPVE(anharmonic)", "Be", "B0", "ae", "De"
                 corresponding to the spectroscopic constants in cm-1
    """

    1.0 / constants.bohr2angstroms
    angstrom_to_meter = 10e-10

    # Make sure the input is valid
    if len(rvals) != len(energies):
        raise ValidationError("The number of energies must match the number of distances")
    npoints = len(rvals)
    if npoints < 5:
        raise ValidationError("At least 5 data points must be provided to compute anharmonicity")
    print("\n\nPerforming a fit to %d data points\n" % npoints)

    # Sort radii and values first from lowest to highest radius
    indices = np.argsort(rvals)
    rvals = np.array(rvals)[indices]
    energies = np.array(energies)[indices]

    # Make sure the molecule the user provided is the active one
    molecule.update_geometry()
    natoms = molecule.natom()
    if natoms != 2:
        raise Exception("The current molecule must be a diatomic for this code to work!")
    m1 = molecule.mass(0)
    m2 = molecule.mass(1)

    # Find rval of the minimum of energies, check number of points left and right
    min_index = np.argmin(energies)
    if min_index < 3 :
        print("\nWarning: fewer than 3 points provided with a r < r(min(E))!\n")
    if min_index >= len(energies) - 3:
        print("\nWarning: fewer than 3 points provided with a r > r(min(E))!\n")

    # Optimize the geometry, refitting the surface around each new geometry
    print("\nOptimizing geometry based on current surface:\n\n")
    re = rvals[min_index]
    maxit = 30
    thres = 1.0e-9
    for i in range(maxit):
        derivs = least_squares_fit_polynomial(rvals,energies,localization_point=re)
        e,g,H = derivs[0:3]
        print("       E = %20.14f, x = %14.7f, grad = %20.14f\n" % (e, re, g))
        if abs(g) < thres:
            break
        re -= g/H;
        if i == maxit-1:
            raise ConvergenceError("diatomic geometry optimization", maxit)
    print(" Final E = %20.14f, x = %14.7f, grad = %20.14f\n" % (e, re, g))
    if re < min(rvals):
        raise Exception("Minimum energy point is outside range of points provided.  Use a lower range of r values.")
    if re > max(rvals):
        raise Exception("Minimum energy point is outside range of points provided.  Use a higher range of r values.")

    # Convert to convenient units, and compute spectroscopic constants
    d0,d1,d2,d3,d4 = derivs*constants.hartree2aJ
    print("\nEquilibrium Energy %20.14f Hartrees\n" % e)
    print("Gradient           %20.14f\n" % g)
    print("Quadratic Force Constant %14.7f MDYNE/A\n" % d2)
    print("Cubic Force Constant     %14.7f MDYNE/A**2\n" % d3)
    print("Quartic Force Constant   %14.7f MDYNE/A**3\n" % d4)

    constants.h / (2.0 * np.pi)
    mu = ((m1*m2)/(m1+m2))*constants.amu2kg
    we = 5.3088375e-11 * np.sqrt(d2/mu)
    wexe = (1.2415491e-6)*(we/d2)**2 * ((5.0*d3*d3)/(3.0*d2)-d4)

    # Rotational constant: Be
    I = ((m1*m2)/(m1+m2)) * constants.amu2kg * (re * angstrom_to_meter)**2
    B = constants.h / (8.0 * np.pi**2 * constants.c * I)

    # alpha_e and quartic centrifugal distortion constant
    ae = -(6.0 * B**2 / we) * ((1.05052209e-3*we*d3)/(np.sqrt(B * d2**3))+1.0)
    de = 4.0*B**3 / we**2

    # B0 and r0 (plus re check using Be)
    B0 = B - ae / 2.0
    r0 = np.sqrt(constants.h / (8.0 * np.pi**2 * mu * constants.c * B0))
    recheck = np.sqrt(constants.h / (8.0 * np.pi**2 * mu * constants.c * B))
    r0 /= angstrom_to_meter
    recheck /= angstrom_to_meter

    # Fundamental frequency nu
    nu = we - 2.0 * wexe
    zpve_nu = 0.5 * we - 0.25 * wexe
    zpve_we = 0.5 * we

    # Generate pretty pictures, if requested
    if(plot_fit):
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            msg = "\n\tPlot not generated; matplotlib is not installed on this machine.\n\n"
            print(msg)

        # Correct the derivatives for the missing factorial prefactors
        dvals = np.zeros(5)
        dvals[0:5] = derivs[0:5]
        dvals[2] /= 2
        dvals[3] /= 6
        dvals[4] /= 24

        # Default plot range, before considering energy levels
        minE = np.min(energies)
        maxE = np.max(energies)
        minR = np.min(rvals)
        maxR = np.max(rvals)

        # Plot vibrational energy levels
        we_au = we / constants.hartree2wavenumbers
        wexe_au = wexe / constants.hartree2wavenumbers
        coefs2 = [ dvals[2], dvals[1], dvals[0] ]
        coefs4 = [ dvals[4], dvals[3], dvals[2], dvals[1], dvals[0] ]
        for n in range(3):
            Eharm = we_au*(n+0.5)
            Evpt2 = Eharm - wexe_au*(n+0.5)**2
            coefs2[-1] = -Eharm
            coefs4[-1] = -Evpt2
            roots2 = np.roots(coefs2)
            roots4 = np.roots(coefs4)
            xvals2 = roots2 + re
            xvals4 = np.choose(np.where(np.isreal(roots4)), roots4)[0].real + re
            Eharm += dvals[0]
            Evpt2 += dvals[0]
            plt.plot(xvals2, [Eharm, Eharm], 'b', linewidth=1)
            plt.plot(xvals4, [Evpt2, Evpt2], 'g', linewidth=1)
            maxE = Eharm
            maxR = np.max([xvals2,xvals4])
            minR = np.min([xvals2,xvals4])

        # Find ranges for the plot
        dE = maxE - minE
        minE -= 0.2*dE
        maxE += 0.4*dE
        dR = maxR - minR
        minR -= 0.2*dR
        maxR += 0.2*dR

        # Generate the fitted PES
        xpts = np.linspace(minR, maxR, 1000)
        xrel = xpts - re
        xpows = xrel[:, None] ** range(5)
        fit2 = np.einsum('xd,d', xpows[:,0:3], dvals[0:3])
        fit4 = np.einsum('xd,d', xpows, dvals)

        # Make / display the plot
        plt.plot(xpts, fit2, 'b', linewidth=2.5, label='Harmonic (quadratic) fit')
        plt.plot(xpts, fit4, 'g', linewidth=2.5, label='Anharmonic (quartic) fit')
        plt.plot([re, re], [minE, maxE], 'b--', linewidth=0.5)
        plt.plot([r0, r0], [minE, maxE], 'g--', linewidth=0.5)
        plt.scatter(rvals, energies, c='Black', linewidth=3, label='Input Data')
        plt.legend()

        plt.xlabel('Bond length (Angstroms)')
        plt.ylabel('Energy (Eh)')
        plt.xlim(minR, maxR)
        plt.ylim(minE, maxE)
        if plot_fit == 'screen':
            plt.show()
        else:
            plt.savefig(plot_fit)
            print("\n\tPES fit saved to %s.\n\n" % plot_fit)

    print("\nre       = %10.6f A  check: %10.6f\n" % (re, recheck))
    print("r0       = %10.6f A\n" % r0)
    print("E at re  = %17.10f Eh\n" % e)
    print("we       = %10.4f cm-1\n" % we)
    print("wexe     = %10.4f cm-1\n" % wexe)
    print("nu       = %10.4f cm-1\n" % nu)
    print("ZPVE(we) = %10.4f cm-1\n" % zpve_we)
    print("ZPVE(nu) = %10.4f cm-1\n" % zpve_nu)
    print("Be       = %10.4f cm-1\n" % B)
    print("B0       = %10.4f cm-1\n" % B0)
    print("ae       = %10.4f cm-1\n" % ae)
    print("De       = %10.7f cm-1\n" % de)
    results = {
               "re"               :  re,
               "r0"               :  r0,
               "we"               :  we,
               "wexe"             :  wexe,
               "nu"               :  nu,
               "E(re)"            :  e,
               "ZPVE(harmonic)"   :  zpve_we,
               "ZPVE(anharmonic)" :  zpve_nu,
               "Be"               :  B,
               "B0"               :  B0,
               "ae"               :  ae,
               "De"               :  de
              }
    return results
