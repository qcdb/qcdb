import math

import qcelemental as qcel


def process_units(molrec):
    """From any (not both None) combination of `units` and
    `input_units_to_au`, returns both quantities validated. The degree
    of checking is unnecessary if coming from a molrec (prevalidated and
    guaranteed to have "units"), but function is general-purpose.

    """
    units = molrec.get('units', None)
    input_units_to_au = molrec.get('input_units_to_au', None)

    b2a = qcel.constants.bohr2angstroms
    a2b = 1. / b2a

    def perturb_check(candidate, reference):
        return (abs(candidate, reference) < 0.05)

    if units is None and input_units_to_au is not None:
        if perturb_check(input_units_to_au, 1.):
            funits = 'Bohr'
            fiutau = input_units_to_au
        elif perturb_check(input_units_to_au, a2b):
            funits = 'Angstrom'
            fiutau = input_units_to_au
        else:
            raise ValidationError("""No big perturbations to physical constants! {} !~= ({} or {})""".format(
                input_units_to_au, 1.0, a2b))

    elif units in ['Angstrom', 'Bohr'] and input_units_to_au is None:
        funits = units

        if funits == 'Bohr':
            fiutau = 1.
        elif funits == 'Angstrom':
            fiutau = a2b

    elif units in ['Angstrom', 'Bohr'] and input_units_to_au is not None:
        expected_iutau = a2b if units == 'Angstrom' else 1.

        if perturb_check(input_units_to_au, expected_iutau):
            funits = units
            fiutau = input_units_to_au
        else:
            raise ValidationError("""No big perturbations to physical constants! {} !~= {}""".format(
                input_units_to_au, expected_iutau))

    else:
        raise ValidationError('Insufficient information: {} & {}'.format(units, input_units_to_au))

    return funits, fiutau


def conv_float2negexp(val: float) -> int:
    """Least restrictive negative exponent of base 10 that achieves the floating point convergence criterium `val`."""
    return -1 * int(math.floor(math.log(val, 10)))


def program_prefix(hint: str, *, return_program: bool = False, return_prefix: bool = False) -> str:
    """Translate between CMS package name and method prefix."""

    pkgprefix = {
        "p4": "Psi4",
        "c4": "CFOUR",
        "d3": "DFTD3",
        "nwc": "NWChem",
        "gms": "GAMESS",
    }

    lookup = {}
    for pfx, prog in pkgprefix.items():
        dashed_pfx = pfx + "-"
        lookup[pfx] = (prog, dashed_pfx, prog)
        lookup[dashed_pfx] = (prog, dashed_pfx, prog)
        lookup[prog.lower()] = (dashed_pfx, dashed_pfx, prog)

    if return_program:
        return lookup[hint.lower()][2]
    elif return_prefix:
        return lookup[hint.lower()][1]
    else:
        return lookup[hint.lower()][0]
