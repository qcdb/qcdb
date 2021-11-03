import re

from ..exceptions import KeywordValidationError


def enum(inputval):
    allowed = inputval.upper().split()

    def closedenum(x, allowed=allowed):
        if x.upper() in allowed:
            return x.upper()
        else:
            raise KeywordValidationError("""Not allowed value: {} not in {}""".format(inputval, allowed))

    return closedenum


def enum_bool(inputval):
    allowed = inputval.upper().split()

    def closedenum(x, allowed=allowed):
        if x.upper() in allowed:
            if x.upper() == "TRUE":
                return "ON"
            elif x.upper() == "FALSE":
                return "OFF"
            else:
                return x.upper()
        else:
            raise KeywordValidationError("""Not allowed value: {} not in {}""".format(inputval, allowed))

    return closedenum


def onoff_boolean(inputval):
    yes = re.compile(r"^(yes|true|on|1)", re.IGNORECASE)
    no = re.compile(r"^(no|false|off|0)", re.IGNORECASE)

    if yes.match(str(inputval)):
        return "on"
    elif no.match(str(inputval)):
        return "off"
    else:
        raise KeywordValidationError("""Can't interpret into boolean: {}""".format(inputval))


def intenum(inputval, nullable=False):
    allowed = [int(x) for x in inputval.split()]
    if nullable:
        allowed.append(None)

    def closedenum(x, allowed=allowed):
        if x in allowed:
            return x
        else:
            raise KeywordValidationError("""Not allowed integer value: {} not in {}""".format(inputval, allowed))

    return closedenum


def mixedenum(inputval, nullable=False):
    allowed = []
    for x in inputval.split():
        try:
            val = int(x)
        except ValueError:
            val = x

        allowed.append(val)

    if nullable:
        allowed.append(None)

    def closedenum(x, allowed=allowed):
        if x in allowed:
            return x
        else:
            raise KeywordValidationError("""Not allowed mixed value: {} not in {}""".format(inputval, allowed))

    return closedenum


def casesensitive_enum(inputval):
    allowed = inputval.split()

    def closedenum(x, allowed=allowed):
        if x in allowed:
            return x
        else:
            raise KeywordValidationError("""Not allowed case-sensitive value: {} not in {}""".format(inputval, allowed))

    return closedenum


def boolean(inputval):
    yes = re.compile(r"^(yes|true|on|1)", re.IGNORECASE)
    no = re.compile(r"^(no|false|off|0)", re.IGNORECASE)

    if yes.match(str(inputval)):
        return True
    elif no.match(str(inputval)):
        return False
    else:
        raise KeywordValidationError("""Can't interpret into boolean: {}""".format(inputval))


def sphcart(inputval):
    sph = re.compile(r"^(yes|true|on|1|sph|spherical)", re.IGNORECASE)
    cart = re.compile(r"^(no|false|off|0|cart|cartesian)", re.IGNORECASE)

    if sph.match(str(inputval)):
        return True
    elif cart.match(str(inputval)):
        return False
    else:
        raise KeywordValidationError("""Can't interpret into boolean True (sph) or False (cart): {}""".format(inputval))


def gridradang(inputval):

    if isinstance(inputval, dict):
        retdict = {}
        for k, v in inputval.items():
            rad, ang = v
            retdict[k] = (positive_integer(rad), positive_integer(ang))
        return dict(sorted(retdict.items()))  # to place '' key first
    else:
        rad, ang = inputval
        return {"": (positive_integer(rad), positive_integer(ang))}


def bool_or_elem_dict(inputval):
    if isinstance(inputval, dict):
        retdict = {k.capitalize(): positive_integer(v) for k, v in inputval.items()}
        return dict(sorted(retdict.items()))
    else:
        return boolean(inputval)


def atompair(inputval):

    if isinstance(inputval, dict):
        retdict = {}
        for k, v in inputval.items():
            atom1, atom2 = v
            retdict[k] = (positive_integer(atom1), positive_integer(atom2))
        return dict(sorted(retdict.items()))  # to place '' key first
    else:
        return {""}


def percentage(inputval):
    if 0.0 <= inputval <= 100.0:
        return float(inputval)
    else:
        raise KeywordValidationError("Percentage should be between 0 and 100: {}".format(inputval))


def nonnegative_float(inputval):
    if 0.0 <= inputval:
        return float(inputval)
    else:
        raise KeywordValidationError("Float should be non-negative: {}".format(inputval))


def positive_integer(inputval):
    if inputval > 0 and float(inputval).is_integer():
        return int(inputval)
    else:
        raise KeywordValidationError("Positive integer, if you please: {}".format(inputval))


def nonnegative_integer(inputval):
    if inputval > -1 and float(inputval).is_integer():
        return int(inputval)
    else:
        raise KeywordValidationError("Non-negative integer number, if you please: {}".format(inputval))


def integer(inputval):
    if float(inputval).is_integer():
        return int(inputval)
    else:
        raise KeywordValidationError("Integer number, if you please: {}".format(inputval))


def parse_convergence(inputval):

    if inputval > 0 and isinstance(inputval, int):
        return pow(10.0, -inputval)
    elif inputval > 0 and inputval < 5:
        return inputval
    else:
        raise KeywordValidationError("wth! you call this a convergence criterion? {}".format(inputval))


def parse_memory(inputval, min_mem_allowed=262144000):
    """Validates expression for total memory allocation. Takes memory value
    `inputval` as type int, float, or str; int and float are taken literally
    as bytes to be set, string taken as a unit-containing value (e.g., 30 mb)
    which is case-insensitive.

    :returns: *memory_amount* (float) Number of bytes of memory

    :raises: ValidationError when <500MiB or disallowed type or misformatted

    :examples:

    >>> # [1] Passing absolute number of bytes
    >>> psi4.set_memory(600000000)
    >>> psi4.get_memory()
    Out[1]: 600000000L

    >>> # [2] Passing memory value as string with units
    >>> psi4.set_memory('30 GB')
    >>> psi4.get_memory()
    Out[2]: 30000000000L

    >>> # Good examples
    >>> psi4.set_memory(800000000)        # 800000000
    >>> psi4.set_memory(2004088624.9)     # 2004088624
    >>> psi4.set_memory(1.0e9)            # 1000000000
    >>> psi4.set_memory('600 mb')         # 600000000
    >>> psi4.set_memory('600.0 MiB')      # 629145600
    >>> psi4.set_memory('.6 Gb')          # 600000000
    >>> psi4.set_memory(' 100000000kB ')  # 100000000000
    >>> psi4.set_memory('2 eb')           # 2000000000000000000

    >>> # Bad examples
    >>> psi4.set_memory({})         # odd type
    >>> psi4.set_memory('')         # no info
    >>> psi4.set_memory("8 dimms")  # unacceptable units
    >>> psi4.set_memory("1e5 gb")   # string w/ exponent
    >>> psi4.set_memory("5e5")      # string w/o units
    >>> psi4.set_memory(2000)       # mem too small
    >>> psi4.set_memory(-5e5)       # negative (and too small)

    """
    # Handle memory given in bytes directly (int or float)
    if isinstance(inputval, (int, float)):
        val = inputval
        units = ""
    # Handle memory given as a string
    elif isinstance(inputval, str):
        memory_string = re.compile(r"^\s*(\d*\.?\d+)\s*([KMGTPBE]i?B)\s*$", re.IGNORECASE)
        matchobj = re.search(memory_string, inputval)
        if matchobj:
            val = float(matchobj.group(1))
            units = matchobj.group(2)
        else:
            raise KeywordValidationError(
                """Invalid memory specification: {}. Try 5e9 or '5 gb'.""".format(repr(inputval))
            )
    else:
        raise KeywordValidationError(
            """Invalid type {} in memory specification: {}. Try 5e9 or '5 gb'.""".format(type(inputval), repr(inputval))
        )

    # Units decimal or binary?
    multiplier = 1000
    if "i" in units.lower():
        multiplier = 1024
        units = units.lower().replace("i", "").upper()

    # Build conversion factor, convert units
    unit_list = ["", "KB", "MB", "GB", "TB", "PB", "EB"]
    mult = 1
    for unit in unit_list:
        if units.upper() == unit:
            break
        mult *= multiplier

    memory_amount = int(val * mult)

    # Check minimum memory requirement
    if memory_amount < min_mem_allowed:
        raise KeywordValidationError(
            """set_memory(): Requested {:.3} MiB ({:.3} MB); minimum 250 MiB (263 MB). Please, sir, I want some more.""".format(
                memory_amount / 1024 ** 2, memory_amount / 1000 ** 2
            )
        )

    return memory_amount


def parse_memory_nomin(inputval):
    return parse_memory(inputval, min_mem_allowed=0)
