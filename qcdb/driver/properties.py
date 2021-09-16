"""Module with a *procedures* dictionary specifying available quantum
chemical methods and functions driving the main quantum chemical
functionality, namely properties calculations

"""
# import copy
import pprint

from ..keywords import register_kwds
from . import pe  # keep this at top of imports
from . import driver_helpers, driver_util
from .proc_table import procedures

pp = pprint.PrettyPrinter(width=120)


@register_kwds(pe.nu_options)
def properties(name, **kwargs):
    r"""Function to compute the single-point electronic properties."""

    from . import load_proc_table

    kwargs = driver_util.kwargs_lower(kwargs)

    if "options" in kwargs:
        driver_helpers.set_options(kwargs.pop("options"))

    # Bounce if name is function
    if hasattr(name, "__call__"):
        return name(properties, kwargs.pop("label", "custom function"), ptype="properties", **kwargs)

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

    # return_wfn = kwargs.pop('return_wfn', False)
    pe.active_qcvars = {}

    package = driver_util.get_package(lowername, kwargs)

    # jobrec = procedures['properties'][package][lowername](lowername, molecule=molecule, options=pe.nu_options, ptype='properties', **kwargs)
