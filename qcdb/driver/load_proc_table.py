from qcelemental.util import which, which_import

from qcengine.programs import register_program

from ..programs.cfour import QcdbCFOURHarness, cfour_gradient_list, cfour_hessian_list, cfour_list, run_cfour
from ..programs.dftd3 import alt_run_dftd3, dftd3_list
from ..programs.gamess import QcdbGAMESSHarness, gamess_gradient_list, gamess_hessian_list, gamess_list, run_gamess
from ..programs.nwchem import (
    QcdbNWChemHarness,
    nwchem_gradient_list,
    nwchem_hessian_list,
    nwchem_list,
    run_nwchem,
)  # nwchem_properties_list, run_nwchem
from ..programs.psi4 import QcdbPsi4Harness, run_psi4
from .proc_table import procedures

register_program(QcdbCFOURHarness(name="QCDB-CFOUR"))
register_program(QcdbGAMESSHarness(name="QCDB-GAMESS"))
register_program(QcdbPsi4Harness(name="QCDB-Psi4"))
register_program(QcdbNWChemHarness(name="QCDB-NWChem"))

procedures["energy"]["psi4"] = {}
procedures["properties"]["psi4"] = {}
procedures["gradient"]["psi4"] = {}
procedures["hessian"]["psi4"] = {}

procedures["energy"]["cfour"] = {}
procedures["gradient"]["cfour"] = {}
procedures["hessian"]["cfour"] = {}

procedures["energy"]["nwchem"] = {}
procedures["gradient"]["nwchem"] = {}
procedures["hessian"]["nwchem"] = {}

procedures["energy"]["gamess"] = {}
procedures["gradient"]["gamess"] = {}
procedures["hessian"]["gamess"] = {}

if which("psi4") and which_import("psi4"):
    import psi4

    print("loaded psi4 from:", psi4.__file__)

    # integrate Psi4 with driver routines
    for mtd in psi4.driver.proc_table.procedures["energy"]:
        procedures["energy"]["psi4"][mtd.lower()] = run_psi4
        procedures["energy"]["psi4"]["p4-" + mtd.lower()] = run_psi4
    procedures["energy"]["psi4"]["p4-mrccsd(t)"] = run_psi4
    procedures["energy"]["psi4"]["p4-mrccsdtq"] = run_psi4

    for mtd in psi4.driver.proc_table.procedures["properties"]:
        procedures["properties"]["psi4"][mtd.lower()] = run_psi4
        procedures["properties"]["psi4"]["p4-" + mtd.lower()] = run_psi4

    for mtd in psi4.driver.proc_table.procedures["gradient"]:
        procedures["gradient"]["psi4"][mtd.lower()] = run_psi4
        procedures["gradient"]["psi4"]["p4-" + mtd.lower()] = run_psi4
    procedures["gradient"]["psi4"]["mrccsdt-1a"] = run_psi4  # fd, debug
    procedures["gradient"]["psi4"]["mrccsdt-1b"] = run_psi4  # fd, debug
    procedures["gradient"]["psi4"]["mrccsdt-2"] = run_psi4  # fd, debug
    procedures["gradient"]["psi4"]["mrccsdt-3"] = run_psi4  # fd, debug

    for mtd in psi4.driver.proc_table.procedures["hessian"]:
        procedures["hessian"]["psi4"][mtd.lower()] = run_psi4
        procedures["hessian"]["psi4"]["p4-" + mtd.lower()] = run_psi4
    procedures["hessian"]["psi4"]["p4-mp2"] = run_psi4  # fd
    procedures["hessian"]["psi4"]["p4-ccsd"] = run_psi4  # fd
    procedures["hessian"]["psi4"]["p4-ccsd(t)"] = run_psi4  # fd
    procedures["hessian"]["psi4"]["p4-ccsdt"] = run_psi4  # fd

if which("xcfour"):

    # integrate CFOUR with driver routines
    for mtd in cfour_list():
        procedures["energy"]["cfour"][mtd.lower()] = run_cfour

    for mtd in cfour_gradient_list():
        procedures["gradient"]["cfour"][mtd.lower()] = run_cfour

    for mtd in cfour_hessian_list():
        procedures["hessian"]["cfour"][mtd.lower()] = run_cfour

if which("dftd3"):

    # integrate DFTD3 with driver routines
    procedures["energy"]["dftd3"] = {}
    for mtd in dftd3_list():
        procedures["energy"]["dftd3"][mtd.lower()] = alt_run_dftd3

if which("nwchem"):

    # integrate NWChem with driver routines
    for mtd in nwchem_list():
        procedures["energy"]["nwchem"][mtd.lower()] = run_nwchem

    for mtd in nwchem_gradient_list():
        procedures["gradient"]["nwchem"][mtd.lower()] = run_nwchem

    for mtd in nwchem_hessian_list():
        procedures["hessian"]["nwchem"][mtd.lower()] = run_nwchem

# procedures['properties']['nwchem'] = {}
# for mtd in nwchem_properties_list():
#     procedures['properties']['nwchem'][mtd.lower()] = run_nwchem

if which("rungms"):

    # integrate GAMESS with driver routines
    for mtd in gamess_list():
        procedures["energy"]["gamess"][mtd.lower()] = run_gamess

    for mtd in gamess_gradient_list():
        procedures["gradient"]["gamess"][mtd.lower()] = run_gamess

    for mtd in gamess_hessian_list():
        procedures["hessian"]["gamess"][mtd.lower()] = run_gamess

# import pprint
# pprint.pprint(procedures, width=200)
