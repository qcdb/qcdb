from qcelemental.util import which, which_import

from .proc_table import procedures
from .. import intf_psi4
from .. import intf_cfour
from .. import intf_dftd3
from .. import intf_nwchem
from .. import intf_gamess

from qcengine.programs import register_program
register_program(intf_cfour.QcdbCFOURHarness(name='QCDB-CFOUR'))
register_program(intf_gamess.QcdbGAMESSHarness(name='QCDB-GAMESS'))
register_program(intf_psi4.QcdbPsi4Harness(name='QCDB-Psi4'))
register_program(intf_nwchem.QcdbNWChemHarness(name='QCDB-NWChem'))

if which('psi4') and which_import('psi4'):
    import psi4
    print('loaded psi4 from:', psi4.__file__)

    # integrate Psi4 with driver routines
    procedures['energy']['psi4'] = {}
    for mtd in psi4.driver.proc_table.procedures['energy']:
        procedures['energy']['psi4'][mtd.lower()] = intf_psi4.run_psi4
        procedures['energy']['psi4']['p4-' + mtd.lower()] = intf_psi4.run_psi4
    procedures['energy']['psi4']['p4-mrccsd(t)'] = intf_psi4.run_psi4
    procedures['energy']['psi4']['p4-mrccsdtq'] = intf_psi4.run_psi4

    procedures['properties']['psi4'] = {}
    for mtd in psi4.driver.proc_table.procedures['properties']:
        procedures['properties']['psi4'][mtd.lower()] = intf_psi4.run_psi4
        procedures['properties']['psi4']['p4-' + mtd.lower()] = intf_psi4.run_psi4

    procedures['gradient']['psi4'] = {}
    for mtd in psi4.driver.proc_table.procedures['gradient']:
        procedures['gradient']['psi4'][mtd.lower()] = intf_psi4.run_psi4
        procedures['gradient']['psi4']['p4-' + mtd.lower()] = intf_psi4.run_psi4

    procedures['hessian']['psi4'] = {}
    for mtd in psi4.driver.proc_table.procedures['hessian']:
        procedures['hessian']['psi4'][mtd.lower()] = intf_psi4.run_psi4
        procedures['hessian']['psi4']['p4-' + mtd.lower()] = intf_psi4.run_psi4

if which('xcfour'):

    # integrate CFOUR with driver routines
    procedures['energy']['cfour'] = {}
    for mtd in intf_cfour.cfour_list():
        procedures['energy']['cfour'][mtd.lower()] = intf_cfour.run_cfour

    procedures['gradient']['cfour'] = {}
    for mtd in intf_cfour.cfour_gradient_list():
       procedures['gradient']['cfour'][mtd.lower()] = intf_cfour.run_cfour

    procedures['hessian']['cfour'] = {}
    for mtd in intf_cfour.cfour_hessian_list():
       procedures['hessian']['cfour'][mtd.lower()] = intf_cfour.run_cfour

if which('dftd3'):

    # integrate DFTD3 with driver routines
    procedures['energy']['dftd3'] = {}
    for mtd in intf_dftd3.dftd3_list():
        procedures['energy']['dftd3'][mtd.lower()] = intf_dftd3.alt_run_dftd3

if which('nwchem'):

    # integrate NWChem with driver routines
    procedures['energy']['nwchem'] = {}
    for mtd in intf_nwchem.nwchem_list():
        procedures['energy']['nwchem'][mtd.lower()] = intf_nwchem.run_nwchem

    procedures['gradient']['nwchem'] = {}
    for mtd in intf_nwchem.nwchem_gradient_list():
       procedures['gradient']['nwchem'][mtd.lower()] = intf_nwchem.run_nwchem

    procedures['hessian']['nwchem'] = {}
    for mtd in intf_nwchem.nwchem_hessian_list():
        procedures['hessian']['nwchem'][mtd.lower()] = intf_nwchem.run_nwchem

if which('rungms'):

    # integrate GAMESS with driver routines
    procedures['energy']['gamess'] = {}
    for mtd in intf_gamess.gamess_list():
        procedures['energy']['gamess'][mtd.lower()] = intf_gamess.run_gamess

    # integrate GAMESS with driver routines
    procedures['gradient']['gamess'] = {}
    for mtd in intf_gamess.gamess_gradient_list():
        procedures['gradient']['gamess'][mtd.lower()] = intf_gamess.run_gamess
