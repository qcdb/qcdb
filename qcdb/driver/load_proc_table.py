from qcelemental.util import which, which_import

from qcengine.programs import register_program

from ..programs.cfour import QcdbCFOURHarness, cfour_gradient_list, cfour_hessian_list, cfour_list, run_cfour
from ..programs.dftd3 import alt_run_dftd3, dftd3_list
from ..programs.gamess import QcdbGAMESSHarness, gamess_gradient_list, gamess_hessian_list, gamess_list, run_gamess
from ..programs.nwchem import QcdbNWChemHarness, nwchem_gradient_list, nwchem_hessian_list, nwchem_list, run_nwchem #nwchem_properties_list, run_nwchem
from ..programs.psi4 import QcdbPsi4Harness, run_psi4
from .proc_table import procedures

register_program(QcdbCFOURHarness(name='QCDB-CFOUR'))
register_program(QcdbGAMESSHarness(name='QCDB-GAMESS'))
register_program(QcdbPsi4Harness(name='QCDB-Psi4'))
register_program(QcdbNWChemHarness(name='QCDB-NWChem'))

if which('psi4') and which_import('psi4'):
    import psi4
    print('loaded psi4 from:', psi4.__file__)

    # integrate Psi4 with driver routines
    procedures['energy']['psi4'] = {}
    for mtd in psi4.driver.proc_table.procedures['energy']:
        procedures['energy']['psi4'][mtd.lower()] = run_psi4
        procedures['energy']['psi4']['p4-' + mtd.lower()] = run_psi4
    procedures['energy']['psi4']['p4-mrccsd(t)'] = run_psi4
    procedures['energy']['psi4']['p4-mrccsdtq'] = run_psi4

    procedures['properties']['psi4'] = {}
    for mtd in psi4.driver.proc_table.procedures['properties']:
        procedures['properties']['psi4'][mtd.lower()] = run_psi4
        procedures['properties']['psi4']['p4-' + mtd.lower()] = run_psi4

    procedures['gradient']['psi4'] = {}
    for mtd in psi4.driver.proc_table.procedures['gradient']:
        procedures['gradient']['psi4'][mtd.lower()] = run_psi4
        procedures['gradient']['psi4']['p4-' + mtd.lower()] = run_psi4

    procedures['hessian']['psi4'] = {}
    for mtd in psi4.driver.proc_table.procedures['hessian']:
        procedures['hessian']['psi4'][mtd.lower()] = run_psi4
        procedures['hessian']['psi4']['p4-' + mtd.lower()] = run_psi4

if which('xcfour'):

    # integrate CFOUR with driver routines
    procedures['energy']['cfour'] = {}
    for mtd in cfour_list():
        procedures['energy']['cfour'][mtd.lower()] = run_cfour

    procedures['gradient']['cfour'] = {}
    for mtd in cfour_gradient_list():
        procedures['gradient']['cfour'][mtd.lower()] = run_cfour

    procedures['hessian']['cfour'] = {}
    for mtd in cfour_hessian_list():
        procedures['hessian']['cfour'][mtd.lower()] = run_cfour

if which('dftd3'):

    # integrate DFTD3 with driver routines
    procedures['energy']['dftd3'] = {}
    for mtd in dftd3_list():
        procedures['energy']['dftd3'][mtd.lower()] = alt_run_dftd3

if which('nwchem'):

    # integrate NWChem with driver routines
    procedures['energy']['nwchem'] = {}
    for mtd in nwchem_list():
        procedures['energy']['nwchem'][mtd.lower()] = run_nwchem

    procedures['gradient']['nwchem'] = {}
    for mtd in nwchem_gradient_list():
        procedures['gradient']['nwchem'][mtd.lower()] = run_nwchem

    procedures['hessian']['nwchem'] = {}
    for mtd in nwchem_hessian_list():
        procedures['hessian']['nwchem'][mtd.lower()] = run_nwchem

   # procedures['properties']['nwchem'] = {}
   # for mtd in nwchem_properties_list():
   #     procedures['properties']['nwchem'][mtd.lower()] = run_nwchem

if which('rungms'):

    # integrate GAMESS with driver routines
    procedures['energy']['gamess'] = {}
    for mtd in gamess_list():
        procedures['energy']['gamess'][mtd.lower()] = run_gamess

    procedures['gradient']['gamess'] = {}
    for mtd in gamess_gradient_list():
        procedures['gradient']['gamess'][mtd.lower()] = run_gamess

    procedures['hessian']['gamess'] = {}
    for mtd in gamess_hessian_list():
        procedures['hessian']['gamess'][mtd.lower()] = run_gamess
