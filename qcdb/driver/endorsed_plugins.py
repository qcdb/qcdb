from .proc_table import procedures
from .. import iface_psi4
from .. import iface_cfour
from .. import iface_dftd3
from .. import iface_nwchem
from .. import intf_gamess

#from .util import pl
#
#try:
#    import psi4
#except ImportError:
#    pass
#
#if addons('psi4'):

#def load_up_psi():
if True:
    import psi4
    print('loaded psi4 from:', psi4.__file__)

    # integrate Psi4 with driver routines
    procedures['energy']['psi4'] = {}
    for mtd in psi4.driver.proc_table.procedures['energy']:
        procedures['energy']['psi4'][mtd.lower()] = iface_psi4.run_psi4
        procedures['energy']['psi4']['p4-' + mtd.lower()] = iface_psi4.run_psi4

    procedures['properties']['psi4'] = {}
    for mtd in psi4.driver.proc_table.procedures['properties']:
        procedures['properties']['psi4'][mtd.lower()] = iface_psi4.run_psi4

    procedures['gradient']['psi4'] = {}
    for mtd in psi4.driver.proc_table.procedures['gradient']:
        procedures['gradient']['psi4'][mtd.lower()] = iface_psi4.run_psi4

    procedures['hessian']['psi4'] = {}
    for mtd in psi4.driver.proc_table.procedures['hessian']:
        procedures['hessian']['psi4'][mtd.lower()] = iface_psi4.run_psi4

if True:

    # integrate CFOUR with driver routines
    procedures['energy']['cfour'] = {}
    for mtd in iface_cfour.cfour_list():
        procedures['energy']['cfour'][mtd.lower()] = iface_cfour.run_cfour

    procedures['gradient']['cfour'] = {}
    for mtd in iface_cfour.cfour_gradient_list():
       procedures['gradient']['cfour'][mtd.lower()] = iface_cfour.run_cfour

    procedures['hessian']['cfour'] = {}
    for mtd in iface_cfour.cfour_hessian_list():
       procedures['hessian']['cfour'][mtd.lower()] = iface_cfour.run_cfour

if True:

    # integrate DFTD3 with driver routines
    procedures['energy']['dftd3'] = {}
    for mtd in iface_dftd3.dftd3_list():
        procedures['energy']['dftd3'][mtd.lower()] = iface_dftd3.alt_run_dftd3

if True:

    # integrate NWChem with driver routines
    procedures['energy']['nwchem'] = {}
    for mtd in iface_nwchem.nwchem_list():
        procedures['energy']['nwchem'][mtd.lower()] = iface_nwchem.run_nwchem

    procedures['gradient']['nwchem'] = {}
    for mtd in iface_nwchem.nwchem_gradient_list():
       procedures['gradient']['nwchem'][mtd.lower()] = iface_nwchem.run_nwchem

if True:

    # integrate GAMESS with driver routines
    procedures['energy']['gamess'] = {}
    for mtd in intf_gamess.gamess_list():
        procedures['energy']['gamess'][mtd.lower()] = intf_gamess.run_gamess
