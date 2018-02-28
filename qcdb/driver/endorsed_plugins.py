from .proc_table import procedures
from .. import iface_psi4
from .. import iface_cfour

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
    
    # Integrate Psi4 with driver routines
    procedures['energy']['psi4'] = {}
    for mtd in psi4.driver.proc_table.procedures['energy']:
        procedures['energy']['psi4'][mtd.lower()] = iface_psi4.run_psi4
    
if True:
    
    # Integrate CFOUR with driver routines
    procedures['energy']['cfour'] = {}
    for mtd in iface_cfour.cfour_list():
        procedures['energy']['cfour'][mtd.lower()] = iface_cfour.run_cfour

    #procedures['gradient']['cfour'] = {}
    #for mtd in iface_cfour.cfour_gradient_list():
    #   procedures['gradient']['cfour'][mtd.lower()] = iface_cfour.run_cfour
