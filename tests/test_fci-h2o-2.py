from addons import using_psi4, using_gamess
from utils import *

import qcdb


_refnuc   =   9.2342185209120
_refscf   = -75.985323665263
_refci    = -76.1210978474779
_refcorr  = _refci - _refscf


def system_water():
    h2o = qcdb.set_molecule("""
   O       .0000000000         .0000000000        -.0742719254
   H       .0000000000       -1.4949589982       -1.0728640373
   H       .0000000000        1.4949589982       -1.0728640373
units bohr
    """)
    
    h2o.update_geometry()
    assert compare_values(_refnuc, h2o.nuclear_repulsion_energy(), 9, "Nuclear repulsion energy")

    return h2o


@using_psi4
def test_fci_rhf_psi4():
    #! 6-31G H2O Test FCI Energy Point

    h2o = system_water()
    qcdb.set_options({
        'basis': '6-31G',
        #'psi4_detci__icore': 0,
    })

    E = qcdb.energy('p4-fci', molecule=h2o)

    assert compare_values(_refnuc, h2o.nuclear_repulsion_energy(), 9, "nre")
    assert compare_values(_refscf, qcdb.get_variable("HF total energy"), 8, "hf total energy")
    assert compare_values(_refci, E, 7, "return E")
    assert compare_values(_refci, qcdb.get_variable("FCI TOTAL ENERGY"), 7, "fci total energy")
    assert compare_values(_refcorr, qcdb.get_variable("FCI CORRELATION ENERGY"), 7, "fci correlation energy")
    assert compare_values(_refci, qcdb.get_variable("CI TOTAL ENERGY"), 7, "ci total energy")
    assert compare_values(_refcorr, qcdb.get_variable("CI CORRELATION ENERGY"), 7, "ci correlation energy")


@using_gamess
def test_fci_rhf_gamess():
    #! 6-31G H2O Test FCI Energy Point

    h2o = system_water()
    qcdb.set_options({
        'basis': '6-31G',
        'gamess_cidet__ncore': 0,
    })

    E = qcdb.energy('gms-fci', molecule=h2o)

    assert compare_values(_refnuc, h2o.nuclear_repulsion_energy(), 9, "nre")
    assert compare_values(_refscf, qcdb.get_variable("HF total energy"), 8, "hf total energy")
    assert compare_values(_refci, E, 7, "return E")
    assert compare_values(_refci, qcdb.get_variable("FCI TOTAL ENERGY"), 7, "fci total energy")
    assert compare_values(_refcorr, qcdb.get_variable("FCI CORRELATION ENERGY"), 7, "fci correlation energy")
    assert compare_values(_refci, qcdb.get_variable("CI TOTAL ENERGY"), 7, "ci total energy")
    assert compare_values(_refcorr, qcdb.get_variable("CI CORRELATION ENERGY"), 7, "ci correlation energy")
