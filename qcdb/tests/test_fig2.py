from collections import defaultdict
import pprint

import pytest

import qcelemental as qcel
import qcengine as qcng
import qcdb

from .utils import *

_abbr = {v: k for k, v in qcdb.driver.driver_util.pkgprefix.items()}

_the_fc_energy = -55.74674563
_the_energy = -55.74896542928988
_the_short_energy = str(_the_energy)[:10]


def _qcskmol_nh2():
    nh2 = qcel.models.Molecule.from_data(
        """
0 2
N  0.000  0.000 -0.146
H  0.000 -1.511  1.014
H  0.000  1.511  1.014
units au
"""
    )

    return nh2


_dnh2 = {
    "geometry": [0.000, 0.000, -0.146, 0.000, -1.511, 1.014, 0.000, 1.511, 1.014],
    "symbols": ["N", "H", "H"],
    "molecular_charge": 0,
    "molecular_multiplicity": 2,
}

_ins = defaultdict(lambda: defaultdict(dict))


_ins["cfour"][
    "a"
] = """
comment
N  0.000  0.000 -0.146
H  0.000 -1.511  1.014
H  0.000  1.511  1.014

*CFOUR(REFERENCE=ROHF,BASIS=AUG-PVDZ
CALC_LEVEL=CCSD,CHARGE=0,UNITS=BOHR
COORDINATES=CARTESIAN,MULTIPLICITY=2)
"""

# memory added
_ins["gamess"][
    "a"
] = """
 $ccinp ncore=0 $end
 $basis gbasis=accd $end
 $contrl cctyp=ccsd coord=prinaxis
  icharg=0 ispher=1 mult=2
  runtyp=energy scftyp=rohf
  units=bohr $end
 $data

 C1
 N 7  0.000  0.000 -0.146
 H 1  0.000 -1.511  1.014
 H 1  0.000  1.511  1.014
 $end
 $system mwords=2000 $end
"""

_ins["nwchem"][
    "a"
] = """
geometry units bohr
N  0.000  0.000 -0.146
H  0.000 -1.511  1.014
H  0.000  1.511  1.014
end
charge 0
basis spherical
  h library aug-cc-pvdz
  n library aug-cc-pvdz
end
scf
  rohf
  nopen 1
end
tce
  ccsd
end
task tce energy
"""

_ins["psi4"][
    "a"
] = """
molecule {
0 2
N  0.000  0.000 -0.146
H  0.000 -1.511  1.014
H  0.000  1.511  1.014
units au
}

set reference rohf

energy('ccsd/aug-cc-pvdz')
"""

#   (A)(A)
#
#   cp infile . && xcfour        > outfile
#                  rungms infile > outfile
#                  nwchem infile > outfile
#                  psi4   infile   outfile


#   (B)(B)
#
#                         'cfour'
# ene = qcng.compute(json, 'gamess')
#                         'nwchem'
#                         'psi4'

#   (C)(C)
#   (D)(D)
#
#                   c4
# ene = qcdb.energy('gms-ccsd/aug-cc-pvdz')
#                   nwc
#                   p4


_ins["cfour"]["b"] = {
    "molecule": _dnh2,
    "driver": "energy",
    "model": {"method": "ccsd", "basis": "aug-pvdz"},
    "keywords": {"reference": "rohf"},
}

_ins["gamess"]["b"] = {
    "molecule": _dnh2,
    "driver": "energy",
    "model": {"method": "ccsd", "basis": "accd"},
    "keywords": {"contrl__ispher": 1, "contrl__scftyp": "rohf", "ccinp__ncore": 0},
}

_ins["nwchem"]["b"] = {
    "molecule": _dnh2,
    "driver": "energy",
    "model": {"method": "ccsd", "basis": "aug-cc-pvdz"},
    "keywords": {"basis__spherical": True, "scf__rohf": True, "qc_module": "tce"},
}

_ins["psi4"]["b"] = {
    "molecule": _dnh2,
    "driver": "energy",
    "model": {"method": "ccsd", "basis": "aug-cc-pvdz"},
    "keywords": {"reference": "rohf"},
}

_ins["cfour"]["c"] = {
    "molecule": _dnh2,
    "driver": "energy",
    "model": {"method": "ccsd", "basis": "aug-cc-pvdz"},
    "keywords": {"cfour_reference": "rohf"},
}

_ins["gamess"]["c"] = {
    "molecule": _dnh2,
    "driver": "energy",
    "model": {"method": "ccsd", "basis": "aug-cc-pvdz"},
    "keywords": {"gamess_contrl__scftyp": "rohf", "gamess_ccinp__ncore": 0},
}

_ins["nwchem"]["c"] = {
    "molecule": _dnh2,
    "driver": "energy",
    "model": {"method": "ccsd", "basis": "aug-cc-pvdz"},
    "keywords": {"nwchem_scf__rohf": True, "qc_module": "tce"},
}

_ins["psi4"]["c"] = {
    "molecule": _dnh2,
    "driver": "energy",
    "model": {"method": "ccsd", "basis": "aug-cc-pvdz"},
    "keywords": {"psi4_reference": "rohf"},
}

_ins["all"]["d"] = {
    "molecule": _dnh2,
    "driver": "energy",
    "model": {"method": "ccsd", "basis": "aug-cc-pvdz"},
    "keywords": {"reference": "rohf", "freeze_core": False},
}


@pytest.mark.parametrize("qcprog", ["cfour", "gamess", "nwchem", "psi4"])
def test_fig2a_txt(qcprog, request):
    sin = _ins[qcprog]["a"]
    from qcelemental.util import which
    from qcengine.util import execute

    if qcprog == "cfour":
        success, dexe = execute([which("xcfour")], {"ZMAT": sin})
    if qcprog == "gamess":
        success, dexe = execute([which("rungms"), "input"], {"input.inp": sin})
    if qcprog == "nwchem":
        success, dexe = execute([which("nwchem"), "input"], {"input": sin})
    if qcprog == "psi4":
        success, dexe = execute([which("psi4"), "input", "-o", "stdout"], {"input": sin})

    pprint.pprint(dexe, width=180)
    assert success
    assert _the_short_energy in dexe["stdout"]


@pytest.mark.parametrize("qcprog", ["cfour", "gamess", "nwchem", "psi4"])
def test_fig2b_json(qcprog, request):
    datin = _ins[qcprog]["b"]

    atres = qcng.compute(datin, qcprog)
    pprint.pprint(atres.dict(), width=180)
    assert compare_values(_the_energy, atres.return_result, atol=1.0e-6, label=request.node.name)


@pytest.mark.parametrize("qcprog", ["cfour", "gamess", "nwchem", "psi4"])
def test_fig2c_api(qcprog, request):
    datin = _ins[qcprog]["c"]

    qcskmol = qcel.models.Molecule(**datin["molecule"])
    qcdbmol = qcdb.Molecule.from_schema(qcskmol.dict())
    qcdb.set_keywords(datin["keywords"])
    driver = {
        "energy": qcdb.energy,
        "gradient": qcdb.gradient,
    }[datin["driver"]]
    call = _abbr[qcprog] + datin["model"]["method"] + "/" + datin["model"]["basis"]

    ene = driver(call, molecule=qcdbmol)
    assert compare_values(_the_energy, ene, atol=1.0e-6, label=request.node.name)


@pytest.mark.parametrize("qcprog", ["cfour", "gamess", "nwchem", "psi4"])
def wanted_test_fig2c_json(qcprog, request):
    datin = _ins[qcprog]["c"]

    atres = qcdb.compute(datin, qcprog)
    assert compare_values(_the_energy, atres.return_result, atol=1.0e-6, label=request.node.name)


@pytest.mark.parametrize("qcprog", ["cfour", "gamess", "nwchem", "psi4"])
def test_fig2d_api(qcprog, request):
    datin = _ins["all"]["d"]

    qcskmol = qcel.models.Molecule(**datin["molecule"])
    qcdbmol = qcdb.Molecule.from_schema(qcskmol.dict())
    qcdb.set_keywords(datin["keywords"])
    driver = {
        "energy": qcdb.energy,
        "gradient": qcdb.gradient,
    }[datin["driver"]]
    call = _abbr[qcprog] + datin["model"]["method"] + "/" + datin["model"]["basis"]

    ene = driver(call, molecule=qcdbmol)
    assert compare_values(_the_energy, ene, atol=1.0e-6, label=request.node.name)
