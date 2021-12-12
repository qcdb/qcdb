import pprint
from collections import defaultdict

import pytest
import qcelemental as qcel
import qcengine as qcng
from qcengine.testing import using

import qcdb

from .utils import *

pp = pprint.PrettyPrinter(width=120)


_the_fig2_fc_energy = -55.74674563
_the_fig2_energy = -55.74896542928988
_the_short_fig2_energy = str(_the_fig2_energy)[:10]


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

_fig2_ins = defaultdict(lambda: defaultdict(dict))


_fig2_ins["cfour"][
    "a"
] = """
comment
N  0.000  0.000 -0.146
H  0.000 -1.511  1.014
H  0.000  1.511  1.014

*CFOUR(REFERENCE=ROHF
CALC_LEVEL=CCSD,BASIS=AUG-PVDZ
CHARGE=0,MULTIPLICITY=2
COORDINATES=CARTESIAN,UNITS=BOHR)
"""

# memory added
_fig2_ins["gamess"][
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

_fig2_ins["nwchem"][
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

_fig2_ins["psi4"][
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


_fig2_ins["cfour"]["b"] = {
    "molecule": _dnh2,
    "driver": "energy",
    "model": {"method": "ccsd", "basis": "aug-pvdz"},
    "keywords": {"reference": "rohf"},
}

_fig2_ins["gamess"]["b"] = {
    "molecule": _dnh2,
    "driver": "energy",
    "model": {"method": "ccsd", "basis": "accd"},
    "keywords": {"contrl__ispher": 1, "contrl__scftyp": "rohf", "ccinp__ncore": 0},
}

_fig2_ins["nwchem"]["b"] = {
    "molecule": _dnh2,
    "driver": "energy",
    "model": {"method": "ccsd", "basis": "aug-cc-pvdz"},
    "keywords": {"basis__spherical": True, "scf__rohf": True, "qc_module": "tce"},
}

_fig2_ins["psi4"]["b"] = {
    "molecule": _dnh2,
    "driver": "energy",
    "model": {"method": "ccsd", "basis": "aug-cc-pvdz"},
    "keywords": {"reference": "rohf"},
}

_fig2_ins["cfour"]["c"] = {
    "molecule": _dnh2,
    "driver": "energy",
    "model": {"method": "ccsd", "basis": "aug-cc-pvdz"},
    "keywords": {"cfour_reference": "rohf"},
}

_fig2_ins["gamess"]["c"] = {
    "molecule": _dnh2,
    "driver": "energy",
    "model": {"method": "ccsd", "basis": "aug-cc-pvdz"},
    "keywords": {"gamess_contrl__scftyp": "rohf", "gamess_ccinp__ncore": 0},
}

_fig2_ins["nwchem"]["c"] = {
    "molecule": _dnh2,
    "driver": "energy",
    "model": {"method": "ccsd", "basis": "aug-cc-pvdz"},
    "keywords": {"nwchem_scf__rohf": True, "qc_module": "tce"},
}

_fig2_ins["psi4"]["c"] = {
    "molecule": _dnh2,
    "driver": "energy",
    "model": {"method": "ccsd", "basis": "aug-cc-pvdz"},
    "keywords": {"psi4_reference": "rohf"},
}

_fig2_ins["all"]["d"] = {
    "molecule": _dnh2,
    "driver": "energy",
    "model": {"method": "ccsd", "basis": "aug-cc-pvdz"},
    "keywords": {"reference": "rohf", "freeze_core": False},
}


@pytest.mark.parametrize(
    "qcprog",
    [
        pytest.param("cfour", marks=using("cfour")),
        pytest.param("gamess", marks=using("gamess")),
        pytest.param("nwchem", marks=using("nwchem")),
        pytest.param("psi4", marks=using("psi4")),
    ],
)
def test_fig2a_txt(qcprog):
    """
    Execution (a; txt)

    cp infile . && xcfour        > outfile
                   rungms infile > outfile
                   nwchem infile > outfile
                   psi4   infile   outfile
    """

    sin = _fig2_ins[qcprog]["a"]
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
    assert _the_short_fig2_energy in dexe["stdout"]


@pytest.mark.parametrize(
    "qcprog",
    [
        pytest.param("cfour", marks=using("cfour")),
        pytest.param("gamess", marks=using("gamess")),
        pytest.param("nwchem", marks=using("nwchem")),
        pytest.param("psi4", marks=using("psi4")),
    ],
)
def test_fig2b_json(qcprog, request):
    """
    Execution (b; json)

                             'cfour'
    ene = qcng.compute(json, 'gamess')
                             'nwchem'
                             'psi4'

    """
    datin = _fig2_ins[qcprog]["b"]

    atres = qcng.compute(datin, qcprog)
    pprint.pprint(atres.dict(), width=180)
    assert compare_values(_the_fig2_energy, atres.return_result, atol=1.0e-6, label=request.node.name)


@pytest.mark.parametrize(
    "qcprog",
    [
        pytest.param("cfour", marks=using("cfour")),
        pytest.param("gamess", marks=using("gamess")),
        pytest.param("nwchem", marks=using("nwchem")),
        pytest.param("psi4", marks=using("psi4")),
    ],
)
@pytest.mark.parametrize("part", ["c", "d"])
def test_fig2cd_api(qcprog, part, request):
    """
    Execution (c, d; api)

                       c4
    ene = qcdb.energy('gms-ccsd/aug-cc-pvdz')
                       nwc
                       p4

    """
    if part == "c":
        datin = _fig2_ins[qcprog]["c"]
        mo = {"mode": "sandwich"}
    elif part == "d":
        datin = _fig2_ins["all"]["d"]
        mo = {"mode": "unified"}

    qcskmol = qcel.models.Molecule(**datin["molecule"])
    qcdbmol = qcdb.Molecule.from_schema(qcskmol.dict())
    qcdb.set_keywords(datin["keywords"])
    driver = {
        "energy": qcdb.energy,
        "gradient": qcdb.gradient,
    }[datin["driver"]]
    call = qcdb.util.program_prefix(qcprog) + datin["model"]["method"] + "/" + datin["model"]["basis"]

    ene = driver(call, molecule=qcdbmol, mode_options=mo)
    assert compare_values(_the_fig2_energy, ene, atol=1.0e-6, label=request.node.name)


@pytest.mark.parametrize(
    "qcprog",
    [
        pytest.param("cfour", marks=using("cfour")),
        pytest.param("gamess", marks=using("gamess")),
        pytest.param("nwchem", marks=using("nwchem")),
        pytest.param("psi4", marks=using("psi4")),
    ],
)
@pytest.mark.parametrize("part", ["c", "d"])
def test_fig2cd_json(qcprog, part, request):
    """
    Execution (c, d; json)

                             'cfour'
    ene = qcdb.compute(json, 'gamess')
                             'nwchem'
                             'psi4'

    """
    if part == "c":
        datin = _fig2_ins[qcprog]["c"]
        mo = None
    elif part == "d":
        datin = _fig2_ins["all"]["d"]
        mo = {"module_fallback": True}

    atres = qcdb.compute(datin, qcprog, mode_options=mo)
    pprint.pprint(atres.dict(), width=200)

    assert compare_values(_the_fig2_energy, atres.return_result, atol=1.0e-6, label=request.node.name)


def test_snippet2():
    snippet1 = """
O 0 0 0
H 2 0 0
--
@22Ne 5 0 0
units bohr
"""

    snippet2 = {
        "atom_labels": ["", "", ""],
        "atomic_numbers": [8, 1, 10],
        "fix_com": False,
        "fix_orientation": False,
        "fragment_charges": [0.0, 0.0],
        "fragment_multiplicities": [2, 1],
        "fragments": [[0, 1], [2]],
        "geometry": [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [5.0, 0.0, 0.0]],
        "mass_numbers": [16, 1, 22],
        "masses": [15.99491462, 1.00782503, 21.99138511],
        "molecular_charge": 0.0,
        "molecular_multiplicity": 2,
        "name": "HNeO",
        "provenance": {"creator": "QCElemental", "routine": "qcelemental.molparse.from_schema", "version": "v0.8.0"},
        "real": [True, True, False],
        "schema_name": "qcschema_molecule",
        "schema_version": 2,
        "symbols": ["O", "H", "Ne"],
        "validated": True,
    }

    mol1 = qcel.models.Molecule.from_data(snippet1)

    # the `qcel.models.Molecule.dict()` method removes default and regeneratable fields for compact
    #   storage in database, so test those separately
    assert compare(snippet2.pop("atom_labels"), mol1.atom_labels, "atom_labels")
    assert compare(snippet2.pop("atomic_numbers"), mol1.atomic_numbers, "atomic_numbers")
    assert compare_recursive(snippet2, mol1.dict(), "snippet1 -> snippet2", forgive=["provenance"])


@using("psi4")
@using("cfour")
@pytest.mark.parametrize(
    "bsse,ans",
    [
        pytest.param(None, -229.34491301195257, id="electronic_ene"),
        pytest.param("cp", -0.00017499, id="counterpoise_interaction_ene"),
    ],
)
def test_snippet3(bsse, ans):
    import qcdb

    nefh = qcdb.set_molecule(
        """Ne
           --
           F 1 R
           H 2 1.0 1 135.0"""
    )
    qcdb.set_options({"e_convergence": 7, "mp2_type": "df"})
    results = {r / 100: None for r in range(200, 400, 10)}
    for intra in results:
        nefh.R = intra
        results[intra] = qcdb.energy("p4-mp2/jun-cc-pvtz")
    rmin = min(results, key=results.get)
    qcdb.set_options({"e_convergence": 9})
    nefh.R = rmin
    model = "p4-mp2/aug-cc-pv[tq]z + d:c4-ccsd(t)/aug-cc-pvtz"
    ene = qcdb.energy(model, bsse_type=bsse)
    print(f"Ne...FH at optimal dist. {rmin} A has IE {ene} E_h.")

    assert compare_values(3.2, rmin, atol=1.0e-6)
    assert compare_values(ans, ene, atol=1.0e-6)


@using("cfour")
def test_snippet4a():
    qcdb.set_molecule(
        """
     H
     H 1 0.74
    """
    )

    qcdb.set_keywords(
        {
            "basis": "6-31g",  # ok, new info
            "cfour_calc_level": "ccsd",  # clash w/"c4-hf" below
        }
    )

    with pytest.raises(qcdb.exceptions.KeywordReconciliationError) as e:
        qcdb.energy("c4-hf")

    assert "Conflicting option requirements btwn user (CCSD) and driver (SCF) for CALC_LEVEL" in str(e)


@using("cfour")
def test_snippet4b():
    qcdb.set_molecule(
        """
     H
     H 1 0.74
    """
    )

    qcdb.set_keywords(
        {
            "basis": "6-31g",  # ok, new info
            "cfour_deriv_level": "first",  # clash w/energy() below (use gradient())
        }
    )

    with pytest.raises(qcdb.exceptions.KeywordReconciliationError) as e:
        qcdb.energy("c4-hf")

    assert "Conflicting option requirements btwn user (FIRST) and driver (ZERO) for DERIV_LEVEL" in str(e)


@using("cfour")
def test_snippet4c():
    qcdb.set_molecule(
        """
     H
     H 1 0.74
    """
    )

    qcdb.set_keywords(
        {
            "basis": "6-31g",  # ok, new info
            "cfour_multiplicity": 3,  # clash w/implicit singlet of mol above
            "cfour_units": "angstrom",  # ok, consistent w/mol above
        }
    )

    with pytest.raises(qcdb.exceptions.KeywordReconciliationError) as e:
        qcdb.energy("c4-hf")

    assert "Conflicting option requirements btwn user (3) and driver (1) for MULTIPLICITY" in str(e)


@using("cfour")
def test_snippet4d():
    qcdb.set_molecule(
        """
     H
     H 1 0.74
    """
    )

    qcdb.set_keywords(
        {
            "basis": "6-31g",  # ok, new info
            "cfour_memory_size": 9000000,  # clash w/1 gib below
        }
    )

    with pytest.raises(qcdb.exceptions.KeywordReconciliationError) as e:
        qcdb.energy("c4-hf", local_options={"memory": 1})

    assert "Conflicting option requirements btwn user (9000000) and driver (134217728) for MEMORY_SIZE" in str(e)


@using("cfour")
@using("gamess")
@using("nwchem")
@using("psi4")
def test_sectionIII():
    """Helper Functions"""

    def gen_rvals(center, stepsize, npoints):
        if npoints % 2 == 0:
            return []
        center_ind = npoints // 2
        rvals = [center + stepsize * (ind - center_ind) for ind in range(npoints)]
        return rvals

    """ QCDB Setup """

    bh = qcdb.set_molecule(
        """
            B
            H 1 R
            units angstrom
            """
    )

    local_options = {"memory": 35, "nnodes": 1, "ncores": 1}
    qcdb.set_options(
        {
            "e_convergence": 1e-11,
            "scf__d_convergence": 1e-9,
            "nwchem_ccsd__maxiter": 100,
            "psi4_mp2_type": "conv",
            "psi4_scf_type": "direct",
            "psi4_df_scf_guess": "false",
        }
    )

    """ Geometry Scan """

    Re = 1.229980860371199  # self consistent with fci optimization
    dR = 0.005
    npoints = 5

    R_arr = gen_rvals(Re, dR, npoints)

    # energy of base method
    E_base = [0.0] * npoints

    # various corrections to base energy
    dE_basis = [0.0] * npoints
    dE_dboc = [0.0] * npoints
    dE_x2c = [0.0] * npoints
    dE_fci = [0.0] * npoints
    dE_ccsdtq = [0.0] * npoints

    for i, R in enumerate(R_arr):
        print(f"<<<\n<<< POINT {i}\n<<<")

        bh.R = R
        print(f"~~~ Calculation at Re={R} Angstrom ({i+1}/{npoints}) ~~~")

        # read from file in case of restart
        #  if i < 1:
        #  if True:
        #  if i > 0:
        if False:
            #  if i < 4:
            with open(f"dist_{i}", "r+") as f:
                lines = f.readlines()
                lines = [float(line.split()[3]) for line in lines]
                E_base[i] = lines[0]
                dE_basis[i] = lines[1]
                dE_dboc[i] = lines[2]
                dE_x2c[i] = lines[3]
                dE_fci[i] = lines[4]
                dE_ccsdtq[i] = lines[5]
                print(lines)
            continue

        # ccsdtq correction: (CCSDTQ - CCSD(T)) / cc-pVDZ
        qcdb.set_options(
            {
                "cfour_dropmo": [1],
            }
        )
        _, jrec = qcdb.energy("c4-ccsd(t)/cc-pVTZ", return_wfn=True, local_options=local_options)
        E_ccsdpt = float(jrec["qcvars"]["CCSD(T) TOTAL ENERGY"].data)
        _, jrec = qcdb.energy("c4-ccsdtq/cc-pVTZ", return_wfn=True, local_options=local_options)
        E_ccsdtq = float(jrec["qcvars"]["CCSDTQ TOTAL ENERGY"].data)
        qcdb.set_options({"cfour_dropmo": None})
        dE_ccsdtq[i] = E_ccsdtq - E_ccsdpt
        print(f"~~~ CCSDTQ Correction={dE_ccsdtq[-1]} Har. ({i+1}/{npoints}) ~~~")

        # base calculation: CCSD(T) / cc-pCV[Q5]Z
        #  qcdb.set_options({'memory': '10 gb'})
        # E, jrec = qcdb.energy('nwc-ccsd(t)/cc-pCV[T,Q]Z', return_wfn=True)
        E, jrec = qcdb.energy("nwc-ccsd(t)/cc-pCVTZ", return_wfn=True)
        #  qcdb.set_options({'memory': '55 gb'})
        E_base[i] = E
        print(f"~~~ Base Energy={E} Har. ({i+1}/{npoints}) ~~~")

        # basis set correction: MP2 / (aug-cc-pCV[56]Z) - cc-pCV[Q5]Z)
        E_small, _ = qcdb.energy("p4-mp2/cc-pCV[T,Q]Z", return_wfn=True)
        E_large, _ = qcdb.energy("p4-mp2/aug-cc-pCV[T,Q]Z", return_wfn=True)
        dE_basis[i] = E_large - E_small
        print(f"~~~ Basis Correction={dE_basis[-1]} Har. ({i+1}/{npoints}) ~~~")

        # relativistic correction: (X2C-CCSD(T) - CCSD(T)) / cc-pCVTZ-DK
        qcdb.set_options({"psi4_relativistic": "x2c"})
        E_x2c_on, jrec = qcdb.energy("p4-ccsd(t)/aug-cc-pCVTZ-DK", return_wfn=True)
        qcdb.set_options({"psi4_relativistic": "no"})
        E_x2c_off, jrec = qcdb.energy("p4-ccsd(t)/aug-cc-pCVTZ-DK", return_wfn=True)
        dE_x2c[i] = E_x2c_on - E_x2c_off
        print(f"~~~ Relativistic Correction={dE_x2c[-1]} Har. ({i+1}/{npoints}) ~~~")

        # fci correction: (FCI - CCSD(T)) / cc-pVDZ
        E_cc, _ = qcdb.energy("gms-ccsd(t)/cc-pVDZ", return_wfn=True, local_options=local_options)
        E_fci, _ = qcdb.energy("gms-fci/cc-pVDZ", return_wfn=True, local_options=local_options)
        dE_fci[i] = E_fci - E_cc
        print(f"~~~ FCI Correction={dE_fci[-1]} Har. ({i+1}/{npoints}) ~~~")

        with open(f"dist_{i}", "a+") as f:
            f.write(f"Re {R_arr[i]} ,E_base     {E_base[i]} ")
            f.write("\n")
            f.write(f"Re {R_arr[i]} ,dE_basis   {dE_basis[i]} ")
            f.write("\n")
            f.write(f"Re {R_arr[i]} ,dE_dboc    {dE_dboc[i]} ")
            f.write("\n")
            f.write(f"Re {R_arr[i]} ,dE_x2c     {dE_x2c[i]} ")
            f.write("\n")
            f.write(f"Re {R_arr[i]} ,dE_fci     {dE_fci[i]} ")
            f.write("\n")
            f.write(f"Re {R_arr[i]} ,dE_ccsdtq  {dE_ccsdtq[i]} ")
            f.write("\n")

        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        print(f"Re {R_arr[i]} ,E_base     {E_base[i]}")
        print(f"Re {R_arr[i]} ,dE_basis   {dE_basis[i]}")
        print(f"Re {R_arr[i]} ,dE_x2c     {dE_x2c[i]}")
        print(f"Re {R_arr[i]} ,dE_fci     {dE_fci[i]}")
        print(f"Re {R_arr[i]} ,dE_ccsdtq  {dE_ccsdtq[i]}")
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

    # Energy with a single correction
    E_basis = []
    E_dboc = []
    E_x2c = []
    E_fci = []
    E_ccsdtq = []
    # Total energies, using all corrections
    E_tot_fci = []
    E_tot_ccsdtq = []
    for i in range(npoints):

        E_basis.append(E_base[i] + dE_basis[i])
        E_dboc.append(E_base[i] + dE_dboc[i])
        E_x2c.append(E_base[i] + dE_x2c[i])
        E_fci.append(E_base[i] + dE_fci[i])
        E_ccsdtq.append(E_base[i] + dE_ccsdtq[i])

        E_tot_fci.append(E_base[i] + dE_basis[i] + dE_dboc[i] + dE_x2c[i] + dE_fci[i])
        E_tot_ccsdtq.append(E_base[i] + dE_basis[i] + dE_dboc[i] + dE_x2c[i] + dE_ccsdtq[i])

    for i in range(len(R_arr)):
        print(
            "Re",
            R_arr[i],
        )
        print("   E_base", E_base[i])
        print("   dE_basis", dE_basis[i])
        print("   dE_x2c", dE_x2c[i])
        print("   dE_fci", dE_fci[i])
        print("   dE_ccsdtq", dE_ccsdtq[i])
        print("   E_tot_fci", E_tot_fci[i])
        print("   E_tot_ccsdtq", E_tot_ccsdtq[i])

    """ Spectroscopic Constants """

    # need a Psi4 molecule to give atomic masses to the diatomic module
    import psi4

    psi4.geometry(
        """
            B
            H 1 1000
            units au
            """
    )

    print("Total correction w/ FCI")
    phys_consts_tot_fci = qcdb.diatomic(R_arr, E_tot_fci, molecule=bh)
    pp.pprint(phys_consts_tot_fci)

    print("Total correction w/ CCSDTQ")
    phys_consts_tot_ccsdtq = qcdb.diatomic(R_arr, E_tot_ccsdtq, molecule=bh)
    pp.pprint(phys_consts_tot_ccsdtq)

    print("Base Energy")
    phys_consts_base = qcdb.diatomic(R_arr, E_base, molecule=bh)
    pp.pprint(phys_consts_base)

    print("Base Energy w/ only basis correction")
    phys_consts_basis = qcdb.diatomic(R_arr, E_basis, molecule=bh)
    pp.pprint(phys_consts_basis)

    print("Base Energy w/ only DBOC")
    phys_consts_dboc = qcdb.diatomic(R_arr, E_dboc, molecule=bh)
    pp.pprint(phys_consts_dboc)

    print("Base Energy w/ only rel. correction")
    phys_consts_x2c = qcdb.diatomic(R_arr, E_x2c, molecule=bh)
    pp.pprint(phys_consts_x2c)

    print("Base Energy w/ only FCI correction")
    phys_consts_fci = qcdb.diatomic(R_arr, E_fci, molecule=bh)
    pp.pprint(phys_consts_fci)

    print("Base Energy w/ only CCSDTQ correction")
    phys_consts_ccsdtq = qcdb.diatomic(R_arr, E_ccsdtq, molecule=bh)
    pp.pprint(phys_consts_ccsdtq)
