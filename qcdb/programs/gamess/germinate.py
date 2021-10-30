import itertools
import math
import re
import uuid
from collections import defaultdict
from typing import Dict, List, Tuple

import qcelemental as qcel
import qcengine as qcng
from qcelemental.util import which

from ...exceptions import ValidationError
from ...molecule import Molecule
from ...util import conv_float2negexp


def _get_symmetry_card(pg: str, full_pg_n: int) -> Tuple[str, str]:
    """Assemble card -2- for GAMESS input from ``pg`` (result of Molecule.full_point_group_with_n()) and ``full_pg_n`` (result of Molecule.full_pg_n())."""

    # PSI: FullPointGroupList = ["ATOM", "C_inf_v", "D_inf_h", "C1", "Cs", "Ci", "Cn", "Cnv", "Cnh", "Sn", "Dn", "Dnd", "Dnh", "Td", "Oh", "Ih"]
    # GMS:                                                      C1    Cs    Ci    Cn    Cnv    Cnh          Dn    Dnd    Dnh    Td    Oh
    # GMS:                        Dnh-2   Cnv-4      Dnh-4                                            S2n
    # GMS:    T, Th, O
    # GAMESS Manual: "For linear molecules, choose either Cnv or Dnh, and enter NAXIS as 4. Enter atoms as Dnh with NAXIS=2."

    if pg == "ATOM":
        return "Dnh", 2
    elif pg == "C_inf_v":
        return "Cnv", 2
    elif pg == "D_inf_h":
        return "Dnh", 4
    elif pg == "Sn":
        return "S2n", int(full_pg_n / 2)
    elif "n" in pg:
        return pg, full_pg_n
    else:
        return pg, ""


def get_master_frame(
    kmol: "qcelemental.models.Molecule", scratch_directory
) -> Tuple["qcelemental.models.Molecule", Dict[str, str]]:
    """Do whatever it takes to figure out the GAMESS master frame by which ``kmol`` can be run with full symmetry."""

    harness = qcng.get_program("gamess")

    # want the full frame-independent symmetry, so allow reorientation to Psi4 master frame
    qmol = Molecule.from_schema(kmol.dict() | {"fix_com": False, "fix_orientation": False})
    pgn, naxis = _get_symmetry_card(qmol.full_point_group_with_n(), qmol.full_pg_n())

    # run exetyp=check asserting full symmetry to extract master frame from GAMESS
    # * fix_*=F so harness returns the internal GAMESS frame, not the naive input frame
    # * uses an arbitrary UHF/6-31G model
    # * most common failure mode is high-symmetry or wrong-quadrant geometry when GAMESS generates too many atoms
    #   * haven't found a reliable programmatic path out that uses full internal symmetry, so fall back to C2v, then C1
    internal_symmetry_card = f"{pgn} {naxis}".strip()

    if not all(kmol.real):
        # TODO is this the best way to handle ghosts?
        # * the C1 early return is to avoid the qcng.compute that uses coord=prinaxis that can't take ghost specification
        # * is the gamess master_frame the input frame for C1? have tentatively assumed so
        data = {
            "unique": list(range(len(kmol.symbols))),
            "symmetry_card": "C1",
        }
        return kmol, data

    for symmetry_card in [internal_symmetry_card, "Cnv 2", "C1"]:
        naive_kmol = kmol.copy(
            update={
                "fix_symmetry": symmetry_card,
                "atom_labels": list(range(len(kmol.symbols))),
                "fix_com": False,
                "fix_orientation": False,
            }
        )

        atin = qcel.models.AtomicInput(
            **{
                "driver": "energy",
                "keywords": {
                    "contrl__exetyp": "check",
                    "contrl__scftyp": "uhf",
                    "basis__ngauss": 6,
                },
                "model": {
                    "method": "hf",
                    "basis": "n31",
                },
                "molecule": naive_kmol,
            }
        )

        try:
            atres = qcng.compute(
                atin, "gamess", local_options={"nnodes": 1, "ncores": 1, "memory": 1}, raise_error=True
            )
        except qcng.exceptions.UnknownError as e:
            mobj = re.search(
                # fmt: off
                r"^\s+" + r"AFTER PRINCIPAL AXIS TRANSFORMATION, THE PROGRAM" + r"\s*" +
                r"^\s+" + r"HAS CHOSEN THE FOLLOWING ATOMS AS BEING UNIQUE:" + r"\s*" +
                r"((?:\s+([A-Za-z]\w*)\s+\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s*\n)+)" +
                r"^\s+" + r"EXECUTION OF GAMESS TERMINATED -ABNORMALLY-",
                # fmt: on
                str(e),
                re.MULTILINE | re.IGNORECASE,
            )
            if mobj:
                pass
            else:
                raise e
        else:
            break

    # `mf_kmol` and `mf_qmol` are now in GAMESS master frame
    # * all atoms present (as natural for Molecule classes), and we don't know which atoms are unique
    mf_kmol, _ = kmol.align(atres.molecule, atoms_map=False, mols_align=True, run_mirror=True, verbose=0)

    mf_qmol = Molecule.from_schema(mf_kmol.dict())
    assert (
        mf_qmol.full_point_group_with_n() == qmol.full_point_group_with_n()
    ), f"{mf_qmol.full_point_group_with_n()} (mf) != {qmol.full_point_group_with_n()} (in)"
    assert mf_qmol.full_pg_n() == qmol.full_pg_n(), f"{mf_qmol.full_pg_n()} (mf) != {qmol.full_pg_n()} (in)"
    assert abs(mf_qmol.nuclear_repulsion_energy() - qmol.nuclear_repulsion_energy()) < 1.0e-3, "NRE"

    # nunique machinery in psi4/qcdb.Molecule class works within Abelian point groups, so start there
    d2h_subgroup = qmol.point_group().symbol()
    d2h_unique = [mf_qmol.unique(i) for i in range(mf_qmol.nunique())]

    if mf_qmol.get_full_point_group() == d2h_subgroup:
        # Abelian! home free
        full_pg_unique = d2h_unique

    else:
        # `possibly_nonunique` are atom indices that could be redundant with post-D2h symmetry operations
        # * formed from D2h unique less the first index for each element
        d2h_unique_by_element = defaultdict(list)
        for i in d2h_unique:
            d2h_unique_by_element[mf_qmol.symbol(i)].append(i)
        possibly_nonunique = []
        for k, v in d2h_unique_by_element.items():
            possibly_nonunique.extend(v[1:])

        # `trials` are a brute-force set of all possible atom indices, one or more of which must be the post-D2h unique list
        trials = []
        for drop in range(len(possibly_nonunique) + 1):
            for aa in itertools.combinations(possibly_nonunique, drop):
                trial = sorted(list(set(d2h_unique) - set(aa)))
                trials.append(trial)

        all_atom_lines = mf_kmol.to_string(dtype="gamess").splitlines()[3:]

        for selected_atoms in trials:
            selected_atom_lines = "\n".join([ln for iln, ln in enumerate(all_atom_lines) if iln in selected_atoms])
            inp = _get_exetype_check_input(
                selected_atoms,
                symmetry_card,
                mf_kmol.molecular_charge,
                mf_kmol.molecular_multiplicity,
                selected_atom_lines,
            )

            # run exetyp=check asserting full symmetry to find a unique list that doesn't generate overlapping atoms
            # * can't use AtomicInput into usual harness b/c partial atom list not expressible in Molecule
            gamessrec = {
                "infiles": {
                    "gamess.inp": inp,
                },
                "command": [which("rungms"), "gamess", "00", "1"],
                "scratch_directory": scratch_directory,
                "scratch_messy": False,
            }
            success, dexe = harness.execute(gamessrec)
            # pprint.pprint(dexe, width=200)

            if "THERE ARE ATOMS LESS THAN   0.100 APART, QUITTING..." not in dexe["stdout"]:
                # "THE NUCLEAR REPULSION ENERGY IS       12.7621426235"
                break

        full_pg_unique = selected_atoms

    data = {
        "unique": full_pg_unique,
        "symmetry_card": symmetry_card,
    }

    return mf_kmol, data


def _get_exetype_check_input(
    comment: str, symmetry_card: str, molecular_charge: float, molecular_multiplicity: int, atom_lines: str
) -> str:
    return f"""
 $basis gbasis=n31 ngauss=6 $end
 $contrl coord=unique exetyp=check icharg={molecular_charge} mult={molecular_multiplicity} runtyp=energy
  units=bohr $end
 $system memddi=0 mwords=100 $end
 $data
trial {comment}
{symmetry_card}

{atom_lines}
 $end

"""


def muster_molecule_and_basisset(
    molrec: Dict, qbs: "BasisSet", ropts: "Keywords", full_pg_unique: List[int], symmetry_card: str, verbose: int = 1
) -> str:
    kwgs = {"accession": uuid.uuid4(), "verbose": verbose}
    units = "Bohr"

    native_puream = qbs.has_puream()
    all_atom_basisset = qbs.print_detail_gamess(return_list=True)

    data_group_cart = qcel.molparse.to_string(
        molrec, dtype="gamess", units=units, atom_format=None, ghost_format=None, width=17, prec=12
    ).splitlines()
    all_atom_lines = data_group_cart[3:]

    data_group_uniq = data_group_cart[:2]  # $data and card -1-
    data_group_uniq.append(f""" {symmetry_card}""")  # card -2-
    if symmetry_card != "C1":
        data_group_uniq.append("")  # empty cards -3- and -4-

    for iat in range(len(molrec["elem"])):
        if iat in full_pg_unique:
            data_group_uniq.append(all_atom_lines[iat])  # card -5U-
            data_group_uniq.extend(all_atom_basisset[iat].splitlines()[1:])  # cards -6U- and -7U-
            data_group_uniq.append("")  # card -8U-

    data_group_uniq.append(""" $end""")

    ropts.require("GAMESS", "contrl__coord", "unique", **kwgs)
    ropts.require("GAMESS", "contrl__units", {"Bohr": "bohr", "Angstrom": "angs"}[units], **kwgs)
    ropts.require("GAMESS", "contrl__icharg", int(molrec["molecular_charge"]), **kwgs)
    ropts.require("GAMESS", "contrl__mult", molrec["molecular_multiplicity"], **kwgs)
    ropts.require("GAMESS", "contrl__ispher", {True: 1, False: -1}[native_puream], **kwgs)

    return "\n".join(data_group_uniq)


def muster_modelchem(name: str, dertype: int, ropts: "Keywords", sysinfo: Dict, verbose: int = 1) -> None:
    lowername = name.lower()
    accession = uuid.uuid4()

    # runtyp = {'energy': 'energy',
    #          'gradient': 'gradient',
    #          'hessian': 'hessian',
    #          'properties': 'prop',
    #         }[driver]

    runtyp = {
        0: "energy",
        1: "gradient",
        2: "hessian",
        #'properties': 'prop',
    }[dertype]

    ropts.require("GAMESS", "contrl__runtyp", runtyp, accession=accession, verbose=verbose)

    if lowername == "gms-gamess":
        pass

    elif lowername in ["gms-scf", "gms-hf"]:
        ropts.require("GAMESS", "contrl__mplevl", 0, accession=accession, verbose=verbose)
        ropts.require("GAMESS", "contrl__cityp", "none", accession=accession, verbose=verbose)
        ropts.require("GAMESS", "contrl__cctyp", "none", accession=accession, verbose=verbose)

    elif lowername == "gms-mp2":
        ropts.require("GAMESS", "contrl__mplevl", 2, accession=accession, verbose=verbose)
        ropts.require("GAMESS", "contrl__cityp", "none", accession=accession, verbose=verbose)
        ropts.require("GAMESS", "contrl__cctyp", "none", accession=accession, verbose=verbose)

    elif lowername == "gms-cisd":
        # todo qc_module CITYP
        ropts.require("GAMESS", "contrl__mplevl", 0, accession=accession, verbose=verbose)
        ropts.require("GAMESS", "contrl__cctyp", "none", accession=accession, verbose=verbose)

        qopt = ropts.scroll["QCDB"]["FREEZE_CORE"]
        # if qopt.disputed():
        if True:
            if qopt.value is True:
                fcae = "fc"
            elif qopt.value is False:
                fcae = "ae"
        nocc = int(math.ceil(sysinfo[fcae]["nels"] / 2))
        nbeta = int(math.floor(sysinfo[fcae]["nels"] / 2))

        if ropts.scroll["QCDB"]["QC_MODULE"].value == "fsoci":
            ropts.require("GAMESS", "contrl__cityp", "fsoci", accession=accession, verbose=verbose)
            ropts.suggest("GAMESS", "cidet__ncore", sysinfo[fcae]["ncore"], accession=accession, verbose=verbose)
            ropts.suggest("GAMESS", "cidet__nact", nocc, accession=accession, verbose=verbose)
            ropts.suggest("GAMESS", "cidet__nels", sysinfo[fcae]["nels"], accession=accession, verbose=verbose)
        elif ropts.scroll["QCDB"]["QC_MODULE"].value == "guga":
            ropts.require("GAMESS", "contrl__cityp", "guga", accession=accession, verbose=verbose)
            ropts.require("GAMESS", "cidrt__iexcit", 2, accession=accession, verbose=verbose)
            ropts.suggest("GAMESS", "cidrt__nfzc", sysinfo[fcae]["ncore"], accession=accession, verbose=verbose)
            ropts.suggest("GAMESS", "cidrt__group", "C1", accession=accession, verbose=verbose)  # TODO tailor to mol
            # ropts.suggest('GAMESS', 'cidrt__ndoc', nocc, accession=accession, verbose=verbose)

            ropts.suggest("GAMESS", "cidrt__ndoc", nbeta, accession=accession, verbose=verbose)
            ropts.suggest("GAMESS", "cidrt__nalp", nocc - nbeta, accession=accession, verbose=verbose)
            ropts.suggest("GAMESS", "cidrt__nval", sysinfo[fcae]["nact"] - nocc, accession=accession, verbose=verbose)
            # ropts.suggest('GAMESS', 'cidrt__intact', True, accession=accession, verbose=verbose)
    # $cidrt  GROUP=C2V IEXCIT=4 NFZC=0 NDOC=5 NVAL=9 $END

    elif lowername == "gms-lccd":
        ropts.require("GAMESS", "contrl__mplevl", 0, accession=accession, verbose=verbose)
        ropts.require("GAMESS", "contrl__cityp", "none", accession=accession, verbose=verbose)
        ropts.require("GAMESS", "contrl__cctyp", "lccd", accession=accession, verbose=verbose)

    elif lowername == "gms-ccd":
        ropts.require("GAMESS", "contrl__mplevl", 0, accession=accession, verbose=verbose)
        ropts.require("GAMESS", "contrl__cityp", "none", accession=accession, verbose=verbose)
        ropts.require("GAMESS", "contrl__cctyp", "ccd", accession=accession, verbose=verbose)

    elif lowername == "gms-ccsd":
        ropts.require("GAMESS", "contrl__mplevl", 0, accession=accession, verbose=verbose)
        ropts.require("GAMESS", "contrl__cityp", "none", accession=accession, verbose=verbose)
        ropts.require("GAMESS", "contrl__cctyp", "ccsd", accession=accession, verbose=verbose)

    elif lowername in ["gms-ccsd+t(ccsd)", "gms-ccsd(t)"]:
        ropts.require("GAMESS", "contrl__mplevl", 0, accession=accession, verbose=verbose)
        ropts.require("GAMESS", "contrl__cityp", "none", accession=accession, verbose=verbose)
        ropts.require("GAMESS", "contrl__cctyp", "ccsd(t)", accession=accession, verbose=verbose)

    elif lowername == "gms-cr-ccl":
        ropts.require("GAMESS", "contrl__mplevl", 0, accession=accession, verbose=verbose)
        ropts.require("GAMESS", "contrl__cityp", "none", accession=accession, verbose=verbose)
        ropts.require("GAMESS", "contrl__cctyp", "cr-ccl", accession=accession, verbose=verbose)

    elif lowername == "gms-ccsd(tq)":
        ropts.require("GAMESS", "contrl__mplevl", 0, accession=accession, verbose=verbose)
        ropts.require("GAMESS", "contrl__cityp", "none", accession=accession, verbose=verbose)
        ropts.require("GAMESS", "contrl__cctyp", "ccsd(tq)", accession=accession, verbose=verbose)

    elif lowername == "gms-fci":
        ropts.require("GAMESS", "contrl__mplevl", 0, accession=accession, verbose=verbose)
        ropts.require("GAMESS", "contrl__cityp", "aldet", accession=accession, verbose=verbose)
        ropts.require("GAMESS", "contrl__cctyp", "none", accession=accession, verbose=verbose)

        qopt = ropts.scroll["QCDB"]["FREEZE_CORE"]
        if True:
            if qopt.value is True:
                fcae = "fc"
            elif qopt.value is False:
                fcae = "ae"

        ropts.suggest("GAMESS", "cidet__ncore", sysinfo[fcae]["ncore"], accession=accession, verbose=verbose)
        ropts.suggest("GAMESS", "cidet__nact", sysinfo[fcae]["nact"], accession=accession, verbose=verbose)
        ropts.suggest("GAMESS", "cidet__nels", sysinfo[fcae]["nels"], accession=accession, verbose=verbose)

    elif lowername == "gms-pbe":
        ropts.require("GAMESS", "contrl__dfttyp", "pbe", accession=accession, verbose=verbose)

    elif lowername == "gms-b3lyp":
        ropts.require("GAMESS", "contrl__dfttyp", "b3lypv1r", accession=accession, verbose=verbose)

    elif lowername == "gms-b3lyp5":
        ropts.require("GAMESS", "contrl__dfttyp", "b3lyp", accession=accession, verbose=verbose)

    # unused from Nuwan

    #    elif lowername == 'gms-dft':
    #        if dertype == 0:
    #            options ['GAMESS']['GAMESS_CONTRL_RUNTYP']['value'] = 'energy'
    #            options ['GAMESS']['GAMESS_CONTRL_DFTTYP']['value'] = 'b3lyp'
    #        elif dertype == 1:
    #            options ['GAMESS']['GAMESS_CONTRL_RUNTYP']['value'] = 'gradient'
    #            options ['GAMESS']['GAMESS_CONTRL_DFTTYP']['value'] = 'b3lyp'
    #        elif dertype == 2:
    #            options ['GAMESS']['GAMESS_CONTRL_RUNTYP']['value'] = 'hessian'
    #            options ['GAMESS']['GAMESS_CONTRL_DFTTYP']['value'] = 'b3lyp'
    #
    #    elif lowername == 'gms-eom-ccsd':
    #        if dertype == 0:
    #            options ['GAMESS']['GAMESS_CONTRL_RUNTYP']['value'] = 'energy'
    #            options ['GAMESS']['GAMESS_CONTRL_CCTYP']['value']  = 'eom-ccsd'
    #        elif dertype == 1:
    #            options ['GAMESS']['GAMESS_CONTRL_RUNTYP']['value'] = 'gradient'
    #            options ['GAMESS']['GAMESS_CONTRL_CCTYP']['value']  = 'eom-ccsd'
    #        elif dertype == 2:
    #            options ['GAMESS']['GAMESS_CONTRL_RUNTYP']['value'] = 'hessian'
    #            options ['GAMESS']['GAMESS_CONTRL_CCTYP']['value']  = 'eom-ccsd'
    #
    #    elif lowername == 'gms-cis':
    #        if dertype == 0:
    #            options ['GAMESS']['GAMESS_CONTRL_RUNTYP']['value'] = 'energy'
    #            options ['GAMESS']['GAMESS_CONTRL_CITYP']['value']  = 'cis'
    #        elif dertype == 1:
    #            options ['GAMESS']['GAMESS_CONTRL_RUNTYP']['value'] = 'gradient'
    #            options ['GAMESS']['GAMESS_CONTRL_CITYP']['value']  = 'cis'
    #        elif dertype == 2:
    #            options ['GAMESS']['GAMESS_CONTRL_RUNTYP']['value'] = 'hessian'
    #            options ['GAMESS']['GAMESS_CONTRL_CITYP']['value']  = 'cis'
    #
    #    elif lowername == 'gms-efp':
    #        if dertype == 0:
    #            options ['GAMESS']['GAMESS_CONTRL_RUNTYP']['value'] = 'energy'
    #            options ['GAMESS']['GAMESS_CONTRL_COORD']['value']  = 'fragonly'
    #        elif dertype == 1:
    #            options ['GAMESS']['GAMESS_CONTRL_COORD']['value']  = 'gradient'
    #            options ['GAMESS']['GAMESS_CONTRL_COORD']['value']  = 'fragonly'
    #        elif dertype == 2:
    #            options ['GAMESS']['GAMESS_CONTRL_RUNTYP']['value'] = 'hessian'
    #            options ['GAMESS']['GAMESS_CONTRL_COORD']['value']  = 'fragonly'

    else:
        raise ValidationError(f"""Requested GAMESS computational methods {lowername} is not available.""")


def muster_inherited_keywords(ropts: "Keywords", sysinfo: Dict, verbose: int = 1) -> None:
    accession = uuid.uuid4()

    kwgs = {"accession": accession, "verbose": verbose}
    do_translate = ropts.scroll["QCDB"]["TRANSLATE_QCDB"].value

    # qcdb/memory [B] --> gamess/system__mwords [M QW]
    qopt = ropts.scroll["QCDB"]["MEMORY"]
    if do_translate or qopt.is_required():
        mem = int(qopt.value / 8e6)
        print("\n\nMEMORY", mem, "\n\n")
        ropts.suggest("GAMESS", "system__mwords", mem, **kwgs)

    # qcdb/reference --> gamess/contrl__scftyp
    # TODO ref or scf__ref?
    qref = ropts.scroll["QCDB"]["SCF__REFERENCE"].value
    if qref in ["RHF", "UHF", "ROHF"]:
        # ref = {'RHF': 'RHF',
        #       'UHF': 'UHF',
        #       'ROHF': 'ROHF'}[ropts.scroll['QCDB']['REFERENCE'].value]
        ropts.suggest("GAMESS", "contrl__scftyp", qref, **kwgs)

    # qcdb/scf__d_convergence --> gamess/scf__conv
    qopt = ropts.scroll["QCDB"]["SCF__D_CONVERGENCE"]
    if qopt.disputed():
        conv = conv_float2negexp(qopt.value)
        ropts.suggest("GAMESS", "scf__conv", conv, **kwgs)

    # --> gamess/ccinp__iconv
    qopt = ropts.scroll["QCDB"]["E_CONVERGENCE"]
    if qopt.disputed():
        conv = conv_float2negexp(qopt.value)
        ropts.suggest("GAMESS", "ccinp__iconv", conv, **kwgs)

    # qcdb/freeze_core --> gamess/
    qopt = ropts.scroll["QCDB"]["FREEZE_CORE"]
    if qopt.disputed():
        if qopt.value is True:
            ncore = sysinfo["fc"]["ncore"]
        elif qopt.value is False:
            ncore = 0
        ropts.suggest("GAMESS", "ccinp__ncore", ncore, accession=accession, verbose=verbose)
        ropts.suggest("GAMESS", "mp2__nacore", ncore, accession=accession, verbose=verbose)
