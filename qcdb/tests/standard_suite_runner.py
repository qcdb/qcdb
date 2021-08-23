import re
import pprint
from typing import Any, Dict

import pytest
import numpy as np
from qcelemental.molutil import compute_scramble
from qcengine.programs.tests.standard_suite_contracts import (
    contractual_hf,
    contractual_mp2,
    contractual_mp2p5,
    contractual_mp3,
    contractual_mp4_prsdq_pr,
    contractual_mp4,
    contractual_cisd,
    contractual_qcisd,
    contractual_qcisd_prt_pr,
    contractual_lccd,
    contractual_lccsd,
    contractual_ccd,
    contractual_ccsd,
    contractual_ccsdpt_prccsd_pr,
    contractual_ccsd_prt_pr,
    contractual_accsd_prt_pr,
    contractual_ccsdt1a,
    contractual_ccsdt1b,
    contractual_ccsdt2,
    contractual_ccsdt3,
    contractual_ccsdt,
    contractual_ccsdt_prq_pr,
    contractual_ccsdtq,
    contractual_dft_current,
    contractual_current,
    query_has_qcvar,
    query_qcvar,
)
from qcengine.programs.tests.standard_suite_ref import answer_hash, std_suite
from qcengine.programs.util import mill_qcvars

from .utils import compare, compare_values

pp = pprint.PrettyPrinter(width=120)


def runner_asserter(inp, ref_subject, method, basis, tnm, scramble, frame):

    qcprog = inp["qc_module"].split("-")[0]
    qc_module_in = inp["qc_module"]  # returns "<qcprog>"|"<qcprog>-<module>"  # input-specified routing
    qc_module_xptd = (
        (qcprog + "-" + inp["xptd"]["qc_module"]) if inp.get("xptd", {}).get("qc_module", None) else None
    )  # expected routing
    driver = inp["driver"]
    reference = inp["reference"]
    fcae = inp["fcae"]

    # <<<  Molecule  >>>

    # 1. ref mol: `ref_subject` nicely oriented mol taken from standard_suite_ref.py
    ref_subject.update_geometry()
    min_nonzero_coords = np.count_nonzero(np.abs(ref_subject.geometry(np_out=True)) > 1.e-10)

    if scramble is None:
        subject = ref_subject
        ref2in_mill = compute_scramble(subject.natom(), do_resort=False, do_shift=False, do_rotate=False, do_mirror=False)  # identity AlignmentMill

    else:
        subject, scramble_data = ref_subject.scramble(**scramble, do_test=False, fix_mode="copy")
        ref2in_mill = scramble_data["mill"]

    # 2. input mol: `subject` now ready for `atin.molecule`. may have been scrambled away from nice ref orientation

    # <<<  Reference Values  >>>

    # ? precedence on next two
    mp2_type = inp.get("corl_type", inp["keywords"].get("mp2_type", "df"))  # hard-code of read_options.cc MP2_TYPE
    mp_type = inp.get("corl_type", inp["keywords"].get("mp_type", "conv"))  # hard-code of read_options.cc MP_TYPE
    ci_type = inp.get("corl_type", inp["keywords"].get("ci_type", "conv"))  # hard-code of read_options.cc CI_TYPE
    cc_type = inp.get("corl_type", inp["keywords"].get("cc_type", "conv"))  # hard-code of read_options.cc CC_TYPE
    corl_natural_values = {
        "hf": "conv",  # dummy to assure df/cd/conv scf_type refs available
        "mp2": mp2_type,
        "mp3": mp_type,
        "mp4(sdq)": mp_type,
        "mp4": mp_type,
        "cisd": ci_type,
        "qcisd": ci_type,
        "qcisd(t)": ci_type,
        "lccd": cc_type,
        "lccsd": cc_type,
        "ccd": cc_type,
        "ccsd": cc_type,
        "ccsd+t(ccsd)": cc_type,
        "ccsd(t)": cc_type,
        "a-ccsd(t)": cc_type,
        "ccsdt-1a": cc_type,
        "ccsdt-1b": cc_type,
        "ccsdt-2": cc_type,
        "ccsdt-3": cc_type,
        "ccsdt": cc_type,
        "ccsdt(q)": cc_type,
        "ccsdtq": cc_type,

        "pbe": "conv",
        "b3lyp": "conv",
        "b3lyp5": "conv",
    }
    corl_type = corl_natural_values[method]

    natural_ref = {"conv": "pk", "df": "df", "cd": "cd"}
    scf_type = inp["keywords"].get("scf_type", natural_ref[corl_type])
    natural_values = {"pk": "pk", "direct": "pk", "df": "df", "mem_df": "df", "disk_df": "df", "cd": "cd"}
    scf_type = natural_values[scf_type]

    is_dft = (method in ["pbe", "b3lyp", "b3lyp5"])

    # * absolute and relative tolerances function approx as `or` operation. see https://numpy.org/doc/stable/reference/generated/numpy.allclose.html
    # * can't go lower on atol_e because hit digit limits accessible for reference values
    # * dz gradients tend to be less accurate than larger basis sets/mols
    # * analytic Hessian very loose to catch gms/nwc HF Hessian
    atol_e, rtol_e = 2.e-7, 1.e-16
    atol_g, rtol_g = 5.e-7, 2.e-5
    atol_h, rtol_h = 1.e-5, 2.e-5
    if is_dft:
        atol_g = 6.e-6
    using_fd = "xptd" in inp and "fd" in inp["xptd"]  # T/F: notate fd vs. anal for docs table
    loose_fd = inp.get("xptd", {}).get("fd", False)  # T/F: relax conv crit for 3-pt internal findif fd
    if loose_fd:
        if basis == "cc-pvdz":
            atol_g = 1.e-4
            atol_h, rtol_h = 1.e-4, 5.e-4
        else:
            atol_g = 2.e-5
            atol_h, rtol_h = 5.e-5, 2.e-4

# VIEW    atol_e, atol_g, atol_h, rtol_e, rtol_g, rtol_h = 1.e-9, 1.e-9, 1.e-9, 1.e-16, 1.e-16, 1.e-16

    chash = answer_hash(
        system=subject.name(),
        basis=basis,
        fcae=fcae,
        scf_type=scf_type,
        reference=reference,
        corl_type=corl_type,
    )
    ref_block = std_suite[chash]

    # check all calcs against conventional reference to looser tolerance
    atol_conv = 1.0e-4
    rtol_conv = 1.0e-3
    chash_conv = answer_hash(
        system=subject.name(),
        basis=basis,
        fcae=fcae,
        reference=reference,
        corl_type="conv",
        scf_type="pk",
    )
    ref_block_conv = std_suite[chash_conv]

    # <<<  Prepare Calculation and Call API  >>>

    import qcdb

    driver_call = {"energy": qcdb.energy, "gradient": qcdb.gradient, "hessian": qcdb.hessian}
    local_options = {"nnodes": 1, "ncores": 1, "scratch_messy": False, "memory": 4}

    qcdb.set_options(
        {
            # "guess": "sad",
            # "e_convergence": 8,
            # "d_convergence": 7,
            # "r_convergence": 7,
            "e_convergence": 10,
            "d_convergence": 9,
            # "r_convergence": 9,
            # "points": 5,
        }
    )
    extra_kwargs = inp["keywords"].pop("function_kwargs", {})
    qcdb.set_options(inp["keywords"])

    if "error" in inp:
        errtype, errmatch, reason = inp["error"]
        with pytest.raises(errtype) as e:
            driver_call[driver](inp["call"], molecule=subject, local_options=local_options, **extra_kwargs)

        assert re.search(errmatch, str(e.value)), f"Not found: {errtype} '{errmatch}' in {e.value}"
        _recorder(qcprog, qc_module_in, driver, method, reference, fcae, scf_type, corl_type, "error", "nyi: " + reason)
        return

    ret, wfn = driver_call[driver](
        inp["call"], molecule=subject, return_wfn=True, local_options=local_options, **extra_kwargs
    )

    print("WFN")
    pp.pprint(wfn)

    qc_module_out = wfn["provenance"]["creator"].lower()
    if "module" in wfn["provenance"]:
        qc_module_out += "-" + wfn["provenance"]["module"]  # returns "<qcprog>-<module>"
    # assert 0, f"{qc_module_xptd=} {qc_module_in=} {qc_module_out=}"  # debug

    # 3. output mol: `wfn.molecule` after calc. orientation for nonscalar quantities may be different from `subject` if fix_=False
    wfn_molecule = qcdb.Molecule.from_schema(wfn["molecule"])

    _, ref2out_mill, _ = ref_subject.B787(wfn_molecule, atoms_map=False, mols_align=True, fix_mode="true", verbose=0)

    if subject.com_fixed() and subject.orientation_fixed():
        assert frame == "fixed"
        with np.printoptions(precision=3, suppress=True):
            assert compare_values(subject.geometry(), wfn_molecule.geometry(), atol=5.e-8), f"coords: atres ({wfn_molecule.geometry(np_out=True)}) != atin ({subject.geometry(np_out=True)})"  # 10 too much
        assert (
            ref_subject.com_fixed()
            and ref_subject.orientation_fixed()
            and subject.com_fixed()
            and subject.orientation_fixed()
            and wfn_molecule.com_fixed()
            and wfn_molecule.orientation_fixed()
        ), f"fixed, so all T: {ref_subject.com_fixed()} {ref_subject.orientation_fixed()} {subject.com_fixed()} {subject.orientation_fixed()} {wfn_molecule.com_fixed()} {wfn_molecule.orientation_fixed()}"

        ref_block = mill_qcvars(ref2in_mill, ref_block)
        ref_block_conv = mill_qcvars(ref2in_mill, ref_block_conv)

    else:
        assert frame == "free" or frame == ""  # "": direct from standard_suite_ref.std_molecules
        with np.printoptions(precision=3, suppress=True):
            assert compare(min_nonzero_coords, np.count_nonzero(np.abs(wfn_molecule.geometry(np_out=True)) > 1.e-10), tnm + " !0 coords wfn"), f"ncoords {wfn_molecule.geometry(np_out=True)} != {min_nonzero_coords}"
        assert (
            (not ref_subject.com_fixed())
            and (not ref_subject.orientation_fixed())
            and (not subject.com_fixed())
            and (not subject.orientation_fixed())
            and (not wfn_molecule.com_fixed())
            and (not wfn_molecule.orientation_fixed())
        ), f"free, so all F: {ref_subject.com_fixed()} {ref_subject.orientation_fixed()} {subject.com_fixed()} {subject.orientation_fixed()} {wfn_molecule.com_fixed()} {wfn_molecule.orientation_fixed()}"

        if scramble is None:
            # wfn exactly matches ref_subject and ref_block
            with np.printoptions(precision=3, suppress=True):
                assert compare_values(ref_subject.geometry(), wfn_molecule.geometry(), atol=5.e-8), f"coords: atres ({wfn_molecule.geometry(np_out=True)}) != atin ({ref_subject.geometry(np_out=True)})"
        else:
            # wfn is "pretty" (max zeros) but likely not exactly ref_block (by axis exchange, phasing, atom shuffling) since Psi4 ref frame is not unique
            ref_block = mill_qcvars(ref2out_mill, ref_block)
            ref_block_conv = mill_qcvars(ref2out_mill, ref_block_conv)


    # <<<  Comparison Tests  >>>

    assert wfn["success"] is True
    assert (
        wfn["provenance"]["creator"].lower() == qcprog
    ), f'ENGINE used ({ wfn["provenance"]["creator"].lower()}) != requested ({qcprog})'

    # qcvars
    contractual_args = [
        qc_module_out,
        driver,
        reference,
        method,
        corl_type,
        fcae,
    ]
    asserter_args = [
        [qcdb, wfn["qcvars"]],
        ref_block,
        [atol_e, atol_g, atol_h],
        [rtol_e, rtol_g, rtol_h],
        ref_block_conv,
        atol_conv,
        rtol_conv,
        tnm,
    ]

    def qcvar_assertions():
        print("BLOCK", chash, contractual_args)
        if method == "hf":
            _asserter(asserter_args, contractual_args, contractual_hf)
        elif method == "mp2":
            _asserter(asserter_args, contractual_args, contractual_mp2)
        elif method == "mp3":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_mp2p5)
            _asserter(asserter_args, contractual_args, contractual_mp3)
        elif method == "mp4(sdq)":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_mp2p5)
            _asserter(asserter_args, contractual_args, contractual_mp3)
            _asserter(asserter_args, contractual_args, contractual_mp4_prsdq_pr)
        elif method == "mp4":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_mp2p5)
            _asserter(asserter_args, contractual_args, contractual_mp3)
            _asserter(asserter_args, contractual_args, contractual_mp4_prsdq_pr)
            _asserter(asserter_args, contractual_args, contractual_mp4)
        elif method == "cisd":
            _asserter(asserter_args, contractual_args, contractual_cisd)
        elif method == "qcisd":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_qcisd)
        elif method == "qcisd(t)":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_qcisd)
            _asserter(asserter_args, contractual_args, contractual_qcisd_prt_pr)
        elif method == "lccd":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_lccd)
        elif method == "lccsd":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_lccsd)
        elif method == "ccd":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_ccd)
        elif method == "ccsd":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_ccsd)
        elif method == "ccsd+t(ccsd)":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_ccsd)
            _asserter(asserter_args, contractual_args, contractual_ccsdpt_prccsd_pr)
        elif method == "ccsd(t)":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_ccsd)
            _asserter(asserter_args, contractual_args, contractual_ccsd_prt_pr)
        elif method == "a-ccsd(t)":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_ccsd)
            _asserter(asserter_args, contractual_args, contractual_accsd_prt_pr)
        elif method == "ccsdt-1a":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_ccsdt1a)
        elif method == "ccsdt-1b":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_ccsdt1b)
        elif method == "ccsdt-2":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_ccsdt2)
        elif method == "ccsdt-3":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_ccsdt3)
        elif method == "ccsdt":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_ccsdt)
        elif method == "ccsdt(q)":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_ccsdt)
            _asserter(asserter_args, contractual_args, contractual_ccsdt_prq_pr)
        elif method == "ccsdtq":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_ccsdtq)
        # separations here for DFT appropriate when qcvars are labeled by functional

    if "wrong" in inp:
        if basis == "cc-pvdz" and contractual_args in [["cfour-ecc", "gradient", "rhf", mtd, "conv", "fc"] for mtd in ["ccsdt-1a", "ccsdt-1b", "ccsdt-2", "ccsdt-3"]]:
            # these four tests have pass/fail too close for dz to "get it right" with general tolerances
            pass
        else:
            errmatch, reason = inp["wrong"]
            with pytest.raises(AssertionError) as e:
                qcvar_assertions()

            assert errmatch in str(e.value), f"Not found: AssertionError '{errmatch}' for '{reason}' in {e.value}"
            _recorder(qcprog, qc_module_out, driver, method, reference, fcae, scf_type, corl_type, "wrong", reason + f" First wrong at `{errmatch}`.")
            pytest.xfail(reason)

    # primary label checks
    qcvar_assertions()

    # routing checks
    if qc_module_in != qcprog:
        assert qc_module_out == qc_module_in, f"QC_MODULE used ({qc_module_out}) != requested ({qc_module_in})"
    if qc_module_xptd:
        assert qc_module_out == qc_module_xptd, f"QC_MODULE used ({qc_module_out}) != expected ({qc_module_xptd})"

    # aliases checks
    if is_dft:
        _asserter(asserter_args, contractual_args, contractual_dft_current)
    else:
        _asserter(asserter_args, contractual_args, contractual_current)

    # returns checks
    if driver == "energy":
        assert compare_values(
            ref_block[f"{method.upper()} TOTAL ENERGY"], wfn["return_result"], tnm + " wfn", atol=atol_e, rtol=rtol_e
        )
        assert compare_values(
            ref_block[f"{method.upper()} TOTAL ENERGY"], wfn["properties"]["return_energy"], tnm + " prop", atol=atol_e, rtol=rtol_e
        )
        assert compare_values(ref_block[f"{method.upper()} TOTAL ENERGY"], ret, tnm + " return")

    elif driver == "gradient":
        assert compare_values(
            ref_block[f"{method.upper()} TOTAL GRADIENT"], wfn["return_result"], tnm + " grad wfn", atol=atol_g, rtol=rtol_g
        )
        assert compare_values(
            ref_block[f"{method.upper()} TOTAL ENERGY"], wfn["properties"]["return_energy"], tnm + " prop", atol=atol_e, rtol=rtol_e
        )
        # assert compare_values(ref_block[f"{method.upper()} TOTAL GRADIENT"], wfn["properties"]["return_gradient"], tnm + " grad prop", atol=atol_g)
        assert compare_values(ref_block[f"{method.upper()} TOTAL GRADIENT"], ret, tnm + " grad return", atol=atol_g, rtol=rtol_g)

    elif driver == "hessian":
        assert compare_values(
            ref_block[f"{method.upper()} TOTAL HESSIAN"], wfn["return_result"], tnm + " hess wfn", atol=atol_h, rtol=rtol_h
        )
        # assert compare_values(ref_block[f"{method.upper()} TOTAL GRADIENT"], wfn["properties"]["return_gradient"], tnm + " grad prop", atol=atol_g)
        assert compare_values(
            ref_block[f"{method.upper()} TOTAL ENERGY"], wfn["properties"]["return_energy"], tnm + " prop", atol=atol_e, rtol=rtol_e
        )
        assert compare_values(ref_block[f"{method.upper()} TOTAL HESSIAN"], ret, tnm + " hess return", atol=atol_h, rtol=rtol_h)

    # generics checks
    # yapf: disable
    assert compare(ref_block["N BASIS FUNCTIONS"], wfn["properties"]["calcinfo_nbasis"], tnm + " nbasis wfn"), f"nbasis {wfn.properties.calcinfo_nbasis} != {ref_block['N BASIS FUNCTIONS']}"
    assert compare(ref_block["N MOLECULAR ORBITALS"], wfn["properties"]["calcinfo_nmo"], tnm + " nmo wfn"), f"nmo {wfn.properties.calcinfo_nmo} != {ref_block['N MOLECULAR ORBITALS']}"
    assert compare(ref_block["N ALPHA ELECTRONS"], wfn["properties"]["calcinfo_nalpha"], tnm + " nalpha wfn"), f"nalpha {wfn.properties.calcinfo_nalpha} != {ref_block['N ALPHA ELECTRONS']}"
    assert compare(ref_block["N BETA ELECTRONS"], wfn["properties"]["calcinfo_nbeta"], tnm + " nbeta wfn"), f"nbeta {wfn.properties.calcinfo_nbeta} != {ref_block['N BETA ELECTRONS']}"
    # yapf: enable

    # record
    _recorder(qcprog, qc_module_out, driver, method, reference, fcae, scf_type, corl_type, "fd" if using_fd else "pass", "")

    # assert 0


def _asserter(asserter_args, contractual_args, contractual_fn):
    """For expectations in `contractual_fn`, check that the QCVars are present in P::e.globals and wfn and match expected ref_block."""

    qcvar_stores, ref_block, atol_egh, rtol_egh, ref_block_conv, atol_conv, rtol_conv, tnm = asserter_args

    for obj in qcvar_stores:
        for rpv, pv, present in contractual_fn(*contractual_args):
            label = tnm + " " + pv
            atol = atol_egh["EGH".index(rpv.split()[-1][0])]
            rtol = rtol_egh["EGH".index(rpv.split()[-1][0])]

            if present:
                # verify exact match to method (may be df) and near match to conventional (non-df) method
                tf, errmsg = compare_values(
                    ref_block[rpv], query_qcvar(obj, pv), label, atol=atol, rtol=rtol, return_message=True, quiet=True
                )
                assert compare_values(ref_block[rpv], query_qcvar(obj, pv), label, atol=atol, rtol=rtol), errmsg
                tf, errmsg = compare_values(
                    ref_block_conv[rpv], query_qcvar(obj, pv), label, atol=atol_conv, rtol=rtol_conv, return_message=True, quiet=True
                )
                assert compare_values(ref_block_conv[rpv], query_qcvar(obj, pv), label, atol=atol_conv, rtol=rtol_conv), errmsg

                # Note that the double compare_values lines are to collect the errmsg in the first for assertion in the second.
                #   If the errmsg isn't present in the assert, the string isn't accessible through `e.value`.
                #   If a plain bool is compared in the assert, the printed message will show booleans and not numbers.
            else:
                # verify and forgive known contract violations
                assert compare(False, query_has_qcvar(obj, pv), label + " SKIP"), f"{label} wrongly present"


def _recorder(engine, module, driver, method, reference, fcae, scf_type, corl_type, status, note):
    with open("stdsuite_qcng.txt", "a") as fp:
        stuff = {"module": module, "driver": driver, "method": method, "reference": reference, "fcae": fcae, "scf_type": scf_type, "corl_type": corl_type, "status": status, "note": note}
        fp.write(f"{stuff!r}\n")
