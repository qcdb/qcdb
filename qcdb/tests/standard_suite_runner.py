import pprint

import pytest
import numpy as np
from qcengine.programs.tests.standard_suite_contracts import (
    contractual_hf,
    contractual_mp2,
    contractual_mp2p5,
    contractual_mp3,
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
    contractual_current,
    query_has_qcvar,
    query_qcvar,
)
from qcengine.programs.tests.standard_suite_ref import answer_hash, std_suite

from .utils import compare, compare_values

pp = pprint.PrettyPrinter(width=120)


def runner_asserter(inp, subject, method, basis, tnm):

    qcprog = inp["qc_module"].split("-")[0]
    qc_module_in = inp["qc_module"]  # returns "<qcprog>"|"<qcprog>-<module>"  # input-specified routing
    qc_module_xptd = (
        (qcprog + "-" + inp["xptd"]["qc_module"]) if inp.get("xptd", {}).get("qc_module", None) else None
    )  # expected routing
    driver = inp["driver"]
    reference = inp["reference"]
    fcae = inp["fcae"]

    # <<<  Reference Values  >>>

    # ? precedence on next two
    mp2_type = inp.get("corl_type", inp["keywords"].get("mp2_type", "df"))  # hard-code of read_options.cc MP2_TYPE
    mp_type = inp.get("corl_type", inp["keywords"].get("mp_type", "conv"))  # hard-code of read_options.cc MP_TYPE
    cc_type = inp.get("corl_type", inp["keywords"].get("cc_type", "conv"))  # hard-code of read_options.cc CC_TYPE
    corl_natural_values = {
        "hf": "conv",  # dummy to assure df/cd/conv scf_type refs available
        "mp2": mp2_type,
        "mp3": mp_type,
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
    }
    corl_type = corl_natural_values[method]

    natural_ref = {"conv": "pk", "df": "df", "cd": "cd"}
    scf_type = inp["keywords"].get("scf_type", natural_ref[corl_type])
    natural_values = {"pk": "pk", "direct": "pk", "df": "df", "mem_df": "df", "disk_df": "df", "cd": "cd"}
    scf_type = natural_values[scf_type]

    atol_e, atol_g, atol_h = 1.e-6, 2.e-6, 5.e-6
    if inp.get("xptd", {}).get("fd", False):
        # relax atol for qcprog internal findif since 3-pt is crude and arrays less precise
        atol_g, atol_h = 2.e-5, 8.e-5
    chash = answer_hash(
        system=subject.name(),
        basis=basis,
        fcae=fcae,
        scf_type=scf_type,
        reference=reference,
        corl_type=corl_type,
    )

    # check all calcs against conventional reference to looser tolerance
    atol_conv = 1.0e-4
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
    local_options = {"nnodes": 1, "ncores": 1, "scratch_messy": False}

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
        errtype, errmsg = inp["error"]
        with pytest.raises(errtype) as e:
            driver_call[driver](inp["call"], molecule=subject, local_options=local_options, **extra_kwargs)

        assert errmsg in str(e.value), f"Not found: {errtype} '{errmsg}' in {e.value}"
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

    # <<<  Comparison Tests  >>>

    assert wfn["success"] is True
    assert (
        wfn["provenance"]["creator"].lower() == qcprog
    ), f'ENGINE used ({ wfn["provenance"]["creator"].lower()}) != requested ({qcprog})'
    if qc_module_in != qcprog:
        assert qc_module_out == qc_module_in, f"QC_MODULE used ({qc_module_out}) != requested ({qc_module_in})"
    if qc_module_xptd:
        assert qc_module_out == qc_module_xptd, f"QC_MODULE used ({qc_module_out}) != expected ({qc_module_xptd})"

    ref_block = std_suite[chash]

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
        ref_block_conv,
        atol_conv,
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
        elif method == "ccd":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_ccd)
        elif method == "ccsd":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_ccsd)
#        elif method == "ccsd(t)":
#            _asserter(asserter_args, contractual_args, contractual_mp2)
#            _asserter(asserter_args, contractual_args, contractual_ccsd)
#            _asserter(asserter_args, contractual_args, contractual_ccsd_prt_pr)
        elif method == "ccsd+t(ccsd)":
            _asserter(asserter_args, contractual_args, contractual_mp2)
            _asserter(asserter_args, contractual_args, contractual_ccsd)
            _asserter(asserter_args, contractual_args, contractual_ccsdpt_prccsd_pr)
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

    if "wrong" in inp:
        errmsg, reason = inp["wrong"]
        with pytest.raises(AssertionError) as e:
            qcvar_assertions()

        # print("WRONG", errmsg, reason, str(e.value), "ENDW")
        assert errmsg in str(e.value)
        pytest.xfail(reason)

    qcvar_assertions()

    # aliases
    _asserter(asserter_args, contractual_args, contractual_current)

    # returns
    if driver == "energy":
        assert compare_values(
            ref_block[f"{method.upper()} TOTAL ENERGY"], wfn["return_result"], tnm + " wfn", atol=atol_e
        )
        assert compare_values(
            ref_block[f"{method.upper()} TOTAL ENERGY"], wfn["properties"]["return_energy"], tnm + " prop", atol=atol_e
        )
        assert compare_values(ref_block[f"{method.upper()} TOTAL ENERGY"], ret, tnm + " return")

    elif driver == "gradient":
        assert compare_values(
            ref_block[f"{method.upper()} TOTAL GRADIENT"], wfn["return_result"], tnm + " grad wfn", atol=atol_g
        )
        assert compare_values(
            ref_block[f"{method.upper()} TOTAL ENERGY"], wfn["properties"]["return_energy"], tnm + " prop", atol=atol_e
        )
        # assert compare_values(ref_block[f"{method.upper()} TOTAL GRADIENT"], wfn["properties"]["return_gradient"], tnm + " grad prop", atol=atol_g)
        assert compare_values(ref_block[f"{method.upper()} TOTAL GRADIENT"], ret, tnm + " grad return", atol=atol_g)

    elif driver == "hessian":
        assert compare_values(
            ref_block[f"{method.upper()} TOTAL HESSIAN"], wfn["return_result"], tnm + " hess wfn", atol=atol_h
        )
        # assert compare_values(ref_block[f"{method.upper()} TOTAL GRADIENT"], wfn["properties"]["return_gradient"], tnm + " grad prop", atol=atol_g)
        assert compare_values(
            ref_block[f"{method.upper()} TOTAL ENERGY"], wfn["properties"]["return_energy"], tnm + " prop", atol=atol_e
        )
        assert compare_values(ref_block[f"{method.upper()} TOTAL HESSIAN"], ret, tnm + " hess return", atol=atol_h)

    # generics
    # yapf: disable
    assert compare(ref_block["N BASIS FUNCTIONS"], wfn["properties"]["calcinfo_nbasis"], tnm + " nbasis wfn"), f"nbasis {wfn.properties.calcinfo_nbasis} != {ref_block['N BASIS FUNCTIONS']}"
    assert compare(ref_block["N MOLECULAR ORBITALS"], wfn["properties"]["calcinfo_nmo"], tnm + " nmo wfn"), f"nmo {wfn.properties.calcinfo_nmo} != {ref_block['N MOLECULAR ORBITALS']}"
    assert compare(ref_block["N ALPHA ELECTRONS"], wfn["properties"]["calcinfo_nalpha"], tnm + " nalpha wfn"), f"nalpha {wfn.properties.calcinfo_nalpha} != {ref_block['N ALPHA ELECTRONS']}"
    assert compare(ref_block["N BETA ELECTRONS"], wfn["properties"]["calcinfo_nbeta"], tnm + " nbeta wfn"), f"nbeta {wfn.properties.calcinfo_nbeta} != {ref_block['N BETA ELECTRONS']}"
    # yapf: enable


def _asserter(asserter_args, contractual_args, contractual_fn):
    """For expectations in `contractual_fn`, check that the QCVars are present in P::e.globals and wfn and match expected ref_block."""

    qcvar_stores, ref_block, atol_egh, ref_block_conv, atol_conv, tnm = asserter_args

    for obj in qcvar_stores:
        for rpv, pv, present in contractual_fn(*contractual_args):
            label = tnm + " " + pv
            atol = atol_egh["EGH".index(rpv.split()[-1][0])]

            if present:
                # verify exact match to method (may be df) and near match to conventional (non-df) method
                tf, errmsg = compare_values(
                    ref_block[rpv], query_qcvar(obj, pv), label, atol=atol, return_message=True, quiet=True
                )
                assert compare_values(ref_block[rpv], query_qcvar(obj, pv), label, atol=atol), errmsg
                tf, errmsg = compare_values(
                    ref_block_conv[rpv], query_qcvar(obj, pv), label, atol=atol_conv, return_message=True, quiet=True
                )
                assert compare_values(ref_block_conv[rpv], query_qcvar(obj, pv), label, atol=atol_conv), errmsg

                # Note that the double compare_values lines are to collect the errmsg in the first for assertion in the second.
                #   If the errmsg isn't present in the assert, the string isn't accessible through `e.value`.
                #   If a plain bool is compared in the assert, the printed message will show booleans and not numbers.
            else:
                # verify and forgive known contract violations
                assert compare(False, query_has_qcvar(obj, pv), label + " SKIP"), f"{label} wrongly present"
