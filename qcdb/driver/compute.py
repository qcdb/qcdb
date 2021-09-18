import sys
import traceback
from typing import Any, Dict, Optional, Union

import qcelemental as qcel
import qcengine as qcng

from ..util import program_prefix
from .driver_util import pkgprefix

prefixpkg = {v: k for k, v in pkgprefix.items()}


def compute(
    input_data: Union[Dict[str, Any], "AtomicInput"],
    program: str,
    # raise_error: bool = False,
    local_options: Optional[Dict[str, Any]] = None,
    mode_options: Optional[Dict[str, str]] = None,
    # return_dict: bool = False,
) -> "AtomicResult":
    """Run an analytic single-point specified in ``input_data`` through program ``qcprog``."""

    try:
        input_model = qcng.util.model_wrapper(input_data, qcel.models.AtomicInput)

        ## Echo the infile on the outfile
        # core.print_out("\n  ==> Input QCSchema <==\n")
        # core.print_out("\n--------------------------------------------------------------------------\n")
        # core.print_out(pp.pformat(json.loads(input_model.json())))
        # core.print_out("\n--------------------------------------------------------------------------\n")

        # keep_wfn = input_model.protocols.wavefunction != 'none'

        import qcdb

        driver_call = {"energy": qcdb.energy, "gradient": qcdb.gradient, "hessian": qcdb.hessian}
        qmol = qcdb.Molecule.from_schema(input_model.molecule.dict())

        mode_config = qcdb.driver.config.get_mode_config(mode_options=mode_options)
        if mode_config.implicit_program == "qcdb":
            keywords = input_model.keywords
            mtd_call = prefixpkg[program.lower()] + input_model.model.method + "/" + input_model.model.basis
        else:
            # check consis mode_options and program?
            keywords = {(program + "_" + k): v for k, v in input_model.keywords.items()}
            mtd_call = program_prefix(program) + input_model.model.method + "/" + input_model.model.basis

        retres, ret = driver_call[input_model.driver](
            mtd_call,
            molecule=qmol,
            options=keywords,
            return_wfn=True,
            local_options=local_options,
            mode_options=mode_options,
        )

        ## qcschema should be copied
        # ret_data = run_json_qcschema(input_model.dict(), clean, False, keep_wfn=keep_wfn)
        # ret_data["provenance"].update({
        #    "creator": "Psi4",
        #    "version": __version__,
        #    "routine": "psi4.schema_runner.run_qcschema"
        # })

        # exit_printing(start_time=start_time, success=True)

        # import pprint
        # pprint.pprint(ret)
        # ret = atres

        ret.pop("qcvars")  # duplicated in extras["qcvars"]
        ret = qcel.models.AtomicResult(**ret)

    except Exception as exc:

        if not isinstance(input_data, dict):
            input_data = input_data.dict()

        input_data = input_data.copy()
        ret = qcel.models.FailedOperation(
            input_data=input_data,
            success=False,
            error={
                "error_type": type(exc).__name__,
                "error_message": "".join(traceback.format_exception(*sys.exc_info())),
            },
        )

    return ret
