import copy
import inspect
import pprint
from typing import Any, Dict, Optional

import qcelemental as qcel
import qcengine as qcng
from qcelemental.models import AtomicInput
from qcengine.programs.psi4 import Psi4Harness
from qcengine.programs.util import PreservingDict

from ... import qcvars
from ...driver.config import get_mode_config
from ...util import print_jobrec, provenance_stamp
from .germinate import muster_inherited_keywords

pp = pprint.PrettyPrinter(width=200)


def run_psi4(name: str, molecule: "Molecule", options: "Keywords", **kwargs) -> Dict:

    local_options = kwargs.get("local_options", None)
    mode_options = get_mode_config(mode_options=kwargs.get("mode_options"))

    resi = AtomicInput(
        **{
            "driver": inspect.stack()[1][3],
            "extras": {
                "qcdb:options": copy.deepcopy(options),
                "qcdb:mode_config": mode_options,
            },
            "model": {
                "method": name,
                "basis": "(auto)",
            },
            "molecule": molecule.to_schema(dtype=2),
            "provenance": provenance_stamp(__name__),
        }
    )

    jobrec = qcng.compute(resi, "qcdb-psi4", local_options=local_options, raise_error=True).dict()
    hold_qcvars = jobrec["extras"].pop("qcdb:qcvars")
    jobrec["qcvars"] = {key: qcel.Datum(**dval) for key, dval in hold_qcvars.items()}

    return jobrec


class QcdbPsi4Harness(Psi4Harness):
    def compute(self, input_model: AtomicInput, config: "JobConfig") -> "AtomicResult":
        self.found(raise_error=True)

        verbose = 1
        print_jobrec(f"[1] {self.name} RESULTINPUT PRE-PLANT", input_model.dict(), verbose >= 3)

        input_data = self.qcdb_build_input(input_model, config)
        input_model = AtomicInput(**input_data)

        print_jobrec(f"[2] {self.name} RESULTINPUT PRE-ENGINE", input_model.dict(), verbose >= 4)

        # 'PATH': (':'.join([os.path.abspath(x) for x in os.environ.get('PSIPATH', '').split(':') if x != '']) +
        #          ':' + os.environ.get('PATH')),# +
        # 'PSI_SCRATCH': tmpdir,
        # 'PYTHONPATH': os.environ.get('PYTHONPATH'),
        # 'LD_LIBRARY_PATH': os.environ.get('LD_LIBRARY_PATH')

        output_model = Psi4Harness.compute(self, input_model=input_model, config=config)

        print_jobrec(f"[3] {self.name} RESULT POST-ENGINE", output_model.dict(), verbose >= 4)

        # ???
        if not output_model.success:
            return output_model

        print_jobrec(f"[4a] {self.name} RESULT POST-HARVEST", output_model.dict(), verbose >= 5)

        output_model = self.qcdb_post_parse_output(input_model, output_model)

        print_jobrec(f"[4] {self.name} RESULT POST-POST-HARVEST", output_model.dict(), verbose >= 2)

        return output_model

    def qcdb_build_input(
        self, input_model: AtomicInput, config: "JobConfig", template: Optional[str] = None
    ) -> Dict[str, Any]:
        input_data = input_model.dict()

        ropts = input_model.extras["qcdb:options"]
        mode_config = input_model.extras["qcdb:mode_config"]

        ropts.require("QCDB", "MEMORY", f"{config.memory} gib", accession="00000000", verbose=False)

        muster_inherited_keywords(ropts, mode_config)
        mtd = input_data["model"]["method"]
        mtd = mtd[3:] if mtd.startswith("p4-") else mtd
        input_data["model"]["method"] = mtd

        # should we put this memory in the JobConfig object? I don't think the units agree
        ropts.scroll["QCDB"].pop("MEMORY")
        # print(config.memory, '!!')
        # config.memory = omem.value #???
        # print(config.memory, '!!')

        input_data["extras"] = {"wfn_qcvars_only": True}
        # input_data['kwargs'] = jobrec['kwargs']
        # input_data['return_output'] = True

        # print("Touched Keywords")  # debug
        # print(ropts.print_changed(history=True))  # debug

        popts = {}
        function_kwargs = {}
        # was recently active
        # for k, v in ropts.scroll["QCDB"].items():
        #     if v.disputed():
        #         popts[k] = v.value

        for k, v in ropts.scroll["PSI4"].items():
            if v.disputed2():
                if k.startswith("FUNCTION_KWARGS_"):
                    function_kwargs[k[16:]] = v.value
                else:
                    popts[k] = v.value
        input_data["keywords"] = popts
        input_data["keywords"]["function_kwargs"] = function_kwargs

        # print("Collected Keywords")  # debug
        # pp.pprint(popts)  # debug

        if "BASIS" in input_data["keywords"]:
            input_data["model"]["basis"] = input_data["keywords"]["BASIS"]

        return input_data

    def qcdb_post_parse_output(self, input_model: AtomicInput, output_model: "AtomicResult") -> "AtomicResult":

        dqcvars = PreservingDict(copy.deepcopy(output_model.extras["qcvars"]))
        for k in list(dqcvars.keys()):
            if k in ["DETCI AVG DVEC NORM", "MCSCF TOTAL ENERGY"]:
                dqcvars.pop(k)
        qcvars.build_out(dqcvars)
        calcinfo = qcvars.certify_and_datumize(dqcvars, plump=True, nat=len(output_model.molecule.symbols))
        output_model.extras["qcdb:qcvars"] = calcinfo

        return output_model
