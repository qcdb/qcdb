import copy
import inspect
import pprint
from typing import Any, Dict, Optional

import qcelemental as qcel
import qcengine as qcng
from qcelemental.models import AtomicInput
from qcelemental.util import which
from qcengine.exceptions import InputError
from qcengine.programs.gamess import GAMESSHarness
from qcengine.programs.gamess.keywords import format_keywords
from qcengine.programs.util import PreservingDict

from ... import qcvars
from ...basisset import BasisSet
from ...molecule import Molecule
from ...util import accession_stamp, print_jobrec, provenance_stamp
from .germinate import get_master_frame, muster_inherited_keywords, muster_modelchem, muster_molecule_and_basisset

pp = pprint.PrettyPrinter(width=120)


def run_gamess(name: str, molecule: "Molecule", options: "Keywords", **kwargs) -> Dict:

    local_options = kwargs.get("local_options", None)

    resi = AtomicInput(
        **{
            "driver": inspect.stack()[1][3],
            "extras": {
                "qcdb:options": copy.deepcopy(options),
            },
            "model": {
                "method": name,
                "basis": "(auto)",
            },
            "molecule": molecule.to_schema(dtype=2) | {"fix_com": True, "fix_orientation": True},
            "provenance": provenance_stamp(__name__),
        }
    )

    jobrec = qcng.compute(resi, "qcdb-gamess", local_options=local_options, raise_error=True).dict()

    hold_qcvars = jobrec["extras"].pop("qcdb:qcvars")
    jobrec["qcvars"] = {key: qcel.Datum(**dval) for key, dval in hold_qcvars.items()}
    jobrec["molecule"]["fix_com"] = molecule.com_fixed()
    jobrec["molecule"]["fix_orientation"] = molecule.orientation_fixed()

    return jobrec


class QcdbGAMESSHarness(GAMESSHarness):
    def compute(self, input_model: AtomicInput, config: "JobConfig") -> "AtomicResult":
        self.found(raise_error=True)

        verbose = 1

        print_jobrec(f"[1] {self.name} RESULTINPUT PRE-PLANT", input_model.dict(), verbose >= 3)

        job_inputs = self.qcdb_build_input(input_model, config)

        print_jobrec(f"[2] {self.name}REC PRE-ENGINE", job_inputs, verbose >= 4)

        success, dexe = self.execute(job_inputs)

        print_jobrec(f"[3] {self.name}REC POST-ENGINE", dexe, verbose >= 4)

        if "INPUT HAS AT LEAST ONE SPELLING OR LOGIC MISTAKE" in dexe["stdout"]:
            raise InputError(dexe["stdout"])

        if not success:
            output_model = input_model
            output_model["error"] = {"error_type": "execution_error", "error_message": dexe["stderr"]}

        dexe["outfiles"]["stdout"] = dexe["stdout"]
        dexe["outfiles"]["stderr"] = dexe["stderr"]
        dexe["outfiles"]["dsl_input"] = job_inputs["infiles"][
            "gamess.inp"
        ]  # full DSL input not available in stdout, so stash the file
        output_model = self.parse_output(dexe["outfiles"], input_model)

        print_jobrec(f"[4a] {self.name} RESULT POST-HARVEST", output_model.dict(), verbose >= 5)

        output_model = self.qcdb_post_parse_output(input_model, output_model)

        print_jobrec(f"[4] {self.name} RESULT POST-POST-HARVEST", output_model.dict(), verbose >= 2)

        return output_model

    def qcdb_build_input(
        self, input_model: AtomicInput, config: "JobConfig", template: Optional[str] = None
    ) -> Dict[str, Any]:
        gamessrec = {
            "infiles": {},
            "scratch_messy": config.scratch_messy,
            "scratch_directory": config.scratch_directory,
        }

        kwgs = {"accession": accession_stamp(), "verbose": 1}
        ropts = input_model.extras["qcdb:options"]

        if not all(input_model.molecule.real):
            raise InputError("GAMESS can't handle ghost atoms yet.")

        mf_mol, mf_data = get_master_frame(input_model.molecule, config.scratch_directory)

        # c1 so _all_ atoms written to BasisSet
        mf_qmol_c1 = Molecule.from_schema(mf_mol.dict() | {"fix_symmetry": "c1"})
        _qcdb_basis = ropts.scroll["QCDB"]["BASIS"].value
        # _gamess_basis = ropts.scroll['GAMESS']['BASIS'].value
        qbs = BasisSet.pyconstruct(mf_qmol_c1, "BASIS", _qcdb_basis)

        sysinfo = {}
        # forcing nfc above. all these need to be reocmputed together for a consistent cidet input group
        # this will be default FC  # TODO change these values when user sets custom FC
        nel = mf_mol.nelectrons()
        nfzc = mf_qmol_c1.n_frozen_core(depth=True)
        nels = nel - 2 * nfzc
        nact = qbs.nbf() - nfzc
        sysinfo["fc"] = {
            "nel": nel,
            "ncore": nfzc,
            "nact": nact,
            "nels": nels,
        }
        nfzc = 0
        nels = nel - 2 * nfzc
        nact = qbs.nbf() - nfzc
        sysinfo["ae"] = {
            "nel": nel,
            "ncore": nfzc,
            "nact": nact,
            "nels": nels,
        }

        # Handle qcdb keywords implying gamess keyword values
        muster_inherited_keywords(ropts, sysinfo)

        molbascmd = muster_molecule_and_basisset(
            mf_qmol_c1.to_dict(), qbs, ropts, mf_data["unique"], mf_data["symmetry_card"]
        )

        # Handle calc type and quantum chemical method
        muster_modelchem(input_model.model.method, input_model.driver.derivative_int(), ropts, sysinfo)

        ropts.require("QCDB", "MEMORY", f"{config.memory} gib", **kwgs)

        # Handle memory
        # * [GiB] --> [M QW]
        # * docs on mwords: "This is given in units of 1,000,000 words (as opposed to 1024*1024 words)"
        # * docs: "the memory required on each processor core for a run using p cores is therefore MEMDDI/p + MWORDS."
        # * int() rounds down
        mwords_total = int(config.memory * (1024 ** 3) / 8e6)

        asdf = "mem trials\n"
        for mem_frac_replicated in (1, 0.5, 0.1, 0.75):
            mwords, memddi = self._partition(mwords_total, mem_frac_replicated, config.ncores)
            asdf += f"loop {mwords_total=} {mem_frac_replicated=} {config.ncores=} -> repl: {mwords=} dist: {memddi=} -> percore={memddi/config.ncores + mwords} tot={memddi + config.ncores * mwords}\n"
            trial_opts = {key: ropt.value for key, ropt in sorted(ropts.scroll["GAMESS"].items()) if ropt.disputed()}
            trial_opts["contrl__exetyp"] = "check"
            trial_opts["system__parall"] = not (config.ncores == 1)
            trial_opts["system__mwords"] = mwords
            trial_opts["system__memddi"] = memddi
            trial_gamessrec = {
                "infiles": {"trial_gamess.inp": format_keywords(trial_opts) + molbascmd},
                "command": [which("rungms"), "trial_gamess", "00", str(config.ncores)],
                "scratch_messy": False,
                "scratch_directory": config.scratch_directory,
            }
            success, dexe = self.execute(trial_gamessrec)

            # this would be a lot cleaner if there was a unique or list of memory error strings
            if (
                ("ERROR: ONLY CCTYP=CCSD OR CCTYP=CCSD(T) CAN RUN IN PARALLEL." in dexe["stdout"])
                or ("ERROR: ROHF'S CCTYP MUST BE CCSD OR CR-CCL, WITH SERIAL EXECUTION" in dexe["stdout"])
                or ("CI PROGRAM CITYP=FSOCI    DOES NOT RUN IN PARALLEL." in dexe["stdout"])
            ):
                config.ncores = 1
            elif "INPUT HAS AT LEAST ONE SPELLING OR LOGIC MISTAKE" in dexe["stdout"]:
                raise InputError(dexe["stdout"])
            elif "EXECUTION OF GAMESS TERMINATED -ABNORMALLY-" in dexe["stdout"]:
                pass
            else:
                ropts.require("GAMESS", "SYSTEM__MWORDS", mwords, accession="12341234", verbose=True)
                ropts.require("GAMESS", "SYSTEM__MEMDDI", mwords, accession="12341234", verbose=True)
                asdf += f"breaking {mwords=} {memddi=}\n"
                break

        ropts.print_changed(history=True)
        # print("Touched Keywords")  # debug
        # print(ropts.print_changed(history=True))  # debug

        # Handle conversion of qcsk keyword structure into program format
        skma_options = {key: ropt.value for key, ropt in sorted(ropts.scroll["GAMESS"].items()) if ropt.disputed()}

        optcmd = format_keywords(skma_options)

        gamessrec["infiles"]["gamess.inp"] = optcmd + molbascmd
        gamessrec["command"] = [
            which("rungms"),
            "gamess",
            "00",
            str(config.ncores),
        ]  # rungms JOB VERNO NCPUS >& JOB.log &

        return gamessrec

    def qcdb_post_parse_output(self, input_model: AtomicInput, output_model: "AtomicResult") -> "AtomicResult":

        dqcvars = PreservingDict(copy.deepcopy(output_model.extras["qcvars"]))
        qcvars.build_out(dqcvars)
        calcinfo = qcvars.certify_and_datumize(dqcvars, plump=True, nat=len(output_model.molecule.symbols))
        output_model.extras["qcdb:qcvars"] = calcinfo

        return output_model
