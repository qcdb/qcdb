import copy
import inspect
import sys
import traceback
from typing import Any, Dict, Optional

import qcelemental as qcel
import qcengine as qcng
from qcelemental.models import AtomicInput
from qcengine.exceptions import InputError
from qcengine.programs.cfour import CFOURHarness
from qcengine.programs.cfour.keywords import format_keyword, format_keywords
from qcengine.programs.util import PreservingDict

from ... import qcvars
from ...basisset import BasisSet
from ...util import accession_stamp, print_jobrec, provenance_stamp
from .germinate import (
    extract_basis_from_genbas,
    muster_basisset,
    muster_inherited_keywords,
    muster_modelchem,
    muster_molecule,
)


def run_cfour(name: str, molecule: "Molecule", options: "Keywords", **kwargs) -> Dict:
    """QCDB API to QCEngine connection for CFOUR."""

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

    jobrec = qcng.compute(resi, "qcdb-cfour", local_options=local_options, raise_error=True).dict()

    hold_qcvars = jobrec["extras"].pop("qcdb:qcvars")
    jobrec["qcvars"] = {key: qcel.Datum(**dval) for key, dval in hold_qcvars.items()}
    jobrec["molecule"]["fix_com"] = molecule.com_fixed()
    jobrec["molecule"]["fix_orientation"] = molecule.orientation_fixed()

    return jobrec


class QcdbCFOURHarness(CFOURHarness):
    def compute(self, input_model: AtomicInput, config: "JobConfig") -> "AtomicResult":
        self.found(raise_error=True)

        verbose = 1

        print_jobrec(f"[1] {self.name} RESULTINPUT PRE-PLANT", input_model.dict(), verbose >= 3)

        job_inputs = self.qcdb_build_input(input_model, config)

        print_jobrec(f"[2] {self.name}REC PRE-ENGINE", job_inputs, verbose >= 4)

        success, dexe = self.execute(job_inputs)

        print_jobrec(f"[3] {self.name}REC POST-ENGINE", dexe, verbose >= 4)

        if not success:
            output_model = input_model
            output_model["error"] = {"error_type": "execution_error", "error_message": dexe["stderr"]}

        dexe["outfiles"]["stdout"] = dexe["stdout"]
        dexe["outfiles"]["stderr"] = dexe["stderr"]
        output_model = self.parse_output(dexe["outfiles"], input_model)

        print_jobrec(f"[4a] {self.name} RESULT POST-HARVEST", output_model.dict(), verbose >= 5)

        output_model = self.qcdb_post_parse_output(input_model, output_model)

        print_jobrec(f"[4] {self.name} RESULT POST-POST-HARVEST", output_model.dict(), verbose >= 2)

        return output_model

    def qcdb_build_input(
        self, input_model: AtomicInput, config: "JobConfig", template: Optional[str] = None
    ) -> Dict[str, Any]:
        cfourrec = {
            "infiles": {},
            "scratch_messy": config.scratch_messy,
            "scratch_directory": config.scratch_directory,
        }

        kwgs = {"accession": accession_stamp(), "verbose": 1}
        molrec = qcel.molparse.from_schema(input_model.molecule.dict())
        molrec["fix_symmetry"] = "c1"  # write _all_ atoms to GENBAS
        ropts = input_model.extras["qcdb:options"]

        ropts.require("QCDB", "MEMORY", f"{config.memory} gib", **kwgs)

        molcmd = muster_molecule(molrec, ropts, verbose=1)

        # Handle qcdb keywords implying cfour keyword values
        muster_inherited_keywords(ropts)

        _qcdb_basis = ropts.scroll["QCDB"]["BASIS"].value
        _cfour_basis = ropts.scroll["CFOUR"]["BASIS"].value
        if (
            input_model.model.basis != "(auto)"
            and not ropts.scroll["QCDB"]["BASIS"].disputed()
            and not ropts.scroll["CFOUR"]["BASIS"].disputed()
        ):
            qbs = BasisSet.pyconstruct(molrec, "BASIS", input_model.model.basis)
            cfourrec["infiles"]["GENBAS"] = qbs.print_detail_cfour()
            bascmd = muster_basisset(molrec, ropts, qbs.has_puream())
        elif _qcdb_basis == "":
            _, cased_basis = format_keyword("CFOUR_BASIS", _cfour_basis)
            cfourrec["infiles"]["GENBAS"] = extract_basis_from_genbas(
                cased_basis, input_model.molecule.symbols, exact=False
            )
            bascmd = ""
        else:
            qbs = BasisSet.pyconstruct(molrec, "BASIS", _qcdb_basis)
            # if qbs.has_ECP(): #    raise ValidationError("""ECPs not hooked up for Cfour""")
            cfourrec["infiles"]["GENBAS"] = qbs.print_detail_cfour()  # qbs.genbas()
            bascmd = muster_basisset(molrec, ropts, qbs.has_puream())

        # Handle calc type and quantum chemical method
        muster_modelchem(input_model.model.method, input_model.driver.derivative_int(), ropts)

        # print(jobrec['options'].print_changed(history=False))
        # Handle driver vs input/default keyword reconciliation

        # print("Touched Keywords")  # debug
        # print(ropts.print_changed(history=True))  # debug

        # Handle conversion of psi4 keyword structure into cfour format
        skma_options = {key: ropt.value for key, ropt in sorted(ropts.scroll["CFOUR"].items()) if ropt.disputed()}
        optcmd = format_keywords(skma_options)

        # Assemble ZMAT pieces
        cfourrec["infiles"]["ZMAT"] = molcmd + optcmd + bascmd
        cfourrec["command"] = ["xcfour"]

        return cfourrec

    def qcdb_post_parse_output(self, input_model: AtomicInput, output_model: "AtomicResult") -> "AtomicResult":

        dqcvars = PreservingDict(copy.deepcopy(output_model.extras["qcvars"]))

        # badly placed
        # Cfour's SCS-MP2 is non adjustible and only valid for UHF
        # ROMP2 doesn't print SS & OS
        #        if "MP2 OPPOSITE-SPIN CORRELATION ENERGY" in dqcvars and "MP2 SAME-SPIN CORRELATION ENERGY" in dqcvars:
        #            oss_opt = input_model.extras['qcdb:options'].scroll['QCDB']['MP2_OS_SCALE']
        #            sss_opt = input_model.extras['qcdb:options'].scroll['QCDB']['MP2_SS_SCALE']
        #            custom_scsmp2_corl = \
        #                Decimal(oss_opt.value) * dqcvars["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] + \
        #                Decimal(sss_opt.value) * dqcvars["MP2 SAME-SPIN CORRELATION ENERGY"]
        #            if "MP2 SINGLES ENERGY" in dqcvars:
        #                custom_scsmp2_corl += dqcvars["MP2 SINGLES ENERGY"]
        #            dqcvars["CUSTOM SCS-MP2 CORRELATION ENERGY"] = custom_scsmp2_corl

        try:
            qcvars.build_out(dqcvars)
        except ValueError:
            raise InputError(
                "STDOUT:\n"
                + output_model.stdout
                + "\nTRACEBACK:\n"
                + "".join(traceback.format_exception(*sys.exc_info()))
            )
        calcinfo = qcvars.certify_and_datumize(dqcvars, plump=True, nat=len(output_model.molecule.symbols))
        output_model.extras["qcdb:qcvars"] = calcinfo

        return output_model


#    # find environment by merging PSIPATH and PATH environment variables
#    # * `path` kwarg gets precedence
#    # * filter out None values as subprocess will fault on them
#    lenv = {
#        'HOME': os.environ.get('HOME'),
#        'PATH': (':'.join([os.path.abspath(x) for x in os.environ.get('PSIPATH', '').split(':') if x != '']) +
#                 ':' + os.environ.get('PATH')),# +
##                 ':' + qcdb.get_datadir() + '/basis'),
##        'GENBAS_PATH': qcdb.get_datadir() + '/basis',
#        'CFOUR_NUM_CORES': os.environ.get('CFOUR_NUM_CORES'),
#        'MKL_NUM_THREADS': os.environ.get('MKL_NUM_THREADS'),
#        'OMP_NUM_THREADS': os.environ.get('OMP_NUM_THREADS'),
#        'LD_LIBRARY_PATH': os.environ.get('LD_LIBRARY_PATH')
#        }
#    if 'executable_path' in cfourrec:
#        lenv['PATH'] = cfourrec['executable_path'] + ':' + lenv['PATH']

# print('HH')
# print(jobrec['options'].print_changed(history=False))

#    tmp_zmat = """UHF-SCF energy calculation
# N
# H 1 R
# H 1 R 2 A
#
# R=1.008
# A=105.0"""
#    tmp_zmat = """auto-generated by qcdb from molecule H2N
# N                     0.000000000000     0.000000000000    -0.145912918812
# H                     0.000000000000    -1.511214304079     1.013682601327
# H                     0.000000000000     1.511214304079     1.013682601327
#
# *ACES2(CALC=HF,BASIS=qz2p
# MULT=2,REF=UHF
# COORDINATES=CARTESIAN
# UNITS=BOHR
# OCCUPATION=3-1-1-0/3-0-1-0
# SCF_CONV=12
# MEMORY=20000000)
#
# """
