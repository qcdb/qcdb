import copy
import pprint
import inspect
from typing import Any, Dict, Optional

import qcelemental as qcel
from qcelemental.models import AtomicInput
from qcelemental.util import which

import qcengine as qcng
from qcengine.exceptions import InputError
from qcengine.programs.gamess import GAMESSHarness
from qcengine.programs.gamess.keywords import format_keywords
from qcengine.programs.util import PreservingDict

from ... import qcvars
from ...basisset import BasisSet
from ...util import print_jobrec, provenance_stamp
from .germinate import muster_inherited_keywords, muster_modelchem, muster_molecule_and_basisset

pp = pprint.PrettyPrinter(width=120)


def run_gamess(name: str, molecule: 'Molecule', options: 'Keywords', **kwargs) -> Dict:

    resi = AtomicInput(
        **{
            'driver': inspect.stack()[1][3],
            'extras': {
                'qcdb:options': copy.deepcopy(options),
            },
            'model': {
                'method': name,
                'basis': '(auto)',
            },
            'molecule': molecule.to_schema(dtype=2),
            'provenance': provenance_stamp(__name__),
        })

    jobrec = qcng.compute(resi, "qcdb-gamess", raise_error=True).dict()
    hold_qcvars = jobrec['extras'].pop('qcdb:qcvars')
    jobrec['qcvars'] = {key: qcel.Datum(**dval) for key, dval in hold_qcvars.items()}
    return jobrec


class QcdbGAMESSHarness(GAMESSHarness):
    def compute(self, input_model: AtomicInput, config: 'JobConfig') -> 'AtomicResult':
        self.found(raise_error=True)

        verbose = 1

        print_jobrec(f'[1] {self.name} RESULTINPUT PRE-PLANT', input_model.dict(), verbose >= 3)

        job_inputs = self.qcdb_build_input(input_model, config)

        print_jobrec(f'[2] {self.name}REC PRE-ENGINE', job_inputs, verbose >= 4)

        success, dexe = self.execute(job_inputs)

        print_jobrec(f'[3] {self.name}REC POST-ENGINE', dexe, verbose >= 4)

        if 'INPUT HAS AT LEAST ONE SPELLING OR LOGIC MISTAKE' in dexe["stdout"]:
            raise InputError(dexe["stdout"])

        if not success:
            output_model = input_model
            output_model["error"] = {"error_type": "execution_error", "error_message": dexe["stderr"]}

        dexe["outfiles"]["stdout"] = dexe["stdout"]
        dexe["outfiles"]["stderr"] = dexe["stderr"]
        output_model = self.parse_output(dexe["outfiles"], input_model)

        print_jobrec(f'[4a] {self.name} RESULT POST-HARVEST', output_model.dict(), verbose >= 5)

        output_model = self.qcdb_post_parse_output(input_model, output_model)

        print_jobrec(f'[4] {self.name} RESULT POST-POST-HARVEST', output_model.dict(), verbose >= 2)

        return output_model

    def qcdb_build_input(self, input_model: AtomicInput, config: 'JobConfig',
                         template: Optional[str] = None) -> Dict[str, Any]:
        gamessrec = {
            'infiles': {},
            'scratch_directory': config.scratch_directory,
        }

        molrec = qcel.molparse.from_schema(input_model.molecule.dict())
        molrecc1 = molrec.copy()
        molrecc1['fix_symmetry'] = 'c1'  # write _all_ atoms to input
        ropts = input_model.extras['qcdb:options']

        # Handle qcdb keywords implying gamess keyword values
        muster_inherited_keywords(ropts)

        _qcdb_basis = ropts.scroll['QCDB']['BASIS'].value
        #_gamess_basis = ropts.scroll['GAMESS']['BASIS'].value
        qbs = BasisSet.pyconstruct(molrecc1, 'BASIS', _qcdb_basis)

        molbascmd = muster_molecule_and_basisset(molrec, ropts, qbs)

        nel = sum([z * int(real) for z, real in zip(molrec['elez'], molrec['real'])]) - molrec['molecular_charge']
        nel = int(nel)
        #nfzc = input_model.molecule.nfrozen_core(depth=True)  # this will be default FC  # TODO change these values when user sets custom FC
        # not sure how to do this???
        nfzc = 0

        # forcing nfc above. all these need to be reocmputed together for a consistent cidet input group
        nels = nel - 2 * nfzc
        nact = qbs.nbf() - nfzc
        sysinfo = {
            'nel': nel,
            'ncore': nfzc,
            'nact': nact,
            'nels': nels,
        }

        # Handle calc type and quantum chemical method
        muster_modelchem(input_model.model.method, input_model.driver.derivative_int(), ropts, sysinfo)

        # Handle conversion of psi4 keyword structure into cfour format
        skma_options = {key: ropt.value for key, ropt in sorted(ropts.scroll['GAMESS'].items()) if ropt.disputed()}

        optcmd = format_keywords(skma_options)

        gamessrec['infiles']['gamess.inp'] = optcmd + molbascmd
        gamessrec['command'] = [which("rungms"), "gamess"]

        return gamessrec

    def qcdb_post_parse_output(self, input_model: AtomicInput, output_model: 'AtomicResult') -> 'AtomicResult':

        dqcvars = PreservingDict(copy.deepcopy(output_model.extras['qcvars']))
        qcvars.build_out(dqcvars)
        calcinfo = qcvars.certify_and_datumize(dqcvars, plump=True, nat=len(output_model.molecule.symbols))
        output_model.extras['qcdb:qcvars'] = calcinfo

        return output_model
