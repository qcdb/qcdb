import copy
import pprint
import inspect
from typing import Any, Dict, Optional
from decimal import Decimal

import qcelemental as qcel
from qcelemental.models import FailedOperation, ResultInput

import qcengine as qcng
from qcengine.exceptions import InputError
from qcengine.programs.nwchem import NWChemHarness
from qcengine.programs.nwchem.keywords import format_keywords
from qcengine.programs.util import PreservingDict

from ... import qcvars
from ...basisset import BasisSet
from ...util import print_jobrec, provenance_stamp
from .harvester import format_modelchem_for_nwchem, muster_inherited_options
from .molbasopt import format_molecule, muster_and_format_basis_for_nwchem

pp = pprint.PrettyPrinter(width=120)





def run_nwchem(name, molecule, options, **kwargs):
    """QCDB API to QCEngine connection for NWChem."""

    resi = ResultInput(
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

    jobrec = qcng.compute(resi, "qcdb-nwchem", raise_error=True).dict()

    hold_qcvars = jobrec['extras'].pop('qcdb:qcvars')
    jobrec['qcvars'] = {key: qcel.Datum(**dval) for key, dval in hold_qcvars.items()}

    return jobrec


class QcdbNWChemHarness(NWChemHarness):
    def compute(self, input_model: 'ResultInput', config: 'JobConfig') -> 'Result':
        self.found(raise_error=True)

        verbose = 1

        print_jobrec(f'[1] {self.name} RESULTINPUT PRE-PLANT', input_model.dict(), verbose >= 3)

        job_inputs = self.qcdb_build_input(input_model, config)

        print_jobrec(f'[2] {self.name}REC PRE-ENGINE', job_inputs, verbose >= 4)

        # 'NWCHEM_OMP_NUM_CORES': os.environ.get('NWCHEM_OMP_NUM_CORES'),

        success, dexe = self.execute(job_inputs)

        dexe["outfiles"]["stdout"] = dexe["stdout"]
        dexe["outfiles"]["stderr"] = dexe["stderr"]

        print_jobrec(f'[3] {self.name}REC POST-ENGINE', dexe, verbose >= 4)

        if 'There is an error in the input file' in dexe["stdout"]:
            raise InputError(dexe["stdout"])  # for nwc, stderr also works

        if success:
            output_model = self.parse_output(dexe["outfiles"], input_model)

            print_jobrec(f'[4a] {self.name} RESULT POST-HARVEST', output_model.dict(), verbose >= 5)

            output_model = self.qcdb_post_parse_output(input_model, output_model)

            print_jobrec(f'[4] {self.name} RESULT POST-POST-HARVEST', output_model.dict(), verbose >= 2)

        else:
            output_model = FailedOperation(success=False,
                                           error={
                                               "error_type": "execution_error",
                                               "error_message": dexe["stderr"],
                                           },
                                           input_data=input_model.dict())

        return output_model

    def qcdb_build_input(self, input_model: 'ResultInput', config: 'JobConfig',
                         template: Optional[str] = None) -> Dict[str, Any]:

        nwchemrec = {
            'infiles': {},
            'scratch_directory': config.scratch_directory,
        }

        molrec = qcel.molparse.from_schema(input_model.molecule.dict())
        molrecc1 = molrec.copy()
        molrecc1['fix_symmetry'] = 'c1'  # write _all_ atoms to input
        ropts = input_model.extras['qcdb:options']

        molcmd = format_molecule(molrec, ropts, verbose=1)

    # Handle memory
    # I don't think memory belongs in jobrec. it goes in pkgrec (for pbs) and possibly duplicated in options (for prog)
#    if 'memory' in jobrec:
#        memcmd, memkw = harvester.muster_memory(jobrec['memory'])
#    else:
#        memcmd, memkw = '', {}
    #mem = int(0.000001 * core.get_memory())
    #if mem == 524:
    #    memcmd, memkw = '', {}
    #else:
    #    memcmd, memkw = qcdb.cfour.muster_memory(mem)

        # Handle qcdb keywords implying nwchem keyword values
        # used to be before molecule formatting -- was that necessary?
        muster_inherited_options(ropts)

        _qcdb_basis = ropts.scroll['QCDB']['BASIS'].value
        # TODO _gamess_basis = ropts.scroll['NWCHEM']['BASIS'].value
        if _qcdb_basis == '':
            raise ValueError('nwchem bas not impl. set with `basis cc-pvdz`, etc. to use qcdb (psi4) basis set library.')
        # create a qcdb.Molecule to reset PG to c1 so all atoms. but print_detail isn't transmitting in a way nwc is picking up on, so per-element for now
        #qmol = Molecule(jobrec['molecule'])
        #qmol.reset_point_group('c1')  # need basis printed for *every* atom
        #qbs = BasisSet.pyconstruct(qmol, 'BASIS', _qcdb_basis)
        qbs = BasisSet.pyconstruct(molrec, 'BASIS', _qcdb_basis)

        #if qbs.has_ECP(): #    raise ValidationError("""ECPs not hooked up for Cfour""")
        bascmd = muster_and_format_basis_for_nwchem(molrec, ropts, qbs, verbose=1)

        # Handle calc type and quantum chemical method
  #      harvester.nu_muster_modelchem(jobrec['method'], jobrec['dertype'], ropts])
        mdccmd = format_modelchem_for_nwchem(input_model.model.method, input_model.driver.derivative_int(), ropts, sysinfo=None)

    #PRprint('Touched Keywords')
    #PRprint(ropts.print_changed(history=False))

        # Handle driver vs input/default keyword reconciliation

        # Handle conversion of qcdb keyword structure into nwchem format
#OLD    optcmd = moptions.prepare_options_for_nwchem(jobrec['options'])
#    resolved_options = {k: v.value for k, v in jobrec['options'].scroll['NWCHEM'].items() if v.disputed()}
        skma_options = {key: ropt.value for key, ropt in sorted(ropts.scroll['NWCHEM'].items()) if ropt.disputed()}
        optcmd = format_keywords(skma_options)

    # Handle text to be passed untouched to cfour
#    litcmd = core.get_global_option('LITERAL_CFOUR')

        # Assemble ZMAT pieces
    #zmat = memcmd + molcmd + optcmd + mdccmd + psicmd + bascmd + litcmd
    #zmat = molcmd + optcmd + bascmd
        nwchemrec['infiles']['nwchem.nw'] = 'echo\n' + molcmd + bascmd + optcmd + mdccmd
#OLD    nwchemrec['nwchem.nw'] = write_input(jobrec['method'], jobrec['dertype'], jobrec['molecule'], jobrec['options']) #molecule)
#    print('<<< ZMAT||\n{}\n||>>>\n'.format(nwchemrec['nwchem.nw']))
        nwchemrec['command'] = ['nwchem']  # subnw?

        return nwchemrec

    def qcdb_post_parse_output(self, input_model: 'ResultInput', output_model: 'Result') -> 'Result':

        progvars = PreservingDict(copy.deepcopy(output_model.extras['qcvars']))
        ropts = input_model.extras['qcdb:options']

        # badly placed
        if "MP2 OPPOSITE-SPIN CORRELATION ENERGY" in progvars and "MP2 SAME-SPIN CORRELATION ENERGY" in progvars:
            oss_opt = ropts.scroll['QCDB']['MP2_OS_SCALE']
            sss_opt = ropts.scroll['QCDB']['MP2_SS_SCALE']
            custom_scsmp2_corl = \
                    Decimal(oss_opt.value) * progvars["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] + \
                    Decimal(sss_opt.value) * progvars["MP2 SAME-SPIN CORRELATION ENERGY"]
            if "MP2 SINGLES ENERGY" in progvars:
                custom_scsmp2_corl += progvars["MP2 SINGLES ENERGY"]
            progvars["CUSTOM SCS-MP2 CORRELATION ENERGY"] = custom_scsmp2_corl

        if "CCSD OPPOSITE-SPIN CORRELATION ENERGY" in progvars and "CCSD SAME-SPIN CORRELATION ENERGY" in progvars:
            oss_opt = ropts.scroll['QCDB']['CCSD_OS_SCALE']
            sss_opt = ropts.scroll['QCDB']['CCSD_SS_SCALE']
            custom_scsmp2_corl = \
                    Decimal(oss_opt.value) * progvars["CCSD OPPOSITE-SPIN CORRELATION ENERGY"] + \
                    Decimal(sss_opt.value) * progvars["CCSD SAME-SPIN CORRELATION ENERGY"]
            # TODO check on validity of CCSD SINGLES ENERGY
            if "CCSD SINGLES ENERGY" in progvars:
                custom_scsccsd_corl += progvars["CCSD SINGLES ENERGY"]
            progvars["CUSTOM SCS-CCSD CORRELATION ENERGY"] = custom_scsccsd_corl

        qcvars.build_out(progvars)
        calcinfo = qcvars.certify_and_datumize(progvars, plump=True, nat=len(output_model.molecule.symbols))
        output_model.extras['qcdb:qcvars'] = calcinfo

        return output_model
