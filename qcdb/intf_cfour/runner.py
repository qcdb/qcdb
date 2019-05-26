import copy
import pprint
pp = pprint.PrettyPrinter(width=120)
import inspect
from typing import Any, Dict, Optional
from decimal import Decimal

import qcelemental as qcel
from qcelemental.models import ResultInput

import qcengine as qcng
from qcengine.programs.util import PreservingDict
from qcengine.programs.cfour import CFOURExecutor
from qcengine.programs.cfour.keywords import format_keywords, format_keyword

from .. import qcvars
#from ..driver.driver_helpers import print_variables
from ..libmintsbasisset import BasisSet
from ..util import provenance_stamp
from . import harvester
from .bas import extract_basis_from_genbas, format_basis_for_cfour, format_molecule


def run_cfour(name, molecule, options, **kwargs):
    """QCDB API to QCEngine connection for CFOUR."""

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

    jobrec = qcng.compute(resi, "cfour", raise_error=True).dict()

    hold_qcvars = jobrec['extras'].pop('qcdb:qcvars')
    jobrec['qcvars'] = {key: qcel.Datum(**dval) for key, dval in hold_qcvars.items()}

    return jobrec


def _print_helper(label, dicary, do_print):
    if do_print:
        print(label + ' <<<')
        pp.pprint(dicary)
        print('>>>')


class QcdbCFOURExecutor(CFOURExecutor):
    def compute(self, input_model: 'ResultInput', config: 'JobConfig') -> 'Result':
        self.found(raise_error=True)

        verbose = 2

        _print_helper(f'[1] {self.name} RESULTINPUT PRE-PLANT', input_model.dict(), verbose >= 3)

        job_inputs = self.qcdb_build_input(input_model, config)

        _print_helper(f'[2] {self.name}REC PRE-ENGINE', job_inputs, verbose >= 4)

        success, dexe = self.execute(job_inputs)

        _print_helper(f'[3] {self.name}REC POST-ENGINE', dexe, verbose >= 4)

        if not success:
            output_model = input_model
            output_model["error"] = {"error_type": "execution_error", "error_message": output["stderr"]}

        dexe["outfiles"]["stdout"] = dexe["stdout"]
        dexe["outfiles"]["stderr"] = dexe["stderr"]
        output_model = self.parse_output(dexe["outfiles"], input_model)

        _print_helper(f'[4a] {self.name} RESULT POST-HARVEST', output_model.dict(), verbose >= 5)

        output_model = self.qcdb_post_parse_output(input_model, output_model)

        _print_helper(f'[4] {self.name} RESULT POST-POST-HARVEST', output_model.dict(), verbose >= 2)

        return output_model

    def qcdb_build_input(self, input_model: 'ResultInput', config: 'JobConfig',
                         template: Optional[str] = None) -> Dict[str, Any]:
        cfourrec = {
            'infiles': {},
            'scratch_location': config.scratch_directory,
        }

        molrec = qcel.molparse.from_schema(input_model.molecule.dict())
        ropts = input_model.extras['qcdb:options']

        molcmd = format_molecule(molrec, ropts, verbose=1)

        # Handle qcdb keywords implying cfour keyword values
        harvester.muster_inherited_options(ropts)

        _qcdb_basis = ropts.scroll['QCDB']['BASIS'].value
        _cfour_basis = ropts.scroll['CFOUR']['BASIS'].value
        if _qcdb_basis == '':
            _, cased_basis = format_keyword('CFOUR_BASIS', _cfour_basis)
            cfourrec['infiles']['GENBAS'] = extract_basis_from_genbas(cased_basis,
                                                                      input_model.molecule.symbols,
                                                                      exact=False)
            bascmd = ''
        else:
            qbs = BasisSet.pyconstruct(molrec, 'BASIS', _qcdb_basis)
            #if qbs.has_ECP(): #    raise ValidationError("""ECPs not hooked up for Cfour""")
            cfourrec['infiles']['GENBAS'] = qbs.print_detail_cfour()  #qbs.genbas()
            bascmd = format_basis_for_cfour(molrec, ropts, qbs.has_puream())

        # Handle calc type and quantum chemical method
        harvester.nu_muster_modelchem(input_model.model.method, input_model.driver.derivative_int(), ropts)

        #print(jobrec['options'].print_changed(history=False))
        # Handle driver vs input/default keyword reconciliation

        # Handle conversion of psi4 keyword structure into cfour format
        skma_options = {key: ropt.value for key, ropt in sorted(ropts.scroll['CFOUR'].items()) if ropt.disputed()}
        optcmd = format_keywords(skma_options)

        # Assemble ZMAT pieces
        cfourrec['infiles']['ZMAT'] = molcmd + optcmd + bascmd
        cfourrec['command'] = ['xcfour']

        return cfourrec

    def qcdb_post_parse_output(self, input_model: 'ResultInput', output_model: 'Result') -> 'Result':

        dqcvars = PreservingDict(copy.deepcopy(output_model.extras['qcvars']))

        # badly placed
        # Cfour's SCS-MP2 is non adjustible and only valid for UHF
        # ROMP2 doesn't print SS & OS
        if "MP2 OPPOSITE-SPIN CORRELATION ENERGY" in dqcvars and "MP2 SAME-SPIN CORRELATION ENERGY" in dqcvars:
            oss_opt = input_model.extras['qcdb:options'].scroll['QCDB']['MP2_OS_SCALE']
            sss_opt = input_model.extras['qcdb:options'].scroll['QCDB']['MP2_SS_SCALE']
            custom_scsmp2_corl = \
                Decimal(oss_opt.value) * dqcvars["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] + \
                Decimal(sss_opt.value) * dqcvars["MP2 SAME-SPIN CORRELATION ENERGY"]
            if "MP2 SINGLES ENERGY" in dqcvars:
                custom_scsmp2_corl += dqcvars["MP2 SINGLES ENERGY"]
            dqcvars["CUSTOM SCS-MP2 CORRELATION ENERGY"] = custom_scsmp2_corl

        qcvars.build_out(dqcvars)
        calcinfo = qcvars.certify(dqcvars, plump=True, nat=len(output_model.molecule.symbols))
        output_model.extras['qcdb:qcvars'] = calcinfo

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

#print('HH')
#print(jobrec['options'].print_changed(history=False))

#    tmp_zmat = """UHF-SCF energy calculation
#N
#H 1 R
#H 1 R 2 A
#
#R=1.008
#A=105.0"""
#    tmp_zmat = """auto-generated by qcdb from molecule H2N
#N                     0.000000000000     0.000000000000    -0.145912918812
#H                     0.000000000000    -1.511214304079     1.013682601327
#H                     0.000000000000     1.511214304079     1.013682601327
#
#*ACES2(CALC=HF,BASIS=qz2p
#MULT=2,REF=UHF
#COORDINATES=CARTESIAN
#UNITS=BOHR
#OCCUPATION=3-1-1-0/3-0-1-0
#SCF_CONV=12
#MEMORY=20000000)
#
#"""

#    # amalgamate output
#    text = cfourrec['stdout']
#    text += '\n  <<<  Cfour {} {} Results  >>>\n\n'.format('', '') #name.lower(), calledby.capitalize()))  # banner()
#
#    for fl, contents in cfourrec['outputs'].items():
#        if contents is not None:
#            text += f'\n  Cfour scratch file {fl} has been read.\n'
#            text += contents
#    text += print_variables(calcinfo)
#    jobrec['raw_output'] = text
