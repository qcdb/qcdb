import copy
import pprint
pp = pprint.PrettyPrinter(width=120)
from typing import Any, Dict, Optional
import inspect

import qcelemental as qcel
from qcelemental.models import ResultInput

import qcengine as qcng
from qcengine.programs.util import PreservingDict
from qcengine.programs.psi4 import Psi4Harness
#from qcengine.programs.cfour.keywords import format_keywords, format_keyword

from ... import qcvars
#from ..driver.driver_helpers import print_variables
#from ..exceptions import *
from ...util import provenance_stamp
#from .worker import psi4_subprocess
from .botanist import muster_inherited_options


#def run_psi4_old(name, molecule, options, **kwargs):
#
#    if kwargs['ptype'] not in ['energy', 'properties', 'gradient', 'hessian']:
#        raise ValidationError("""run_psi4: ptype not regonized: {}""".format(ptype))
#
#    prov = {}
#    prov['creator'] = 'QCDB'
#    prov['version'] = __version__
#    prov['routine'] = sys._getframe().f_code.co_name
#
#    jobrec = {}
#    jobrec['error'] = ''
#    jobrec['success'] = None
#    jobrec['return_output'] = True
#    jobrec['provenance'] = [prov]
#    jobrec['molecule'] = molecule.to_dict(np_out=False)
#    jobrec['method'] = name
#    jobrec['driver'] = kwargs['ptype']
#    jobrec['kwargs'] = kwargs
#    jobrec['options'] = copy.deepcopy(options)
#    jobrec['hooks'] = kwargs.get('hooks', {})
#
#    try:
#        jobrec['molecule']
#        jobrec['method']
#    except KeyError as err:
#        jobrec['error'] += repr(err) + 'Required fields missing from ({})'.format(jobrec.keys())
#        return jobrec
#
#
#    psi4rec = {}
#    psi4rec['json'] = {}
#
#    opts = jobrec['options']
#    # NOTE TODO very limited OPTIONS HANDSHAKE
#    muster_inherited_options(opts)
#
#    psi4rec['json']['schema_name'] = 'qcschema_input'
#    psi4rec['json']['schema_version'] = 1
#    omem = opts.scroll['QCDB'].pop('MEMORY')
#    psi4rec['json']['memory'] = omem.value
#    psi4rec['json']['molecule'] = qcel.molparse.to_schema(jobrec['molecule'], dtype=2)
#    psi4rec['json']['driver'] = jobrec['driver']
#    mtd = jobrec['method']
#    psi4rec['json']['model'] = {
#        'method': mtd[3:] if mtd.startswith('p4-') else mtd,
#        'basis': '(auto)',
#    }
#    psi4rec['json']['extras'] = {'wfn_qcvars_only': True}
#    psi4rec['json']['kwargs'] = jobrec['kwargs']
#    psi4rec['json']['return_output'] = True
#
#    popts = {}
#    for k, v in opts.scroll['QCDB'].items():
#        if v.disputed():
#            popts[k] = v.value
#
#    for k, v in opts.scroll['PSI4'].items():
#        if v.disputed():
#            popts[k] = v.value
#    psi4rec['json']['keywords'] = popts
#    if 'BASIS' in psi4rec['json']['keywords']:
#        psi4rec['json']['model']['basis'] =  psi4rec['json']['keywords']['BASIS']
#
#    psi4rec['command'] = ['psi4', '--json', '--nthread', '6']  # TODO
#
#    psi4_subprocess(psi4rec)  # updates psi4rec
#
#    psi4rec = psi4rec['json']  # TODO NOT how this should be done figure out 1-tier/2-tier
#
#    try:
#        pass
#    except KeyError as err:
#        raise KeyError(
#            'Required fields missing from ({})'.format(jobrec.keys())) from err
#
#    try:
#        psi4rec['raw_output']
#    except KeyError as err:
#        raise KeyError('Required fields missing from ({})'.format(
#            psi4rec.keys())) from err
#
#    if not psi4rec['success']:
#        raise RuntimeError(psi4rec['error'])
#
#    for fl in psi4rec.keys():
#        if fl.startswith('outfile_'):
#            jobrec[fl] = psi4rec[fl]
#
#    # Absorb results into qcdb data structures
#    progvars = PreservingDict(psi4rec['extras']['qcvars'])
#    for k in list(progvars.keys()):
#        if k in ['DETCI AVG DVEC NORM', 'MCSCF TOTAL ENERGY']:
#            progvars.pop(k)
#
#    qcvars.build_out(progvars)
#    calcinfo = qcvars.certify_and_datumize(progvars, plump=True, nat=len(jobrec['molecule']['mass']))
#
#    jobrec['raw_output'] = psi4rec['raw_output']
#    jobrec['qcvars'] = calcinfo
#
#    jobrec['provenance'].append(psi4rec['provenance'])
#
#    return jobrec

def run_psi4(name, molecule, options, **kwargs):

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

    jobrec = qcng.compute(resi, "qcdb-psi4", raise_error=True).dict()
    hold_qcvars = jobrec['extras'].pop('qcdb:qcvars')
    jobrec['qcvars'] = {key: qcel.Datum(**dval) for key, dval in hold_qcvars.items()}
    #pp.pprint(jobrec)
    #print(jobrec.keys())
    #print(jobrec['success'])
    return jobrec

def _print_helper(label, dicary, do_print):
    if do_print:
        print(label + ' <<<')
        pp.pprint(dicary)
        print('>>>')




class QcdbPsi4Harness(Psi4Harness):

    def compute(self, input_model: 'ResultInput', config: 'JobConfig') -> 'Result':
        self.found(raise_error=True)

        verbose = 1
        _print_helper(f'[1] {self.name} RESULTINPUT PRE-PLANT', input_model.dict(), verbose >= 3)

        input_data = self.qcdb_build_input(input_model, config)
        input_model = ResultInput(**input_data)

        _print_helper(f'[2] {self.name} RESULTINPUT PRE-ENGINE', input_model.dict(), verbose >= 4)

        # 'PATH': (':'.join([os.path.abspath(x) for x in os.environ.get('PSIPATH', '').split(':') if x != '']) +
        #          ':' + os.environ.get('PATH')),# +
        # 'PSI_SCRATCH': tmpdir,
        # 'PYTHONPATH': os.environ.get('PYTHONPATH'),
        # 'LD_LIBRARY_PATH': os.environ.get('LD_LIBRARY_PATH')

        output_model = Psi4Harness.compute(self, input_model=input_model, config=config)

        _print_helper(f'[3] {self.name} RESULT POST-ENGINE', output_model.dict(), verbose >= 4)

        # ???
        if not output_model.success:
            return output_model

        _print_helper(f'[4a] {self.name} RESULT POST-HARVEST', output_model.dict(), verbose >= 5)

        output_model = self.qcdb_post_parse_output(input_model, output_model)

        _print_helper(f'[4] {self.name} RESULT POST-POST-HARVEST', output_model.dict(), verbose >= 2)

        return output_model

    def qcdb_build_input(self, input_model: 'ResultInput', config: 'JobConfig',
                         template: Optional[str] = None) -> Dict[str, Any]:
        input_data = input_model.dict()


        molrec = qcel.molparse.from_schema(input_model.molecule.dict())
        ropts = input_model.extras['qcdb:options']

        muster_inherited_options(ropts)
        mtd = input_data['model']['method']
        mtd =  mtd[3:] if mtd.startswith('p4-') else mtd
        input_data['model']['method'] = mtd

        # should we put this memory in the JobConfig object? I don't think the units agree
        omem = ropts.scroll['QCDB'].pop('MEMORY')
        #print(config.memory, '!!')
        #config.memory = omem.value #???
        #print(config.memory, '!!')

        input_data['extras'] = {'wfn_qcvars_only' : True}
        #input_data['kwargs'] = jobrec['kwargs']
        #input_data['return_output'] = True

        popts = {}
        for k, v in ropts.scroll['QCDB'].items():
            if v.disputed():
                popts[k] = v.value

        for k, v in ropts.scroll['PSI4'].items():
            if v.disputed():
                popts[k] = v.value
        input_data['keywords'] = popts

        if 'BASIS' in input_data['keywords']:
            input_data['model']['basis'] =  input_data['keywords']['BASIS']

        return input_data

    def qcdb_post_parse_output(self, input_model: 'ResultInput', output_model: 'Result') -> 'Result':

        dqcvars = PreservingDict(copy.deepcopy(output_model.extras['qcvars']))
        for k in list(dqcvars.keys()):
            if k in ['DETCI AVG DVEC NORM', 'MCSCF TOTAL ENERGY']:
                dqcvars.pop(k)
        qcvars.build_out(dqcvars)
        calcinfo = qcvars.certify_and_datumize(dqcvars, plump=True, nat=len(output_model.molecule.symbols))
        output_model.extras['qcdb:qcvars'] = calcinfo

        return output_model

