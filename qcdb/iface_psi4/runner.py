import sys
import copy
import pprint
pp = pprint.PrettyPrinter(width=120)
import inspect

import numpy as np

from .. import __version__
from .. import qcvars
from ..driver.driver_helpers import print_variables
#from ..datastructures import *
from ..exceptions import *
from ..molecule import Molecule
from ..pdict import PreservingDict
from .worker import psi4_subprocess


def run_psi4(name, molecule, options, **kwargs):
    print('\nhit run_psi4', name, kwargs)
    print('INTO PSI QP', options.print_changed())

    #calledby = inspect.stack()
    #print('CALLEDBY')
    #for cur in calledby:
    #    print('CUR', cur[3])
    if kwargs['ptype'] not in ['energy', 'properties', 'gradient', 'hessian']:
        raise ValidationError("""run_psi4: ptype not regonized: {}""".format(ptype))
    print('\n\n KWARGS', kwargs)


    jobrec = {}
    jobrec['error'] = ''
    jobrec['success'] = None
    jobrec['return_output'] = True
    prov = {}
    prov['creator'] = 'QCDB'
    prov['version'] = __version__
    prov['routine'] = sys._getframe().f_code.co_name
    jobrec['provenance'] = [prov]

    jobrec['molecule'] = molecule.to_dict(np_out=False)
    jobrec['method'] = name
    #jobrec['dertype'] = ['energy', 'gradient', 'hessian'].index(inspect.stack()[1][3])
    #jobrec['driver'] = 'energy'
    jobrec['driver'] = kwargs['ptype']
    jobrec['kwargs'] = kwargs
    jobrec['options'] = copy.deepcopy(options)

    jobrec = psi4_driver(jobrec)
    return jobrec


def psi4_driver(jobrec):
    import json

    try:
        jobrec['molecule']
        jobrec['method']
    except KeyError as err:
        #raise KeyError(
        #    'Required fields missing from ({})'.format(jobrec.keys())) from err
        jobrec['error'] += repr(err) + 'Required fields missing from ({})'.format(jobrec.keys())
        return jobrec

    print('PSI4_DRIVER jr', jobrec)
    # this is what the psi4 program needs, not what the job needs
    # *

    psi4rec = psi4_plant(jobrec)
    #psi4rec['scratch_messy'] = True

    # test json roundtrip
    jpsi4rec = json.dumps(psi4rec)
    psi4rec = json.loads(jpsi4rec)

    print('PSI4REC')
    pp.pprint(psi4rec)
    psi4_subprocess(psi4rec)  # updates psi4rec

    psi4_harvest(jobrec, psi4rec)  # ? updates jobrec ?

    print('QC', jobrec['qcvars'])
    print('PSI4 JOBREC')
    pp.pprint(jobrec)
    return jobrec


def psi4_plant(jobrec):  # jobrec@i -> psi4@i
    psi4rec = {}
    psi4rec['json'] = {}

    opts = jobrec['options']

    omem = opts.scroll['QCDB'].pop('MEMORY')
    psi4rec['json']['memory'] = omem.value  #pts.scroll['QCDB']['MEMORY'].value
    psi4rec['json']['molecule'] = {'qm': jobrec['molecule']}
    psi4rec['json']['driver'] = jobrec['driver']
    psi4rec['json']['method'] = jobrec['method']
    #psi4rec['json']['args'] = 
    psi4rec['json']['kwargs'] = jobrec['kwargs']
    #psi4rec['json']['options'] = 
    #psi4rec['json']['scratch_location'] = 
    psi4rec['json']['return_output'] = True

    # NOTE TODO skipped the OPTIONS HANDSHAKE
    popts = {}
    for k, v in opts.scroll['QCDB'].items():
        if v.disputed():
            popts[k] = v.value

    for k, v in opts.scroll['PSI4'].items():
        if 'E_CONVERGENCE' in k:
            print('OO', k, v)
        if v.disputed():
            popts[k] = v.value
    psi4rec['json']['options'] = popts
    print('INTO JSON RUNNER')


    print('in psi4_plant JR')
    pp.pprint(jobrec)


    # Handle qcdb keywords implying cfour keyword values
#    if core.get_option('CFOUR', 'TRANSLATE_PSI4'):
#    harvester.muster_inherited_options(jobrec['options'])

#    opts = jobrec['options']


    print('HH')
    print(jobrec['options'].print_changed())

    # Handle driver vs input/default keyword reconciliation
    #userkw = query_options_defaults_from_psi() #prepare_options_for_modules()

#    userkw = jobrec['options']
    print('\nDEFAULT OPTIONS')
    #print(userkw)
    #pprint.pprint(userkw['CFOUR'])

    # Handle conversion of qcdb keyword structure into psi4 format
#    optcmd = moptions.prepare_options_for_cfour(jobrec['options'])

#    print('<<< ZMAT||{}||>>>\n'.format(zmat))
    psi4rec['command'] = ['psi4', '--json']

    return psi4rec




def psi4_harvest(jobrec, psi4rec):  # jobrec@i, psi4rec@io -> jobrec@io
    """Processes raw results from read-only `psi4rec` into QCAspect fields in returned `jobrec`."""
    print('UK', psi4rec.keys())
    print('LK', psi4rec['json'].keys())
    psi4rec = psi4rec['json']  # TODO figure out 1-tier/2-tier
#
#UK dict_keys(['command', 'stdout', 'json'])
#LK dict_keys(['options', 'success', 'memory', 'variables', 'raw_output', 'return_value', 'error', 'driver', 'provenance', 'psivars', 'method', 'molecule'])

    try:
        pass
        #jobrec['molecule']['real']
        #jobrec['do_gradient']
    except KeyError as err:
        raise KeyError(
            'Required fields missing from ({})'.format(jobrec.keys())) from err

    try:
        psi4rec['raw_output']
        #if jobrec['do_gradient'] is True:
        #    dftd3rec['dftd3_gradient']
    except KeyError as err:
        raise KeyError('Required fields missing from ({})'.format(
            psi4rec.keys())) from err

    if psi4rec['error']:
        raise RuntimeError(psi4rec['error'])
    print('PSIREC POST <<<')
    pp.pprint(psi4rec)
    print('>>>')

    # Absorb results into qcdb data structures
    progvars = PreservingDict(psi4rec['psivars'])
    import psi4
    progarrs = {k: np.array(psi4.core.Matrix.from_serial(v)) for k, v in psi4rec['psiarrays'].items()}
    print('PRGARR', progarrs)
    progvars.update(progarrs)
#    if psi4rec['driver'] in ['energy', 'properties']:
#        progvars['CURRENT ENERGY'] = psi4rec['return_value']
#    elif psi4rec['driver'] in ['hessian']:
#        progvars['CURRENT ENERGY'] = psi4rec['return_value']
#        progvars['CURRENT GRADIENT'] = psi4rec['return_value']
#        progvars['CURRENT HESSIAN'] = psi4rec['return_value']
    qcvars.build_out(progvars)
    calcinfo = qcvars.certify(progvars)

    #jobrec['raw_output'] = text
    jobrec['qcvars'] = calcinfo

    #prov = {}
    #prov['creator'] = 'Psi4'
    #prov['routine'] = sys._getframe().f_code.co_name
    #prov['version'] = version
    jobrec['provenance'].append(psi4rec['provenance'])

    return jobrec






    """
    Required Input Fields
    ---------------------

    Optional Input Fields
    ---------------------

    Output Fields
    -------------

    """


