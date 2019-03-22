import sys
import copy
import pprint
pp = pprint.PrettyPrinter(width=120)
import inspect

import numpy as np
import qcelemental as qcel

from .. import __version__
from .. import qcvars
from ..driver.driver_helpers import print_variables
from ..exceptions import *
from ..molecule import Molecule
from ..pdict import PreservingDict
from .worker import psi4_subprocess
from .botanist import muster_inherited_options


def run_psi4(name, molecule, options, **kwargs):
    #print('run_psi4 options <<<\n', options.print_changed(), '\n>>>')

    #calledby = inspect.stack()
    #print('CALLEDBY')
    #for cur in calledby:
    #    print('CUR', cur[3])
    if kwargs['ptype'] not in ['energy', 'properties', 'gradient', 'hessian']:
        raise ValidationError("""run_psi4: ptype not regonized: {}""".format(ptype))

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
    jobrec['driver'] = kwargs['ptype']
    jobrec['kwargs'] = kwargs
    jobrec['options'] = copy.deepcopy(options)
    jobrec['hooks'] = kwargs.get('hooks', {})

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

    #print('[1] PSI4 JOBREC PRE-PLANT (j@i) <<<')
    #pp.pprint(jobrec)
    #print('>>>')

    psi4rec = psi4_plant(jobrec)

#    # test json roundtrip
#    jpsi4rec = json.dumps(psi4rec)
#    psi4rec = json.loads(jpsi4rec)

    #print('[2] PSI4REC PRE-SUBPROCESS (x@i) <<<')
    #pp.pprint(psi4rec)
    #print('>>>\n')

    psi4_subprocess(psi4rec)  # updates psi4rec

    #print('[3] PSI4REC POST-SUBPROCESS (x@io) <<<')
    #pp.pprint(psi4rec)
    #print('>>>\n')

    psi4_harvest(jobrec, psi4rec)  # updates jobrec

    #print('[4] PSI4 JOBREC POST-HARVEST (j@io) <<<')
    #pp.pprint(jobrec)
    #print('>>>')

    return jobrec


def psi4_plant(jobrec):  # jobrec@i -> psi4@i
    psi4rec = {}
    psi4rec['json'] = {}

    opts = jobrec['options']
    # NOTE TODO very limited OPTIONS HANDSHAKE
    muster_inherited_options(opts)

    psi4rec['json']['schema_name'] = 'qcschema_input'
    psi4rec['json']['schema_version'] = 1
    omem = opts.scroll['QCDB'].pop('MEMORY')
    psi4rec['json']['memory'] = omem.value
    psi4rec['json']['molecule'] = qcel.molparse.to_schema(jobrec['molecule'], dtype=2)
    psi4rec['json']['driver'] = jobrec['driver']
    mtd = jobrec['method']
    psi4rec['json']['model'] = {
        'method': mtd[3:] if mtd.startswith('p4-') else mtd,
        'basis': '(auto)',
    }
    psi4rec['json']['extras'] = {'wfn_qcvars_only': True}
    #psi4rec['json']['args'] = 
    psi4rec['json']['kwargs'] = jobrec['kwargs']
    #psi4rec['json']['scratch_location'] = 
    psi4rec['json']['return_output'] = True

    #for hookkey, hookfunc in jobrec['hooks']['pre'].items():
    #    psi4rec['json']['in_' + hookkey] = hookfunc()

    # March 2019 temporarily suspending GRIDDAT
    #if opts.scroll['PSI4']['GRIDDAT'].value != '':
    #    psi4rec['json']['infile_' + 'grid.dat'] = opts.scroll['PSI4']['GRIDDAT'].value
    
    popts = {}
    for k, v in opts.scroll['QCDB'].items():
        if v.disputed():
            popts[k] = v.value

    for k, v in opts.scroll['PSI4'].items():
        if v.disputed():
            popts[k] = v.value
    psi4rec['json']['keywords'] = popts
    if 'BASIS' in psi4rec['json']['keywords']:
        psi4rec['json']['model']['basis'] =  psi4rec['json']['keywords']['BASIS']

    # Handle qcdb keywords implying cfour keyword values
#    if core.get_option('CFOUR', 'TRANSLATE_PSI4'):
#    harvester.muster_inherited_options(jobrec['options'])

    # Handle conversion of qcdb keyword structure into psi4 format
    # * psi wants python anyways, so no action needed

    #psi4rec['command'] = ['psi4', '--json']
    psi4rec['command'] = ['psi4', '--json', '--nthread', '6']  # TODO

    return psi4rec




def psi4_harvest(jobrec, psi4rec):  # jobrec@i, psi4rec@io -> jobrec@io
    """Processes raw results from read-only `psi4rec` into QCAspect fields in returned `jobrec`."""

    psi4rec = psi4rec['json']  # TODO NOT how this should be done figure out 1-tier/2-tier

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

    if not psi4rec['success']:
        raise RuntimeError(psi4rec['error'])

    #c4files = {}
    for fl in psi4rec.keys():
        if fl.startswith('outfile_'):
            jobrec[fl] = psi4rec[fl]
    #for fl in ['GRD', 'FCMFINAL', 'DIPOL']:
    #    field = 'output_' + fl.lower()
    #    if field in cfourrec:
    #        text += '  Cfour scratch file {} has been read\n'.format(fl)
    #        text += cfourrec[field]
    #        c4files[fl] = cfourrec[field]

    # Absorb results into qcdb data structures
    progvars = PreservingDict(psi4rec['extras']['qcvars'])
    for k in list(progvars.keys()):
        if k in ['DETCI AVG DVEC NORM', 'MCSCF TOTAL ENERGY']:
            progvars.pop(k)

    qcvars.build_out(progvars)
    calcinfo = qcvars.certify(progvars, plump=True, nat=len(jobrec['molecule']['mass']))

    jobrec['raw_output'] = psi4rec['raw_output']
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

