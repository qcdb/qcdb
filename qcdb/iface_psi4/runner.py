import sys

from .. import __version__
from .. import qcvars
from ..exceptions import *
from ..pdict import PreservingDict


def run_psi4_deferred(name, molecule, options, **kwargs):
    print('\nhit run_psi4_deferred', name, kwargs)
    print(options)

    jobrec = {}
    prov = {}
    prov
    prov['creator'] = 'QCDB'
    prov['version'] = __version__
    prov['routine'] = sys._getframe().f_code.co_name
    jobrec['provenance'] = [prov]
    jobrec['molecule'] = {'qm': molecule.to_dict(np_out=False)}

    try:
        mem = options.scroll['QCDB'].pop('MEMORY')
    except KeyError:
        pass
    else:
#    if options.scroll['QCDB']['MEMORY'].has_changed:
#        mem = options.scroll['QCDB']['MEMORY'].value
#        del options.scroll['QCDB']['MEMORY']
#    #if 'MEMORY' in options['GLOBALS']:
#    #    mem = options['GLOBALS'].pop('MEMORY')
        print('MEM', mem)
        jobrec['memory'] = mem.value #mem['value']

#    try:
#        dftd3_driver(jobrec)
#        jobrec['success'] = True
#    except Exception as err:
#        jobrec['success'] = False
#        #json_data["error"] += repr(error)

    jobrec['error'] = ''
    jobrec['success'] = False
    #jobrec['stdout']
    jobrec['return_output'] = True

    popts = {}
    for k, v in options.scroll['QCDB'].items():
#    for k, v in options['GLOBALS'].items():
#        print('P4 opt', k, v)
#   #     psi4.core.set_global_option(k.upper(), v['value'])
#        popts[k] = v['value']
        print('WWWW', k, v.value, v.is_default())
        if not v.is_default():
            popts[k] = v.value
            
    jobrec['options'] = popts

    jobrec['driver'] = 'energy'
    jobrec['method'] = name
    jobrec['return_output'] = False #True
    import pprint
    print('JOBREC PREE <<<')
    pprint.pprint(jobrec)
    print('>>>')
    import psi4
    psi4.core.clean()
    #try:
    psi4.json_wrapper.run_json(jobrec)
    #    jobrec['success'] = True
    #except Exception as err:
    #    jobrec['success'] = False
    #    print ('STOP TROUBLE')
    if jobrec['error']:
        raise RuntimeError(jobrec['error'])
    print('JOBREC POST <<<')
    pprint.pprint(jobrec)
    print('>>>')

    progvars = PreservingDict(jobrec['psivars'])
    #morevars = qcvars.fill_in(progvars)
    #calcinfo = qcvars.certify_qcvars(morevars)
    qcvars.fill_in(progvars)
    calcinfo = qcvars.certify_qcvars(progvars)
    jobrec['qcvars'] = calcinfo

    #jobrec['qcvars'] = {info.lbl: info for info in calcinfo}
#    jobrec['wfn'] = wfn
    return jobrec

#run_psi4 = run_psi4_realtime
run_psi4 = run_psi4_deferred

def write_job(name, dertype, molecule, options):
    pass

    # Handle memory
    memopt = options['GLOBALS']['MEMORY']
    if memopt.is_default():
        memcmd, memkw = '', {}
    else:
        memcmd, memkw = muster_memory(memopt.value)

#    mem = int(0.000001 * core.get_memory())
#    if mem == 524:
#        memcmd, memkw = '', {}
#    else:
#        memcmd, memkw = qcdb.cfour.muster_memory(mem)


def muster_memory(mem):
    """Transform input `mem` in bytes options for psi4.

    """
    text = ''

    # prepare memory keywords to be set as c-side keywords
    options = defaultdict(lambda: defaultdict(dict))
    options['PSI']['MEMORY']['value'] = int(mem)

    for item in options['PSI']:
        options['PSI'][item]['clobber'] = True
    return text, options


def run_psi4_realtime(name, molecule, options, **kwargs):
    print('\nhit run_psi4_realtime', name, options.keys(), options, kwargs)
    import psi4

    if 'MEMORY' in options['GLOBALS']:
        mem = options['GLOBALS'].pop('MEMORY')
        print('MEM', mem)
        psi4.set_memory(mem['value'])

    for k, v in options['GLOBALS'].items():
        print('P4 opt', k, v)
        psi4.core.set_global_option(k.upper(), v['value'])

    pmol = psi4.core.Molecule.from_dict(molecule.to_dict())
    _, wfn = psi4.energy(name, molecule=pmol, return_wfn=True, **kwargs)

    calcinfo = []
    for pv, var in wfn.variables().items():
        if pv in qcvardefs.keys():
            calcinfo.append(QCAspect(pv, qcvardefs[pv]['units'], var, ''))
        else:
            raise ValidationError('Undefined QCvar!: {}'.format(pv))

    jobrec = {}
    jobrec['qcvars'] = {info.lbl: info for info in calcinfo}
    jobrec['wfn'] = wfn
    return jobrec


