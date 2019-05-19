import sys
import copy
import pprint
pp = pprint.PrettyPrinter(width=120)
import inspect
#from decimal import Decimal

import qcengine as qcng

from .. import __version__
#from .. import moptions
from .. import qcvars
#from ..driver.driver_helpers import print_variables
#from ..exceptions import *
#from ..iface_psi4.options import query_options_defaults_from_psi
from ..libmintsbasisset import BasisSet
from ..molecule import Molecule
from ..pdict import PreservingDict
from .molbasopt import muster_and_format_molecule_and_basis_for_gamess, format_options_for_gamess
from .harvester import muster_modelchem, muster_inherited_options


def run_gamess(name, molecule, options, **kwargs):
#    print('\nhit run_gamess', name, kwargs)

#    #calledby = inspect.stack()
#    #print('CALLEDBY')
#    #for cur in calledby:
#    #    print('CUR', cur[3])

    jobrec = {}
    jobrec['error'] = ''
    jobrec['success'] = None
    prov = {}
    prov
    prov['creator'] = 'QCDB'
    prov['version'] = __version__
    prov['routine'] = sys._getframe().f_code.co_name
    jobrec['provenance'] = [prov]

    jobrec['molecule'] = molecule.to_dict(np_out=False)
    jobrec['method'] = name
    jobrec['dertype'] = ['energy', 'gradient', 'hessian'].index(inspect.stack()[1][3])

#    #if 'MEMORY' in options['GLOBALS']:
#    #    mem = options['GLOBALS'].pop('MEMORY')
#    #    print('MEM', mem)
#    #    jobrec['memory'] = mem['value']
#
##    popts = collections.defaultdict(lambda: collections.defaultdict(dict))
##    #for k, v in options['GLOBALS'].items():
##    for k, v in options['CFOUR'].items():
##        print('C4 opt', k, v)
##   #     psi4.core.set_global_option(k.upper(), v['value'])
##        popts[k] = v['value']
##    jobrec['options'] = popts

    qcschema_input = {
        'schema_name': 'qcschema_input',
        'schema_version': 1,
        'driver': inspect.stack()[1][3],
        'model': {
            'method': name,
            'basis': '(auto)',
        },
        'molecule': molecule.to_schema(dtype=2),
        'extras': {},
    }

    jobrec['options'] = copy.deepcopy(options)
    jobrec['qcschema_input'] = qcschema_input


#    print('comin in')
#    #print(jobrec['options'])
#PR    print(jobrec['options'].print_changed())

    jobrec = gamess_driver(jobrec)

    return jobrec


def gamess_driver(jobrec):
#    import json
#
#    try:
#        #jobrec['dashlevel']
#        #jobrec['dashparams']
#        #jobrec['functional']
#        jobrec['molecule']
#        #jobrec['do_gradient']
#    except KeyError as err:
#        #raise KeyError(
#        #    'Required fields missing from ({})'.format(jobrec.keys())) from err
#        jobrec['error'] += repr(err) + 'Required fields missing from ({})'.format(jobrec.keys())
#        return jobrec

#    print('GAMESS_DRIVER jr', jobrec)

    gamessrec = gamess_plant(jobrec)
#    #cfourrec = cfour_plant(jobrec)
#    ##cfourrec['scratch_messy'] = True
#
#    ## test json roundtrip
#    #jcfourrec = json.dumps(cfourrec)
#    #cfourrec = json.loads(jcfourrec)

#    print('GAMESSREC')
#    pp.pprint(gamessrec)

    qcschema_input = jobrec['qcschema_input']
    qcschema_input['extras']['gamess.inp'] = gamessrec['gamess.inp']

    ret = qcng.compute(qcschema_input, 'gamess').dict()
#    print('[5] {} JOBREC POST-HARVEST (j@io) <<<'.format('CFOUR'))
#    pp.pprint(ret)
#    print('>>>')
    jobrec.pop('qcschema_input')
    progvars = PreservingDict(ret['extras']['qcvars'])
    qcvars.build_out(progvars)
    calcinfo = qcvars.certify(progvars, plump=True, nat=len(jobrec['molecule']['mass']))
    jobrec['raw_output'] = ret['stdout']
    jobrec['qcvars'] = calcinfo
    jobrec['success'] = True

#    jobrec['qcvars'] = ret['extras']['qcvars']
#    pp.pprint(jobrec)
#    print('END [6]')




#    #cfour_subprocess(cfourrec)  # updates cfourrec
#
#    #cfour_harvest(jobrec, cfourrec)  # ? updates jobrec ?
#
#    #print('QC', jobrec['qcvars'])
#    #return jobrec
#
#    verbose = 1
#
#    if verbose >= 3:
#        print('[1] {} JOBREC PRE-PLANT (j@i) <<<'.format('CFOUR'))
#        pp.pprint(jobrec)
#        print('>>>')
#
#    cfourrec = cfour_plant(jobrec)
#
#    # test json roundtrip
#    jcfourrec = json.dumps(cfourrec)
#    cfourrec = json.loads(jcfourrec)
#
#    if verbose >= 4:
#        print('[2] {}REC PRE-SUBPROCESS (m@i) <<<'.format('CFOUR'))
#        pp.pprint(cfourrec)
#        print('>>>\n')
#
#    cfour_subprocess(cfourrec)  # updates cfourrec
#
#    if verbose >= 4:
#        print('[3] {}REC POST-SUBPROCESS (m@io) <<<'.format('CFOUR'))
#        pp.pprint(cfourrec)
#        print('>>>\n')
#
#    cfour_harvest(jobrec, cfourrec)  # updates jobrec
#
#    if verbose >= 2:
#        print('[4] {} JOBREC POST-HARVEST (j@io) <<<'.format('CFOUR'))
#        pp.pprint(jobrec)
#        print('>>>')

    return jobrec


#def psi4_plant(jobrec):  # jobrec@i -> psi4@i
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
#    #psi4rec['json']['args'] = 
#    psi4rec['json']['kwargs'] = jobrec['kwargs']
#    #psi4rec['json']['scratch_location'] = 
#    psi4rec['json']['return_output'] = True
#
#    #for hookkey, hookfunc in jobrec['hooks']['pre'].items():
#    #    psi4rec['json']['in_' + hookkey] = hookfunc()
#    if opts.scroll['PSI4']['GRIDDAT'].value != '':
#        psi4rec['json']['infile_' + 'grid.dat'] = opts.scroll['PSI4']['GRIDDAT'].value
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


def gamess_plant(jobrec):  # jobrec@i -> gamess@i
    gamessrec = {}
#
#    #return cfourrec
#
##def write_zmat(name, dertype, molecule):  # -> zmat
#    """Returns string with contents of Cfour ZMAT file as gathered from
#    active molecule, current keyword settings, and cfour {...} block.
#
#    """
#    import qcdb
#    # Handle memory
#    # I don't think memory belongs in jobrec. it goes in pkgrec (for pbs) and possibly duplicated in options (for prog)
##    if 'memory' in jobrec:
##        memcmd, memkw = harvester.muster_memory(jobrec['memory'])
##    else:
##        memcmd, memkw = '', {}
#    #mem = int(0.000001 * core.get_memory())
#    #if mem == 524:
#    #    memcmd, memkw = '', {}
#    #else:
#    #    memcmd, memkw = qcdb.cfour.muster_memory(mem)

#    print('in gamess_plant')
#    pp.pprint(jobrec)

#    molcmd = format_molecule_and_basis_for_gamess(jobrec['molecule'], jobrec['options'], verbose=1)

    # Handle qcdb keywords implying gamess keyword values
#    # if core.get_option('CFOUR', 'TRANSLATE_PSI4'):
    muster_inherited_options(jobrec['options'])

    qmol = Molecule(jobrec['molecule'])
    _qcdb_basis = jobrec['options'].scroll['QCDB']['BASIS'].value
    _gamess_basis = jobrec['options'].scroll['GAMESS']['BASIS__GBASIS'].value
#    #if core.get_global_option('BASIS') == '':
    if _qcdb_basis == '':
        raise ValueError('gamess bas not impl')
#        _, cased_basis = moptions.format_option_for_cfour('CFOUR_BASIS', _cfour_basis)
#        cfourrec['genbas'] = extract_basis_from_genbas(cased_basis, jobrec['molecule']['elem'], exact=False)
#        bascmd = ''
    #else:
    qbs = BasisSet.pyconstruct(jobrec['molecule'], 'BASIS', _qcdb_basis)

#        #if qbs.has_ECP(): #    raise ValidationError("""ECPs not hooked up for Cfour""")
    molbascmd = muster_and_format_molecule_and_basis_for_gamess(jobrec['molecule'], jobrec['options'], qbs, verbose=1)

#    print(molbascmd)

    # Handle calc type and quantum chemical method
    nel = sum([z * int(real) for z, real in zip(jobrec['molecule']['elez'], jobrec['molecule']['real'])]) - jobrec['molecule']['molecular_charge']
    nel = int(nel)
    nfzc = qmol.nfrozen_core(depth=True)  # this will be default FC  # TODO change these values when user sets custom FC
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

    muster_modelchem(jobrec['method'], jobrec['dertype'], jobrec['options'], sysinfo)

    #print('HH')
    #print(jobrec['options'].print_changed(history=False))

#    # Handle driver vs input/default keyword reconciliation
#
    # Handle conversion of qcdb keyword structure into gamess format
    resolved_options = {k: v.value for k, v in jobrec['options'].scroll['GAMESS'].items() if v.disputed()}
    optcmd = format_options_for_gamess(resolved_options)

#    # Handle text to be passed untouched to cfour
#    litcmd = core.get_global_option('LITERAL_CFOUR')

    # Assemble input pieces
    gamessrec['gamess.inp'] = optcmd + molbascmd

    return gamessrec


#def cfour_harvest(jobrec, cfourrec):  # jobrec@i, cfourrec@io -> jobrec@io
#    """Processes raw results from read-only `cfourrec` into QCAspect fields in returned `jobrec`."""
#
#    try:
#        pass
#        #jobrec['molecule']['real']
#        #jobrec['do_gradient']
#    except KeyError as err:
#        raise KeyError(
#            'Required fields missing from ({})'.format(jobrec.keys())) from err
#
#    try:
#        cfourrec['stdout']
#        #if jobrec['do_gradient'] is True:
#        #    dftd3rec['dftd3_gradient']
#    except KeyError as err:
#        raise KeyError('Required fields missing from ({})'.format(
#            cfourrec.keys())) from err
#
#    # amalgamate output
#    text = cfourrec['stdout']
#    text += '\n  <<<  Cfour {} {} Results  >>>\n\n'.format('', '') #name.lower(), calledby.capitalize()))  # banner()
#
#    c4files = {}
#    for fl in ['GRD', 'FCMFINAL', 'DIPOL']:
#        field = 'output_' + fl.lower()
#        if field in cfourrec:
#            text += '  Cfour scratch file {} has been read\n'.format(fl)
#            text += cfourrec[field]
#            c4files[fl] = cfourrec[field]
#
#
#    #if molecule.name() == 'blank_molecule_psi4_yo':
#    #    qcdbmolecule = None
#    #else:
#    #    molecule.update_geometry()
#    #    qcdbmolecule = qcdb.Molecule(molecule.create_psi4_string_from_molecule())
#    #    qcdbmolecule.update_geometry()
#    qmol = Molecule(jobrec['molecule'])
#
#    # c4mol, if it exists, is dinky, just a clue to geometry of cfour results
#    psivar, c4hess, c4grad, c4mol, version, errorTMP = harvester.harvest(qmol, cfourrec['stdout'], **c4files)
#
#    jobrec['error'] += errorTMP
#    # Absorb results into psi4 data structures
#    #for key in psivar.keys():
#    #    core.set_variable(key.upper(), float(psivar[key]))
#    #calcinfo = qcvars.certify_qcvars(psivar)
#    #jobrec['qcvars'] = {info.lbl: info for info in calcinfo}
#    progvars = PreservingDict(psivar)
#
#    #if qcdbmolecule is None and c4mol is not None:
#    #    molecule = geometry(c4mol.create_psi4_string_from_molecule(), name='blank_molecule_psi4_yo')
#    #    molecule.update_geometry()
#    #    # This case arises when no Molecule going into calc (cfour {} block) but want
#    #    #   to know the orientation at which grad, properties, etc. are returned (c4mol).
#    #    #   c4mol is dinky, w/o chg, mult, dummies and retains name
#    #    #   blank_molecule_psi4_yo so as to not interfere with future cfour {} blocks
#
#    if c4grad is not None:
#        progvars['CURRENT GRADIENT'] = c4grad
#        #mat = core.Matrix.from_list(c4grad)
#        #core.set_gradient(mat)
#
#    if c4hess is not None:
#        progvars['CURRENT HESSIAN'] = c4hess
#
#    # badly placed
#    # Cfour's SCS-MP2 is non adjustible and only valid for UHF
#    # ROMP2 doesn't print SS & OS
#    if "MP2 OPPOSITE-SPIN CORRELATION ENERGY" in progvars and "MP2 SAME-SPIN CORRELATION ENERGY" in progvars:
#        oss_opt = jobrec['options'].scroll['QCDB']['MP2_OS_SCALE']
#        sss_opt = jobrec['options'].scroll['QCDB']['MP2_SS_SCALE']
#        custom_scsmp2_corl = \
#            Decimal(oss_opt.value) * progvars["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] + \
#            Decimal(sss_opt.value) * progvars["MP2 SAME-SPIN CORRELATION ENERGY"]
#        if "MP2 SINGLES ENERGY" in progvars:
#            custom_scsmp2_corl += progvars["MP2 SINGLES ENERGY"]
#        progvars["CUSTOM SCS-MP2 CORRELATION ENERGY"] = custom_scsmp2_corl
#
#    qcvars.build_out(progvars)
#    calcinfo = qcvars.certify(progvars)
#    text += print_variables(calcinfo)
#
#    jobrec['raw_output'] = text
#    jobrec['qcvars'] = calcinfo
#
#    prov = {}
#    prov['creator'] = 'Cfour'
#    prov['routine'] = sys._getframe().f_code.co_name
#    prov['version'] = version
#    jobrec['provenance'].append(prov)
#
#    return jobrec
