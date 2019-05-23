import sys
import copy
import pprint
pp = pprint.PrettyPrinter(width=120)
import inspect
from decimal import Decimal

import qcengine as qcng

from .. import __version__
from .. import moptions
from .. import qcvars
from ..driver.driver_helpers import print_variables
from ..exceptions import *
from ..intf_psi4.options import query_options_defaults_from_psi
from ..libmintsbasisset import BasisSet
from ..molecule import Molecule
from ..pdict import PreservingDict
from . import harvester
from .worker import cfour_subprocess
from .bas import extract_basis_from_genbas, format_basis_for_cfour, format_molecule_for_cfour


def run_cfour(name, molecule, options, **kwargs):
    #print('\nhit run_cfour', name, kwargs)

    #calledby = inspect.stack()
    #print('CALLEDBY')
    #for cur in calledby:
    #    print('CUR', cur[3])

    jobrec = {}
    jobrec['error'] = ''
    jobrec['success'] = None
    prov = {}
    prov
    prov['creator'] = 'QCDB'
    prov['version'] = __version__
    prov['routine'] = sys._getframe().f_code.co_name
#    jobrec['provenance'] = [prov]

    #jobrec['molecule'] = {'qm': molecule.to_dict(np_out=False)}
    jobrec['molecule'] = molecule.to_dict(np_out=False)
    jobrec['method'] = name
    jobrec['dertype'] = ['energy', 'gradient', 'hessian'].index(inspect.stack()[1][3])

    #if 'MEMORY' in options['GLOBALS']:
    #    mem = options['GLOBALS'].pop('MEMORY')
    #    print('MEM', mem)
    #    jobrec['memory'] = mem['value']

#    popts = collections.defaultdict(lambda: collections.defaultdict(dict))
#    #for k, v in options['GLOBALS'].items():
#    for k, v in options['CFOUR'].items():
#        print('C4 opt', k, v)
#   #     psi4.core.set_global_option(k.upper(), v['value'])
#        popts[k] = v['value']
#    jobrec['options'] = popts

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

    #print('comin in')
    #print(jobrec['options'])
    #PRprint(jobrec['options'].print_changed())

    #try:
    jobrec = cfour_driver(jobrec)

    pp.pprint(jobrec)
    return jobrec


def cfour_driver(jobrec):
    import json

    try:
        #jobrec['dashlevel']
        #jobrec['dashparams']
        #jobrec['functional']
        jobrec['molecule']
        #jobrec['do_gradient']
    except KeyError as err:
        #raise KeyError(
        #    'Required fields missing from ({})'.format(jobrec.keys())) from err
        jobrec['error'] += repr(err) + 'Required fields missing from ({})'.format(jobrec.keys())
        return jobrec

    #print('CFOUR_DRIVER jr', jobrec)
    ## this is what the cfour program needs, not what the job needs
    ## *

    #cfourrec = cfour_plant(jobrec)
    ##cfourrec['scratch_messy'] = True

    ## test json roundtrip
    #jcfourrec = json.dumps(cfourrec)
    #cfourrec = json.loads(jcfourrec)

    #print('CFOURREC')
    #pp.pprint(cfourrec)
    #cfour_subprocess(cfourrec)  # updates cfourrec

    #cfour_harvest(jobrec, cfourrec)  # ? updates jobrec ?

    #print('QC', jobrec['qcvars'])
    #return jobrec

    verbose = 1

    if verbose >= 3:
        print('[1] {} JOBREC PRE-PLANT (j@i) <<<'.format('CFOUR'))
        pp.pprint(jobrec)
        print('>>>')

    cfourrec = cfour_plant(jobrec)

    # test json roundtrip
    jcfourrec = json.dumps(cfourrec)
    cfourrec = json.loads(jcfourrec)

    if verbose >= 4:
        print('[2] {}REC PRE-SUBPROCESS (m@i) <<<'.format('CFOUR'))
        pp.pprint(cfourrec)
        print('>>>\n')

#    cfour_subprocess(cfourrec)  # updates cfourrec


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


    qcschema_input = jobrec['qcschema_input']
    qcschema_input['extras']['infiles'] = {
        'ZMAT': cfourrec['zmat'],
        'GENBAS': cfourrec['genbas'],
    }
    ret = qcng.compute(qcschema_input, 'cfour').dict()
    jobrec.pop('qcschema_input')
    if verbose >= 4:
        print('[3] {}REC POST-ENGINE (m@io) <<<'.format('CFOUR'))
        pp.pprint(ret)
        print('>>>\n')

#    cfourrec['stdout'] = ret['stdout']
    jobrec['raw_output'] = ret['stdout']
    jobrec['provenance'] = ret['provenance']
    cfourrec['outfiles'] = ret['extras']['outfiles']
    cfourrec['progvars'] = ret['extras']['qcvars']
    #qcvars.build_out(progvars)

    if verbose >= 5:
        print('[3] {}REC POST-SUBPROCESS (m@io) <<<'.format('CFOUR'))
        pp.pprint(cfourrec)
        print('>>>\n')

    cfour_harvest(jobrec, cfourrec)  # updates jobrec

    if verbose >= 2:
        print('[4] {} JOBREC POST-HARVEST (j@io) <<<'.format('CFOUR'))
        pp.pprint(jobrec)
        print('>>>')

    return jobrec


    #ret = qcng.compute(qcschema_input, 'gamess').dict()
    #jobrec.pop('qcschema_input')
    #progvars = PreservingDict(ret['extras']['qcvars'])
    #qcvars.build_out(progvars)
    #calcinfo = qcvars.certify(progvars, plump=True, nat=len(jobrec['molecule']['mass']))
    #jobrec['raw_output'] = ret['stdout']
    #jobrec['qcvars'] = calcinfo
    #jobrec['success'] = True



def cfour_plant(jobrec):  # jobrec@i -> cfour@i
    cfourrec = {}

    #return cfourrec

#def write_zmat(name, dertype, molecule):  # -> zmat
    """Returns string with contents of Cfour ZMAT file as gathered from
    active molecule, current keyword settings, and cfour {...} block.

    """
    import qcdb
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

    #print('in cfour_plant')
    #pp.pprint(jobrec)

    molcmd = format_molecule_for_cfour(jobrec['molecule'], jobrec['options'], verbose=1)

    # Handle qcdb keywords implying cfour keyword values
    # if core.get_option('CFOUR', 'TRANSLATE_PSI4'):
    harvester.muster_inherited_options(jobrec['options'])

    _qcdb_basis = jobrec['options'].scroll['QCDB']['BASIS'].value
    _cfour_basis = jobrec['options'].scroll['CFOUR']['BASIS'].value
    #if core.get_global_option('BASIS') == '':
    if _qcdb_basis == '':
        _, cased_basis = moptions.format_option_for_cfour('CFOUR_BASIS', _cfour_basis)
        cfourrec['genbas'] = extract_basis_from_genbas(cased_basis, jobrec['molecule']['elem'], exact=False)
        bascmd = ''
    else:
        qbs = BasisSet.pyconstruct(jobrec['molecule'], 'BASIS', _qcdb_basis)
        #if qbs.has_ECP(): #    raise ValidationError("""ECPs not hooked up for Cfour""")
        cfourrec['genbas'] = qbs.print_detail_cfour() #qbs.genbas()
        bascmd = format_basis_for_cfour(jobrec['molecule'], jobrec['options'], qbs.has_puream())

    # Handle calc type and quantum chemical method
    harvester.nu_muster_modelchem(jobrec['method'], jobrec['dertype'], jobrec['options'])

    #print('HH')
    #print(jobrec['options'].print_changed(history=False))

    # Handle driver vs input/default keyword reconciliation

    # Handle conversion of psi4 keyword structure into cfour format
    optcmd = moptions.prepare_options_for_cfour(jobrec['options'])

    # Handle text to be passed untouched to cfour
#    litcmd = core.get_global_option('LITERAL_CFOUR')


#    popts = {}
#    for k, v in options.scroll['QCDB'].items():
#        if not v.is_default():
#            print('QQQQ', k, v.value, v.is_default())
#            popts[k] = v.value
#
#    for k, v in options.scroll['PSI4'].items():
#        if not v.is_default():
#            print('PPPP', k, v.value, v.is_default())
#            popts[k] = v.value
#    jobrec['options'] = popts

    # Assemble ZMAT pieces
    #zmat = memcmd + molcmd + optcmd + mdccmd + psicmd + bascmd + litcmd
    zmat = molcmd + optcmd + bascmd
    cfourrec['zmat'] = zmat
    #print('<<< ZMAT||{}||>>>\n'.format(zmat))
    cfourrec['command'] = ['xcfour']

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
#    cfourrec['zmat'] = tmp_zmat
#    print('<<< TMP_ZMAT\n', tmp_zmat, '\n>>>\n')

#    if len(re.findall(r'^\*(ACES2|CFOUR|CRAPS)\(', zmat, re.MULTILINE)) != 1:
#        core.print_out('\n  Faulty ZMAT constructed:\n%s' % (zmat))
#        raise ValidationError("""
#Multiple *CFOUR(...) blocks in input. This usually arises
#because molecule or options are specified both the psi4 way through
#molecule {...} and set ... and the cfour way through cfour {...}.""")

    #return zmat
    return cfourrec


def cfour_harvest(jobrec, cfourrec):  # jobrec@i, cfourrec@io -> jobrec@io

    progvars = PreservingDict(cfourrec['progvars'])

    # badly placed
    # Cfour's SCS-MP2 is non adjustible and only valid for UHF
    # ROMP2 doesn't print SS & OS
    if "MP2 OPPOSITE-SPIN CORRELATION ENERGY" in progvars and "MP2 SAME-SPIN CORRELATION ENERGY" in progvars:
        oss_opt = jobrec['options'].scroll['QCDB']['MP2_OS_SCALE']
        sss_opt = jobrec['options'].scroll['QCDB']['MP2_SS_SCALE']
        custom_scsmp2_corl = \
            Decimal(oss_opt.value) * progvars["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] + \
            Decimal(sss_opt.value) * progvars["MP2 SAME-SPIN CORRELATION ENERGY"]
        if "MP2 SINGLES ENERGY" in progvars:
            custom_scsmp2_corl += progvars["MP2 SINGLES ENERGY"]
        progvars["CUSTOM SCS-MP2 CORRELATION ENERGY"] = custom_scsmp2_corl

    qcvars.build_out(progvars)
    calcinfo = qcvars.certify(progvars, plump=True, nat=len(jobrec['molecule']['mass']))
    jobrec['qcvars'] = calcinfo

    jobrec['success'] = True

    return jobrec

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

