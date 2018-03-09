import sys
import copy
import pprint
pp = pprint.PrettyPrinter(width=120)
import inspect

from .. import __version__
from .. import molparse
from .. import moptions
from .. import qcvars
from ..driver.driver_helpers import print_variables
from ..exceptions import *
from ..iface_psi4.options import query_options_defaults_from_psi
from ..libmintsbasisset import BasisSet
from ..molecule import Molecule
from ..pdict import PreservingDict
from . import harvester
from .worker import cfour_subprocess
from .bas import extract_basis_from_genbas, format_basis_for_cfour, format_molecule_for_cfour


def run_cfour(name, molecule, options, **kwargs):
    print('\nhit run_cfour', name, kwargs)

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
    jobrec['provenance'] = [prov]

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

    jobrec['options'] = copy.deepcopy(options)


    print('comin in')
    #print(jobrec['options'])
    print(jobrec['options'].print_changed())

    #try:
    jobrec = cfour_driver(jobrec)

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

    print('CFOUR_DRIVER jr', jobrec)
    # this is what the cfour program needs, not what the job needs
    # *

    cfourrec = cfour_plant(jobrec)
    #cfourrec['scratch_messy'] = True

    # test json roundtrip
    jcfourrec = json.dumps(cfourrec)
    cfourrec = json.loads(jcfourrec)

    print('CFOURREC')
    pp.pprint(cfourrec)
    cfour_subprocess(cfourrec)  # updates cfourrec

    cfour_harvest(jobrec, cfourrec)  # ? updates jobrec ?

    print('QC', jobrec['qcvars'])
    return jobrec


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

    print('in cfour_plant')
    pp.pprint(jobrec)

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

    print('HH')
    print(jobrec['options'].print_changed())

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
    print('<<< ZMAT||{}||>>>\n'.format(zmat))
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
    """Processes raw results from read-only `cfourrec` into QCAspect fields in returned `jobrec`."""

    try:
        pass
        #jobrec['molecule']['real']
        #jobrec['do_gradient']
    except KeyError as err:
        raise KeyError(
            'Required fields missing from ({})'.format(jobrec.keys())) from err

    try:
        cfourrec['stdout']
        #if jobrec['do_gradient'] is True:
        #    dftd3rec['dftd3_gradient']
    except KeyError as err:
        raise KeyError('Required fields missing from ({})'.format(
            cfourrec.keys())) from err

    # amalgamate output
    text = cfourrec['stdout']
    text += '\n  <<<  Cfour {} {} Results  >>>\n\n'.format('', '') #name.lower(), calledby.capitalize()))  # banner()

    c4files = {}
    for fl in ['GRD', 'FCMFINAL', 'DIPOL']:
        field = 'output_' + fl.lower()
        if field in cfourrec:
            text += '  Cfour scratch file {} has been read\n'.format(fl)
            text += cfourrec[field]
            c4files[fl] = cfourrec[field]


    #if molecule.name() == 'blank_molecule_psi4_yo':
    #    qcdbmolecule = None
    #else:
    #    molecule.update_geometry()
    #    qcdbmolecule = qcdb.Molecule(molecule.create_psi4_string_from_molecule())
    #    qcdbmolecule.update_geometry()
    qmol = Molecule(jobrec['molecule'])

    # c4mol, if it exists, is dinky, just a clue to geometry of cfour results
    psivar, c4hess, c4grad, c4mol, version, errorTMP = harvester.harvest(qmol, cfourrec['stdout'], **c4files)

    print ('errorTMP', errorTMP)
    jobrec['error'] += errorTMP
    # Absorb results into psi4 data structures
    #for key in psivar.keys():
    #    core.set_variable(key.upper(), float(psivar[key]))
    #calcinfo = qcvars.certify_qcvars(psivar)
    #jobrec['qcvars'] = {info.lbl: info for info in calcinfo}
    progvars = PreservingDict(psivar)

    #if qcdbmolecule is None and c4mol is not None:
    #    molecule = geometry(c4mol.create_psi4_string_from_molecule(), name='blank_molecule_psi4_yo')
    #    molecule.update_geometry()
    #    # This case arises when no Molecule going into calc (cfour {} block) but want
    #    #   to know the orientation at which grad, properties, etc. are returned (c4mol).
    #    #   c4mol is dinky, w/o chg, mult, dummies and retains name
    #    #   blank_molecule_psi4_yo so as to not interfere with future cfour {} blocks

    if c4grad is not None:
        progvars['CURRENT GRADIENT'] = c4grad
        #mat = core.Matrix.from_list(c4grad)
        #core.set_gradient(mat)

    if c4hess is not None:
        progvars['CURRENT HESSIAN'] = c4hess

    qcvars.build_out(progvars)
    calcinfo = qcvars.certify(progvars)
    text += print_variables(calcinfo)

    jobrec['raw_output'] = text
    jobrec['qcvars'] = calcinfo

    prov = {}
    prov['creator'] = 'Cfour'
    prov['routine'] = sys._getframe().f_code.co_name
    prov['version'] = version
    jobrec['provenance'].append(prov)

    return jobrec



    """
    Required Input Fields
    ---------------------

    Optional Input Fields
    ---------------------

    Output Fields
    -------------

    """


