import sys
import copy
import pprint
pp = pprint.PrettyPrinter(width=120)
import inspect

from .. import __version__
from .. import moptions
from .. import qcvars
from ..driver.driver_helpers import print_variables
from ..exceptions import *
#from ..iface_psi4.options import query_options_defaults_from_psi
from ..libmintsbasisset import BasisSet
from ..molecule import Molecule
from ..pdict import PreservingDict
from . import harvester
from .worker import nwchem_subprocess
#from .molbas import * #format_molecule_for_nwchem, format_basis_for_nwchem, format_basis_for_nwchem_puream, local_prepare_options_for_modules
from .molbasopt import muster_and_format_molecule_for_nwchem, muster_and_format_basis_for_nwchem, format_options_for_nwchem
from .harvester import muster_inherited_options, format_modelchem_for_nwchem
from ..moptions.options import reconcile_options


def run_nwchem(name, molecule, options, **kwargs):
    print('\nhit run_nwchem', name, kwargs)

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

    jobrec['molecule'] = molecule.to_dict(np_out=False)
    jobrec['method'] = name
    jobrec['dertype'] = ['energy', 'gradient', 'hessian'].index(inspect.stack()[1][3])

    jobrec['options'] = copy.deepcopy(options)

    #print('comin in')
    #print(jobrec['options'])
    #PRprint(jobrec['options'].print_changed())

    #try:
    jobrec = nwchem_driver(jobrec)

    return jobrec


def nwchem_driver(jobrec):
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

#    print('[1] NWCHEM JOBREC PRE-PLANT (j@i) <<<')
#    pp.pprint(jobrec)
#    print('>>>')

    nwchemrec = nwchem_plant(jobrec)

    # test json roundtrip
    jnwchemrec = json.dumps(nwchemrec)
    nwchemrec = json.loads(jnwchemrec)

#    print('[2] NWCHEMREC PRE-SUBPROCESS (e@i) <<<')
#    pp.pprint(nwchemrec)
#    print('>>>\n')

    nwchem_subprocess(nwchemrec)  # updates nwchemrec

#    print('[3] NWCHEMREC POST-SUBPROCESS (e@io) <<<')
#    pp.pprint(nwchemrec)
#    print('>>>\n')

    nwchem_harvest(jobrec, nwchemrec)  # updates jobrec

#    print('[4] NWCHEM JOBREC POST-HARVEST (j@io) <<<')
#    pp.pprint(jobrec)
#    print('>>>')

    return jobrec


def nwchem_plant(jobrec):  # jobrec@i -> engine@i
    nwchemrec = {}

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

    print('in nwchem_plant')

    # Handle qcdb keywords implying nwchem keyword values
    muster_inherited_options(jobrec['options'])

    molcmd = muster_and_format_molecule_for_nwchem(jobrec['molecule'], jobrec['options'], verbose=1)

    _qcdb_basis = jobrec['options'].scroll['QCDB']['BASIS'].value
    _gamess_basis = jobrec['options'].scroll['NWCHEM']['BASIS'].value
    if _qcdb_basis == '':
        raise ValueError('nwchem bas not impl. set with `basis cc-pvdz`, etc. to use qcdb (psi4) basis set library.')
    # create a qcdb.Molecule to reset PG to c1 so all atoms. but print_detail isn't transmitting in a way nwc is picking up on, so per-element for now
    #qmol = Molecule(jobrec['molecule'])
    #qmol.reset_point_group('c1')  # need basis printed for *every* atom
    #qbs = BasisSet.pyconstruct(qmol, 'BASIS', _qcdb_basis)
    qbs = BasisSet.pyconstruct(jobrec['molecule'], 'BASIS', _qcdb_basis)

    #if qbs.has_ECP(): #    raise ValidationError("""ECPs not hooked up for Cfour""")
    bascmd = muster_and_format_basis_for_nwchem(jobrec['molecule'], jobrec['options'], qbs, verbose=1)

    # Handle calc type and quantum chemical method
  #  harvester.nu_muster_modelchem(jobrec['method'], jobrec['dertype'], jobrec['options'])
    mdccmd = format_modelchem_for_nwchem(jobrec['method'], jobrec['dertype'], jobrec['options'], sysinfo=None)

    #PRprint('Touched Keywords')
    #PRprint(jobrec['options'].print_changed(history=False))

    # Handle driver vs input/default keyword reconciliation

    # Handle conversion of qcdb keyword structure into nwchem format
#OLD    optcmd = moptions.prepare_options_for_nwchem(jobrec['options'])
    resolved_options = {k: v.value for k, v in jobrec['options'].scroll['NWCHEM'].items() if v.disputed()}
    optcmd = format_options_for_nwchem(resolved_options)

    # Handle text to be passed untouched to cfour
#    litcmd = core.get_global_option('LITERAL_CFOUR')

    # Assemble ZMAT pieces
    #zmat = memcmd + molcmd + optcmd + mdccmd + psicmd + bascmd + litcmd
    #zmat = molcmd + optcmd + bascmd
    nwchemrec['nwchem.nw'] = 'echo\n' + molcmd + bascmd + optcmd + mdccmd
#OLD    nwchemrec['nwchem.nw'] = write_input(jobrec['method'], jobrec['dertype'], jobrec['molecule'], jobrec['options']) #molecule)
    #print('<<< ZMAT||\n{}\n||>>>\n'.format(nwchemrec['nwchem.nw']))
    nwchemrec['command'] = ['nwchem']  # subnw?

    return nwchemrec




def nwchem_harvest(jobrec, nwchemrec):  # jobrec@i, enginerec@io -> jobrec@io
    """Processes raw results from read-only `nwchemrec` into QCAspect fields in returned `jobrec`."""

    try:
        pass
        #jobrec['molecule']['real']
        #jobrec['do_gradient']
    except KeyError as err:
        raise KeyError(
            'Required fields missing from ({})'.format(jobrec.keys())) from err

    try:
        nwchemrec['stdout']
        #if jobrec['do_gradient'] is True:
        #    dftd3rec['dftd3_gradient']
    except KeyError as err:
        raise KeyError('Required fields missing from ({})'.format(
            nwchemrec.keys())) from err

    # amalgamate output
    text = nwchemrec['stdout']
    text += '\n  <<<  NWChem {} {} Results  >>>\n\n'.format('', '') #name.lower(), calledby.capitalize()))  # banner()

    nwfiles = {}
    #for fl in ['GRD', 'FCMFINAL', 'DIPOL']:
    #    field = 'output_' + fl.lower()
    #    if field in nwchemrec:
    #        text += '  NWChem scratch file {} has been read\n'.format(fl)
    #        text += nwchemrec[field]
    #        nwfiles[fl] = nwchemrec[field]


    #if molecule.name() == 'blank_molecule_psi4_yo':
    #    qcdbmolecule = None
    #else:
    #    molecule.update_geometry()
    #    qcdbmolecule = qcdb.Molecule(molecule.create_psi4_string_from_molecule())
    #    qcdbmolecule.update_geometry()
    qmol = Molecule(jobrec['molecule'])

    # nwmol, if it exists, is dinky, just a clue to geometry of nwchem results
    psivar, nwhess, nwgrad, nwmol, version, errorTMP = harvester.harvest(qmol, nwchemrec['stdout'], **nwfiles)

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

    if nwgrad is not None:
        progvars['CURRENT GRADIENT'] = nwgrad

    if nwhess is not None:
        progvars['CURRENT HESSIAN'] = nwhess

    qcvars.build_out(progvars)
    calcinfo = qcvars.certify(progvars)
    text += print_variables(calcinfo)

    jobrec['raw_output'] = text
    jobrec['qcvars'] = calcinfo

    prov = {}
    prov['creator'] = 'NWChem'
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




def write_input(name, dertype, molecule, ropts):
    """Returns string with contents of NWChem INPUT file as gathered from
    active molecule, current keyword settings, and nwchem {...} block.
    """
    molecule = Molecule(molecule)
    qcdbmolecule = molecule

    mem = int(0.000001 * ropts.scroll['QCDB']['MEMORY'].value)
    if mem == 524:
        memcmd, memkw = '', {}

    else:
        memcmd, memkw = harvester.muster_memory(mem) #printing memory command here always 
        #print ("memcmd: ", memcmd) 
        #print ("memkw: ", memkw) 

    # Handle molecule and basis set
    if molecule.name() == 'blank_molecule_psi4_yo':
        molcmd, molkw = '', {}
        bascmd, baskw = '', {}
        bascmd2, baskw2 = '', {}
        #core.set_local_option('NWCHEM', 'TRANSLATE_PSI4', False)
    else:
        molecule.update_geometry()
#        print(molecule.create_psi4_string_from_molecule())
#        qcdbmolecule = Molecule(molecule.create_psi4_string_from_molecule())
#        qcdbmolecule.tagline = molecule.name()
        molcmd, molkw = qcdbmolecule.format_molecule_for_nwchem()
        #print ("molcmd: ", molcmd) #form: geometry units " unit(angstrom, bohr ...) "
        #print ("molkw: ", molkw)  #NWCHEM_CHARGE, 

        #if core.get_global_option('BASIS') == '':
        if ropts.scroll['QCDB']['BASIS'].value == '':
            bascmd, baskw = '', {}
            bascmd2, baskw2 = '', {}
            #print ("bascmd blank: ", bascmd)
            #print ("baskw blank: ", baskw)

        else:
            user_pg = molecule.schoenflies_symbol()
            molecule.reset_point_group('c1')  # need basis printed for *every* atom
            molecule.reset_point_group(user_pg)
            molecule.update_geometry()

#            qcdbmolecule = Molecule(jobrec['molecule'])
            _qcdb_basis = ropts.scroll['QCDB']['BASIS'].value
            qbs = BasisSet.pyconstruct(qcdbmolecule, 'BASIS', _qcdb_basis)
        ### New way to call basis sets
            tmp_gobas = _qcdb_basis if _qcdb_basis != '' else 'sto-3g'
            #tmp_gobas = core.get_global_option('BASIS') if core.get_global_option('BASIS') else 'sto-3g'
            tmp_basis = BasisSet.pyconstruct(qcdbmolecule, "BASIS", tmp_gobas)
            bascmd, baskw = qcdbmolecule.format_basis_for_nwchem_puream(qbs.has_puream())
            bascmd2, baskw2 = qcdbmolecule.format_basis_for_nwchem(tmp_basis)


#    _cfour_basis = jobrec['options'].scroll['CFOUR']['BASIS'].value
#    #if core.get_global_option('BASIS') == '':
#    if _qcdb_basis == '':
#        _, cased_basis = moptions.format_option_for_cfour('CFOUR_BASIS', _cfour_basis)
#        cfourrec['genbas'] = extract_basis_from_genbas(cased_basis, jobrec['molecule']['elem'], exact=False)
#        bascmd = ''
#    else:
#        qbs = BasisSet.pyconstruct(jobrec['molecule'], 'BASIS', _qcdb_basis)
#        #if qbs.has_ECP(): #    raise ValidationError("""ECPs not hooked up for Cfour""")
#        cfourrec['genbas'] = qbs.print_detail_cfour() #qbs.genbas()
#        bascmd = format_basis_for_cfour(jobrec['molecule'], jobrec['options'], qbs.has_puream())



        ### Old way (translating basis set name only into nwchem-speak)
            # bascmd2, baskw2 = qcdbmolecule.format_basis_for_nwchem(core.get_global_option('BASIS')) 
        ################

            #print ("bascmd: ", bascmd)
            #print ("bascmd2: ", bascmd2)
            #print ("baskw2: ", baskw2)
       
    # Handle psi4 keywords implying nwchem keyword values
    if ropts.scroll['QCDB']['TRANSLATE_QCDB'].value:
    #if core.get_option('NWCHEM', 'TRANSLATE_PSI4'):
        psicmd, psikw = harvester.muster_psi4options(local_prepare_options_for_modules(changedOnly=True), qcdbmolecule)
        #print ("psicmd: ", psicmd)
        #print ("psikw: ", psikw)
    else:
        psicmd, psikw = '', {}

    # Handle calc type and quantum chemical method
    mdccmd, mdckw = harvester.muster_modelchem(name, dertype)
    #print ("mdccmd: ", mdccmd)
    #print ("mdckw: ", mdckw)
    # Handle calc type and quantum chemical method
    mdccmd, mdckw = harvester.muster_modelchem(name, dertype)

    # Handle driver vs input/default keyword reconciliation
    userkw = local_prepare_options_for_modules() 
    userkw = reconcile_options(userkw, memkw) 
    userkw = reconcile_options(userkw, molkw)
    userkw = reconcile_options(userkw, baskw) 
    userkw = reconcile_options(userkw, baskw2)
    userkw = reconcile_options(userkw, psikw)
    userkw = reconcile_options(userkw, mdckw)

    # Handle conversion of psi4 keyword structure into nwchem format
    optcmd = prepare_options_for_nwchem(userkw)
    optcmd2 = prepare_options_for_task_nwchem(userkw)
    #print ("printing optcmd below")
    #print (optcmd) 

    # Handle text to be passed untouched to nwchem
 #   litcmd = core.get_global_option('LITERAL_NWCHEM')
    litcmd = ''

    # Assemble NWChem input pieces
    # bascmd = puream true or not [form: basis cartesian/spherical]
    # bascmd2 = basis sets for atoms
    # optcmd2 = task "theory" "dertype", This cmd has to be the last string all the time

    nwinput = 'echo\n' + memcmd + molcmd + bascmd + bascmd2 + psicmd + optcmd + litcmd + mdccmd +optcmd2 
   # print ("#########NWCHEM INPUT PRINT##########")
   # print (nwinput)
   # print ("#########NWCHEM INPUT END##########")

    return nwinput

    if "MP2 OPPOSITE-SPIN CORRELATION ENERGY" in progvars and "MP2 SAME-SPIN CORRELATION ENERGY" in progvars:
        oss_opt = jobrec['options'].scroll['QCDB']['MP2_OS_SCALE']
        sss_opt = jobrec['options'].scroll['QCDB']['MP2_SS_SCALE'] 
        custom_scsmp2_corl = \
                Decimal(oss_opt.value) * progvars["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] + \
                Decimal(sss_opt.value) * progvars["MP2 SAME-SPIN CORRELATION ENERGY"]
        if "MP2 SINGLES ENERGY" in progvars:
            custom_scsmp2_corl += progvars["MP2 SINGLES ENERGY"]
        progvars["CUSTOM SCS-MP2 CORRELATION ENERGY"] = custom_scsmp2_corl
