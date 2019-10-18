from typing import Any, Dict, Optional
import sys
import copy
import pprint
pp = pprint.PrettyPrinter(width=120)
import inspect
#from decimal import Decimal

import qcelemental as qcel
from qcelemental.models import ResultInput
from qcelemental.util import which

import qcengine as qcng
from qcengine.exceptions import InputError
from qcengine.programs.util import PreservingDict
from qcengine.programs.gamess import GAMESSHarness
from qcengine.programs.gamess.keywords import format_keywords

#from .. import moptions
from ... import qcvars
#from ..driver.driver_helpers import print_variables
#from ..exceptions import *
#from ..iface_psi4.options import query_options_defaults_from_psi
from ...libmintsbasisset import BasisSet
from ...util import provenance_stamp
from ...molecule import Molecule
from ...pdict import PreservingDict
from .molbasopt import muster_and_format_molecule_and_basis_for_gamess
from .harvester import muster_modelchem, muster_inherited_options


#def run_gamess_old(name, molecule, options, **kwargs):
#
#    prov = {}
#    prov['creator'] = 'QCDB'
#    prov['version'] = __version__
#    prov['routine'] = sys._getframe().f_code.co_name
#    
#
#    qcschema_input = {
#        'schema_name': 'qcschema_input',
#        'schema_version': 1,
#        'driver': inspect.stack()[1][3],
#        'model': {
#            'method': name,
#            'basis': '(auto)',
#        },
#        'molecule': molecule.to_schema(dtype=2),
#        'extras': {},
#    }
#
#    jobrec = {}
#    jobrec['provenance'] = [prov]
#    jobrec['error'] = ''
#    jobrec['success'] = None
#    jobrec['molecule'] = molecule.to_dict(np_out=False)
#    jobrec['method'] = name
#    jobrec['dertype'] = ['energy', 'gradient', 'hessian'].index(inspect.stack()[1][3])
#
#
#    jobrec['options'] = copy.deepcopy(options)
#    jobrec['qcschema_input'] = qcschema_input
#
#    gamessrec = {}
#    muster_inherited_options(jobrec['options'])
#
#    qmol = Molecule(jobrec['molecule'])
#    _qcdb_basis = jobrec['options'].scroll['QCDB']['BASIS'].value
#    _gamess_basis = jobrec['options'].scroll['GAMESS']['BASIS__GBASIS'].value
#    if _qcdb_basis == '':
#        raise ValueError('gamess bas not impl')
#    qbs = BasisSet.pyconstruct(jobrec['molecule'], 'BASIS', _qcdb_basis)
#
#    molbascmd = muster_and_format_molecule_and_basis_for_gamess(jobrec['molecule'], jobrec['options'], qbs, verbose=1)
#
#
#    # Handle calc type and quantum chemical method
#    nel = sum([z * int(real) for z, real in zip(jobrec['molecule']['elez'], jobrec['molecule']['real'])]) - jobrec['molecule']['molecular_charge']
#    nel = int(nel)
#    nfzc = qmol.nfrozen_core(depth=True)  # this will be default FC  # TODO change these values when user sets custom FC
#    nfzc = 0
#
#    # forcing nfc above. all these need to be reocmputed together for a consistent cidet input group
#    nels = nel - 2 * nfzc
#    nact = qbs.nbf() - nfzc
#    sysinfo = {
#        'nel': nel,
#        'ncore': nfzc,
#        'nact': nact,
#        'nels': nels,
#    }
#
#    muster_modelchem(jobrec['method'], jobrec['dertype'], jobrec['options'], sysinfo)
#
#    # Handle conversion of qcdb keyword structure into gamess format
#    resolved_options = {k: v.value for k, v in jobrec['options'].scroll['GAMESS'].items() if v.disputed()}
#    optcmd = format_options_for_gamess(resolved_options)
#
#    # Assemble input pieces
#    gamessrec['gamess.inp'] = optcmd + molbascmd
#
#    qcschema_input = jobrec['qcschema_input']
#    qcschema_input['extras']['gamess.inp'] = gamessrec['gamess.inp']
#
#
#    ret = qcng.compute(qcschema_input, 'gamess').dict()
#    jobrec.pop('qcschema_input')
#    progvars = PreservingDict(ret['extras']['qcvars'])
#    qcvars.build_out(progvars)
#    calcinfo = qcvars.certify(progvars, plump=True, nat=len(jobrec['molecule']['mass']))
#    jobrec['raw_output'] = ret['stdout']
#    jobrec['qcvars'] = calcinfo
#    jobrec['success'] = True
#
#    return jobrec


def run_gamess(name, molecule, options, **kwargs):
#def run_gamess2(name, molecule, options, **kwargs):

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

    jobrec = qcng.compute(resi, "qcdb-gamess", raise_error=True).dict()
    hold_qcvars = jobrec['extras'].pop('qcdb:qcvars')
    jobrec['qcvars'] = {key: qcel.Datum(**dval) for key, dval in hold_qcvars.items()}
    return jobrec

def _print_helper(label, dicary, do_print):
    if do_print:
        print(label + ' <<<')
        pp.pprint(dicary)
        print('>>>')


class QcdbGAMESSHarness(GAMESSHarness):
    def compute(self, input_model: 'ResultInput', config: 'JobConfig') -> 'Result':
        self.found(raise_error=True)

        verbose = 1

        _print_helper(f'[1] {self.name} RESULTINPUT PRE-PLANT', input_model.dict(), verbose >= 3)

        job_inputs = self.qcdb_build_input(input_model, config)

        _print_helper(f'[2] {self.name}REC PRE-ENGINE', job_inputs, verbose >= 4)

        success, dexe = self.execute(job_inputs)

        _print_helper(f'[3] {self.name}REC POST-ENGINE', dexe, verbose >= 4)

        if 'INPUT HAS AT LEAST ONE SPELLING OR LOGIC MISTAKE' in dexe["stdout"]:
            raise InputError(dexe["stdout"])

        if not success:
            output_model = input_model
            output_model["error"] = {"error_type": "execution_error", "error_message": dexe["stderr"]}

        dexe["outfiles"]["stdout"] = dexe["stdout"]
        dexe["outfiles"]["stderr"] = dexe["stderr"]
        output_model = self.parse_output(dexe["outfiles"], input_model)

        _print_helper(f'[4a] {self.name} RESULT POST-HARVEST', output_model.dict(), verbose >= 5)

        output_model = self.qcdb_post_parse_output(input_model, output_model)

        _print_helper(f'[4] {self.name} RESULT POST-POST-HARVEST', output_model.dict(), verbose >= 2)

        return output_model

    def qcdb_build_input(self, input_model: 'ResultInput', config: 'JobConfig',
                         template: Optional[str] = None) -> Dict[str, Any]:
        gamessrec = {
            'infiles': {},
            'scratch_directory': config.scratch_directory,
        }

        molrec = qcel.molparse.from_schema(input_model.molecule.dict())
        molrecc1 = molrec.copy()
        molrecc1['fix_symmetry'] = 'c1'  # write _all_ atoms to input
        ropts = input_model.extras['qcdb:options']

        # Handle qcdb keywords implying gamess keyword values
        muster_inherited_options(ropts)

        _qcdb_basis = ropts.scroll['QCDB']['BASIS'].value
        #_gamess_basis = ropts.scroll['GAMESS']['BASIS'].value
        qbs = BasisSet.pyconstruct(molrecc1, 'BASIS', _qcdb_basis)

        molbascmd = muster_and_format_molecule_and_basis_for_gamess(molrec, ropts, qbs)

        nel = sum([z * int(real) for z, real in zip(molrec['elez'], molrec['real'])]) - molrec['molecular_charge']
        nel = int(nel)
        #nfzc = input_model.molecule.nfrozen_core(depth=True)  # this will be default FC  # TODO change these values when user sets custom FC
        # not sure how to do this???
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

        # Handle calc type and quantum chemical method
        muster_modelchem(input_model.model.method, input_model.driver.derivative_int(), ropts, sysinfo)


        # Handle conversion of psi4 keyword structure into cfour format
        skma_options = {key: ropt.value for key, ropt in sorted(ropts.scroll['GAMESS'].items()) if ropt.disputed()}

        optcmd = format_keywords(skma_options)

        gamessrec['infiles']['gamess.inp'] = optcmd + molbascmd
        gamessrec['command'] = [which("rungms"), "gamess"]

        return gamessrec

    def qcdb_post_parse_output(self, input_model: 'ResultInput', output_model: 'Result') -> 'Result':

        dqcvars = PreservingDict(copy.deepcopy(output_model.extras['qcvars']))
        qcvars.build_out(dqcvars)
        calcinfo = qcvars.certify(dqcvars, plump=True, nat=len(output_model.molecule.symbols))
        output_model.extras['qcdb:qcvars'] = calcinfo

        return output_model



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
