import re
import sys
import uuid
import struct
from collections import defaultdict
from decimal import Decimal

from ..pdict import PreservingDict
#from ..periodictable import *
#from ..physconst import *
from ..exceptions import *
from ..molecule import Molecule
from ..orient import OrientMols
from decimal import Decimal
from ..moptions.options import conv_float2negexp
#from . import nwc_movecs


def muster_inherited_options(ropts, verbose=1):
    accession = sys._getframe().f_code.co_name + '_' + str(uuid.uuid4())
    kwgs = {'accession': accession, 'verbose': verbose}

    do_translate = ropts.scroll['QCDB']['TRANSLATE_QCDB'].value

# REVISIT MEMORY
    # qcdb/memory [B] --> nwchem/total_memory [MB]
    qopt = ropts.scroll['QCDB']['MEMORY']
    if do_translate or qopt.is_required():
        mem = str(int(0.000001 * qopt.value)) + ' mb'
        ropts.suggest('NWCHEM', 'memory', mem, **kwgs)

#    # qcdb/puream --> cfour/spherical
#    ropts.suggest('CFOUR', 'SPHERICAL', ropts.scroll['QCDB']['PUREAM'].value, **kwgs)

    # qcdb/reference --> nwchem/scf_[r|u|ro]hf
    # TODO ref or scf__ref?
    qref = ropts.scroll['QCDB']['SCF__REFERENCE'].value
    print('<<<< QREF {} >>>'.format(qref))
    if qref in ['RHF', 'UHF', 'ROHF']:
        rhf = (qref == 'RHF')
        uhf = (qref == 'UHF')
        rohf = (qref == 'ROHF')
        ropts.suggest('NWCHEM', 'scf__rhf', rhf, **kwgs)
        ropts.suggest('NWCHEM', 'scf__uhf', uhf, **kwgs)
        ropts.suggest('NWCHEM', 'scf__rohf', rohf, **kwgs)

    # qcdb/scf__d_convergence --> nwchem/scf__thresh
    # TODO should be e_conv?
    qopt = ropts.scroll['QCDB']['SCF__D_CONVERGENCE']
    if qopt.disputed():
        conv = conv_float2negexp(qopt.value)
        ropts.suggest('NWCHEM', 'scf__thresh', conv, **kwgs)

#    # qcdb/scf__maxiter --> cfour/scf_maxcyc
#    ropts.suggest('CFOUR', 'SCF_MAXCYC', ropts.scroll['QCDB']['SCF__MAXITER'].value, **kwgs)
#
#    # qcdb/scf__damping_percentage --> cfour/scf_damping
#    damp = int(10 * ropts.scroll['QCDB']['SCF__DAMPING_PERCENTAGE'].value)
#    ropts.suggest('CFOUR', 'SCF_DAMPING', damp, **kwgs)

    # qcdb/e_convergence --> nwchem/ccsd(t)__thresh or dft__convergence__energy
    qopt = ropts.scroll['QCDB']['E_CONVERGENCE']
    if qopt.disputed():
        conv = conv_float2negexp(qopt.value)
        ropts.suggest('NWCHEM', 'ccsd__thresh', conv, **kwgs)
        ropts.suggest('NWCHEM', 'dft__convergence__energy', conv, **kwgs)


def muster_memory(mem):
    """Transform input *mem* in MB into nwchem-speak 
    
    """
    text = ''

    options = defaultdict(lambda: defaultdict(dict))
    options['NWCHEM']['NWCHEM_MEMORY']['value'] = [int(mem), 'mb']  #Total memory set in mb

    for item in options['NWCHEM']:
        options['NWCHEM'][item]['clobber'] = True

    return text, options


def muster_psi4options(opt, mol):
    """Translate psi4 keywords *opt* that have been explicitly set into
    their NWChem counterparts. Since explicitly set NWChem module keyword
    values will always be used preferentially to these inferred from
    psi4, the 'clobber' property is set to False.
    *mol* = qcdb.Molecule(molecule.create_psi4_string_from_molecule())   

    """
    text = ''
    options = defaultdict(lambda: defaultdict(dict))

    mult = mol.multiplicity()
    #print ('mult: ', mult)

    if 'SCF' in opt:
        if 'REFERENCE' in opt['SCF']:
            if opt['SCF']['REFERENCE']['value'] == 'UKS':
                options['NWCHEM']['NWCHEM_DFT']['value'] = 'ODFT'

            elif opt['SCF']['REFERENCE']['value'] == 'RKS':
                if mult != 1:
                    options['NWCHEM']['NWCHEM_DFT']['value'] = 'RODFT'
                else:
                    pass
            else:
                options['NWCHEM']['NWCHEM_SCF']['value'] = \
                        {'RHF': 'RHF',
                         'UHF': 'UHF',
                         'ROHF': 'ROHF'}[opt['SCF']['REFERENCE']['value']]

        if 'SCF_TYPE' in opt['SCF']:
            if opt['SCF']['SCF_TYPE']['value'] == 'DIRECT':
                options['NWCHEM']['NWCHEM_DFT_DIRECT']['value'] = 'TRUE'
                options['NWCHEM']['NWCHEM_SCF_DIRECT']['value'] = 'TRUE'

        if 'D_CONVERGENCE' in opt['SCF']:
            options['NWCHEM']['NWCHEM_DFT_CONVERGENCE']['value'] = \
                ['density', opt['SCF']['D_CONVERGENCE']['value']]
            #This is not correct so has to be fixed
            options['NWCHEM']['NWCHEM_SCF_THRESH']['value'] = \
                opt['SCF']['D_CONVERGENCE']['value']

        if 'E_CONVERGENCE' in opt['SCF']:
            options['NWCHEM']['NWCHEM_DFT_CONVERGENCE']['value'] += \
                ['energy', opt['SCF']['E_CONVERGENCE']['value']]
            options['NWCHEM']['NWCHEM_SCF_THRESH']['value'] = \
                opt['SCF']['E_CONVERGENCE']['value']

        if 'MAXITER' in opt['SCF']:
            options['NWCHEM']['NWCHEM_DFT_ITERATIONS']['value'] = \
               opt['SCF']['MAXITER']['value']
            options['NWCHEM']['NWCHEM_SCF_MAXITER']['value'] = \
               opt['SCF']['MAXITER']['value']

        if 'DFT_FUNCTIONAL' in opt['SCF']:
            dft_dispersion_list = \
                ['B3LYP-D3','BLYP-D3','BP86-D3','M05-D3','M05-2X-D3', \
                'PBE0-D3']
            if opt['SCF']['DFT_FUNCTIONAL']['value'] in dft_dispersion_list:
                options['NWCHEM']['NWCHEM_DFT_XC']['value'], options['NWCHEM']['NWCHEM_DFT_DISP']['value'] = \
            muster_dft_w_dispersion(opt)

            else:
                options['NWCHEM']['NWCHEM_DFT_XC']['value'] = muster_dft_functionals(opt)

        if 'DFT_RADIAL_POINTS' and 'DFT_SPHERICAL_POINTS' in opt['SCF']:
            rad_val = opt['SCF']['DFT_RADIAL_POINTS']['value']
            sph_val = opt['SCF']['DFT_SPHERICAL_POINTS']['value']

            sph_val_list = \
                [38,50,74,86,110,146,170,194,230,266,302,350,434,590,770,974,1202,1454,1730,2030,2354, \
                2702,3074,3470,3890,4334,4802,5294,5810]

            iangquad = sph_val_list.index(sph_val) + 1

            grid_val = []
            grid_val.append(rad_val)
            grid_val.append(iangquad)

            #print ('grid_val? ', grid_val)

            options['NWCHEM']['NWCHEM_DFT_GRID']['value'] = grid_val

    for item in options['NWCHEM']:
        options['NWCHEM'][item]['clobber'] = False
    return text, options


def muster_dft_functionals(opt):
    """Translate psi4 dft functionals w/o dispersion into NWChem counterparts. 
    
    """

    text = ''
    options = defaultdict(lambda: defaultdict(dict))
    dft_functional = ''

    dft_same_name_list = ['B3LYP','PBE0','M05','M05-2X',\
      'FT97','HCTH','HCTH120', 'DLDF', \
      'HCTHP14','HCTH407P','B1B95','M06','M06-2X','M06-HF',\
      'M08-HX','M08-SO','M11','M11-L','MPW1B95','MPW1K','MPWB1K',\
      'PW6B95','PWB6K', 'BOP']

    val = opt['SCF']['DFT_FUNCTIONAL']['value']
    if val in dft_same_name_list:
        options['NWCHEM']['DFT_XC']['value'] = \
                    {'B3LYP': 'B3LYP',
                     'PBE0' : 'PBE0',
                       'M05': 'M05',
                    'M05-2X':'M05-2X',
                     'FT97' : 'FT97',
                     'HCTH' : 'HCTH',
                   'HCTH120': 'HCTH120',
                     'DLDF' : 'DLDF',
                   'HCTHP14': 'HCTHP14',
                  'HCTH407P': 'HCTH407P',
                     'B1B95': 'B1B95',
                       'M06': 'M06',
                    'M06-2X': 'M06-2X',
                    'M06-HF': 'M06-HF',
                    'M08-HX': 'M08-HX',
                    'M08-SO': 'M08-SO',
                       'M11': 'M11',
                     'M11-L': 'M11-L',
                   'MPW1B95': 'MPW1B95',
                    'MPW1K' : 'MPW1K',
                    'MPWB1K': 'MPWB1K',
                    'PW6B95': 'PW6B95',
                     'PWB6K': 'PWB6K',
                      'BOP' : 'BOP'
                    }[opt['SCF']['DFT_FUNCTIONAL']['value']]

    elif (val == 'HF'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
            ['HFexch', 1.0]

    elif (val == 'B2PLYP'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
            ['HFexch', 0.53, 'becke88', 0.47, 'lyp', 0.73, 'mp2', 0.27]

    elif (val == 'BLYP'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
            ['becke88', 'lyp']

    elif (val == 'BP86'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
            ['becke88', 'perdew86']

    elif (val == 'PBE'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
            ['xpbe96', 'cpbe96']

    elif (val == 'PW91'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
            ['xperdew91', 'perdew91']

    elif (val == 'B97-1'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['becke97-1']
    elif (val == 'B97-2'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['becke97-2']
    elif (val == 'B97-GGA1'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['becke97gga1']

    elif (val == 'TPSSH'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['xctpssh']

    elif (val == 'BHANDH'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['beckehandh']

    elif (val == 'CAM-B3LYP'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
        ['xcamb88','1.00','lyp','0.81','vwn_5','0.19','hfexch','1.00\n','cam','0.33',\
        'cam_alpha','0.19','cam_beta','0.46']

    elif (val == 'B1LYP'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
            ['HFexch', 0.25,'becke88', 0.75, 'lyp']

    elif (val == 'B1PW91'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
            ['HFexch', 0.25,'becke88', 0.75, 'perdew91']

    elif (val == 'B3LYP5'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
            ['vwn_5', 0.19, 'lyp', 0.81, 'HFexch', 0.20, 'slater',0.8,'becke88','nonlocal',0.72]

    elif (val == 'B86BPBE'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
            ['becke86b', 'cpbe96']

    elif (val == 'B97') or (val == 'B97-0'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['becke97', 'HFexch', 0.1943]
    elif (val == 'B97-1P'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['becke97gga1', 'HFexch', 0.1500]
    elif (val == 'B97-3'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['becke97-3', 'HFexch', 0.2693]

    elif (val == 'BHANDHLYP') or (val == 'BHHLYP'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['becke88', 0.500, 'HFexch', 0.500, 'lyp']

    elif (val == 'LRC-WPBE'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['xwpbe',1.00,'cpbe96',1.0,'hfexch','1.00\n',\
          'cam',0.3,'cam_alpha',0.00,'cam_beta',1.00]
    elif (val == 'LRC-WPBEH'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['xwpbe',0.80,'cpbe96',1.0,'hfexch','1.00\n',\
          'cam',0.2,'cam_alpha',0.20,'cam_beta',0.80]

    elif (val == 'MPW1PW'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['mpw91', 0.75, 'HFexch', 0.25, 'perdew91']
    elif (val == 'MPWLYP1M'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['mpw91', 0.95, 'HFexch', 0.05, 'lyp']
    elif (val == 'MPWLYP1W'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['mpw91', 'vwn_5', 0.12, 'lyp', 0.88]

    elif (val == 'PBE0-13'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['xpbe96', 0.6667, 'HFexch', 0.3333, 'cpbe96']
    elif (val == 'PBEH'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['xpbe96', 0.75, 'HFexch', 0.25, 'cpbe96']
    elif (val == 'PBELYP1W'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['xpbe96', 'vwn_5', 0.26, 'lyp', 0.74]

    elif (val == 'PW86PBE'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['xperdew86', 'cpbe96']

    elif (val == 'TPSSLYP1W'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['xtpss03', 'vwn_5', 0.26, 'lyp', 0.74]

    elif (val == 'XLYP'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
        ['slater','-0.0690','becke88',0.722,'xperdew91',0.347,'lyp']

    elif (val == 'WPBE'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
            ['xwpbe', 1.0, 'cpbe96', 1.0, 'HFexch', 1.0,\
             '\ncam', 0.40, 'cam_alpha', 0.0, 'cam_beta', 1.0]

    elif (val == 'SB98-1A'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
            ['becke98','HFexch',0.229]

    elif (val == 'PBEH3C'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['xpbe96', 0.58, 'HFexch', 0.42, 'cpbe96']
    elif (val == 'PBE1W'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['xpbe96', 'vwn_5', 0.26, 'cpbe96', 0.74]
    elif (val == 'PBE0-2'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['xpbe96', 0.2063, 'HFexch', 0.7937, 'cpbe96', 0.5, 'mp2', 0.5]
    elif (val == 'O3LYP'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = [
            'slater', 0.0706, 'optx', 0.8133, 'HFexch', 0.1161, 'vwn_5', 0.19, 'lyp', 0.81
        ]
    elif (val == 'HSE03'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['xpbe96',1.0,'xcampbe96','-0.25','cpbe96','1.0\n',\
          'cam',0.1061,'cam_alpha',0.0,'cam_beta',0.25]
    elif (val == 'HSE06'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['xpbe96',1.0,'xcampbe96','-0.25','cpbe96',1.0,'srhfexch','0.25\n',\
          'cam',0.11,'cam_alpha',0.0,'cam_beta',1.0]
    elif (val == 'WPBE0'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
            ['xwpbe', 1.0, 'cpbe96', 1.0, 'HFexch', 1.0, \
            '\ncam', 0.30, 'cam_alpha', 0.25, 'cam_beta', 0.75]

    else:
        val = str(val)
        raise ValidationError(
            """Requested functional %s is not available. Check NWChem manual and try nwchem_dft_xc keyword instead. """
            % (val))

    dft_functional = options['NWCHEM']['NWCHEM_DFT_XC']['value']

    return dft_functional


def muster_dft_w_dispersion(opt):
    """Translate psi4 dft functionals w/dispersion into NWChem counterparts. 
    
    """

    text = ''
    options = defaultdict(lambda: defaultdict(dict))
    dft_functional = ''
    dft_dispersion = ''

    val = opt['SCF']['DFT_FUNCTIONAL']['value']

    if (val == 'B3LYP-D3'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['b3lyp']
        options['NWCHEM']['NWCHEM_DFT_DISP']['value'] = ['vdw', 3]

    elif (val == 'BLYP-D3'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
            ['becke88', 'lyp']
        options['NWCHEM']['NWCHEM_DFT_DISP']['value'] = ['vdw', 3]

    elif (val == 'BP86-D3'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
            ['becke88', 'perdew86']
        options['NWCHEM']['NWCHEM_DFT_DISP']['value'] = ['vdw', 3]

    elif (val == 'M05-D3'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['m05']
        options['NWCHEM']['NWCHEM_DFT_DISP']['value'] = ['vdw', 3]

    elif (val == 'M05-2X-D3'):
        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['m05-2x']
        options['NWCHEM']['NWCHEM_DFT_DISP']['value'] = ['vdw', 3]

    elif (val == 'PBE0-D3'):
        options['NWCHEM']['NWCHEM__DFT_XC']['value'] = ['pbe0']
        options['NWCHEM']['NWCHEM_DFT_DISP']['value'] = ['vdw', 3]

    else:
        val = str(val)
        raise ValidationError(
            """Requested functional %s is not available. Check NWChem manual and try nwchem_dft_xc keyword instead. """
            % (val))

    dft_functional = options['NWCHEM']['NWCHEM_DFT_XC']['value']
    dft_dispersion = options['NWCHEM']['NWCHEM_DFT_DISP']['value']

    return dft_functional, dft_dispersion


def format_modelchem_for_nwchem(name, dertype, ropts, sysinfo, verbose=1):
    accession = sys._getframe().f_code.co_name + '_' + str(uuid.uuid4())
    kwgs = {'accession': accession, 'verbose': verbose}

    lowername = name.lower()

    #runtyp = {'energy': 'energy',
    #          'gradient': 'gradient',
    #          'hessian': 'hessian',
    #          'properties': 'prop',
    #         }[driver]

    runtyp = {0: 'energy',
              1: 'gradient',
              2: 'hessian',
              #'properties': 'prop',
             }[dertype]

    #theory = {0: 'scf',
    #          1: 'dft',
    #          }[dertype] #temp, may need to change only two options

    if lowername == 'nwc-nwchem':
        pass

    elif lowername in ['nwc-scf', 'nwc-hf']:
        #ropts.require('NWCHEM', 'task__scf', runtyp, **kwgs)
        mdccmd = f'task scf {runtyp}\n\n'
    
    #property
    #elif lowername == 'nwc-property':
    #    mdccmd = f'task {theory} property \n\n'
        #ropts.suggest('NWCHEM', task__scf', theory, **kwgs)

    #MPn options
    elif lowername == 'nwc-mp2':
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            mdccmd = f'task tce {runtyp}\n\n'
            ropts.require('NWCHEM', 'tce__mp2', True, **kwgs)
        else:
            mdccmd = f'task mp2 {runtyp} \n\n'
    elif lowername == 'nwc-mp3':
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            mdccmd = f'task tce {runtyp}\n\n'
            ropts.require('NWCHEM', 'tce__mp3', True, **kwgs)
    elif lowername == 'nwc-mp4':    
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            mdccmd = f'task tce {runtyp}\n\n'
            ropts.require('NWCHEM', 'tce__mp4', True, **kwgs)
    #CC options
    elif lowername == 'nwc-ccd':
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            mdccmd = f'task tce {runtyp} \n\n'
            ropts.require('NWCHEM', 'tce__ccd', True, **kwgs)
    elif lowername == 'nwc-ccsd':
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            mdccmd = f'task tce {runtyp}\n\n'
            ropts.require('NWCHEM', 'tce__ccsd', True, **kwgs)
        else:
             mdccmd = f'task ccsd {runtyp}\n\n'
    elif lowername == 'nwc-ccsdt':
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            mdccmd = f'task tce {runtyp} \n\n'
            ropts.require('NWCHEM', 'tce__ccsdt', True, **kwgs)
        else:
            mdccmd = f'task ccsdt {runtyp}\n\n'
    elif lowername == 'nwc-ccsd(t)':
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            mdccmd = f'task tce {runtyp} \n\n'
            ropts.require('NWCHEM', 'tce__ccsd(t)', True, **kwgs)
        else:
            mdccmd = f'task ccsd(t) {runtyp}\n\n'
    elif lowername == 'nwc-ccsdtq':
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            mdccmd = f'task tce {runtyp} \n\n'
            ropts.require('NWCHEM', 'tce__ccsdtq', True, **kwgs)
    elif lowername =='nwc-eom-ccsd':
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            mdccmd = f'task tce {runtyp} \n\n'
            ropts.require('NWCHEM', 'tce__ccsd', True, **kwgs)
            ropts.suggest('NWCHEM', 'tce__nroots', 4, **kwgs)
    elif lowername == 'nwc-eom-ccsdt':
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            ropts.require('NWCHEM', 'tce__ccsdt', True, ** kwgs)
            ropts.suggest('NWCHEM', 'tce__nroots', 4, **kwgs)

   # elif lowername == 'nwc-ccsd_act':
    
    #TCE only opts
    elif lowername == 'nwc-cisd':
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            mdccmd = f'task tce {runtyp} \n\n'
            ropts.require('NWCHEM', 'tce__cisd', True, **kwgs)
    elif lowername == 'nwc-cisdt':
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            mdccmd = f'task tce {runtyp} \n\n'
            ropts.require('NWCHEM', 'tce__cisdt', True, **kwgs)
    elif lowername == 'nwc-cisdtq':
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            mdccmd = f'task tce {runtyp} \n\n'
            ropts.require('NWCHEM', 'tce__cisdtq', True, **kwgs)
    elif lowername == 'nwc-qcisd':
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            mdccmd = f'task tce {runtyp} \n\n'
            ropts.require('NWCHEM', 'tce__qcisd', True, **kwgs)

    elif lowername == 'nwc-lccd':
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            mdccmd = f'task tce {runtyp} \n\n'
            ropts.require('NWCHEM', 'tce__lccd', True, **kwgs)
    elif lowername == 'nwc-lccsd':
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            mdccmd = f'task tce {runtyp} \n\n'
            ropts.require('NWCHEM', 'tce__lccsd', True, **kwgs)
    elif lowername == 'nwc-cc2':
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            mdccmd = f'task tce {runtyp} \n\n'
            ropts.require('NWCHEM', 'tce__cc2', True, **kwgs)
    elif lowername == 'nwc-lr-ccsd':
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            mdccmd = f'task tce {runtyp} \n\n'
            ropts.require('NWCHEM', 'tce__lr-ccsd', True, **kwgs)
    elif lowername == 'nwc-ccsdta':
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            mdccmd = f'task tce {runtyp} \n\n'
            ropts.require('NWCHEM', 'tce__ccsdta', True, **kwgs)
    elif lowername == 'nwc-eaccsd':
        if ropts.scroll['QCDB']['QC_MODULE'].value == 'tce':
            mdccmd = f'task tce {runtyp} \n\n'
            ropts.require('NWCHEM', 'tce__eaccsd', True, **kwgs)
            ropts.suggest('NWCHEM', 'tce__nroots', 4, **kwgs)

    #DFT xc functionals
    elif lowername == 'nwc-pbe0':
        ropts.require('NWCHEM', 'xc', 'pbe0', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-b3lyp':
        ropts.require('NWCHEM', 'xc', 'b3lyp', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-b1b95':
        ropts.require('NWCHEM', 'xc', 'b1b95', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-b97-1':
        ropts.require('NWCHEM', 'xc', 'becke97-1', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-b97-2':
        ropts.require('NWCHEM', 'xc', 'becke97-2', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-b97gga1':
        ropts.require('NWCHEM', 'xc', 'becke97gga1', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-bhandh':
        ropts.require('NWCHEM', 'xc', 'beckehandh', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-bop':
        ropts.require('NWCHEM', 'xc', 'bop', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-d1df':
        ropts.require('NWCHEM', 'xc', 'd1df', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-ft97':
        ropts.require('NWCHEM', 'xc', 'ft97', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-hcth':
        ropts.require('NWCHEM', 'xc', 'hcth', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-hcth120':
        ropts.require('NWCHEM', 'xc', 'htch120', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-hcth407p':
        ropts.require('NWCHEM', 'xc', 'hcth407p', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-hcthp14':
        ropts.require('NWCHEM', 'xc', 'hcthp14', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-m05':
        ropts.require('NWCHEM', 'xc', 'm05', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-m05-2x':
        ropts.require('NWCHEM', 'xc', 'm05-2x', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-m06':
        ropts.require('NWCHEM', 'xc', 'm06', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-m06-2x':
        ropts.require('NWCHEM', 'xc', 'm06-2x', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-m06-hf':
        ropts.require('NWCHEM', 'xc', 'm06-hf', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-m06-hf':
        ropts.require('NWCHEM', 'xc', 'm06-hf', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-m08-hx':
        ropts.require('NWCHEM', 'xc', 'm08-hx', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-m08-so':
        ropts.require('NWCHEM', 'xc', 'm08-so', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-m11':
        ropts.require('NWCHEM', 'xc', 'm11', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-m11-l':
        ropts.require('NWCHEM', 'xc', 'm11-l', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-mpw1b95':
        ropts.require('NWCHEM', 'xc', 'mpw1b95', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-mpw1k':
        ropts.require('NWCHEM', 'xc', 'mpw1k', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-mpwb1k':
        ropts.require('NWCHEM', 'xc', 'mpwb1k', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'pw6b95':
        ropts.require('NWCHEM', 'xc', 'pw6b95', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-pwb6k':
        ropts.require('NWCHEM', 'xc', 'pwb6k', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-tpssh':
        ropts.require('NWCHEM', 'xc', 'xctpssh', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    #DFT functionals potential issues - multiple options or req conditions
    elif lowername == 'nwc-b86bpbe':
        ropts.require('NWCHEM', 'xc', 'becke86b', **kwgs)
        ropts.require('NWCHEM', 'xc', 'cpbe96', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
    elif lowername == 'nwc-blyp':
        ropts.require('NWCHEM', 'xc', 'becke88', **kwgs)
        ropts.require('NWCHEM', 'xc', 'lyp', **kwgs)
        mdccmd = f'task dft {runtyp} \n\n'
     
    elif lowername == 'nwc-tddft':
        mdccmd = f'task tddft {runtyp} \n\n'


    else:
        raise ValidationError(f"""Requested NWChem computational method {lowername} is not available.""")

    return mdccmd


def nwchem_list():
    """Return an array of NWChem methods with energies. Appended
    to procedures[qcdb.energy()].

    """
    val = []
    val.append('nwchem')
    val.append('nwc-scf')
    val.append('nwc-hf')
    val.append('nwc-mp2')
    val.append('nwc-mp3')
    val.append('nwc-mp4')
    val.append('nwc-dft')
    val.append('nwc-cc2')
    val.append('nwc-ccsd')
    val.append('nwc-ccsdt')
    val.append('nwc-ccsdtq')
    val.append('nwc-ccsd[t]')
    val.append('nwc-ccsd(t)')
    val.append('nwc-eaccsd')
    val.append('nwc-eom-ccsd')
    val.append('nwc-eom-ccsdt')
    val.append('nwc-eom-ccsdtq')
    val.append('nwc-cisd')
    val.append('nwc-cisdt')
    val.append('nwc-cisdtq')
    val.append('nwc-lccd')
    val.append('nwc-lccsd')
    val.append('nwc-ccd')
    val.append('nwc-eom-ccsd')
    val.append('nwc-eom-ccsdt')
    val.append('nwc-eom-ccsdtq')
    val.extend(dft_functionals_list())
    val.append('nwc-tddft')

    return val


def nwchem_gradient_list():
    """Return an array of NWChem methods with energies. Appended
    to procedures['gradient'].

    """
    val = []
#    val.append('nwchem')
    val.append('nwc-scf')
    val.append('nwc-hf')
    val.append('nwc-mp2')
    val.append('nwc-dft')
    val.append('nwc-ccsd')
    val.append('nwc-ccsdt')
    val.append('nwc-ccsdtq')
    val.append('nwc-ccsd(t)')
#    val.append('nwc-lccd')
#    val.append('nwc-eom-ccsd')
#    val.append('nwc-eom-ccsdt')
#    val.append('nwc-eom-ccsdtq')
#    val.extend(dft_functionals_list())
    val.append('nwc-tddft')

    return val


def nwchem_hessian_list():
    return []


def dft_functionals_list():
    """Return an array of available method of dft functionals through nwchem interface. 

    """
    val = []
    val.append('nwc-b3lyp')
    val.append('nwc-pbe0')
    val.append('nwc-m05')
    val.append('nwc-m05-2x')
    val.append('nwc-ft97')
    val.append('nwc-hcth')
    val.append('nwc-hcth120')
    val.append('nwc-dldf')
    val.append('nwc-b2plyp')
    val.append('nwc-blyp')
    val.append('nwc-bp86')
    val.append('nwc-pbe')
    val.append('nwc-b97-0')
    val.append('nwc-b97-1')
    val.append('nwc-b97-2')
    val.append('nwc-hcthp14')
    val.append('nwc-hcth407p')
    val.append('nwc-b1b95')
    val.append('nwc-m06')
    val.append('nwc-m06-2x')
    val.append('nwc-m06-hf')
    val.append('nwc-m08-hx')
    val.append('nwc-m08-so')
    val.append('nwc-m11')
    val.append('nwc-m11-l')
    val.append('nwc-mpw1b95')
    val.append('nwc-mpw1k')
    val.append('nwc-mpwb1k')
    val.append('nwc-pw6b95')
    val.append('nwc-pwb6k')
    val.append('nwc-b97-gga1')
    val.append('nwc-b97')
    val.append('nwc-tpssh')
    val.append('nwc-bhandh')
    val.append('nwc-cam-b3lyp')
    val.append('nwc-bop')
    val.append('nwc-hf')
    val.append('nwc-pw91')
    val.append('nwc-b1lyp')
    val.append('nwc-b1pw91')
    val.append('nwc-b3lyp5')
    val.append('nwc-b86bpbe')
    val.append('nwc-b97-1p')
    val.append('nwc-b97-3')
    val.append('nwc-bhandhlyp')
    val.append('nwc-bhhlyp')
    val.append('nwc-lrc-wpbe')
    val.append('nwc-lrc-wpbeh')
    val.append('nwc-mpw1pw')
    val.append('nwc-mpwlyp1m')
    val.append('nwc-mpwlyp1w')
    val.append('nwc-pbe0-13')
    val.append('nwc-pbeh')
    val.append('nwc-pbelyp1w')
    val.append('nwc-pw86pbe')
    val.append('nwc-tpsslyp1w')
    val.append('nwc-xlyp')
    val.append('nwc-wpbe')

    return val
