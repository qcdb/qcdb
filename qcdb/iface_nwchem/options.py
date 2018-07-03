import collections

from ..moptions.read_options2 import RottenOption
import copy
import uuid

from ..exceptions import *
from . import parsers
#Need to check if nwchem or nwc is option?
### Example of options 
#options.add('nwchem', RottenOption(
    #keyword = name within NWCHEM documentation
    #default = psi4 options from common-driver-psi4
    #validator = either a parser or lambda transforming it to uppercase or something similar- check the parsers.py file
    #glossary = description

#NWCHEM only working to translate with Psi4 at the moment via QCDB
def load_nwchem_defaults(options):
    #Haven't added NCHEM_OMP_NUM_THREADS
    #options.add('nwchem', RottenOption(
        #keyword='omp_num_threads
        #default= 1
        #validator=
        #glossary='Sets OMP_NUM_THREADS environment variable before calling NWCHEM'
    options.add('nwchem', RottenOption(
        keyword='translate_psi4',
        default= True,
        validator= parsers.boolean,
        glossary='Will translate Psi4 options to NWCHEM counterparts'))

    options.add('nwchem', RottenOption(
        keyword='memory',
        default='total 400 mb',
        validator=parsers.parse_memory,
        glossary='Total memory allocation, which can be specified furth to stack, heap, and global.'))

    options.add('nwchem', RottenOption(
        keyword='charge',
        default= 0,
        validator= parsers.parse_charge,
        glossary='charge from active molecule'))
#SCF block
    options.add('nwchem', RottenOption(
        keyword='scf',
        default='RHF',
        validator=parsers.enum,
        glossary='Reference wave function: RHF, UHF, ROHF'))

    options.add('nwchem', RottenOption(
        keyword='scf_nopen',
        default= 0,
        validator= lambda x: float(x),
        glossary='Specifies the number of open shells in wave function.'))

    options.add('nwchem', RottenOption(
        keyword='scf_thresh',
        default=1.e-4,
        validator=parsers.parse_convergence,
        glossary='SCF Convergence threshold'))

    options.add('nwchem', RottenOption(
        keyword='scf_maxiter',
        default= 20,
        validator= parsers.postive_integer,
        glossary='Max SCF iteration'))
        
    options.add('nwchem', RottenOption(
        keyword='scf_diis',
        default= False,
        validator= parsers.boolean,
        glossary='DIIS Convergence on or off, default is off'))

    options.add('nwchem', RottenOption(
        keyword='scf_direct',
        default= False,
        validator= parsers.boolean,
        glossary='SCF Direct calculation on/off'))

#    options.add('nwchem', RottenOption(
 #       keyword='scf_semidirect',
  #      default= '',
   #     validator= , #validator type?
    #    glossary='''Semidirection calculation which defines amount of disk space and cache memory size. 
     #   Filesize, memsize, and filename can be specified in array type. 
      #  Check NWCHEM manual for details.'''))

    options.add('nwchem', RottenOption(
        keyword='scf_sym',
        default= True,
        validator= parsers.boolean,
        glossary='Symmetry specification in SCF for Fock matrix construction on/off'))

    options.add('nwchem', RottenOption(
        keyword='scf_adapt',
        default= True,
        validator= parsers.boolean,
        glossary='Force symmetry adaption of the molecular orbitals in SCF'))

    options.add('nwchem', RottenOption(
        keyword='scf_tol2e',
        default= 1.e-7,
        validator= parsers.parse_convergence,
        glossary='Integral screening threshold for the evaluation of the energy and related Fock-like matrices'))
### Specifies the source and destination of the molecular orbital vectors.
#        Input type has to be array type. 
#         Refer NWChem manual for available options.
#         e.g) 
    #          set nwchem_scf_vectors [input, try1.movecs, output, try2.movecs, ...] -*/
     #     options.add("NWCHEM_SCF_VECTORS", new ArrayType())

#    options.add('nwchem', RottenOption(
 #       keyword='scf_vectors'
  #      default=''
   #     validator=''
    #    glossary='Specify source and destination of molecular orbital vectors.'))

    options.add('nwchem', RottenOption(
        keyword='scf_profile',
        default= False,
        validator= parsers.boolean,
        glossary='SCF performance file true/false. Default is false.'))

    options.add('nwchem', RottenOption(
        keyword='scf_nr',
        default= 0.1,
        validator= lambda x: float(x),
        glossary='Control Netwon-Raphson value.'))

#TODO #Array block in SCF
#    options.add('nwchem', RottenOption(
 #       keyword='scf_level',
  #      default='',
   #     validator='',
    #    glossary=''

   # options.add('nwchem', RottenOption(
    #    keyword='scf_print',
     #   default='',
      #  validator='',
       # glossary=''

#    options.add('nwchem', RottenOption(
 #       keyword='scf_noprint',
  #      default='',
   #     validator='',
    #    glossary=''

#MP2 block
    options.add('nwchem', RottenOption(
        keyword='mp2_tight',
        default= False,
        validator= parsers.boolean,
        glossary='''Increase precision of MP2 energy and gradients. Will also change SCF and CPHF precision.
        Tightens thresholds for AO and MO integrals within MP2 code. Default is off'''))

#   options.add('nwchem', RottenOption(
#        keyword='mp2_freeze',
 #       default='',
  #      validator='',
   #     glossary=''#TODO #Another array

    options.add('nwchem', RottenOption(
        keyword='mp2_scs',
        default= True,
        validator=parsers.boolean,
        glossary='')) #TODO  #Description needed

    options.add('nwchem', RottenOption(        
        keyword='mp2_fss',
        default= 1.2,
        validator=lambda x: float(x),
        glossary='Scaling factor for same spin'))

    options.add('nwchem', RottenOption(
        keyword='mp2_fos',
        default= 0.3,
        validator=lambda x: float(x),
        glossary='Scaling factor for opposite spin'))

#DFT block
    options.add('nwchem', RottenOption(
        keyword='dft',
        default='',
        validator= parsers.enum,
        glossary='Defining DFT wavefunction: RDFT, RODFT, UDFT, ODFT (Open shell, singlet).'))

#TODO #Array block for DFT - dft_xc, dft_grid, dft_convergence

    options.add('nwchem', RottenOption(
        keyword='dft_iterations',
        default= 30,
        validator= parsers.postive_integer,
        glossary='Specify DFT iterations'))

    options.add('nwchem', RottenOption(
        keyword='dft_mult',
        default= 1,
        validator= lambda x: float(x),
        glossary='DFT Mulitiplicity'))

#    options.add('nwchem', RottenOption(
 #       keyword='dft_vectors',
  #      default='',
   #     validator='',
    #    glossary=''
    options.add('nwchem', RottenOption( 
        keyword='dft_max_ovl', 
        default= False, 
        validator= parsers.boolean, 
        glossary='Lock the ordering of orbitals on/off. Default is off.'))
    
    options.add('nwchem', RottenOption( 
        keyword='dft_smear', 
        default= 0.001, 
        validator= lambda x: float (x), 
        glossary='Smear keyword allows fractional occuption of the MOs.')) 
    
    options.add('nwchem', RottenOption(
        keyword='dft_mulliken',
        default= False,
        validator= parsers.boolean,
        glossary='Mulliken analysis of charge distribution: on/off. Default is off.'))

    options.add('nwchem', RottenOption(
        keyword='dft_direct',
        default= False,
        validator= parsers.boolean,
        glossary='Direct calculation of DFT: on/off. Default is off.'))
#DFT_Semidirect is an array #TODO 
#Add here
    options.add('nwchem', RottenOption(
        keyword='dft_cgmin',
        default= False,
        validator= parsers.boolean,
        glossary='DFT quadratic convergence algorithm: true/false. Default is false.'))

    options.add('nwchem', RottenOption(
        keyword='dft_fukui',
        default= False,
        validator= parsers.boolean,
        glossary='Fukui indices analysis: on/off. Default is off(false).'))
#DFT Arrays- dft_disp, dft_print, dft_noprint #TODO

#TCE block
    options.add('nwchem', RottenOption(
        keyword='tce_dft',
        default=False,
        validator= parsers.boolean,
        glossary='Use DFT as TCE reference wave function. If not specified, default is SCF(HF).'))

    options.add('nwchem', RottenOption(
        keyword='tce',
        default='',
        validator= parsers.enum,
        glossary='''Specify TCE correlation models. Options include:
        LCCD, CCD, LCCSD, CCSD, CCSD_ACT, LR-CCSD, EACCSD, IPCCSD, CC2, CCSDT, CCSDTA, CCSDTQ, CCSD(T), CCSD[T]
        CR-CCSD[T], CR-CCSD(T), CCSD(2)_T, CCSD(2)_TQ, CCSDT(2)_Q, LR-CCSD(T), LR-CCSD(TQ)-1, CREOMSD(T),
        CREOM(T)AC, QCISD, CISD, CISDT, CISDTQ, MBPT2, MBPT3, MBPT4.
        MBP2= MP2, MBPT3= MP3, MBPT4= MP4.'''))


    options.add('nwchem', RottenOption(
        keyword='tce_thresh',
        default=1.e-4,
        validator= parsers.parse_convergence,
        glossary='TCE convergence threshold'))

    options.add('nwchem', RottenOption(
        keyword='tce_maxiter',
        default= 100,
        validator= parsers.positive_integer,
        glosssary='TCE max iterations'))

    options.add('nwchem', RottenOption(
        keyword='tce_io',
        default='',
        validator= lambda x: x.upper(),
        glossary='''Parallel I/O scheme specification. Available:
        FORTRAN, EAF, GA, SF, REPLICATED, DRA, GAEAF'''))

    options.add('nwchem', RottenOption(
        keyword='tce_diis',
        default=5,
        validator= parsers.postive_integer,
        glossary='''Number of iterations for a DIIS extrapolation to be performed. Will accelerate excitation
        amplitude convergence. Default is 5.'''))

    #options.add('nwchem', RottenOption(
     #   keyword='tce_freeze',
      #Array TODO

    options.add('nwchem', RottenOption(
        keyword='tce_nroots',
        default= 1,
        validator= parsers.postive_integer,
        glossary='Number of excited states. Default is 1.'))

    options.add('nwchem', RottenOption(
        keyword='tce_target',
        default= 1,
        validator= parsers.postive_integer,
        glossary='TCE target root. Default is 1.'))

    options.add('nwchem', RottenOption(
        keyword='tce_targetsym',
        default='',
        validator= lambda x: x.upper(),
        glossary='TCE target symmetry. Default is None.'))

    options.add('nwchem', RottenOption(
        keyword='tce_2eorb',
        default= False,
        validator= parsers.boolean,
        glossary='''Economical option of storing two-electron integrals used in coupled cluster calculations,
        taking the difference of the RHF and ROHF values: on/off. Default is off.'''))

#TASK block
    options.add('nwchem', RottenOption(
        keyword='task_scf',
        default='energy',
        validator= lambda x: x.upper(),
        glossary='Specify scf task between: energy, gradient, and hessian. Default is energy.'))

    options.add('nwchem', RottenOption(
        keyword='task_mp2',
        default='energy',
        validator= lambda x: x.upper(),
        glossary='Specify MP2 [semi-direct] task between: energy, gradient, and hessian. Default is energy.'))

    options.add('nwchem', RottenOption(
        keyword='task_direct_mp2',
        default='energy',
        validator=lambda x: x.upper(),
        glossary='Specify MP2 [direct] task between: energy, gradient, and hessian. Default is energy.'))

    options.add('nwchem', RottenOption(
        keyword='task_rimp2',
        default='energy',
        validator= lambda x: x.upper(),
        glossary='Specify RIMP2 task between: energy, gradient, and hessian. Default is energy.'))

    options.add('nwchem', RottenOption(
        keyword='task_dft',
        default='energy',
        validator= lambda x: x.upper(),
        glossary='Specify DFT task between: energy, gradient, and hessian. Default is energy.'))

    options.add('nwchem', RottenOption(
        keyword='task_tce',
        default='energy',
        validator= lambda x: x.upper(),
        glossary='Specify TCE task between: energy, gradient, and hessian. Default is energy.'))

    options.add('nwchem', RottenOption(
        keyword='task_ccsd(t)',
        default='energy',
        validator=lambda x: x.upper(),
        glossary='Specify CCSD(T) task between energy, gradient, and hessian. Default is energy.'))
