import collections

from ..moptions.read_options2 import RottenOption
import copy
import uuid

from ..exceptions import *
from . import parsers
#Need to check if nwchem or nwc is option?
### Example of options 
    #options.add('nwchem', RottenOption( 
        #keyword = name within NWCHEM documentation, 
        #default = psi4 options from common-driver-psi4, 
        #validator = either a parser or lambda transforming it to uppercase or something similar- check the parsers.py file,
        #glossary = description))

#NWCHEM only working to translate with Psi4 at the moment via QCDB
def load_nwchem_defaults(options):
    #Haven't added NCHEM_OMP_NUM_THREADS
    #options.add('nwchem', RottenOption(
        #keyword='omp_num_threads
        #default= 1
        #validator=
        #glossary='Sets OMP_NUM_THREADS environment variable before calling NWCHEM'
   #option for NWChem may be 'nwchem' or 'nwc'
    options.add('nwchem', RottenOption(
        keyword='translate_psi4',
        default= True,
        validator= parsers.boolean,
        glossary='Will translate Psi4 options to NWCHEM counterparts'))
#Memory specifications
    options.add('nwchem', RottenOption(
        keyword='total_memory',
        default='total 400 mb',
        validator=parsers.parse_memory,
        glossary='Total memory allocation, which can be specified furth to stack, heap, and global.'))
    options.add('nwchem', RottenOption(
        keyword='stack_memory',
        default='',
        validator=parsers.parse_memory,
        glossary= 'Stack memory allocation that can be specified beyond the total memory.'))
    options.add('nwchem', RottenOption(
        keyword='heap_memory',
        default='',
        validator=parsers.parse_memory,
        glossary= 'Heap memory allocation that can be specified beyond the total memory allocation.'))
    options.add('nwchem', RottenOption(
        keyword='global_memory',
        default='',
        validator=parsers.parse_memory,
        glossary='Global memory allocation.'

#SCF block
    options.add('nwchem', RottenOption(
        keyword='scf',
        default='RHF',
        validator=parsers.enum("RHF UHF ROHF"),
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
    
    #SCF semidirect options- Array okay? Check with LB
    options.add('nwchem', RottenOption(
        keyword='scfsd__filesize',
        default= 0,
        validator= parsers.postive_integer,
        glossary= '''SCF Semidirect options: File size allows the user to specify the amount of disk space used per process for storing the integrals in 64-bit words. 
        Default of semidirect leads to directive DIRECT.'''))
    options.add('nwchem', RottenOption(
        keyword='scfsd__memsize',
        default= 0,
        validator= parsers.positive_integer,
        glossary= '''SCF Semidirect options: Memory size allows user to specify the number of 64-bit words to be used per process for caching integrals in memory.
        If the amount of storage space is not available, the code cuts the value in half and checks again, and will continue to do so until request is satisfied.'''))
    options.add('nwchem', RottenOption(
        keyword= 'scfsd__filename',
        default='',
        validator= lambda x: x.lower(), #need to specify file extensions?
        glossary= '''SCF Semidirect options: 
        File name default for integral files are placed into scratch directory (see NWChem manual for details). 
        User-specified string in File name keyword will override default.'''))

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

    options.add('nwchem', RottenOption(
        keyword= 'scf_print__level',
        default= 'default',
        validator= parsers.enum('default debug low high never'),
        glossary='Options of output in SCF module controlled by PRINT directive'))

    options.add('nwchem', RottenOption(
        keyword= 'scf_print__debugname',
        default= '',
        validator= parsers.enum('"atomic guess density" "atomic scf" "initial vectors" "intermediate vectors" "final vectors" "intermediate evals" "screening statistics" "symmetry" "geombas"'), #white space?
        glossary='debug options'))

    options.add('nwchem', RottenOption(
        keyword= 'scf_print__defaultname',
        default= ''
        validator= parsers.enum('"mo guess" "final vectors analysis" "final evals" "vectors i/o" parameters convergence timing '), #Not sure if this will be parsed out correctly
        glossary= 'default options'))

    options.add('nwchem', RottenOption(
        keyword= 'scf_print__lowname',
        default= '',
        validator= parsers.enum('information'),
        glossary= 'will print information if choose low print level'))

    options.add('nwchem', RottenOption(
        keyword= 'scf_print__highname',
        default= '',
        validator= parsers.enum('schwarz geometry basis'),
        glossary= 'options to print in the high level'))

    options.add('nwchem', RottenOption(
        keyword= 'scf_print__nevername',
        default= '',
        validator= parsers.enum(' "initial vectors analysis" "mulliken ao" ')
        glossary= 'options in the never level'))

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
        glossary='Spin Component Scaled MP2 calculations. Default is on.')) #TODO  #Description needed

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
        validator= parsers.enum("RDFT RODFT UDFT ODFT"),
        glossary='Defining DFT wavefunction: RDFT, RODFT, UDFT, ODFT (Open shell, singlet).'))

#TODO #Array block for DFT - dft_xc, dft_grid, dft_convergence
#DFT Convergence block
    options.add('nwchem', RottenOption(
        keyword='dft_convergence',
        default= '',
        validator= parsers.enum('energy density, gradient hl_tol dampon dampoff ncydp ncyds ncysh damp diison diisoff diis levlon levloff lshift rabuck'),
        glossary='DFT Convergence options to specify.'))
    options.add('nwchem', RottenOption(
        keyword='convergence__energy',
        default= 1.e-6,
        validator= parsers.parse_convergence,
        glossary= 'total energy convergence DFT'))
    options.add('nwchem', RottenOption(
        keyword='convergence__density',
        default= 1.e-5,
        validator= parsers.parse_convergence,
        glossary= 'DFT convergence of total density.'))
    options.add('nwchem', RottenOption(
        keyword= 'convergence__gradient',
        default= 5.e-4,
        validator= parsers.parse_convergence,
        glossary=''))
    options.add('nwchem', RottenOption(
        keyword='convergence__hltol',
        default= 0.1
        validator= parsers.postive_integers,
        glossary=''))
    options.add('nwchem', RottenOption(
        keyword='convergence__dampon' 
        default=0.0
        validator= lambda x: float(x),
        glossary='Turns on damping when reaching user-specified energy level.'))
    options.add('nwchem', RottenOption(
        keyword='convergence__dampoff',
        default= 0.0,
        validator= lambda x: float(x),
        glossary= 'Turns off damping when reaching user-specified energy level.'))
    options.add('nwchem', RottenOption(
        keyword='convergence__damp',
        default= 0,
        validator=parsers.positive_integers,
        glossary='Percent of previous iterations mixed with current iterations density.'))
    options.add('nwchem', RottenOption(
        keyword='convergence__ncydp',
        default= 2,
        validator= parsers.positive_integers,
        glossary= 'Specifies number of damping cycles. Default is 2.'))
    options.add('nwchem', RottenOption(
        keyword='convergence__diison'
        default= 0.0,
        validator= lambda x: float(x),
        glossary='')) 
    options.add('nwchem', RottenOption(
        keyword='convergence__diisoff',
        default=  0.0,
        validator= lambda x: float(x),
        glossary= ''))
    options.add('nwchem', RottenOption(
        keyword='convergence__diis',
        default= 10,
        validator= lambda x: float(x),
        glossary= 'Number of Fock matrices used in direct inversion of iterative subspace [DIIS] extrapolation'))
    options.add('nwchem', RottenOption(
        keyword='convergence__ncyds',
        default= 30,
        validator= parsers.positive_integers,
        glossary='Specifies number of DIIS [Direct inversion of iterative subspace] cycles needed. Default is 30.'))
    options.add('nwchem', RottenOption(
        keyword='convergence__levlon',
        default=0.0,
        validator=lambda x: float (x),
        glossary=''))
    options.add('nwchem', RottenOption(
        keyword='convergence__levloff',
        default= 0.0,
        validator= lambda x: float(x),
        glossary=''))
    options.add('nwchem', RottenOption(
        keyword='convergence__lshift',
        default= 0.5,
        validator= lambda x: float(x),
        glossary= ''))
    options.add('nwchem', RottenOption(
        keyword='convergence__ncysh',
        default= 0,
        validator= parsers.positive_integers,
        glossary= 'Specifies the number of level-shifting cycles are used in the input. Default is 0.'))
    options.add('nwchem', RottenOption(
        keyword='convergence__rabuck'
        default= 25,
        validator= parsers.positive_integers,
        glossary=''))

#DFT block [continued]
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
#DFT_Semidirect #Check ATL 
    options.add('nwchem', RottenOption(
        keyword='dft__semidirect',
        default= '',
        validator= parsers.enum('filesize memsize filename'),
        glossary= '''DFT Semidirect options: File size allows the user to specify the amount of disk space used per process for storing the integrals in 64-bit words. 
        Default of semidirect leads to directive DIRECT.'''))
   options.add('nwchem', RottenOption(
       keyword='semidirect__filesize',
       default= '', #default is disksize
       validator=parsers.positive_integers,
       glossary=''))
   options.add('nwchem', RottenOption(
       keyword='semidirect__memsize',
       default='',
       validator=parsers.parse_memory,
       glossary=''))
   options.add('nwchem', RottenOption(
       keyword='semidirect__filename',
       default= ''
       validator= lambda x
       glosary=''))

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
        validator= parsers.enum("LCCD CCD LCCSD CCSD CCSD_ACT LR-CCSD EACCSD IPCCSD CC2 CCSDT CCSDTA CCSDTQ CCSDTQ CCSD(T) CCSD[T] CR-CCSD[T] CR-CCSD(T) CCSD(2)_T CCSD(2)_TQ CCSDT(2)_Q LR-CCSD(T) LR-CCSD(TQ)-1 CREOMSD(T) CREOM(T)AC QCISD CISD CISDT CISDTQ MBPT2 MBPT3 MBPT4 MP2 MP3 MP4"),
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
