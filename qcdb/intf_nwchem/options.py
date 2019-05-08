import collections

from ..moptions.read_options2 import RottenOption
from ..moptions import parsers

#from ..exceptions import *
#from . import parsers


def load_nwchem_defaults(options):
    options.add(
        'nwchem',
        RottenOption(
            keyword='translate_psi4',
            default=True,
            validator=parsers.boolean,
            glossary='Will translate Psi4 options to NWCHEM counterparts'))

    #Memory specifications
    options.add(
        'nwchem',
        RottenOption(
            keyword='memory',
            default='400 mb',
            validator=parsers.parse_memory,
            glossary='Total memory allocation, which can be specified further to stack, heap, and global.'))

    options.add(
        'nwchem',
        RottenOption(
            keyword='total_memory',
            default='400 mb',
            validator=parsers.parse_memory,
            glossary='Total memory allocation, which can be specified further to stack, heap, and global.'))

    options.add(
        'nwchem',
        RottenOption(
            keyword='stack_memory',
            default='100 mb',
            validator=parsers.parse_memory_nomin,
            glossary='Stack memory allocation that can be specified beyond the total memory.'))

    options.add(
        'nwchem',
        RottenOption(
            keyword='heap_memory',
            default='100 mb',
            validator=parsers.parse_memory_nomin,
            glossary='Heap memory allocation that can be specified beyond the total memory allocation.'))

    options.add(
        'nwchem',
        RottenOption(
            keyword='global_memory',
            default='200 mb',
            validator=parsers.parse_memory_nomin,
            glossary='Global memory allocation.'))

    #Geometry options
    options.add(
        'nwchem',
        RottenOption(
            keyword='geometry_center',
            default=True,
            validator=parsers.boolean,
            glossary=
            'Enables or disables the translation of the center of nuclear charge to the origin. Default is move the center of nuclear charge to the origin (center or True).'
        ))

# superseded by molecule spec
#    options.add(
#        'nwchem',
#        RottenOption(
#            keyword='geometry__units',
#            default='Angstroms',
#            validator=parsers.enum('ANGSTROMS AN AU ATOMIC BOHR NM NANOMETERS PM PICOMETERS'),
#            glossary='''keyword specifying the <units> string variable. Default for geometry input is Angstroms. 
#        Other options include atomic units (au, atomic or bohr are acceptable keywords), nanometers (nm), or picometers (pm)
#        '''))

    options.add(
        'nwchem',
        RottenOption(
            keyword='geometry_autosym',
            default=True,
            validator=parsers.boolean,
            glossary='The option to automatically determine the symmetry of the geometric system. Default is on.'))

    options.add(
        'nwchem',
        RottenOption(
            keyword='geometry_autoz',
            default=True,
            validator=parsers.boolean,
            glossary=
            'NWChem generates redundant internal coordinates from user input Cartesian coordinates. Default is on.'))

#    #Top level options
#    options.add(
#        'nwchem',
#        RottenOption(
#            keyword='symmetry',
#            default='c1',
#            validator=parsers.enum(' c1 c2v c3v d2h '),
#            glossary='Optional symmetry directive that specifies the group of the molecule, \
#        which can be automatically determined via the geometry autosym function.'))

    options.add(
        'nwchem',
        RottenOption(
            keyword='charge',
            default=0,
            validator=lambda x: float(x),
            glossary='Charge of the molecule. Nuclear charges plus this should equal an integral number of electrons.'))

    options.add(
        'nwchem',
        RottenOption(
            keyword='print',
            default='medium',
            validator=parsers.enum('debug high medium low none'),
            glossary=
            'Top level print options which can be divided from most specific to least: debug, high, low, and none. Default is medium.'
        ))

    #Basis block
    ###Need to figure out how to delineate per atom basis attribution #TODO
    options.add(
        'nwchem',
        RottenOption(
            keyword='basis',
            default='',
            validator=lambda x: x.lower(),
            glossary='This is just a dummy so can throw an error. Only basis-passing implemented at present.'))


    #SCF block
    options.add(
        'nwchem',
        RottenOption(
            keyword='scf__rhf',
            default=True,
            validator=parsers.boolean,  # rather have one kw as enum but that's not the nwc way ("RHF UHF ROHF"),
            glossary='Reference wave function: RHF, UHF, ROHF. RHF is closed-shell default'))

    options.add(
        'nwchem',
        RottenOption(
            keyword='scf__uhf',
            default=False,
            validator=parsers.boolean,  # rather have one kw as enum but that's not the nwc way ("RHF UHF ROHF"),
            glossary='Reference wave function: RHF, UHF, ROHF. UHF is never default'))

    options.add(
        'nwchem',
        RottenOption(
            keyword='scf__rohf',
            default=True,
            validator=parsers.boolean,  # rather have one kw as enum but that's not the nwc way ("RHF UHF ROHF"),
            glossary='Reference wave function: RHF, UHF, ROHF. ROHF is closed-shell default'))

    options.add(
        'nwchem',
        RottenOption(
            keyword='scf__nopen',
            default=1,
            validator=parsers.nonnegative_integer,
            glossary='Specifies the number of open shells in wave function.'))

    options.add(
        'nwchem',
        RottenOption(
            keyword='scf__thresh',
            default=1.e-4,
            validator=parsers.parse_convergence,
            glossary='SCF Convergence threshold'))

    options.add(
        'nwchem',
        RottenOption(
            keyword='scf__maxiter',
            default=20,
            validator=parsers.nonnegative_integer,
            glossary='Max SCF iteration'))

    options.add(
        'nwchem',
        RottenOption(
            keyword='scf__diis',
            default=False,
            validator=parsers.boolean,
            glossary='DIIS Convergence on or off, default is off'))

    options.add(
        'nwchem',
        RottenOption(
            keyword='scf__direct', 
            default=False, 
            validator=parsers.boolean, 
            glossary='SCF Direct calculation on/off'))

    #SCF semidirect options- Array okay? Check with LB
    options.add(
        'nwchem',
        RottenOption(
            keyword='dft__semidirect__filesize',
            default=0,
            validator=parsers.nonnegative_integer,
            glossary=
            '''SCF Semidirect options: File size allows the user to specify the amount of disk space used per process for storing the integrals in 64-bit words. 
        Default of semidirect leads to directive DIRECT.'''))

    options.add(
        'nwchem',
        RottenOption(
            keyword='dft__semidirect__memsize',
            default=0,
            validator=parsers.nonnegative_integer,
            glossary=
            '''SCF Semidirect options: Memory size allows user to specify the number of 64-bit words to be used per process for caching integrals in memory.
        If the amount of storage space is not available, the code cuts the value in half and checks again, and will continue to do so until request is satisfied.'''
        ))

    options.add(
        'nwchem',
        RottenOption(
            keyword='scfsd__filename',
            default='',
            validator=lambda x: x.lower(),  #need to specify file extensions?
            glossary='''SCF Semidirect options: 
        File name default for integral files are placed into scratch directory (see NWChem manual for details). 
        User-specified string in File name keyword will override default.'''))

    options.add(
        'nwchem',
        RottenOption(
            keyword='scf__sym',
            default=True,
            validator=parsers.boolean,
            glossary='Symmetry specification in SCF for Fock matrix construction on/off'))

    options.add(
        'nwchem',
        RottenOption(
            keyword='scf__adapt',
            default=True,
            validator=parsers.boolean,
            glossary='Force symmetry adaption of the molecular orbitals in SCF'))

    options.add(
        'nwchem',
        RottenOption(
            keyword='scf__tol2e',
            default=1.e-7,
            validator=parsers.parse_convergence,
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

    options.add(
        'nwchem',
        RottenOption(
            keyword='scf__profile',
            default=False,
            validator=parsers.boolean,
            glossary='SCF performance file true/false. Default is false.'))

    options.add(
        'nwchem',
        RottenOption(
            keyword='scf__nr', 
            default=0.1, 
            validator=lambda x: float(x), 
            glossary='Control Netwon-Raphson value.'))

    options.add(
        'nwchem',
        RottenOption(
            keyword='scf__print',
            default='medium',
            validator=parsers.enum('none low medium high debug'),
            glossary='Print options within the SCF block. Default is medium.'))

    options.add(
        'nwchem',
        RottenOption(
            keyword='scf__noprint',
            default='none',
            validator=parsers.enum('none low medium high debug'),
            glossary='Options to not print into output file within the SCF block. Default is none.'))

    #MCSCF block
    #required options
    options.add(
        'nwchem',
        RottenOption(
            keyword='mcscf__active',
            default=0,
            validator=parsers.nonnegative_integer,
            glossary='Number of orbitals in the complete active space self consistent theory (CASSCF).'))

    options.add(
        'nwchem',
        RottenOption(
            keyword='mcscf__actelec',
            default=0,
            validator=parsers.nonnegative_integer,
            glossary='Number of electrons in CASSCF active space. Error will occur if discrepancy is spotted.'))

    options.add(
        'nwchem',
        RottenOption(
            keyword='mcscf__multiplicity',
            default=1,
            validator=parsers.nonnegative_integer,
            glossary='Spin multiplicity in CASSCF/MCSCF block, must be specified for MCSCF to work.'))
    #alternative to mcscf_multiplicity & mcscf_symmetry can use mcscf_state
#    options.add(
#        'nwchem',
#        RottenOption(
#            keyword='mcscf_state',
#            default='',
#            validator=parsers.enum(' 1a1 3b1'),
#            glossary=
#            'Defines the spatial symmetry and multiplicity. Format is [multiplicity][state], e.g. 3b2 for triplet in B2.'))
    options.add(
            'nwchem',
            RottenOption(
                keyword='mcscf_level',
                default=0.1,
                validator=parsers.nonnegative_float,
                glossary= 'The Hessian used in the MCSCF optimization by default level shifted by 0.1 until the orbital gradient normal falls below 0.01 at which point the level shift is reduced to zero. The initial value of 0.1 can be changed as increasing the level shift may make convergence more stable in some instances.'))

    #TDDFT block
    options.add(
        'nwchem',
        RottenOption(
            keyword='tddft__nroots',
            default=1,
            validator=lambda x: float(x),
            glossary='The number of excited state roots in a TDDFT caclulation.'))
    options.add(
        'nwchem',
        RottenOption(
            keyword='tddft__singlet',
            default=True,
            validator=parsers.boolean,
            glossary=
            'Default is on. Requests the TDDFT calculatioin of singlet excited states when reference wave function is closed.'
        ))
    options.add(
        'nwchem',
        RottenOption(
            keyword='tddft__triplet',
            default=True,
            validator=parsers.boolean,
            glossary='Default is on. Request the calculation of triplet excited states when reference wave function is closed.')
    )
    options.add(
        'nwchem',
        RottenOption(
            keyword='tddft__thresh',
            default=1.0e-4,
            validator=lambda x: float(x),
            glossary=
            'The convergence threshold of Davidson\'s iterative alogorithm to solve a matrix eigenvalue problem. Default is 1e-4.'
        ))
    options.add(
        'nwchem',
        RottenOption(
            keyword='tddft__maxiter',
            default=100,
            validator=lambda x: float(x),
            glossary='Max iterations. Default is 100.'))

    #MP2 block
    options.add(
        'nwchem',
        RottenOption(
            keyword='mp2__tight',
            default=False,
            validator=parsers.boolean,
            glossary='''Increase precision of MP2 energy and gradients. Will also change SCF and CPHF precision.
        Tightens thresholds for AO and MO integrals within MP2 code. Default is off'''))

    #   options.add('nwchem', RottenOption(
    #        keyword='mp2_freeze',
    #       default='',
    #      validator='',
    #     glossary=''#TODO #Another array

    options.add(
        'nwchem',
        RottenOption(
            keyword='mp2__scs',
            default=True,
            validator=parsers.boolean,
            glossary='Spin Component Scaled MP2 calculations. Default is on.'))

    options.add(
        'nwchem',
        RottenOption(
            keyword='mp2__fss', 
            default=1.2, 
            validator=lambda x: float(x), 
            glossary='Scaling factor for same spin'))

    options.add(
        'nwchem',
        RottenOption(
            keyword='mp2__fos', 
            default=0.3, 
            validator=lambda x: float(x), 
            glossary='Scaling factor for opposite spin'))

    #DFT block
    #options.add(
    #    'nwchem',
    #    RottenOption(
    #        keyword='dft',
    #        default='rdft',
    #        validator=parsers.enum("RDFT RODFT UDFT ODFT"),
    #        glossary='Defining DFT wavefunction: RDFT, RODFT, UDFT, ODFT (Open shell, singlet).'))
    options.add(
            'nwchem',
            RottenOption(
                keyword= 'dft__rdft',
                default= True,
                validator= parsers.boolean, #switching over to boolean instead of enum for driver simplicity
                glossary= 'Defining DFT wavefunction, using restricted DFT. Default is on.'))

    options.add(
            'nwchem',
            RottenOption(
                keyword= 'dft__rodft',
                default= False,
                validator= parsers.boolean,
                glossary= 'Defining DFT wavefunction, using restricted open shell DFT. Default is off(False).'))
    
    options.add(
            'nwchem',
            RottenOption(
                keyword= 'dft__udft',
                default= False,
                validator= parsers.boolean,
                glossary= 'Defining DFT wavefunction, using unrestricted DFT. Default is off (False).'))
    
    options.add(
            'nwchem',
            RottenOption(
                keyword= 'dft__odft',
                default= False,
                validator= parsers.boolean,
                glossary= 'Defining DFT wavefunction, using open shell, singlet DFT wavefunction. Default is off (False).'))

    #TODO #Array block for dft_grid, dft_vectors
    #DFT Convergence
    options.add(
        'nwchem',
        RottenOption(
            keyword='dft__convergence__energy',
            default=1.e-6,
            #validator=lambda x: float(x),
            validator=parsers.parse_convergence,
            glossary='total energy convergence within the DFT block'))

    options.add(
        'nwchem',
        RottenOption(
            keyword='dft__convergence__density',
            default=1.e-5,
            validator=lambda x: float(x),
            glossary='Total density convergence within DFT block that has RMS difference N and N-1 iterations'))

    options.add(
        'nwchem',
        RottenOption(
            keyword='dft__convergence__gradient',
            default=5.e-4,
            validator=lambda x: float(x),
            glossary=
            'Convergence of the orbital gradient, defined as the DIIS error vector becomes less than a certain value. Default is 5.e-4.'
        ))

    options.add(
        'nwchem',
        RottenOption(
            keyword='dft__convergence__hltol',
            default=0.1,
            validator=lambda x: float(x),
            glossary='HUMO LUMO gap tolerance. Default 0.1 au.'))

    options.add(
        'nwchem',
        RottenOption(
            keyword='dft__convergence__dampon',
            default=0.0,
            validator=lambda x: float(x),
            glossary='Turns on damping when reaching user-specified energy level.'))

    options.add(
        'nwchem',
        RottenOption(
            keyword='dft__convergence__dampoff',
            default=0.0,
            validator=lambda x: float(x),
            glossary='Turns off damping when reaching user-specified energy level.'))

    options.add(
        'nwchem',
        RottenOption(
            keyword='dft__convergence__damp',
            default=0,
            validator=parsers.percentage,
            glossary="Percent of previous iterations mixed with current iteration's density."))

    options.add(
        'nwchem',
        RottenOption(
            keyword='dft__convergence__ncydp',
            default=2,
            validator=parsers.positive_integer,
            glossary='Specifies number of damping cycles. Default is 2.'))

    options.add(
        'nwchem',
        RottenOption(
            keyword='dft__convergence__diison',
            default=0.0,
            validator=lambda x: float(x),
            glossary='Direct inversion of the iterative space can turned on at a user-specified energy.'))

    options.add(
        'nwchem',
        RottenOption(
            keyword='dft__convergence__diisoff',
            default=0.0,
            validator=lambda x: float(x),
            glossary='Direct inversion of the iterative space can be turned off at a user-specified energy. '))

    options.add(
        'nwchem',
        RottenOption(
            keyword='dft__convergence__diis',
            default=10,
            validator=lambda x: float(x),
            glossary='Number of Fock matrices used in direct inversion of iterative subspace [DIIS] extrapolation'))

    options.add(
        'nwchem',
        RottenOption(
            keyword='dft__convergence__ncyds',
            default=30,
            validator=parsers.positive_integer,
            glossary='Specifies number of DIIS [Direct inversion of iterative subspace] cycles needed. Default is 30.')
    )

    options.add(
        'nwchem',
        RottenOption(
            keyword='dft__convergence__levlon',
            default=0.0,
            validator=lambda x: float(x),
            glossary=
            'Turning on the level shifting, which is the amount of shift applied to diagonal elements of the unoccupied block in the Fock matrix.'
        ))

    options.add(
        'nwchem',
        RottenOption(
            keyword='dft__convergence__levloff',
            default=0.0,
            validator=lambda x: float(x),
            glossary='Turning off the level shifting function'))

    options.add(
        'nwchem',
        RottenOption(
            keyword='dft__convergence__lshift',
            default=0.5,
            validator=lambda x: float(x),
            glossary=
            'Specify the amount of shift applied to diagonal elements of the unoccupied block in the Fock matrix. Default is 0.5 a.u.'
        ))

    options.add(
        'nwchem',
        RottenOption(
            keyword='dft__convergence__ncysh',
            default=0,
            validator=lambda x: float(x),
            glossary='Specifies the number of level-shifting cycles are used in the input. Default is 0.'))

    options.add(
        'nwchem',
        RottenOption(
            keyword='dft__convergence__rabuck',
            default=25,
            validator=parsers.positive_integer,
            glossary='''The Rabuck method can be implemented when the initial guess is poor. 
        Will use fractional occupation of the orbital levels during the initial cycles of SCF convergence (A.D. Rabuck and G. Scuseria, J. Chem. Phys 110, 695(1999). 
        This option specifies the number of cycles the Rabuck method is active.'''))

    options.add(
        'nwchem',
        RottenOption(
            keyword='dft__xc',
            default='none',
            validator=parsers.enum(
                '''none acm b3lyp beckehandh pbe0 becke97 becke97-1 becke97-2 becke97-3 becke98 hcth hcth120 hcth 147 hcth407
        becke97ggal hcth407p optx hcthp14 mpw91 mpwlk xft97 cft97 ft97 op bop pbeop HFexch becke88 xperdew91 xpbe96 gill96 lyp perdew81
        perdew86 perdew 91 cpbe96 pw91lda slater vwn_1 vwn_2 vwn_3 vwn_4 vwn_5 vwn_1_rpa'''),
            glossary=
            'Can specify the exchange-correlation functionals in the DFT Module. See NWChem manual for the full list of options.'
        ))

    #DFT block [continued]
    options.add(
        'nwchem',
        RottenOption(
            keyword='dft__iterations',
            default=30, 
            validator=parsers.positive_integer,
            glossary='Specify DFT iterations'))

    options.add(
        'nwchem',
        RottenOption(
            keyword='dft__mult',
            default=1,
            validator=parsers.positive_integer,
            glossary='DFT Mulitiplicity'))

    options.add(
        'nwchem',
        RottenOption(
            keyword='dft__max_ovl',
            default=False,
            validator=parsers.boolean,
            glossary='Lock the ordering of orbitals on/off. Default is off.'))

    options.add(
        'nwchem',
        RottenOption(
            keyword='dft__smear',
            default=0.001,
            validator=lambda x: float(x),
            glossary='Smear keyword allows fractional occuption of the MOs.'))

    options.add(
        'nwchem',
        RottenOption(
            keyword='dft__mulliken',
            default=False,
            validator=parsers.boolean,
            glossary='Mulliken analysis of charge distribution: on/off. Default is off.'))

    options.add(
        'nwchem',
        RottenOption(
            keyword='dft__direct',
            default=False,
            validator=parsers.boolean,
            glossary='Direct calculation of DFT: on/off. Default is off.'))

#    options.add(
#        'nwchem',
#        RottenOption(
#            keyword='dft__semidirect__filesize',
#            default='',  #default is disksize
#            validator=parsers.positive_integer,
#            glossary='The semidirect options control caching of integrals; option defines the file size, default is disk size.'))

#    options.add(
#        'nwchem',
#        RottenOption(keyword='dft__semidirect__memsize', default='', validator=parsers.parse_memory, glossary='memory size in the semidirect options.'))
#
#    options.add(
#        'nwchem',
#        RottenOption(keyword='dft__semidirect__filename', default='', validator=lambda x: x.lower(), glosary='name of file in semidirect options.'))

    options.add(
        'nwchem',
        RottenOption(
            keyword='dft__cgmin',
            default=False,
            validator=parsers.boolean,
            glossary='DFT quadratic convergence algorithm: true/false. Default is false.'))

    options.add(
        'nwchem',
        RottenOption(
            keyword='dft__fukui',
            default=False,
            validator=parsers.boolean,
            glossary='Fukui indices analysis: on/off. Default is off(false).'))
    #DFT Arrays- dft_disp
    options.add(
        'nwchem',
        RottenOption(
            keyword='dft__print',
            default='medium',
            validator=parsers.enum('none low medium high debug'),
            glossary='''Print level options within the DFT block.
        Default is medium, which prints out convergence, final vector symmetries, multipole, parameters,
        and semidirect information.'''))
    
    options.add(
        'nwchem',
        RottenOption(
            keyword='dft__noprint',
            default='none',
            validator=parsers.enum('none low medium high debug'),
            glossary='No print options for the DFT block. Default is none.'))

    #CCSD block
    options.add(
            'nwchem',
            RottenOption(
                keyword= 'ccsd__maxiter',
                default= 20,
                validator= parsers.positive_integer,
                glossary= 'Maximum numbers of iterations; iterations default is 20.'))

    options.add(
            'nwchem',
            RottenOption(
                keyword='ccsd__thresh',
                default= 1.0e-6,
                validator=parsers.parse_convergence,
                glossary= 'Convergence threshold for the iterative part of the calculation.'))
    
    options.add(
            'nwchem',
            RottenOption(
                keyword='ccsd__freeze',
                default=0,
                validator=parsers.nonnegative_integer,
                glossary= 'Number of orbitals to freeze'))  # expand to core/atomic/virtual

    #TCE block
    options.add('nwchem', RottenOption(
        keyword='tce__dft',
        default=False,
        validator=parsers.boolean,
        glossary='Use DFT as TCE reference wave function. If not specified, default is SCF(HF).'))

    options.add('nwchem', RottenOption(
        keyword='tce__scf',
        default=True,
        validator=parsers.boolean,
        glossary='The default TCE reference wavefunction. Set to True. Use False and tce__dft: True to turn on dft reference function.'))
    
    options.add('nwchem', RottenOption(
        keyword='tce',
        default=True,
        validator=parsers.boolean,
        glossary='The switch for turning on the Tensor Contraction Engine (TCE). Not necessarily needed for couple cluster methods of singles and doubles (CCSD), but necessary for couple cluster theory for singles, doubles, and triples (CCSDT) and couple cluster theory for singles, doubles, triples, and quadruples (CCSDTQ). Default is on.'
            ))
    
    options.add('nwchem', RottenOption(
        keyword='tce__ccsd',
        default= False,
        validator= parsers.boolean,
        glossary= 'TCE module option of coupled cluster singles and doubles (CCSD).'))
   
    options.add('nwchem', RottenOption(
        keyword= 'tce__ccsd_act',
        default= False,
        validator= parsers.boolean,
        glossary= 'TCE module option of coupled-cluster singles and active doubles. Can also signify active-space EOMCCSD'))

    options.add('nwchem',RottenOption(
        keyword= 'tce__eaccsd',
        default= False,
        validator= parsers.boolean,
        glossary= 'TCE module option of electron affinity EOMCCSD'))

    options.add('nwchem', RottenOption(
        keyword= 'tce__ipccsd',
        default= False,
        validator= parsers.boolean,
        glossary= 'TCE module option of ionization potential EOMCCSD.'))

    options.add('nwchem', RottenOption(
        keyword='tce__ccsd_pr_t',
        default= False,
        validator= parsers.boolean,
        glossary= 'TCE module option of coupled cluster singles, doubles with perturbative connected triples'))

    options.add('nwchem', RottenOption(
        keyword= 'tce__ccsdt',
        default= False,
        validator= parsers.boolean,
        glossary= 'TCE module option of coupled cluster singles, doubles, and triples. Also activates EOM-CCSDT if needed.'))

    options.add('nwchem', RottenOption(
        keyword= 'tce__ccsd_br_t',
        default= False,
        validator= parsers.boolean,
        glossary= 'TCE module option of coupled cluster singles, doubles, and perturbative connected triples.'))

    options.add('nwchem', RottenOption(
        keyword= 'tce__ccsdta',
        default= False,
        validator= parsers.boolean,
        glossary= '''TCE module option of coupled-cluster singles, doubles, and active triples. Also call for EOM-CCSDT but there are three (3) variants:
        
        1) T3A_LVL_1 uses largest set of triply excited amplitudes defined by at least one occupied and one unoccupied active spinorbital label
        2) T3A_LVL_2 uses triply excited amplitudes that carry at least two occupied and unoccupied active spinorbital label
        3) T3A_LVL_3 uses triply excited amplitudes that are defined by active indices only. 
        
        All require defining relevant set of occupied active alpha and beta spiorbitals (ACTIVE_OA and ACTIVE_OB) and active unoccupied alpha and beta spinorbitals (ACTIVE_VA and ACTIVE_VB).'''))

    options.add('nwchem', RottenOption(
        keyword= 'tce__ccsdtq',
        default= False,
        validator= parsers.boolean,
        glossary= 'TCE module option of coupled cluster singles, doubles, triples, and quadruples. Also, option for EOM-CCSDTQ')) 

    options.add('nwchem', RottenOption(
        keyword= 'tce__ccsd_pr_2_t', #ccsd(2)_t
        default= False,
        validator= parsers.boolean,
        glossary= 'TCE module option of CCSD and perturbative (T)_t correction.'))

    options.add('nwchem', RottenOption(
        keyword= 'tce__ccsd_pr_2_tq', #ccsd(2)_tq
        default= False,
        validator= parsers.boolean,
        glossary= 'TCE module option of CCSD and perturbative CCSD(2) correction.'))

    options.add('nwchem', RottenOption(
        keyword= 'tce__ccsdt_pr_2_q', #ccsdt(2)_q
        default= False,
        validator= parsers.boolean,
        glossary= 'TCE module option of CCSDT and perturbative (CCSDT_Q) quadruples correction.')) #format from site weird

    options.add('nwchem', RottenOption(
        keyword= 'tce__lccd',
        default= False,
        validator= parsers.boolean,
        glossary= 'TCE module option of linearized coupled-cluster doubles (LCCD).'))

    options.add('nwchem', RottenOption(
        keyword= 'tce__lccsd',
        default= False,
        validator= parsers.boolean,
        glossary= 'TCE module option of linearized coupled-cluster singles and doubles.'))

    options.add('nwchem', RottenOption(
        keyword= 'tce__lrccsd',
        default= False,
        validator= parsers.boolean,
        glossary= 'TCE module option of locally renormalized EOMCCSD.'))

    options.add('nwchem', RottenOption(
        keyword= 'tce__lrccsd_pr_t',
        default= False,
        validator= parsers.boolean,
        glossary= 'TCE module option of CCSD with perturbative locally renomalized CCSD(T) correction.'))

    options.add('nwchem', RottenOption(
        keyword= 'tce__lrccsd_pr_tq1',
        default= False,
        validator= parsers.boolean,
        glossary= 'TCE module option of CCSD with perturbative locally renomalized CCSD(TQ)(LR-CCSD(TQ)-1) correction.'))

    options.add('nwchem', RottenOption(
        keyword= 'tce__creomsd_pr_t',
        default= False,
        validator= parsers.boolean,
        glossary= '''TCE module option of EOM CCSD energies and completely renormalized EOMCCSD(T)(IA) correction. NWChem will print out two components:

        1- total energy of the k-th state
        2- the delta-corrected EOMCCSD excitation energy'''))

    options.add('nwchem',RottenOption(
        keyword= 'tce__creom_pr_t_ac',
        default= False,
        validator= parsers.boolean,
        glossary='TCE module option of active space CR-EOMCCSD(T) approach.'))

    options.add('nwchem',RottenOption(
        keyword= 'tce__qcisd',
        default= False,
        validator= parsers.boolean,
        glossary= 'TCE module option of quadratic configuration interaction of singles and doubles.'))

    options.add('nwchem', RottenOption(
        keyword= 'tce__cisd',
        default= False,
        validator= parsers.boolean,
        glossary= 'TCE module option of configuration interaction singles and doubles.'))

    options.add('nwchem', RottenOption(
        keyword='tce__cisdt',
        default= False,
        validator= parsers.boolean,
        glossary= 'TCE module option of configuration interaction singles, doubles, and triples.'))

    options.add('nwchem', RottenOption(
        keyword= 'tce__cisdtq',
        default= False,
        validator= parsers.boolean,
        glossary= 'TCE module option of configuration interaction singles, doubles, triples, and quadruples.'))
    options.add('nwchem', RottenOption(
        keyword='tce__ccd',
        default= False,
        validator= parsers.boolean,
        glossary= 'TCE module option of couple cluster doubles.'))

    options.add('nwchem', RottenOption(
        keyword= 'tce__mp2',
        default= False,
        validator= parsers.boolean,
        glossary= 'TCE module option of Moller-Plesset perturbation theory to the second order (MP2).'))

    options.add('nwchem', RottenOption(
        keyword= 'tce__mp3',
        default= False,
        validator= parsers.boolean,
        glossary= 'TCE module option of Moller-Plesset perturbation theory to the third order (MP3).'))

    options.add('nwchem', RottenOption(
        keyword= 'tce__mp4',
        default= False,
        validator= parsers.boolean,
        glossary= 'TCE module option of Moller-Plesset perturbation theory to the fourth order (MP4).'))

    options.add('nwchem', RottenOption(
        keyword='tce__thresh',
        default=1.e-4,
        validator=parsers.parse_convergence,
        glossary='TCE convergence threshold'))
    
    options.add('nwchem',RottenOption(
        keyword='tce__maxiter', 
        default=100, 
        validator=parsers.positive_integer, 
        glossary='TCE max iterations'))
    
    options.add('nwchem', RottenOption(
        keyword='tce__io',
        default='',
        validator=lambda x: x.upper(),
        glossary='''Parallel I/O scheme specification. Available:
            FORTRAN, EAF, GA, SF, REPLICATED, DRA, GAEAF'''))
    
    options.add('nwchem', RottenOption(
        keyword='tce__diis',
        default=5,
        validator=parsers.positive_integer,
        glossary='''Number of iterations for a DIIS extrapolation to be performed. Will accelerate excitation
            amplitude convergence. Default is 5.'''))
    
        #options.add('nwchem', RottenOption(
        #   keyword='tce_freeze',
        #Array TODO
    
    options.add('nwchem',RottenOption(
        keyword='tce__nroots',
        default=0,
        validator=parsers.nonnegative_integer,
        glossary='Number of excited states. Default is 0.'))
    
    options.add('nwchem', RottenOption(
        keyword='tce__target',
        default=1,
        validator=parsers.positive_integer,
        glossary='TCE target root. Default is 1.'))
    
    options.add('nwchem', RottenOption(
        keyword='tce__targetsym',
        default='',
        validator=lambda x: x.upper(),
        glossary='TCE target symmetry. Default is None.'))
    
    options.add('nwchem', RottenOption(
        keyword='tce__2eorb',
        default=False,
        validator=parsers.boolean,
        glossary='''Economical option of storing two-electron integrals used in coupled cluster calculations,
        taking the difference of the RHF and ROHF values: on/off. Default is off.'''))
    
    options.add('nwchem',RottenOption(
        keyword='tce__2emet', 
        default=1, 
        validator=parsers.positive_integer, 
        glossary='Default is 1.'))
    
    options.add('nwchem', RottenOption(
        keyword='tce__active__oa',
        default=1,  # ?
        validator=parsers.positive_integer,
        glossary='Specify the number of occupied alpha spin-orbitals.'))
    
    options.add('nwchem',RottenOption(
        keyword='tce__active__ob',
        default=1,
        validator=parsers.positive_integer,
        glossary='Specify the number of occupied beta spin-orbitals.'))
    
    options.add('nwchem',RottenOption(
        keyword='tce__active__va',
        default=1,
        validator=parsers.positive_integer,
        glossary='Specify the number of unoccupied alpha spin-orbitals.'))
    
    options.add('nwchem', RottenOption(
        keyword='tce__active__vb',
        default=1,
        validator=parsers.positive_integer,
        glossary='Specify the number of unoccupied beta spin-orbitals.'))
    
    options.add('nwchem', RottenOption(
        keyword='tce__tilesize', 
        default=1, 
        validator=parsers.positive_integer, 
        glossary='Tile size in Tensor Contraction Engine (TCE).'))
    
    #TASK block- do we need? pytests ensure what action we're implementing into qcdb
    options.add('nwchem', RottenOption(
        keyword='task__hf',
        default='energy',
        validator=lambda x: x.upper(),
        glossary='Specify HF (via SCF) task between: energy, gradient, and hessian. Default is energy.'))
    
    options.add('nwchem', RottenOption(
        keyword='task__scf',
        default='',
        validator=lambda x: x.upper(),
        glossary=
        'Specify Self-consistent theory task between: energy, gradient, and hessian. Default is no specified option and will run energy.'
            ))
    
    options.add('nwchem',  RottenOption(
        keyword='task__mcscf',
        default='',
        validator=lambda x: x.upper(),
        glossary=
        'Specify Multiconfiguration self-consistent (MCSCF) theory task between: energy, gradient, and hessian. Default is no specified option and will run energy.'
            ))
    
    options.add('nwchem', RottenOption(
        keyword='task__mp2',
        default='energy',
        validator=lambda x: x.upper(),
        glossary='Specify MP2 [semi-direct] task between: energy, gradient, and hessian. Default is energy.'))
    
    options.add('nwchem', RottenOption(
        keyword='task__direct__mp2',
        default='energy',
        validator=lambda x: x.upper(),
        glossary='Specify MP2 [direct] task between: energy, gradient, and hessian. Default is energy.'))
    
   # options.add('nwchem', RottenOption(
       # keyword='task__rimp2',
      #  default='energy',
     #   validator=lambda x: x.upper(),
    #    glossary='Specify RIMP2 task between: energy, gradient, and hessian. Default is energy.'))
    
   # options.add('nwchem', RottenOption(
       # keyword='task__dft',
       # default='energy',
        #validator=lambda x: x.upper(),
        #glossary='Specify DFT task between: energy, gradient, and hessian. Default is energy.'))
    
   # options.add('nwchem',RottenOption(
    #    keyword='task__sodft',
     #   default='energy',
      #  validator=lambda x: x.upper(),
       # glossary='Specify SODFT task between: energy, gradient, and hessian. Default is energy.'))
    
    #options.add('nwchem', RottenOption(
     #   keyword='task__tce',
      #  default='energy',
       # validator=lambda x: x.upper(),
        #glossary='Specify TCE task between: energy, gradient, and hessian. Default is energy.'))
    
    options.add('nwchem',RottenOption(
                keyword='task__tce__mp2',  #not sure if need here but need to distinguish between MP2/MP3/MP4 and TCE MBn theory
                default='energy',
                validator=lambda x: x.upper(),
                glossary='Specify TCE MP2 task between: energy, gradient, and hessian. Default is energy.'))
    
    options.add('nwchem', RottenOption(
         keyword='task__ccsd',
         default='energy',
         validator=lambda x: x.upper(),
         glossary='Specify CCSD task between energy, gradient, and hessian. Default is energy.'))
    
    options.add('nwchem', RottenOption(
         keyword='task__ccsd(t)',
         default='energy',
         validator=lambda x: x.upper(),
         glossary='Specify CCSD(T) task between energy, gradient, and hessian. Default is energy.'))
    
    options.add('nwchem', RottenOption(
        keyword='task__ccsdt',
        default='energy',
        validator=lambda x: x.upper(),
        glossary='Specify CCSDT task between energy, gradient, and hessian. Default is energy.'))
