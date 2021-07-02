from ...keywords import AliasKeyword, Keyword, Keywords, parsers


def load_nwchem_keywords(options: Keywords) -> None:
    options.add(
        'nwchem',
        Keyword(keyword='translate_psi4',
                default=True,
                validator=parsers.boolean,
                glossary='Will translate Psi4 options to NWCHEM counterparts'))

    #Memory specifications
    options.add(
        'nwchem',
        Keyword(keyword='memory',
                default='400 mb',
                validator=parsers.parse_memory,
                glossary='Total memory allocation, which can be specified further to stack, heap, and global.'))

    options.add(
        'nwchem',
        Keyword(keyword='memory__total',
                default='400 mb',
                validator=parsers.parse_memory,
                glossary='Total memory allocation, which can be specified further to stack, heap, and global.'))

    options.add(
        'nwchem',
        Keyword(keyword='memory__stack',
                default='100 mb',
                validator=parsers.parse_memory_nomin,
                glossary='Stack memory allocation that can be specified beyond the total memory.'))

    options.add(
        'nwchem',
        Keyword(keyword='memory__heap',
                default='100 mb',
                validator=parsers.parse_memory_nomin,
                glossary='Heap memory allocation that can be specified beyond the total memory allocation.'))

    options.add(
        'nwchem',
        Keyword(keyword='memory__heap',
                default='200 mb',
                validator=parsers.parse_memory_nomin,
                glossary='Global memory allocation.'))

    #Geometry options
    options.add(
        'nwchem',
        Keyword(
            keyword='geometry__center',
            default=True,
            validator=parsers.boolean,
            glossary=
            'Enables or disables the translation of the center of nuclear charge to the origin. Default is move the center of nuclear charge to the origin (center or True).'
        ))

    options.add(
        'nwchem',
        Keyword(keyword='geometry__autosym',
                default=True,
                validator=parsers.boolean,
                glossary='The option to automatically determine the symmetry of the geometric system. Default is on.'))

    options.add(
        'nwchem',
        Keyword(
            keyword='geometry__autoz',
            default=True,
            validator=parsers.boolean,
            glossary=
            'NWChem generates redundant internal coordinates from user input Cartesian coordinates. Default is on.'))

    #Top level options
    options.add(
        'nwchem',
        Keyword(
            keyword='charge',
            default=0,
            validator=lambda x: float(x),
            glossary='Charge of the molecule. Nuclear charges plus this should equal an integral number of electrons.')
    )

    options.add(
        'nwchem',
        Keyword(
            keyword='print',
            default='medium',
            validator=parsers.enum('debug high medium low none'),
            glossary=
            'Top level print options which can be divided from most specific to least: debug, high, low, and none. Default is medium.'
        ))

    #Property block
    options.add(
        'nwchem',
        Keyword(
            keyword='property__all',
            default=False,
            validator=parsers.boolean,
            glossary=
            'All keyword under property generates information about all currently available properties in NWCHEM.'))

    options.add(
        'nwchem',
        Keyword(
            keyword='property__center',
            default='coc',
            validator=parsers.enum('com coc origin arb'),  #arb need cartesian coordinates
            glossary=
            'Center of expansion options for dipole, quadrople, and octupole calculations. com = center of mass; coc = center of charge; origin; arb is any arbitrary point that must be accompanied by cartesian coordinates.'
        ))

    options.add(
        'nwchem',
        Keyword(
            keyword='property__nbofile',
            default=False,
            validator=parsers.boolean,
            glossary=
            'Creates an input file to be used as input to the stand-alone NBO code. All other properties are calculated upon request.'
        ))

    options.add(
        'nwchem',
        Keyword(keyword='property__dipole',
                default=False,
                validator=parsers.boolean,
                glossary='Determine dipole moment of molecule.'))

    options.add(
        'nwchem',
        Keyword(keyword='property__quadrupole',
                default=False,
                validator=parsers.boolean,
                glossary='Determine quadrople moment'))

    options.add(
        'nwchem',
        Keyword(keyword='property__octupole',
                default=False,
                validator=parsers.boolean,
                glossary='Determine octupole moment'))

    options.add(
        'nwchem',
        Keyword(keyword='property__mulliken',
                default=False,
                validator=parsers.boolean,
                glossary='Mulliken population analysis and bond order analysis'))

    options.add(
        'nwchem',
        Keyword(keyword='property__esp',
                default=False,
                validator=parsers.boolean,
                glossary='Electrostatic potential (diamagnetic shielding) at nuclei'))

    options.add(
        'nwchem',
        Keyword(keyword='property__efield',
                default=False,
                validator=parsers.boolean,
                glossary='Electric field at nuclei'))

    options.add(
        'nwchem',
        Keyword(keyword='property__efieldgrad',
                default=False,
                validator=parsers.boolean,
                glossary='Electric field gradient at nuclei'))

    options.add(
        'nwchem',
        Keyword(
            keyword='property__efieldgradz4',  #nwc 6.8 option
            default=False,
            validator=parsers.boolean,
            glossary='Electric field gradients with relativistic effects.'))

    options.add('nwchem',
                Keyword(keyword='property__gshift', default=False, validator=parsers.boolean,
                        glossary='Gshift'))  #nwc 6.8 option

    options.add(
        'nwchem',
        Keyword(keyword='property__electrondensity',
                default=False,
                validator=parsers.boolean,
                glossary='Electron and spin density at nuclei'))

    options.add(
        'nwchem',
        Keyword(keyword='property__hyperfine',
                default=False,
                validator=parsers.boolean,
                glossary='NMR hyperfine coupling (Fermi-Contact and Spin-Dipole expectation values'))
    #discrepancy in site & nwc 6.8 gh wiki but hyperfine also takes multiple integer input needed

    options.add('nwchem',
                Keyword(keyword='property__aimfile',
                        default=False,
                        validator=parsers.boolean,
                        glossary='Generates AIM Wavefunction files (.wfn | .wfx) which can be used with various codes.'
                        ))  #nwc 6.8 option

    options.add(
        'nwchem',
        Keyword(keyword='property__moldenfile',
                default=False,
                validator=parsers.boolean,
                glossary='Generates files using Molden format (.molden), additional compatible files can be created'
                ))  #nwc 6.8 option

    options.add(
        'nwchem',
        Keyword(
            keyword='property__molden_norm',
            default='none',
            validator=parsers.enum('janpa nwchem none'),
            glossary=
            'Generates files using molden format. Has option to normalize to JANPA convention. Highly recommend using spherical basis set'
        ))
    
    options.add(
        'nwchem',
        Keyword(
            keyword= 'property__shielding',
            default= '',
            validator= parsers.atompair,
            glossary= 'Shielding'))

   # options.add(
        #'nwchem',
        #Keyword(
        #    keyword= 'property__spinspin',
        #    default= (1, 1, 2),
        #    validator=,
        #    glossary= ))

            
    #Property adds: shielding, spinspin both use multiple integer input
    #Example: spinspin 1 1 2 reads for 1 coupling pair between atoms 1 and 2

    #Property__reponse: not as nested as other opts. Additional opts on same tier
    #Following property opts need response which req two numeric- one integer for response order; other frequency (float)
    options.add(
        'nwchem',
        Keyword(keyword='property__velocity',
                default=False,
                validator=parsers.boolean,
                glossary='Use modified velocity gauge for electric dipole'))

    options.add(
        'nwchem',
        Keyword(
            keyword='property__orbeta',
            default=False,
            validator=parsers.boolean,
            glossary=' calculation optical rotation \'beta\' directly. See J. Autschbah, Comp. Lett. 3, (2007), 131.'))

    options.add(
        'nwchem',
        Keyword(keyword='property__giao',
                default=False,
                validator=parsers.boolean,
                glossary='Gauge-including atomic orbital optical rotation, will force orbeta'))

    options.add(
        'nwchem',
        Keyword(
            keyword='property__bdtensor',
            default=False,
            validator=parsers.boolean,
            glossary=
            'Fully origin-invariant optical rotation tensor is calculated [the B-tilde]. See J. Autschbach, ChemPhysChem 12 (2011), 3224 and B. Moore II, M. Srebro, J. Autschbach, J. Chem Theory Comput. 8 (2012), 4336'
        ))

    options.add(
        'nwchem',
        Keyword(
            keyword='property__analysis',
            default=False,
            validator=parsers.boolean,
            glossary=
            'Analysis of response tensors in terms of molecular orbitals. if \'pmlocalization\' then the analysis performed in terms of Pipek-Mezey localized MOs'
        ))

    options.add(
        'nwchem',
        Keyword(
            keyword='property__damping',
            default=0.007,
            validator=lambda x: float(x),
            glossary=
            'Complex response functions. See M. Krykunov, M. D. Kundrat, J. Autschbach, J Chem Phy 125 (2006) 194110'))

    options.add(
        'nwchem',
        Keyword(keyword='property__pmlocalization',
                default=False,
                validator=parsers.boolean,
                glossary='Pipek-Mezey localized MOs'))

    options.add(
        'nwchem',
        Keyword(keyword='property__convergence',
                default=1.0e-4,
                validator=parsers.parse_convergence,
                glossary='set CPKS convergence criterion.'))

    #Relativistic block- electronic approximations

    options.add(
        'nwchem',
        Keyword(
            keyword='relativistic__zora',
            default='on',
            validator=parsers.enum_bool('on off true false'),
            glossary=
            'Zeroth Order regular approximation (ZORA) is the spin-free and spin-orbit one-eectron zeroth-order regular approximation. Default is ON.'
        ))

    options.add(
        'nwchem',
        Keyword(
            keyword='relativistic__douglas_kroll',
            default='dkh',
            validator=parsers.enum_bool('on off true false fpp dkh dkfull dk4 k4full'),
            glossary=
            '''Douglas-Kroll approximation based on certain operators. FPP refers to free-particle projection operators; DKH
        and DKFULL are based on external-field projection operators. DKFULL is considerably better approimxations than the former since it includes certain cross-product integral terms ignored in the DKH approach. 
        DK3 refers to third-order Douglas-Kroll approximation without cross-product integral terms; DK3FULL with cross-product integral terms.'''
        ))

    options.add(
        'nwchem',
        Keyword(
            keyword='relativistic__dyall_mod_dirac',
            default='on',
            validator=parsers.enum_bool('on off true false nesc1e nesc2e'),  #NESC2E options too nested
            glossary=
            '''Dyall's modified Dirac Hamiltonian. Default is ON and will default to one-electron approximation.'''))

    #Raman block (used for polarizability calculations; req. property block)
    options.add(
        'nwchem',
        Keyword(keyword='raman__normal',
                default=True,
                validator=parsers.boolean,
                glossary='Normal is default for Raman calculations in NWChem.'))

    options.add(
        'nwchem',
        Keyword(keyword='raman__resonance',
                default=False,
                validator=parsers.boolean,
                glossary='Raman calculations using resonance. Default is normal.'))

    options.add(
        'nwchem',
        Keyword(keyword='raman__lorentzian',
                default=True,
                validator=parsers.boolean,
                glossary='Raman calculations using Lorentzian function. Is default.'))

    options.add(
        'nwchem',
        Keyword(keyword='raman__gaussian',
                default=False,
                validator=parsers.boolean,
                glossary='Raman calculations using the Gaussian function. Default is Lorentzian function.'))

    options.add(
        'nwchem',
        Keyword(
            keyword='raman__low',
            default=0.0,
            validator=lambda x: float(x),
            glossary=
            'Range in that generates Raman spectrum plot in cm^-1. Option LOW and HIGH modify the frequency range'))

    options.add(
        'nwchem',
        Keyword(
            keyword='raman__high',
            default=0.0,
            validator=lambda x: float(x),
            glossary=
            'Range in that generates Raman spectrum plot in cm^-1. Option LOW and HIGH modify the frequency range'))

    options.add(
        'nwchem',
        Keyword(
            keyword='raman__first',
            default=7,
            validator=parsers.nonnegative_integer,
            glossary=
            'Range of indices of normal modes used in the plot. Default is 7. Option FIRST and LAST modify the range of indices.'
        ))

    options.add(
        'nwchem',
        Keyword(
            keyword='raman__last',
            default=8,  #see what would be appropriate range? Just placeholder
            validator=parsers.nonnegative_integer,
            glossary=
            'Range of indices of normal modes used in the plot. Default is 7. Option FIRST and LAST modify the range of indices.'
        ))

    options.add(
        'nwchem',
        Keyword(keyword='raman__width',
                default=20.0,
                validator=lambda x: float(x),
                glossary='Controls width of smoothed peaks in Raman plot. Default is 20.0'))

    options.add(
        'nwchem',
        Keyword(
            keyword='raman__dq',
            default=0.01,
            validator=lambda x: float(x),
            glossary=
            'Size of the steps along the normal modes. Ir related to step size dR in numerical evaultions of the polarizability derivatives.'
        ))

    #SCF block
    options.add(
        'nwchem',
        Keyword(keyword='scf__rhf',
                default=True,
                validator=parsers.boolean,
                glossary='Default reference wave function for SCF calculations. Other options include UHF, ROHF.'))

    options.add(
        'nwchem',
        Keyword(
            keyword='scf__uhf',
            default=False,
            validator=parsers.boolean,  # rather have one kw as enum but that's not the nwc way ("RHF UHF ROHF"),
            glossary='Reference wave function: RHF, UHF, ROHF. UHF is never default'))

    options.add(
        'nwchem',
        Keyword(
            keyword='scf__rohf',
            default=True,
            validator=parsers.boolean,  # rather have one kw as enum but that's not the nwc way ("RHF UHF ROHF"),
            glossary='Reference wave function: RHF, UHF, ROHF. ROHF is closed-shell default'))

    options.add(
        'nwchem',
        Keyword(
            keyword='scf__singlet',
            default=True,
            validator=parsers.boolean,
            glossary=
            'Keyword to specify number of singly occupied orbitals for calculation. Singlet refers to closed shell and is the default.'
        ))
    options.add(
        'nwchem',
        Keyword(keyword='scf__doublet',
                default=False,
                validator=parsers.boolean,
                glossary='Specifies one singly occupied orbital.'))

    options.add(
        'nwchem',
        Keyword(keyword='scf__triplet',
                default=False,
                validator=parsers.boolean,
                glossary='Specifies two singly occupied orbitals.'))

    options.add(
        'nwchem',
        Keyword(keyword='scf__quartet',
                default=False,
                validator=parsers.boolean,
                glossary='Specifies three singly occupied orbitals.'))

    options.add(
        'nwchem',
        Keyword(keyword='scf__quintet',
                default=False,
                validator=parsers.boolean,
                glossary='Specifies four singly occupied orbitals.'))

    options.add(
        'nwchem',
        Keyword(keyword='scf__sextet',
                default=False,
                validator=parsers.boolean,
                glossary='Specifies five singly occupied orbitals.'))

    options.add(
        'nwchem',
        Keyword(keyword='scf__septet',
                default=False,
                validator=parsers.boolean,
                glossary='Specifies six singly occupied orbitals.'))

    options.add(
        'nwchem',
        Keyword(keyword='scf__octet',
                default=False,
                validator=parsers.boolean,
                glossary='Specifies seven singly occupied orbitals.'))

    options.add(
        'nwchem',
        Keyword(
            keyword='scf__nopen',
            default=1,
            validator=parsers.nonnegative_integer,
            glossary=
            'Specifies the number of open shells in wave function, must be used if the number of singly occupied orbitals is more than seven.'
        ))

    options.add(
        'nwchem',
        Keyword(keyword='scf__thresh',
                default=1.e-4,
                validator=parsers.parse_convergence,
                glossary='SCF Convergence threshold'))

    options.add(
        'nwchem',
        Keyword(keyword='scf__maxiter',
                default=20,
                validator=parsers.nonnegative_integer,
                glossary='Max SCF iteration'))

    options.add(
        'nwchem',
        Keyword(keyword='scf__diis',
                default=False,
                validator=parsers.boolean,
                glossary='DIIS Convergence on or off, default is off'))

    options.add(
        'nwchem',
        Keyword(keyword='scf__direct',
                default=False,
                validator=parsers.boolean,
                glossary='SCF Direct calculation on/off'))

    #SCF semidirect options- Array okay? Check with LB
    options.add(
        'nwchem',
        Keyword(
            keyword='dft__semidirect__filesize',
            default=0,
            validator=parsers.nonnegative_integer,
            glossary=
            '''SCF Semidirect options: File size allows the user to specify the amount of disk space used per process for storing the integrals in 64-bit words. 
        Default of semidirect leads to directive DIRECT.'''))

    options.add(
        'nwchem',
        Keyword(
            keyword='dft__semidirect__memsize',
            default=0,
            validator=parsers.nonnegative_integer,
            glossary=
            '''SCF Semidirect options: Memory size allows user to specify the number of 64-bit words to be used per process for caching integrals in memory.
        If the amount of storage space is not available, the code cuts the value in half and checks again, and will continue to do so until request is satisfied.'''
        ))

    options.add(
        'nwchem',
        Keyword(
            keyword='scfsd__filename',
            default='',
            validator=lambda x: x.lower(),  #need to specify file extensions?
            glossary='''SCF Semidirect options: 
        File name default for integral files are placed into scratch directory (see NWChem manual for details). 
        User-specified string in File name keyword will override default.'''))

    options.add(
        'nwchem',
        Keyword(keyword='scf__sym',
                default=True,
                validator=parsers.boolean,
                glossary='Symmetry specification in SCF for Fock matrix construction on/off'))

    options.add(
        'nwchem',
        Keyword(keyword='scf__adapt',
                default=True,
                validator=parsers.boolean,
                glossary='Force symmetry adaption of the molecular orbitals in SCF'))

    options.add(
        'nwchem',
        Keyword(
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

    #    options.add('nwchem', Keyword(
    #       keyword='scf_vectors'
    #      default=''
    #     validator=''
    #    glossary='Specify source and destination of molecular orbital vectors.'))

    options.add(
        'nwchem',
        Keyword(keyword='scf__profile',
                default=False,
                validator=parsers.boolean,
                glossary='SCF performance file true/false. Default is false.'))

    options.add(
        'nwchem',
        Keyword(keyword='scf__nr', default=0.1, validator=lambda x: float(x),
                glossary='Control Newton-Raphson value.'))

    options.add(
        'nwchem',
        Keyword(keyword='scf__print',
                default='medium',
                validator=parsers.enum('none low medium high debug'),
                glossary='Print options within the SCF block. Default is medium.'))

    options.add(
        'nwchem',
        Keyword(keyword='scf__noprint',
                default='none',
                validator=parsers.enum('none low medium high debug'),
                glossary='Options to not print into output file within the SCF block. Default is none.'))

    #MCSCF block
    #required options
    options.add(
        'nwchem',
        Keyword(keyword='mcscf__active',
                default=0,
                validator=parsers.nonnegative_integer,
                glossary='Number of orbitals in the complete active space self consistent theory (CASSCF).'))

    options.add(
        'nwchem',
        Keyword(keyword='mcscf__actelec',
                default=0,
                validator=parsers.nonnegative_integer,
                glossary='Number of electrons in CASSCF active space. Error will occur if discrepancy is spotted.'))

    options.add(
        'nwchem',
        Keyword(keyword='mcscf__multiplicity',
                default=1,
                validator=parsers.nonnegative_integer,
                glossary='Spin multiplicity in CASSCF/MCSCF block, must be specified for MCSCF to work.'))

    options.add('nwchem', Keyword(
        keyword='mcscf__state',
        default='',
        validator= lambda x: x.lower(),
        glossary='Defines the spatial symmetry and multiplicity format is [multiplicity][state], e.g. 3b2 for triplet in B2.'))
    #alternative to mcscf_multiplicity & mcscf_symmetry can use mcscf_state
    #    options.add(
    #        'nwchem',
    #        Keyword(
    #            keyword='mcscf_state',
    #            default='',
    #            validator=parsers.enum(' 1a1 3b1'),
    #            glossary=
    #            'Defines the spatial symmetry and multiplicity. Format is [multiplicity][state], e.g. 3b2 for triplet in B2.'))
    options.add(
        'nwchem',
        Keyword(
            keyword='mcscf_level',
            default=0.1,
            validator=parsers.nonnegative_float,
            glossary=
            'The Hessian used in the MCSCF optimization by default level shifted by 0.1 until the orbital gradient normal falls below 0.01 at which point the level shift is reduced to zero. The initial value of 0.1 can be changed as increasing the level shift may make convergence more stable in some instances.'
        ))

    #Hessians | Vibrational Frequencies
    #Freq block
    options.add(
        'nwchem',
        Keyword(
            keyword='freq__reuse',
            default='',
            validator=lambda x: x.lower(),  #potentially make parsers to identify a proper string with *.hess
            glossary='Allows resuse of previously computed hessian'))
    options.add(
        'nwchem',
        Keyword(
            keyword='freq__animate',
            default=False,
            validator=parsers.boolean,
            glossary='Generates mode animation input files in standard xyz file format for graphics packages like XMol.'
        ))
    #Could also control step size in full nwc
    #other opts: mass redefinition, temp set

    #Driver block
    options.add(
        'nwchem',
        Keyword(
            keyword='driver__loose',
            default=False,
            validator=parsers.boolean,
            glossary=
            'Convergence directive to specify standard sets of values. Can also toggle individual options gmax, grms, xmax, xrms.'
        ))
    options.add(
        'nwchem',
        Keyword(
            keyword='driver__tight',
            default=False,
            validator=parsers.boolean,
            glossary=
            'Convergence directive to specify standard sets of values. Can also toggle individual options gmax, grms, xmax, xrms.'
        ))
    options.add(
        'nwchem',
        Keyword(keyword='driver__gmax',
                default=0,
                validator=lambda x: float(x),
                glossary='Toggling the maximum gradient in the coordinates being used.'))
    options.add(
        'nwchem',
        Keyword(keyword='driver__grms',
                default=0,
                validator=lambda x: float(x),
                glossary='Toggling the root mean square gradient in the coordinates being used.'))
    options.add(
        'nwchem',
        Keyword(keyword='driver__xmax',
                default=0,
                validator=lambda x: float(x),
                glossary='Toggling the maximum of the Cartesian step.'))
    options.add(
        'nwchem',
        Keyword(keyword='driver__gmax',
                default=0,
                validator=lambda x: float(x),
                glossary='Toggling the root mean square of the Cartesian step.'))
    options.add(
        'nwchem',
        Keyword(keyword='driver__trust',
                default=0.3,
                validator=lambda x: float(x),
                glossary='Fixed trust radius to control the step. Default is 0.3.'))
    options.add(
        'nwchem',
        Keyword(
            keyword='driver__sadstp',
            default=0.1,
            validator=lambda x: float(x),
            glossary='Trust radius used for the mode being maximized during a saddle-point search. Default is 0.1'))
    options.add(
        'nwchem',
        Keyword(keyword='driver__maxiter',
                default=20,
                validator=parsers.nonnegative_integer,
                glossary='Geometry optimization iterations'))
    options.add(
        'nwchem',
        Keyword(
            keyword='driver__redoautoz',
            default=False,
            validator=parsers.boolean,
            glossary=
            'Deletes hessian data and regenerates interal coordinates at current geometry. Useful if large change in geometry.'
        ))
    # options.add('nwchem', Keyword(
    #     keyword= 'driver__xyz',
    #     default= '',
    #     validator= lambda x, #string for prefix of *.xyz
    #     glossary= 'Print geometry at each step.'))
    options.add(
        'nwchem',
        Keyword(keyword='driver__noxyz',
                default=False,
                validator=parsers.boolean,
                glossary='No printing of xyz geometries.'))

    #Stepper block
    options.add(
        'nwchem',
        Keyword(keyword='stepper__min',
                default=False,
                validator=parsers.boolean,
                glossary='Specifies minimum energy search.'))
    options.add(
        'nwchem',
        Keyword(keyword='stepper__ts',
                default=False,
                validator=parsers.boolean,
                glossary='Specifies transition state search.'))
    options.add(
        'nwchem',
        Keyword(keyword='stepper__track',
                default=1,
                validator=parsers.nonnegative_integer,
                glossary='Tracks a specific mode during an optimization for transition state search.'))
    options.add(
        'nwchem',
        Keyword(keyword='stepper__maxiter',
                default=20,
                validator=parsers.nonnegative_integer,
                glossary='Stepper iterations maximum'))
    options.add(
        'nwchem',
        Keyword(keyword='stepper__trust',
                default=0.1,
                validator=lambda x: float(x),
                glossary='Size of steps stepper takes are controlled via trust radius (the max step allowed).'))
    options.add(
        'nwchem',
        Keyword(keyword='stepper__convggm',
                default=8.0e-04,
                validator=parsers.parse_convergence,
                glossary='Convergence tolerance for the largest component of the gradient.'))
    options.add(
        'nwchem',
        Keyword(keyword='stepper__convgg',
                default=1.0e-02,
                validator=parsers.parse_convergence,
                glossary='Covergence tolerance for the gradient norm for all degrees of freedom.'))
    options.add(
        'nwchem',
        Keyword(keyword='stepper__convge',
                default=1.0e-04,
                validator=parsers.parse_convergence,
                glossary='Convergence tolerance for the energy differen in the stepper.'))

    #TDDFT block
    options.add(
        'nwchem',
        Keyword(
            keyword='tddft__cis',
            default=False,
            validator=parsers.boolean,
            glossary=
            'Option toggles the Tamm-Dancoff approximation of configuration interaction singles for TDDFT calculation. Default is RPA.'
        ))
    options.add(
        'nwchem',
        Keyword(keyword='tddft__rpa', default=True, validator=parsers.boolean, glossary='Option defaults to RPA.'))
    options.add(
        'nwchem',
        Keyword(keyword='tddft__nroots',
                default=1,
                validator= parsers.positive_integer,
                glossary='The number of excited state roots in a TDDFT caclulation.'))
    options.add(
        'nwchem',
        Keyword(
            keyword='tddft__maxvecs',
            default=1000,
            validator=parsers.nonnegative_integer,
            glossary=
            'Limit the subspace size of Davidson\'s algorithm; in other words, it is the max number of trial vectors that a calculation is allowed to hold. Typically, 10 to 20 trial vectors are needed for each excited state root to be converged. However, it need not exceed the product of the number of occupied orbitals and the number of virtual orbitals. Default is 1000.'
        ))

    options.add(
        'nwchem',
        Keyword(
            keyword='tddft__singlet',
            default=True,
            validator=parsers.boolean,
            glossary=
            'Default is on. Requests the TDDFT calculation of singlet excited states when reference wave function is closed.'
        ))
    options.add(
        'nwchem',
        Keyword(
            keyword='tddft__nosinglet',
            default=False,
            validator=parsers.boolean,
            glossary=
            'Suppresses the calculation of singlet excited states when the reerence wave function is closed shell. Default is SINGLET.'
        ))
    options.add(
        'nwchem',
        Keyword(
            keyword='tddft__triplet',
            default=True,
            validator=parsers.boolean,
            glossary=
            'Default is on. Request the calculation of triplet excited states when reference wave function is closed.')
    )
    options.add(
        'nwchem',
        Keyword(
            keyword='tddft__notriplet',
            default=False,
            validator=parsers.boolean,
            glossary=
            'Option suppresses the calculation of triplet excited states when the reference wave function is closed shell. The default option is TRIPLET.'
        ))
    options.add(
        'nwchem',
        Keyword(
            keyword='tddft__thresh',
            default=1.0e-4,
            validator=lambda x: float(x),
            glossary=
            'The convergence threshold of Davidson\'s iterative alogorithm to solve a matrix eigenvalue problem. Default is 1e-4.'
        ))
    options.add(
        'nwchem',
        Keyword(keyword='tddft__maxiter',
                default=100,
                validator=lambda x: float(x),
                glossary='Max iterations. Default is 100.'))

    options.add(
        'nwchem',
        Keyword(keyword='tddft__target',
                default=1,
                validator=parsers.positive_integer,
                glossary=
                'Option specifies which excited state root is being used for the geometrical derivative calculation.'))

    options.add(
        'nwchem',
        Keyword(
            keyword='tddft__targetsym',
            default='None',  #Need parser for symmetry of interest?
            validator=parsers.enum('None'),
            glossary=
            'Required option when examining excited state geometry optimization. It is common that the order of excited states changes due to the geometry changes over the course of optimization. For frequency calculations, TARGETSYM must be set to none.'
        ))
    options.add(
        'nwchem',
        Keyword(
            keyword='tddft__symmetry',
            default=False,
            validator=parsers.boolean,
            glossary=
            'Module will generate initial guess vectors transforming as the same irreducible representation as TARGETSYM. This causes the final excited state roots to be exclusively dominated by those with the specificed irreducible representation. May be useful for those interested in just optically allowed transitions, or in the geometry optimization of an excited state root with a particular irreducible representation. By default, SYMMETRY not set. Requires TARGETSYM.'
        ))
    options.add(
        'nwchem',
        Keyword(
            keyword='tddft__ecut',
            default=0.0,
            validator=lambda x: float(x),
            glossary=
            'Option enables restricted excitation window TDDFT (REW-TDDFT). This is an approach best suited for core excitations. By specifying this keyword only excitations from occupied states below the energy cutoff will be considered.'
        ))
    #Need EWIN, Alpha, Beta opts that need two integer values
    #Not implementing CIVECS yet
    #GRAD block within TDDFT is too nested

    #MP2 block
    options.add(
        'nwchem',
        Keyword(keyword='mp2__tight',
                default=False,
                validator=parsers.boolean,
                glossary='''Increase precision of MP2 energy and gradients. Will also change SCF and CPHF precision.
        Tightens thresholds for AO and MO integrals within MP2 code. Default is off'''))

    options.add(
        'nwchem',
        Keyword(keyword='mp2__freeze__core',
                default=0,
                validator=parsers.nonnegative_integer,
                glossary='Number of core orbitals to freeze'))

    options.add_alias('nwchem', AliasKeyword(alias='mp2__freeze', target='mp2__freeze__core'))

    options.add(
        'nwchem',
        Keyword(keyword='mp2__freeze__virtual',
                default=0,
                validator=parsers.nonnegative_integer,
                glossary='Number of virtual orbitals to freeze'))

    options.add(
        'nwchem',
        Keyword(keyword='mp2__freeze__atomic',
                default=False,
                validator=parsers.bool_or_elem_dict,
                glossary=
                """Freeze orbitals by element standard lookup (``True``) or specify in dict ``{'O': 1, 'Si': 3}``"""))

    options.add_alias('nwchem', AliasKeyword(alias='mp2__freeze__core__atomic', target='mp2__freeze__atomic'))

    # freeze atomic                                         mp2__freeze__atomic = True
    # freeze atomic O 1 Si 3                                mp2__freeze__atomic = {'O': 1, 'Si': 3}
    # freeze core atomic             freeze atomic          mp2__freeze__core = "atomic"
    # freeze 10                      freeze core 10         mp2__freeze = 10
    # freeze core 10                                        mp2__freeze__core = 10
    # freeze virtual 5                                      mp2__freeze__virtual = 5

    options.add(
        'nwchem',
        Keyword(keyword='mp2__scs',
                default=True,
                validator=parsers.boolean,
                glossary='Spin Component Scaled MP2 calculations. Default is on.'))

    options.add(
        'nwchem',
        Keyword(keyword='mp2__fss', default=1.2, validator=lambda x: float(x),
                glossary='Scaling factor for same spin'))

    options.add(
        'nwchem',
        Keyword(keyword='mp2__fos',
                default=0.3,
                validator=lambda x: float(x),
                glossary='Scaling factor for opposite spin'))

    #DFT block
    options.add(
        'nwchem',
        Keyword(
            keyword='dft__odft',
            default=False,
            validator=parsers.boolean,
            glossary='Defining DFT wavefunction, using open shell, singlet DFT wavefunction. Default is off (False).'))

    #TODO #Array block for dft_grid, dft_vectors
    #DFT Convergence
    options.add(
        'nwchem',
        Keyword(
            keyword='dft__convergence__energy',
            default=1.e-6,
            #validator=lambda x: float(x),
            validator=parsers.parse_convergence,
            glossary='total energy convergence within the DFT block'))

    options.add(
        'nwchem',
        Keyword(keyword='dft__convergence__density',
                default=1.e-5,
                validator=lambda x: float(x),
                glossary='Total density convergence within DFT block that has RMS difference N and N-1 iterations'))

    options.add(
        'nwchem',
        Keyword(
            keyword='dft__convergence__gradient',
            default=5.e-4,
            validator=lambda x: float(x),
            glossary=
            'Convergence of the orbital gradient, defined as the DIIS error vector becomes less than a certain value. Default is 5.e-4.'
        ))

    options.add(
        'nwchem',
        Keyword(keyword='dft__convergence__hltol',
                default=0.1,
                validator=lambda x: float(x),
                glossary='HUMO LUMO gap tolerance. Default 0.1 au.'))

    options.add(
        'nwchem',
        Keyword(keyword='dft__convergence__dampon',
                default=0.0,
                validator=lambda x: float(x),
                glossary='Turns on damping when reaching user-specified energy level.'))

    options.add(
        'nwchem',
        Keyword(keyword='dft__convergence__dampoff',
                default=0.0,
                validator=lambda x: float(x),
                glossary='Turns off damping when reaching user-specified energy level.'))

    options.add(
        'nwchem',
        Keyword(keyword='dft__convergence__damp',
                default=0,
                validator=parsers.percentage,
                glossary="Percent of previous iterations mixed with current iteration's density."))

    options.add(
        'nwchem',
        Keyword(keyword='dft__convergence__ncydp',
                default=2,
                validator=parsers.positive_integer,
                glossary='Specifies number of damping cycles. Default is 2.'))

    options.add(
        'nwchem',
        Keyword(keyword='dft__convergence__diison',
                default=0.0,
                validator=lambda x: float(x),
                glossary='Direct inversion of the iterative space can turned on at a user-specified energy.'))

    options.add(
        'nwchem',
        Keyword(keyword='dft__convergence__diisoff',
                default=0.0,
                validator=lambda x: float(x),
                glossary='Direct inversion of the iterative space can be turned off at a user-specified energy. '))

    options.add(
        'nwchem',
        Keyword(
            keyword='dft__convergence__diis',
            default=10,
            validator=lambda x: float(x),
            glossary='Number of Fock matrices used in direct inversion of iterative subspace [DIIS] extrapolation'))

    options.add(
        'nwchem',
        Keyword(
            keyword='dft__convergence__ncyds',
            default=30,
            validator=parsers.positive_integer,
            glossary='Specifies number of DIIS [Direct inversion of iterative subspace] cycles needed. Default is 30.')
    )

    options.add(
        'nwchem',
        Keyword(
            keyword='dft__convergence__levlon',
            default=0.0,
            validator=lambda x: float(x),
            glossary=
            'Turning on the level shifting, which is the amount of shift applied to diagonal elements of the unoccupied block in the Fock matrix.'
        ))

    options.add(
        'nwchem',
        Keyword(keyword='dft__convergence__levloff',
                default=0.0,
                validator=lambda x: float(x),
                glossary='Turning off the level shifting function'))

    options.add(
        'nwchem',
        Keyword(
            keyword='dft__convergence__lshift',
            default=0.5,
            validator=lambda x: float(x),
            glossary=
            'Specify the amount of shift applied to diagonal elements of the unoccupied block in the Fock matrix. Default is 0.5 a.u.'
        ))

    options.add(
        'nwchem',
        Keyword(keyword='dft__convergence__ncysh',
                default=0,
                validator=lambda x: float(x),
                glossary='Specifies the number of level-shifting cycles are used in the input. Default is 0.'))

    options.add(
        'nwchem',
        Keyword(keyword='dft__convergence__rabuck',
                default=25,
                validator=parsers.positive_integer,
                glossary='''The Rabuck method can be implemented when the initial guess is poor. 
        Will use fractional occupation of the orbital levels during the initial cycles of SCF convergence (A.D. Rabuck and G. Scuseria, J. Chem. Phys 110, 695(1999). 
        This option specifies the number of cycles the Rabuck method is active.'''))

    options.add(
        'nwchem',
        Keyword(keyword='dft__iterations',
                default=30,
                validator=parsers.positive_integer,
                glossary='Specify DFT iterations'))

    options.add(
        'nwchem',
        Keyword(keyword='dft__mult', 
                default=1, 
                validator=parsers.positive_integer, 
                glossary='DFT Mulitiplicity'))

    options.add(
        'nwchem',
        Keyword(keyword='dft__max_ovl',
                default=False,
                validator=parsers.boolean,
                glossary='Lock the ordering of orbitals on/off. Default is off.'))

    options.add(
        'nwchem',
        Keyword(keyword='dft__smear',
                default=0.001,
                validator=lambda x: float(x),
                glossary='Smear keyword allows fractional occuption of the MOs.'))

    options.add(
        'nwchem',
        Keyword(keyword='dft__mulliken',
                default=False,
                validator=parsers.boolean,
                glossary='Mulliken analysis of charge distribution: on/off. Default is off.'))

    options.add(
        'nwchem',
        Keyword(keyword='dft__direct',
                default=False,
                validator=parsers.boolean,
                glossary='Direct calculation of DFT: on/off. Default is off.'))

    #    options.add(
    #        'nwchem',
    #        Keyword(
    #            keyword='dft__semidirect__filesize',
    #            default='',  #default is disksize
    #            validator=parsers.positive_integer,
    #            glossary='The semidirect options control caching of integrals; option defines the file size, default is disk size.'))

    #    options.add(
    #        'nwchem',
    #        Keyword(keyword='dft__semidirect__memsize', default='', validator=parsers.parse_memory, glossary='memory size in the semidirect options.'))
    #
    #    options.add(
    #        'nwchem',
    #        Keyword(keyword='dft__semidirect__filename', default='', validator=lambda x: x.lower(), glosary='name of file in semidirect options.'))

    options.add(
        'nwchem',
        Keyword(keyword='dft__cgmin',
                default=False,
                validator=parsers.boolean,
                glossary='DFT quadratic convergence algorithm: true/false. Default is false.'))

    options.add(
        'nwchem',
        Keyword(keyword='dft__fukui',
                default=False,
                validator=parsers.boolean,
                glossary='Fukui indices analysis: on/off. Default is off(false).'))

    options.add(
        'nwchem',
        Keyword(keyword='dft__print',
                default='medium',
                validator=parsers.enum('none low medium high debug'),
                glossary='''Print level options within the DFT block.
        Default is medium, which prints out convergence, final vector symmetries, multipole, parameters,
        and semidirect information.'''))

    options.add(
        'nwchem',
        Keyword(keyword='dft__noprint',
                default='none',
                validator=parsers.enum('none low medium high debug'),
                glossary='No print options for the DFT block. Default is none.'))
    #DFT XC [Functionals]
    options.add('nwchem',
                Keyword(keyword='dft__xc', 
                        default='', 
                        validator=lambda x: x.lower(), 
                        glossary='DFT functional'))
    #Need to resolve conflicts with nesting for multiple dft functionals in input #TODO

    options.add(
        'nwchem',
        Keyword(keyword='dft__grid__gausleg',
                default=(50, 10),
                validator=parsers.gridradang,
                glossary="""Gauss-Legendre angular grid. Tuple of ``(radial, angular)`` points for
            general setting or dictionary for per-atom with string keys (empty string for general
            case) and tuple values.  Density-Functional-Theory-for-Molecules#gauss-legendre-angular-grid"""))

    options.add(
        'nwchem',
        Keyword(
            keyword='dft__grid__lebedev',
            default=(50, 8),  # totally made up to roughly match gausleg
            validator=parsers.gridradang,
            glossary="""Lebedev angular grid. Tuple of ``(radial, iangquad)`` points for
            general setting or dictionary for per-atom with string keys (empty string for general
            case) and tuple values. See table at link for decoding iangquad index into angular points.
            Density-Functional-Theory-for-Molecules#lebedev-angular-grid"""))

    #CCSD block
    options.add(
        'nwchem',
        Keyword(keyword='ccsd__maxiter',
                default=20,
                validator=parsers.positive_integer,
                glossary='Maximum numbers of iterations; iterations default is 20.'))

    options.add(
        'nwchem',
        Keyword(keyword='ccsd__thresh',
                default=1.0e-6,
                validator=parsers.parse_convergence,
                glossary='Convergence threshold for the iterative part of the calculation.'))

    # freeze 10 == freeze core 10
    # freeze virtual 5
    # freeze atomic
    # freeze atomic O 1 Si 3

    options.add(
        'nwchem',
        Keyword(keyword='ccsd__freeze__core',
                default=0,
                validator=parsers.nonnegative_integer,
                glossary='Number of core orbitals to freeze'))

    options.add_alias('nwchem', AliasKeyword(alias='ccsd__freeze', target='ccsd__freeze__core'))

    options.add(
        'nwchem',
        Keyword(keyword='ccsd__freeze__virtual',
                default=0,
                validator=parsers.nonnegative_integer,
                glossary='Number of virtual orbitals to freeze'))

    options.add(
        'nwchem',
        Keyword(keyword='ccsd__freeze__atomic',
                default=False,
                validator=parsers.bool_or_elem_dict,
                glossary=
                """Freeze orbitals by element standard lookup (``True``) or specify in dict ``{'O': 1, 'Si': 3}``"""))

    options.add_alias('nwchem', AliasKeyword(alias='ccsd__freeze__core__atomic', target='ccsd__freeze__atomic'))

    #TCE block
    options.add(
        'nwchem',
        Keyword(keyword='tce__dft',
                default=False,
                validator=parsers.boolean,
                glossary='Use DFT as TCE reference wave function. If not specified, default is SCF(HF).'))

    options.add(
        'nwchem',
        Keyword(
            keyword='tce__scf',
            default=True,
            validator=parsers.boolean,
            glossary=
            'The default TCE reference wavefunction. Set to True. Use False and tce__dft: True to turn on dft reference function.'
        ))

    options.add(
        'nwchem',
        Keyword(
            keyword='tce',
            default=True,
            validator=parsers.boolean,
            glossary=
            'The switch for turning on the Tensor Contraction Engine (TCE). Not necessarily needed for couple cluster methods of singles and doubles (CCSD), but necessary for couple cluster theory for singles, doubles, and triples (CCSDT) and couple cluster theory for singles, doubles, triples, and quadruples (CCSDTQ). Default is on.'
        ))

    options.add(
        'nwchem',
        Keyword(
            keyword='tce__ccsd',
            default=False,
            validator=parsers.boolean,
            glossary=
            'TCE module option of coupled cluster singles and doubles (CCSD). Can also activate EOM-CCSD with the tce__nroots option set.'
        ))

    options.add(
        'nwchem',
        Keyword(
            keyword='tce__ccsd_act',
            default=False,
            validator=parsers.boolean,
            glossary=
            'TCE module option of coupled-cluster singles and active doubles. Can also signify active-space EOMCCSD'))

    options.add(
        'nwchem',
        Keyword(keyword='tce__eaccsd',
                default=False,
                validator=parsers.boolean,
                glossary='TCE module option of electron affinity EOMCCSD'))

    options.add(
        'nwchem',
        Keyword(keyword='tce__ipccsd',
                default=False,
                validator=parsers.boolean,
                glossary='TCE module option of ionization potential EOMCCSD.'))

    options.add(
        'nwchem',
        Keyword(keyword='tce__ccsd(t)',
                default=False,
                validator=parsers.boolean,
                glossary='TCE module option of coupled cluster singles, doubles with perturbative connected triples'))

    options.add(
        'nwchem',
        Keyword(
            keyword='tce__ccsdt',
            default=False,
            validator=parsers.boolean,
            glossary=
            'TCE module option of coupled cluster singles, doubles, and triples. Also activates EOM-CCSDT if needed.'))

    options.add(
        'nwchem',
        Keyword(keyword='tce__ccsd[t]',
                default=False,
                validator=parsers.boolean,
                glossary='TCE module option of coupled cluster singles, doubles, and perturbative connected triples.'))

    options.add(
        'nwchem',
        Keyword(
            keyword='tce__ccsdta',
            default=False,
            validator=parsers.boolean,
            glossary=
            '''TCE module option of coupled-cluster singles, doubles, and active triples. Also call for EOM-CCSDT but there are three (3) variants:
        
        1) T3A_LVL_1 uses largest set of triply excited amplitudes defined by at least one occupied and one unoccupied active spinorbital label
        2) T3A_LVL_2 uses triply excited amplitudes that carry at least two occupied and unoccupied active spinorbital label
        3) T3A_LVL_3 uses triply excited amplitudes that are defined by active indices only. 
        
        All require defining relevant set of occupied active alpha and beta spiorbitals (ACTIVE_OA and ACTIVE_OB) and active unoccupied alpha and beta spinorbitals (ACTIVE_VA and ACTIVE_VB).'''
        ))

    options.add(
        'nwchem',
        Keyword(
            keyword='tce__ccsdtq',
            default=False,
            validator=parsers.boolean,
            glossary=
            'TCE module option of coupled cluster singles, doubles, triples, and quadruples. Also, option for EOM-CCSDTQ'
        ))

    options.add(
        'nwchem',
        Keyword(
            keyword='tce__ccsd(2)_t',  #ccsd(2)_t
            default=False,
            validator=parsers.boolean,
            glossary='TCE module option of CCSD and perturbative (T)_t correction.'))

    options.add(
        'nwchem',
        Keyword(
            keyword='tce__ccsd(2)_tq',  #ccsd(2)_tq
            default=False,
            validator=parsers.boolean,
            glossary='TCE module option of CCSD and perturbative CCSD(2) correction.'))

    options.add(
        'nwchem',
        Keyword(
            keyword='tce__ccsdt(2)_q',  #ccsdt(2)_q
            default=False,
            validator=parsers.boolean,
            glossary='TCE module option of CCSDT and perturbative (CCSDT_Q) quadruples correction.')
    )  #format from site weird

    options.add(
        'nwchem',
        Keyword(keyword='tce__lccd',
                default=False,
                validator=parsers.boolean,
                glossary='TCE module option of linearized coupled-cluster doubles (LCCD).'))

    options.add(
        'nwchem',
        Keyword(keyword='tce__lccsd',
                default=False,
                validator=parsers.boolean,
                glossary='TCE module option of linearized coupled-cluster singles and doubles.'))

    options.add(
        'nwchem',
        Keyword(keyword='tce__lrccsd',
                default=False,
                validator=parsers.boolean,
                glossary='TCE module option of locally renormalized EOMCCSD.'))

    options.add(
        'nwchem',
        Keyword(keyword='tce__lrccsd(t)',
                default=False,
                validator=parsers.boolean,
                glossary='TCE module option of CCSD with perturbative locally renomalized CCSD(T) correction.'))

    options.add(
        'nwchem',
        Keyword(keyword='tce__lrccsd(tq)1',
                default=False,
                validator=parsers.boolean,
                glossary=
                'TCE module option of CCSD with perturbative locally renomalized CCSD(TQ)(LR-CCSD(TQ)-1) correction.'))

    options.add(
        'nwchem',
        Keyword(
            keyword='tce__creom(t)ac',
                default=False,
                validator=parsers.boolean,
                glossary='TCE module option of active space CR-EOMCCSD(T) approach.'))

    options.add(
        'nwchem',
        Keyword(
            keyword='tce__crccsd[t]',
                default=False,
                validator=parsers.boolean,
                glossary='TCE module option of completely renormalized CCSD[T] method.'))

    options.add(
        'nwchem',
        Keyword(
            keyword='tce__crccsd(t)',
                default=False,
                validator=parsers.boolean,
                glossary='TCE module option of completely renormalized CCSD(T) method.'))

    options.add(
        'nwchem',
        Keyword(keyword='tce__qcisd',
                default=False,
                validator=parsers.boolean,
                glossary='TCE module option of quadratic configuration interaction of singles and doubles.'))

    options.add(
        'nwchem',
        Keyword(keyword='tce__cisd',
                default=False,
                validator=parsers.boolean,
                glossary='TCE module option of configuration interaction singles and doubles.'))

    options.add(
        'nwchem',
        Keyword(keyword='tce__cisdt',
                default=False,
                validator=parsers.boolean,
                glossary='TCE module option of configuration interaction singles, doubles, and triples.'))

    options.add(
        'nwchem',
        Keyword(keyword='tce__cisdtq',
                default=False,
                validator=parsers.boolean,
                glossary='TCE module option of configuration interaction singles, doubles, triples, and quadruples.'))
    options.add(
        'nwchem',
        Keyword(keyword='tce__ccd',
                default=False,
                validator=parsers.boolean,
                glossary='TCE module option of couple cluster doubles.'))

    options.add(
        'nwchem',
        Keyword(keyword='tce__mp2',
                default=False,
                validator=parsers.boolean,
                glossary='TCE module option of Moller-Plesset perturbation theory to the second order (MP2).'))

    options.add(
        'nwchem',
        Keyword(keyword='tce__mp3',
                default=False,
                validator=parsers.boolean,
                glossary='TCE module option of Moller-Plesset perturbation theory to the third order (MP3).'))

    options.add(
        'nwchem',
        Keyword(keyword='tce__mp4',
                default=False,
                validator=parsers.boolean,
                glossary='TCE module option of Moller-Plesset perturbation theory to the fourth order (MP4).'))

    options.add(
        'nwchem',
        Keyword(keyword='tce__dipole',
                default=False,
                validator=parsers.boolean,
                glossary='Dipole moment calculation  built into TCE for both ground- and excited-states.'))

    options.add(
        'nwchem',
        Keyword(keyword='tce__thresh',
                default=1.e-4,
                validator=parsers.parse_convergence,
                glossary='TCE convergence threshold'))

    options.add(
        'nwchem',
        Keyword(keyword='tce__maxiter', default=100, validator=parsers.positive_integer,
                glossary='TCE max iterations'))

    options.add(
        'nwchem',
        Keyword(keyword='tce__io',
                default='',
                validator=lambda x: x.upper(),
                glossary='''Parallel I/O scheme specification. Available:
            FORTRAN, EAF, GA, SF, REPLICATED, DRA, GAEAF'''))

    options.add(
        'nwchem',
        Keyword(keyword='tce__diis',
                default=5,
                validator=parsers.positive_integer,
                glossary='''Number of iterations for a DIIS extrapolation to be performed. Will accelerate excitation
            amplitude convergence. Default is 5.'''))

    options.add(
        'nwchem',
        Keyword(keyword='tce__dipole',
                default=False,
                validator=parsers.boolean,
                glossary='''When set, will undergo another interative step for delta equation to find dipole moments 
                    and one-particle density matrix. Can do for both ground and excited states.'''))

    options.add(
        'nwchem',
        Keyword(keyword='tce__freeze__core',
                default=0,
                validator=parsers.nonnegative_integer,
                glossary='Number of core orbitals to freeze'))

    options.add_alias('nwchem', AliasKeyword(alias='tce__freeze', target='tce__freeze__core'))

    options.add(
        'nwchem',
        Keyword(keyword='tce__freeze__virtual',
                default=0,
                validator=parsers.nonnegative_integer,
                glossary='Number of virtual orbitals to freeze'))

    options.add(
        'nwchem',
        Keyword(keyword='tce__freeze__atomic',
                default=False,
                validator=parsers.bool_or_elem_dict,
                glossary=
                """Freeze orbitals by element standard lookup (``True``) or specify in dict ``{'O': 1, 'Si': 3}``"""))

    options.add_alias('nwchem', AliasKeyword(alias='tce__freeze__core__atomic', target='tce__freeze__atomic'))
    #need to incorporate virtual/core distinction
    #Array TODO

    options.add(
        'nwchem',
        Keyword(keyword='tce__nroots',
                default=0,
                validator=parsers.nonnegative_integer,
                glossary='Number of excited states. Default is 0.'))

    options.add(
        'nwchem',
        Keyword(keyword='tce__target',
                default=1,
                validator=parsers.positive_integer,
                glossary='TCE target root. Default is 1.'))

    options.add(
        'nwchem',
        Keyword(keyword='tce__targetsym',
                default='',
                validator=lambda x: x.upper(),
                glossary='TCE target symmetry. Default is None.'))

    options.add(
        'nwchem',
        Keyword(keyword='tce__2eorb',
                default=False,
                validator=parsers.boolean,
                glossary='''Economical option of storing two-electron integrals used in coupled cluster calculations,
        taking the difference of the RHF and ROHF values: on/off. Default is off.'''))

    options.add('nwchem',
                Keyword(keyword='tce__2emet', default=1, validator=parsers.positive_integer, glossary='Default is 1.'))

    options.add(
        'nwchem',
        Keyword(keyword='tce__active_oa',
                default=1,
                validator=parsers.positive_integer,
                glossary='Specify the number of occupied alpha spin-orbitals.'))

    options.add(
        'nwchem',
        Keyword(keyword='tce__active_ob',
                default=1,
                validator=parsers.positive_integer,
                glossary='Specify the number of occupied beta spin-orbitals.'))

    options.add(
        'nwchem',
        Keyword(keyword='tce__active_va',
                default=1,
                validator=parsers.positive_integer,
                glossary='Specify the number of unoccupied alpha spin-orbitals.'))

    options.add(
        'nwchem',
        Keyword(keyword='tce__active_vb',
                default=1,
                validator=parsers.positive_integer,
                glossary='Specify the number of unoccupied beta spin-orbitals.'))

    options.add(
        'nwchem',
        Keyword(
            keyword='tce__t3a_lvl',
            default=1,
            validator=parsers.positive_integer,
            glossary=
            'Set level in TCE CCSDTA. 1 uses largest set of triply excited amplitudes defined by at least one occupied and one unoccupied active spinorbitals labels. Level 2 uses triply excited amplitudes that carry at least two occupied and unoccupied active spinorbital labeles. Level 3 uses triply excited amplitudes that are defined by active indices only. will need ACTIVE_OA, ACTIVE_OB, ACTIVE_VA, ACTIVE_VB options as well.'
        ))

    options.add(
        'nwchem',
        Keyword(keyword='tce__tilesize',
                default=1,
                validator=parsers.positive_integer,
                glossary='Tile size in Tensor Contraction Engine (TCE).'))

    options.add(
        'nwchem',
        Keyword(
            keyword='tce__print',
            default='medium',
            validator=parsers.enum('debug high medium low none'),
            glossary=
            'TCE level print options which can be divided from most specific to least: debug, high, low, and none. Default is medium.'
        ))

    # Hessian
    options.add(
        "nwchem",
        Keyword(keyword="hessian__thresh",
                default=1.e-6,
                validator=parsers.parse_convergence,
                glossary=
                "Wavefunction convergence for analytic Hessians."
        ))
