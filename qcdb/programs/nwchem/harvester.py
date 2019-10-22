#def muster_psi4options(opt, mol):
#    """Translate psi4 keywords *opt* that have been explicitly set into
#    their NWChem counterparts. Since explicitly set NWChem module keyword
#    values will always be used preferentially to these inferred from
#    psi4, the 'clobber' property is set to False.
#    *mol* = qcdb.Molecule(molecule.create_psi4_string_from_molecule())
#
#    """
#    text = ''
#    options = defaultdict(lambda: defaultdict(dict))
#
#    mult = mol.multiplicity()
#    #print ('mult: ', mult)
#
#    if 'SCF' in opt:
#        if 'REFERENCE' in opt['SCF']:
#            if opt['SCF']['REFERENCE']['value'] == 'UKS':
#                options['NWCHEM']['NWCHEM_DFT']['value'] = 'ODFT'
#
#            elif opt['SCF']['REFERENCE']['value'] == 'RKS':
#                if mult != 1:
#                    options['NWCHEM']['NWCHEM_DFT']['value'] = 'RODFT'
#                else:
#                    pass
#            else:
#                options['NWCHEM']['NWCHEM_SCF']['value'] = \
#                        {'RHF': 'RHF',
#                         'UHF': 'UHF',
#                         'ROHF': 'ROHF'}[opt['SCF']['REFERENCE']['value']]
#
#        if 'SCF_TYPE' in opt['SCF']:
#            if opt['SCF']['SCF_TYPE']['value'] == 'DIRECT':
#                options['NWCHEM']['NWCHEM_DFT_DIRECT']['value'] = 'TRUE'
#                options['NWCHEM']['NWCHEM_SCF_DIRECT']['value'] = 'TRUE'
#
#        if 'D_CONVERGENCE' in opt['SCF']:
#            options['NWCHEM']['NWCHEM_DFT_CONVERGENCE']['value'] = \
#                ['density', opt['SCF']['D_CONVERGENCE']['value']]
#            #This is not correct so has to be fixed
#            options['NWCHEM']['NWCHEM_SCF_THRESH']['value'] = \
#                opt['SCF']['D_CONVERGENCE']['value']
#
#        if 'E_CONVERGENCE' in opt['SCF']:
#            options['NWCHEM']['NWCHEM_DFT_CONVERGENCE']['value'] += \
#                ['energy', opt['SCF']['E_CONVERGENCE']['value']]
#            options['NWCHEM']['NWCHEM_SCF_THRESH']['value'] = \
#                opt['SCF']['E_CONVERGENCE']['value']
#
#        if 'MAXITER' in opt['SCF']:
#            options['NWCHEM']['NWCHEM_DFT_ITERATIONS']['value'] = \
#               opt['SCF']['MAXITER']['value']
#            options['NWCHEM']['NWCHEM_SCF_MAXITER']['value'] = \
#               opt['SCF']['MAXITER']['value']
#
#        if 'DFT_FUNCTIONAL' in opt['SCF']:
#            dft_dispersion_list = \
#                ['B3LYP-D3','BLYP-D3','BP86-D3','M05-D3','M05-2X-D3', \
#                'PBE0-D3']
#            if opt['SCF']['DFT_FUNCTIONAL']['value'] in dft_dispersion_list:
#                options['NWCHEM']['NWCHEM_DFT_XC']['value'], options['NWCHEM']['NWCHEM_DFT_DISP']['value'] = \
#            muster_dft_w_dispersion(opt)
#
#            else:
#                options['NWCHEM']['NWCHEM_DFT_XC']['value'] = muster_dft_functionals(opt)
#
#        if 'DFT_RADIAL_POINTS' and 'DFT_SPHERICAL_POINTS' in opt['SCF']:
#            rad_val = opt['SCF']['DFT_RADIAL_POINTS']['value']
#            sph_val = opt['SCF']['DFT_SPHERICAL_POINTS']['value']
#
#            sph_val_list = \
#                [38,50,74,86,110,146,170,194,230,266,302,350,434,590,770,974,1202,1454,1730,2030,2354, \
#                2702,3074,3470,3890,4334,4802,5294,5810]
#
#            iangquad = sph_val_list.index(sph_val) + 1
#
#            grid_val = []
#            grid_val.append(rad_val)
#            grid_val.append(iangquad)
#
#            #print ('grid_val? ', grid_val)
#
#            options['NWCHEM']['NWCHEM_DFT_GRID']['value'] = grid_val
#
#    for item in options['NWCHEM']:
#        options['NWCHEM'][item]['clobber'] = False
#    return text, options
#
#
#def muster_dft_functionals(opt):
#    """Translate psi4 dft functionals w/o dispersion into NWChem counterparts.
#
#    """
#
#    text = ''
#    options = defaultdict(lambda: defaultdict(dict))
#    dft_functional = ''
#
#    dft_same_name_list = ['B3LYP','PBE0','M05','M05-2X',\
#      'FT97','HCTH','HCTH120', 'DLDF', \
#      'HCTHP14','HCTH407P','B1B95','M06','M06-2X','M06-HF',\
#      'M08-HX','M08-SO','M11','M11-L','MPW1B95','MPW1K','MPWB1K',\
#      'PW6B95','PWB6K', 'BOP']
#
#    val = opt['SCF']['DFT_FUNCTIONAL']['value']
#    if val in dft_same_name_list:
#        options['NWCHEM']['DFT_XC']['value'] = \
#                    {'B3LYP': 'B3LYP',
#                     'PBE0' : 'PBE0',
#                       'M05': 'M05',
#                    'M05-2X':'M05-2X',
#                     'FT97' : 'FT97',
#                     'HCTH' : 'HCTH',
#                   'HCTH120': 'HCTH120',
#                     'DLDF' : 'DLDF',
#                   'HCTHP14': 'HCTHP14',
#                  'HCTH407P': 'HCTH407P',
#                     'B1B95': 'B1B95',
#                       'M06': 'M06',
#                    'M06-2X': 'M06-2X',
#                    'M06-HF': 'M06-HF',
#                    'M08-HX': 'M08-HX',
#                    'M08-SO': 'M08-SO',
#                       'M11': 'M11',
#                     'M11-L': 'M11-L',
#                   'MPW1B95': 'MPW1B95',
#                    'MPW1K' : 'MPW1K',
#                    'MPWB1K': 'MPWB1K',
#                    'PW6B95': 'PW6B95',
#                     'PWB6K': 'PWB6K',
#                      'BOP' : 'BOP'
#                    }[opt['SCF']['DFT_FUNCTIONAL']['value']]
#
#    elif (val == 'HF'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
#            ['HFexch', 1.0]
#
#    elif (val == 'B2PLYP'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
#            ['HFexch', 0.53, 'becke88', 0.47, 'lyp', 0.73, 'mp2', 0.27]
#
#    elif (val == 'BLYP'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
#            ['becke88', 'lyp']
#
#    elif (val == 'BP86'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
#            ['becke88', 'perdew86']
#
#    elif (val == 'PBE'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
#            ['xpbe96', 'cpbe96']
#
#    elif (val == 'PW91'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
#            ['xperdew91', 'perdew91']
#
#    elif (val == 'B97-1'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['becke97-1']
#    elif (val == 'B97-2'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['becke97-2']
#    elif (val == 'B97-GGA1'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['becke97gga1']
#
#    elif (val == 'TPSSH'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['xctpssh']
#
#    elif (val == 'BHANDH'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['beckehandh']
#
#    elif (val == 'CAM-B3LYP'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
#        ['xcamb88','1.00','lyp','0.81','vwn_5','0.19','hfexch','1.00\n','cam','0.33',\
#        'cam_alpha','0.19','cam_beta','0.46']
#
#    elif (val == 'B1LYP'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
#            ['HFexch', 0.25,'becke88', 0.75, 'lyp']
#
#    elif (val == 'B1PW91'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
#            ['HFexch', 0.25,'becke88', 0.75, 'perdew91']
#
#    elif (val == 'B3LYP5'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
#            ['vwn_5', 0.19, 'lyp', 0.81, 'HFexch', 0.20, 'slater',0.8,'becke88','nonlocal',0.72]
#
#    elif (val == 'B86BPBE'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
#            ['becke86b', 'cpbe96']
#
#    elif (val == 'B97') or (val == 'B97-0'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['becke97', 'HFexch', 0.1943]
#    elif (val == 'B97-1P'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['becke97gga1', 'HFexch', 0.1500]
#    elif (val == 'B97-3'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['becke97-3', 'HFexch', 0.2693]
#
#    elif (val == 'BHANDHLYP') or (val == 'BHHLYP'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['becke88', 0.500, 'HFexch', 0.500, 'lyp']
#
#    elif (val == 'LRC-WPBE'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['xwpbe',1.00,'cpbe96',1.0,'hfexch','1.00\n',\
#          'cam',0.3,'cam_alpha',0.00,'cam_beta',1.00]
#    elif (val == 'LRC-WPBEH'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['xwpbe',0.80,'cpbe96',1.0,'hfexch','1.00\n',\
#          'cam',0.2,'cam_alpha',0.20,'cam_beta',0.80]
#
#    elif (val == 'MPW1PW'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['mpw91', 0.75, 'HFexch', 0.25, 'perdew91']
#    elif (val == 'MPWLYP1M'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['mpw91', 0.95, 'HFexch', 0.05, 'lyp']
#    elif (val == 'MPWLYP1W'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['mpw91', 'vwn_5', 0.12, 'lyp', 0.88]
#
#    elif (val == 'PBE0-13'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['xpbe96', 0.6667, 'HFexch', 0.3333, 'cpbe96']
#    elif (val == 'PBEH'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['xpbe96', 0.75, 'HFexch', 0.25, 'cpbe96']
#    elif (val == 'PBELYP1W'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['xpbe96', 'vwn_5', 0.26, 'lyp', 0.74]
#
#    elif (val == 'PW86PBE'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['xperdew86', 'cpbe96']
#
#    elif (val == 'TPSSLYP1W'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['xtpss03', 'vwn_5', 0.26, 'lyp', 0.74]
#
#    elif (val == 'XLYP'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
#        ['slater','-0.0690','becke88',0.722,'xperdew91',0.347,'lyp']
#
#    elif (val == 'WPBE'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
#            ['xwpbe', 1.0, 'cpbe96', 1.0, 'HFexch', 1.0,\
#             '\ncam', 0.40, 'cam_alpha', 0.0, 'cam_beta', 1.0]
#
#    elif (val == 'SB98-1A'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
#            ['becke98','HFexch',0.229]
#
#    elif (val == 'PBEH3C'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['xpbe96', 0.58, 'HFexch', 0.42, 'cpbe96']
#    elif (val == 'PBE1W'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['xpbe96', 'vwn_5', 0.26, 'cpbe96', 0.74]
#    elif (val == 'PBE0-2'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['xpbe96', 0.2063, 'HFexch', 0.7937, 'cpbe96', 0.5, 'mp2', 0.5]
#    elif (val == 'O3LYP'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = [
#            'slater', 0.0706, 'optx', 0.8133, 'HFexch', 0.1161, 'vwn_5', 0.19, 'lyp', 0.81
#        ]
#    elif (val == 'HSE03'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['xpbe96',1.0,'xcampbe96','-0.25','cpbe96','1.0\n',\
#          'cam',0.1061,'cam_alpha',0.0,'cam_beta',0.25]
#    elif (val == 'HSE06'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['xpbe96',1.0,'xcampbe96','-0.25','cpbe96',1.0,'srhfexch','0.25\n',\
#          'cam',0.11,'cam_alpha',0.0,'cam_beta',1.0]
#    elif (val == 'WPBE0'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
#            ['xwpbe', 1.0, 'cpbe96', 1.0, 'HFexch', 1.0, \
#            '\ncam', 0.30, 'cam_alpha', 0.25, 'cam_beta', 0.75]
#
#    else:
#        val = str(val)
#        raise ValidationError(
#            """Requested functional %s is not available. Check NWChem manual and try nwchem_dft_xc keyword instead. """
#            % (val))
#
#    dft_functional = options['NWCHEM']['NWCHEM_DFT_XC']['value']
#
#    return dft_functional
#
#
#def muster_dft_w_dispersion(opt):
#    """Translate psi4 dft functionals w/dispersion into NWChem counterparts.
#
#    """
#
#    text = ''
#    options = defaultdict(lambda: defaultdict(dict))
#    dft_functional = ''
#    dft_dispersion = ''
#
#    val = opt['SCF']['DFT_FUNCTIONAL']['value']
#
#    if (val == 'B3LYP-D3'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['b3lyp']
#        options['NWCHEM']['NWCHEM_DFT_DISP']['value'] = ['vdw', 3]
#
#    elif (val == 'BLYP-D3'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
#            ['becke88', 'lyp']
#        options['NWCHEM']['NWCHEM_DFT_DISP']['value'] = ['vdw', 3]
#
#    elif (val == 'BP86-D3'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = \
#            ['becke88', 'perdew86']
#        options['NWCHEM']['NWCHEM_DFT_DISP']['value'] = ['vdw', 3]
#
#    elif (val == 'M05-D3'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['m05']
#        options['NWCHEM']['NWCHEM_DFT_DISP']['value'] = ['vdw', 3]
#
#    elif (val == 'M05-2X-D3'):
#        options['NWCHEM']['NWCHEM_DFT_XC']['value'] = ['m05-2x']
#        options['NWCHEM']['NWCHEM_DFT_DISP']['value'] = ['vdw', 3]
#
#    elif (val == 'PBE0-D3'):
#        options['NWCHEM']['NWCHEM__DFT_XC']['value'] = ['pbe0']
#        options['NWCHEM']['NWCHEM_DFT_DISP']['value'] = ['vdw', 3]
#
#    else:
#        val = str(val)
#        raise ValidationError(
#            """Requested functional %s is not available. Check NWChem manual and try nwchem_dft_xc keyword instead. """
#            % (val))
#
#    dft_functional = options['NWCHEM']['NWCHEM_DFT_XC']['value']
#    dft_dispersion = options['NWCHEM']['NWCHEM_DFT_DISP']['value']
#
#    return dft_functional, dft_dispersion

#def local_prepare_options_for_modules(changedOnly=False, commandsInsteadDict=False):
#    from ..intf_psi4.options import query_options_defaults_from_psi
#    return query_options_defaults_from_psi(changedOnly=changedOnly)
#
#
#def prepare_options_for_nwchem(options):
#    """Function to take the full snapshot of the liboptions object
#    encoded in dictionary *options*, find the options directable toward
#    NWChem (options['NWCHEM']['NWCHEM_**']) that aren't default, then write
#    a NWCHEM deck with those options.
#
#    """
#    text = ''
#    scf_block_text = ''
#    mp2_block_text = ''
#    dft_block_text = ''
#    tce_block_text = ''
#    ccsd_block_text = ''
#    tddft_block_text = ''
#
#    for opt, val in options['NWCHEM'].items():
#        if opt.startswith('NWCHEM_'):
#            if opt.startswith('NWCHEM_TASK'):
#                pass
#
#            elif opt.startswith('NWCHEM_SCF'):
#                if (opt == 'NWCHEM_SCF_PRINT') or (opt == 'NWCHEM_SCF_NOPRINT'):
#                    if val['has_changed']:
#                        scf_block_text += """%s \"%s\" \n""" % (format_option_for_theory_block_nwchem(
#                            opt[11:], val['value']))
#                else:
#                    if val['has_changed']:
#                        scf_block_text += """%s %s \n""" % (format_option_for_theory_block_nwchem(
#                            opt[11:], val['value']))
#
#            elif opt.startswith('NWCHEM_MP2'):
#                if val['has_changed']:
#                    mp2_block_text += """%s %s \n""" % (format_option_for_theory_block_nwchem(opt[11:], val['value']))
#
#            elif opt.startswith('NWCHEM_DFT'):
#                if (opt == 'NWCHEM_DFT_GRID'):
#                    if val['has_changed']:
#                        dft_block_text += """%s lebedev %s \n""" % (format_option_for_theory_block_nwchem(
#                            opt[11:], val['value']))
#
#                elif (opt == 'NWCHEM_DFT_PRINT') or (opt == 'NWCHEM_DFT_NOPRINT'):
#                    if val['has_changed']:
#                        dft_block_text += """%s \"%s\" \n""" % (format_option_for_theory_block_nwchem(
#                            opt[11:], val['value']))
#                else:
#                    if val['has_changed']:
#                        dft_block_text += """%s %s \n""" % (format_option_for_theory_block_nwchem(
#                            opt[11:], val['value']))
#
#            elif opt.startswith('NWCHEM_CCSD'):
#                if val['has_changed']:
#                    ccsd_block_text += """%s %s \n""" % (format_option_for_theory_block_nwchem(opt[12:], val['value']))
#
#            elif opt.startswith('NWCHEM_TCE'):
#                if (opt == 'NWCHEM_TCE_MODULE'):
#                    if val['has_changed']:
#                        tce_block_text += """%s \n""" % (val['value'])
#                elif (opt == 'NWCHEM_TCE'):
#                    pass
#                else:
#                    if val['has_changed']:
#                        tce_block_text += """%s %s \n""" % (format_option_for_theory_block_nwchem(
#                            opt[11:], val['value']))
#
#            elif opt.startswith('NWCHEM_TDDFT'):
#                if val['has_changed']:
#                    tddft_block_text += """%s %s \n""" % (format_option_for_theory_block_nwchem(
#                        opt[13:], val['value']))
#            else:
#                if val['has_changed']:
#                    text += """%s %s \n""" % (format_option_for_nwchem(opt, val['value']))
#
#    if options['NWCHEM']['NWCHEM_TASK_DFT']['has_changed']:
#        if dft_block_text:
#            text += "dft\n" + dft_block_text + "end\n\n"
#
#    elif options['NWCHEM']['NWCHEM_TASK_TDDFT']['has_changed']:
#        if dft_block_text:
#            text += "dft\n" + dft_block_text + "end\n\n"
#
#        if tddft_block_text:
#            text += "tddft\n" + tddft_block_text + "end\n\n"
#
#    elif options['NWCHEM']['NWCHEM_TCE_DFT']['has_changed']:  # DFT is used as TCE reference wavefunction
#        if dft_block_text:
#            text += "dft\n" + dft_block_text + "end\n\n"
#        if tce_block_text:
#            text += "tce\n" + tce_block_text + "end\n\n"
#
#    else:  #SCF, MP2 block, TCE block w/ scf reference wavefunction,
#        if scf_block_text:
#            text += "scf\n" + scf_block_text + "end\n\n"
#
#        if mp2_block_text:
#            text += "mp2\n" + mp2_block_text + "end\n\n"
#
#        if tce_block_text:
#            text += "tce\n" + tce_block_text + "end\n\n"
#
#        if ccsd_block_text:
#            text += "ccsd\n" + ccsd_block_text + "end\n\n"
#
#    if text:
#        text = text[:-1] + '\n'
#
#    return text

#    elif str(val) == 'RODFT':
#        text += str(val) + '\ncgmin'
