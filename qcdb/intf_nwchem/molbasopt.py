import os
import re
import sys
import uuid
import collections

import qcelemental as qcel

from ..driver import pe


def muster_and_format_molecule_for_nwchem(molrec, ropts, verbose=1):
    accession = sys._getframe().f_code.co_name + '_' + str(uuid.uuid4())
    kwgs = {'accession': accession, 'verbose': verbose}

    units = 'Bohr'
    molcmd = qcel.molparse.to_string(molrec, dtype='nwchem', units=units)

    ropts.require('NWCHEM', 'charge', int(molrec['molecular_charge']), **kwgs)
    if molrec['molecular_multiplicity'] != 1:
        ropts.require('NWCHEM', 'scf__nopen', molrec['molecular_multiplicity'] - 1, **kwgs)
        ropts.require('NWCHEM', 'dft__mult', molrec['molecular_multiplicity'], **kwgs)

    return molcmd


def hide_format_basis_for_cfour(molrec, ropts, native_puream, verbose=1):  #puream):
    """Function to print the BASIS=SPECIAL block for Cfour according
    to the active atoms in Molecule. Special short basis names
    are used by qcdb GENBAS-writer in accordance with
    Cfour constraints.

    """
    accession = uuid.uuid4()

    text = []
    for iat, elem in enumerate(molrec['elem']):
        text.append("""{}:CD_{}""".format(elem.upper(), iat + 1))
    text.append('')
    text.append('')
    text = '\n'.join(text)

    ropts.require('CFOUR', 'BASIS', 'SPECIAL', accession=accession, verbose=verbose)
    #ropts.suggest('CFOUR', 'SPHERICAL', native_puream, accession=accession, verbose=verbose)
    ropts.suggest('QCDB', 'PUREAM', native_puream, accession=accession, verbose=verbose)
    # cfour or qcdb keywords here?
    # req or sugg for puream?

    #options['CFOUR']['CFOUR_BASIS']['value'] = 'SPECIAL'
    #options['CFOUR']['CFOUR_SPHERICAL']['value'] = puream
    #options['CFOUR']['CFOUR_BASIS']['clobber'] = True
    #options['CFOUR']['CFOUR_SPHERICAL']['clobber'] = True
    #options['CFOUR']['CFOUR_BASIS']['superclobber'] = True
    #options['CFOUR']['CFOUR_SPHERICAL']['superclobber'] = True

    return text


def kill_old_format_basis_for_cfour(molrec, puream):
    """Function to print the BASIS=SPECIAL block for Cfour according
    to the active atoms in Molecule. Special short basis names
    are used by qcdb GENBAS-writer in accordance with
    Cfour constraints.

    """
    text = []
    for iat, elem in enumerate(molrec['elem']):
        text.append("""{}:CD_{}""".format(elem.upper(), iat + 1))
    text.append('')
    text.append('')
    text = '\n'.join(text)

    options = collections.defaultdict(lambda: collections.defaultdict(dict))
    options['CFOUR']['CFOUR_BASIS']['value'] = 'SPECIAL'
    options['CFOUR']['CFOUR_SPHERICAL']['value'] = puream

    options['CFOUR']['CFOUR_BASIS']['clobber'] = True
    options['CFOUR']['CFOUR_SPHERICAL']['clobber'] = True

    options['CFOUR']['CFOUR_BASIS']['superclobber'] = True
    options['CFOUR']['CFOUR_SPHERICAL']['superclobber'] = True

    return text, options


def format_basis_for_nwchem(self, basopt):
    """Function to print NWChem-style basis sets into [basis block] according
    to the active atoms in Molecule. Basis sets are loaded from Psi4 basis sets library. 

    """
    text = ''
    for fr in range(self.nfragments()):
        if self.fragment_types[fr] == 'Absent':
            pass
        else:
            text += basopt.print_detail_nwchem()

    text += 'end \n'

    options = defaultdict(lambda: defaultdict(dict))

    return text, options


def format_basis_for_nwchem_puream(self, puream):
    """Function to recognize puream for NWChem 

    """
    text = ''

    if (puream == True):
        #       print ('I am spherical')
        text += "basis spherical \n"
    else:
        #       print ('I am cartesian')
        text += "basis \n"  #Default : cartesian

    options = defaultdict(lambda: defaultdict(dict))

    return text, options


def muster_and_format_basis_for_nwchem(molrec, ropts, qbs, verbose=1):
    kwgs = {'accession': uuid.uuid4(), 'verbose': verbose}
    units = 'Bohr'

    # this is bad b/c user can't reset puream. adjust after figuring out anonymous options better
    native_puream = qbs.has_puream()
    nwc_puream = {True: 'spherical', False: 'cartesian'}[native_puream]
#    nwc_puream = 'cartesian'
    atom_basisset = qbs.print_detail_gamess(return_list=True)

    bascmd = f"""basis {nwc_puream}\n"""  # nwc wants role, not basis name, I guess: f"""basis "{qbs.name}" {nwc_puream}\n"""
    bascmd += qbs.print_detail_nwchem()
    bascmd += "\nend\n"

    #ropts.require('NWCHEM', 'basis__puream', {True: 'spherical', False: 'cartesian'}[native_puream], accession=accession, verbose=verbose)

    return bascmd


def local_prepare_options_for_modules(changedOnly=False, commandsInsteadDict=False):
    from ..intf_psi4.options import query_options_defaults_from_psi
    return query_options_defaults_from_psi(changedOnly=changedOnly)


def fulllocal_prepare_options_for_modules(changedOnly=False, commandsInsteadDict=False):
    """Function to return a string of commands to replicate the
    current state of user-modified options. Used to capture C++
    options information for distributed (sow/reap) input files.

    .. caution:: Some features are not yet implemented. Buy a developer a coffee.

       - Need some option to get either all or changed

       - Need some option to either get dict or set string or psimod command list

       - command return doesn't revoke has_changed setting for unchanged with changedOnly=False

    """
    import collections

    modules = [
        # PSI4 Modules
        "ADC",
        "CCENERGY",
        "CCEOM",
        "CCDENSITY",
        "CCLAMBDA",
        "CCHBAR",
        "CCRESPONSE",
        "CCSORT",
        "CCTRIPLES",
        "CLAG",
        "CPHF",
        "CIS",
        "DCFT",
        "DETCI",
        "DFMP2",
        "DFTSAPT",
        "FINDIF",
        "FNOCC",
        "LMP2",
        "MCSCF",
        "MINTS",
        "MRCC",
        "OCC",
        "OPTKING",
        "PSIMRCC",
        "RESPONSE",
        "SAPT",
        "SCF",
        "STABILITY",
        "THERMO",
        "TRANSQT",
        "TRANSQT2",
        # External Modules
        "CFOUR",
    ]

    options = collections.defaultdict(dict)
    commands = ''
    for opt in core.get_global_option_list():
        if core.has_global_option_changed(opt) or not changedOnly:
            if opt in ['DFT_CUSTOM_FUNCTIONAL', 'EXTERN']:  # Feb 2017 hack
                continue
            val = core.get_global_option(opt)
            options['GLOBALS'][opt] = {'value': val, 'has_changed': core.has_global_option_changed(opt)}
            if isinstance(val, str):
                commands += """core.set_global_option('%s', '%s')\n""" % (opt, val)
            else:
                commands += """core.set_global_option('%s', %s)\n""" % (opt, val)
            #if changedOnly:
            #    print('Appending module %s option %s value %s has_changed %s.' % \
            #        ('GLOBALS', opt, core.get_global_option(opt), core.has_global_option_changed(opt)))
        for module in modules:
            if core.option_exists_in_module(module, opt):
                hoc = core.has_option_changed(module, opt)
                if hoc or not changedOnly:
                    val = core.get_option(module, opt)
                    options[module][opt] = {'value': val, 'has_changed': hoc}
                    if isinstance(val, str):
                        commands += """core.set_local_option('%s', '%s', '%s')\n""" % (module, opt, val)
                    else:
                        commands += """core.set_local_option('%s', '%s', %s)\n""" % (module, opt, val)
                    #if changedOnly:
                    #    print('Appending module %s option %s value %s has_changed %s.' % \
                    #        (module, opt, core.get_option(module, opt), hoc))

    if commandsInsteadDict:
        return commands
    else:
        return options


def format_options_for_nwchem(options):
    """From NWCHEM-directed, non-default options dictionary `options`, write a NWCHEM deck."""

    grouped_options = collections.defaultdict(dict)
    for group_key, val in options.items():
        nesting = group_key.split('__')
        if len(nesting) == 1:
            group, key = 'aaaglobal', nesting[0]
        elif len(nesting) == 2:
            group, key = nesting
        else:
            print(nesting)
            raise ValueError('Nesting!')
#        group, key = group_key.split('__')
        grouped_options[group][key] = val

    grouped_lines = {}
    for group, opts in sorted(grouped_options.items()):
        line = []
#        line.append(f'${group.lower()}')
        for key, val in grouped_options[group].items():
            line.append(' '.join(format_option_for_nwchem(key, val, lop_off=False)))
#        line.append('end\n')
        if group == 'aaaglobal':
            grouped_lines[group] = '\n'.join(line) + '\n'
        else:
            grouped_lines[group] = f'{group.lower()}\n  ' + '\n  '.join(line) + '\nend\n'
#        grouped_lines[group] = textwrap.fill(' '.join(line), initial_indent=' ', subsequent_indent='  ')

    return '\n'.join(grouped_lines.values()) + '\n'


def prepare_options_for_nwchem(options):
    """Function to take the full snapshot of the liboptions object
    encoded in dictionary *options*, find the options directable toward
    NWChem (options['NWCHEM']['NWCHEM_**']) that aren't default, then write
    a NWCHEM deck with those options.

    """
    text = ''
    scf_block_text = ''
    mp2_block_text = ''
    dft_block_text = ''
    tce_block_text = ''
    ccsd_block_text = ''
    tddft_block_text = ''

    for opt, val in options['NWCHEM'].items():
        if opt.startswith('NWCHEM_'):
            if opt.startswith('NWCHEM_TASK'):
                pass

            elif opt.startswith('NWCHEM_SCF'):
                if (opt == 'NWCHEM_SCF_PRINT') or (opt == 'NWCHEM_SCF_NOPRINT'):
                    if val['has_changed']:
                        scf_block_text += """%s \"%s\" \n""" % (format_option_for_theory_block_nwchem(
                            opt[11:], val['value']))
                else:
                    if val['has_changed']:
                        scf_block_text += """%s %s \n""" % (format_option_for_theory_block_nwchem(
                            opt[11:], val['value']))

            elif opt.startswith('NWCHEM_MP2'):
                if val['has_changed']:
                    mp2_block_text += """%s %s \n""" % (format_option_for_theory_block_nwchem(opt[11:], val['value']))

            elif opt.startswith('NWCHEM_DFT'):
                if (opt == 'NWCHEM_DFT_GRID'):
                    if val['has_changed']:
                        dft_block_text += """%s lebedev %s \n""" % (format_option_for_theory_block_nwchem(
                            opt[11:], val['value']))

                elif (opt == 'NWCHEM_DFT_PRINT') or (opt == 'NWCHEM_DFT_NOPRINT'):
                    if val['has_changed']:
                        dft_block_text += """%s \"%s\" \n""" % (format_option_for_theory_block_nwchem(
                            opt[11:], val['value']))
                else:
                    if val['has_changed']:
                        dft_block_text += """%s %s \n""" % (format_option_for_theory_block_nwchem(
                            opt[11:], val['value']))

            elif opt.startswith('NWCHEM_CCSD'):
                if val['has_changed']:
                    ccsd_block_text += """%s %s \n""" % (format_option_for_theory_block_nwchem(opt[12:], val['value']))

            elif opt.startswith('NWCHEM_TCE'):
                if (opt == 'NWCHEM_TCE_MODULE'):
                    if val['has_changed']:
                        tce_block_text += """%s \n""" % (val['value'])
                elif (opt == 'NWCHEM_TCE'):
                    pass
                else:
                    if val['has_changed']:
                        tce_block_text += """%s %s \n""" % (format_option_for_theory_block_nwchem(
                            opt[11:], val['value']))

            elif opt.startswith('NWCHEM_TDDFT'):
                if val['has_changed']:
                    tddft_block_text += """%s %s \n""" % (format_option_for_theory_block_nwchem(
                        opt[13:], val['value']))
            else:
                if val['has_changed']:
                    text += """%s %s \n""" % (format_option_for_nwchem(opt, val['value']))

    if options['NWCHEM']['NWCHEM_TASK_DFT']['has_changed']:
        if dft_block_text:
            text += "dft\n" + dft_block_text + "end\n\n"

    elif options['NWCHEM']['NWCHEM_TASK_TDDFT']['has_changed']:
        if dft_block_text:
            text += "dft\n" + dft_block_text + "end\n\n"

        if tddft_block_text:
            text += "tddft\n" + tddft_block_text + "end\n\n"

    elif options['NWCHEM']['NWCHEM_TCE_DFT']['has_changed']:  # DFT is used as TCE reference wavefunction
        if dft_block_text:
            text += "dft\n" + dft_block_text + "end\n\n"
        if tce_block_text:
            text += "tce\n" + tce_block_text + "end\n\n"

    else:  #SCF, MP2 block, TCE block w/ scf reference wavefunction,
        if scf_block_text:
            text += "scf\n" + scf_block_text + "end\n\n"

        if mp2_block_text:
            text += "mp2\n" + mp2_block_text + "end\n\n"

        if tce_block_text:
            text += "tce\n" + tce_block_text + "end\n\n"

        if ccsd_block_text:
            text += "ccsd\n" + ccsd_block_text + "end\n\n"

    if text:
        text = text[:-1] + '\n'


#    print ("I'm prepare_options_for_nwchem")
#    print (text)

    return text


def prepare_options_for_task_nwchem(options):
    """Function to take the full snapshot of the liboptions object
   encoded in dictionary *options*, find the options directable toward
   NWChem (options['NWCHEM']['NWCHEM_TASK_**']) that aren't default, then write
   a NWCHEM deck with those options. If TCE is set ON, this funcion also write 
   appropriate NWCHEM_TCE_MODULE.  

    """
    task_value = ''
    text = ''

    if options['NWCHEM']['NWCHEM_TCE']['value'] == 'ON':
        for opt, val in options['NWCHEM'].items():
            if opt.startswith('NWCHEM_TASK'):
                if val['has_changed']:
                    if opt.startswith('NWCHEM_TASK_CCSD'):  #task ccsd, ccsd(t) -> task tce
                        text += """tce; %s ; end\n""" % (opt[12:].lower())  # NWCHEM_TCE_MODULE =ccsd or ccsd(t)
                        text += """task tce %s \n""" % (val['value'])
            else:
                pass
    else:
        for opt, val in options['NWCHEM'].items():
            if opt.startswith('NWCHEM_TASK'):
                if val['has_changed']:
                    text += """task %s %s \n""" % (opt[12:].lower(), val['value'])
            else:
                pass

    if text:
        text = text[:-1] + '\n'

    return text


def format_option_for_nwchem(opt, val, lop_off=True):
    """Function to reformat value *val* for option *opt* from python
    into nwchem-speak. 
        
    """
    text = ''
    spaces = ' '

    # Transform string booleans into " "
    if val is True:
        #text += spaces + opt
        return opt.lower(), ''
    elif val is False:
        return '', ''
    #if str(val) == 'TRUE':
    #    text += '%s' % (spaces)
    #elif str(val) == 'FALSE':
    #    pass

    # complete hack
    if opt.upper() == 'MEMORY':
        return opt.lower(), f'{val} byte'


    elif isinstance(val, list):
        for n in range(len(val)):
            text += str(val[n])
            if n < (len(val) - 1):
                text += '%s' % (spaces)

    else:
        text += str(val)

    if lop_off:
        return opt[7:].lower(), text
    else:
        return opt.lower(), text


def format_option_for_theory_block_nwchem(opt, val):
    """Function to reformat value *val* for option *opt* BLOCK from python
    into nwchem-speak. 
        
    """
    text = ''
    spaces = '  '

    # Transform string booleans into " "
    if str(val) == 'TRUE':
        text += '%s' % (spaces)
    elif str(val) == 'FALSE':
        pass

    elif isinstance(val, list):
        for n in range(len(val)):
            text += str(val[n])
            if n < (len(val) - 1):
                text += '%s' % (spaces)

    elif str(val) == 'RODFT':
        text += str(val) + '\ncgmin'

    else:
        text += str(val)

    return opt.lower(), text
