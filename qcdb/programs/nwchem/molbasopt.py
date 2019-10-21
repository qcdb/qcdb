import uuid
from typing import Dict

import qcelemental as qcel


def format_molecule(molrec: Dict, ropts: 'Keywords', verbose: int = 1) -> str:
    kwgs = {'accession': uuid.uuid4(), 'verbose': verbose}

    molcmd, moldata = qcel.molparse.to_string(molrec, dtype='nwchem', units='Bohr', return_data=True)

    for key, val in moldata['keywords'].items():
        ropts.require('NWCHEM', key, val, **kwgs)

    return molcmd


def muster_and_format_basis_for_nwchem(molrec, ropts, qbs, verbose=1):

    # this is bad b/c user can't reset puream. adjust after figuring out anonymous options better
    native_puream = qbs.has_puream()
    nwc_puream = {True: 'spherical', False: 'cartesian'}[native_puream]
    #    nwc_puream = 'cartesian'

    bascmd = f"""basis {nwc_puream}\n"""  # nwc wants role, not basis name, I guess: f"""basis "{qbs.name}" {nwc_puream}\n"""
    bascmd += qbs.print_detail_nwchem()
    bascmd += "\nend\n"

    #ropts.require('NWCHEM', 'basis__puream', {True: 'spherical', False: 'cartesian'}[native_puream], accession=accession, verbose=verbose)

    return bascmd


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
