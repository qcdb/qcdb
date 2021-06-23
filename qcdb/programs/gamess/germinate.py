import uuid
from typing import Dict

import qcelemental as qcel

from ...molecule import Molecule
from ...util import conv_float2negexp


def muster_molecule_and_basisset(molrec: Dict, ropts: 'Keywords', qbs: 'BasisSet', verbose: int = 1) -> str:
    kwgs = {'accession': uuid.uuid4(), 'verbose': verbose}
    units = 'Bohr'

    native_puream = qbs.has_puream()
    atom_basisset = qbs.print_detail_gamess(return_list=True)

    gamess_data_record_cart = qcel.molparse.to_string(molrec,
                                                      dtype='gamess',
                                                      units=units,
                                                      atom_format=None,
                                                      ghost_format=None,
                                                      width=17,
                                                      prec=12)
    all_atom_lines = gamess_data_record_cart.splitlines()[3:]

    qmol = Molecule(molrec)
    qmol.update_geometry()

    # PSI: FullPointGroupList = ["ATOM", "C_inf_v", "D_inf_h", "C1", "Cs", "Ci", "Cn", "Cnv", "Cnh", "Sn", "Dn", "Dnd", "Dnh", "Td", "Oh", "Ih"]
    # GMS:                                                      C1    Cs    Ci    Cn    Cnv    Cnh          Dn    Dnd    Dnh    Td    Oh
    # GMS:                        Dnh-2   Cnv-4      Dnh-4                                            S2n
    # GMS:    T, Th, O
    # GAMESS Manual: "For linear molecules, choose either Cnv or Dnh, and enter NAXIS as 4. Enter atoms as Dnh with NAXIS=2."

    pg = qmol.full_point_group_with_n()
    if pg == 'ATOM':
        pgn, naxis = 'Dnh', 2
    elif pg == 'C_inf_v':
        pgn, naxis = 'Cnv', 2
    elif pg == 'D_inf_h':
        pgn, naxis = 'Dnh', 4
    elif pg == 'Sn':
        pgn, naxis = 'S2n', qmol.full_pg_n() / 2  # n/2n never tested
    else:
        pgn, naxis = pg, qmol.full_pg_n()

    uniq_atombas_lines = gamess_data_record_cart.splitlines()[:2]  # $data and card -1-
    uniq_atombas_lines.append(f""" {pgn} {naxis}""")  # card -2-
    uniq_atombas_lines.append('')  # empty cards -3- and -4-

    for iat in range(qmol.natom()):
        if iat == qmol.unique(qmol.atom_to_unique(iat)):
            uniq_atombas_lines.append(all_atom_lines[iat])  # card -5U-
            uniq_atombas_lines.extend(atom_basisset[iat].splitlines()[1:])  # cards -6U- and -7U-
            uniq_atombas_lines.append('')  # card -8U-

    uniq_atombas_lines.append(""" $end""")

    ropts.require('GAMESS', 'contrl__coord', 'unique', **kwgs)
    ropts.require('GAMESS', 'contrl__units', {'Bohr': 'bohr', 'Angstrom': 'angs'}[units], **kwgs)
    ropts.require('GAMESS', 'contrl__icharg', int(molrec['molecular_charge']), **kwgs)
    ropts.require('GAMESS', 'contrl__mult', molrec['molecular_multiplicity'], **kwgs)
    ropts.require('GAMESS', 'contrl__ispher', {True: 1, False: -1}[native_puream], **kwgs)

    return '\n'.join(uniq_atombas_lines)


def muster_modelchem(name: str, dertype: int, ropts: 'Keywords', sysinfo: Dict, verbose: int = 1) -> None:
    lowername = name.lower()
    accession = uuid.uuid4()

    #runtyp = {'energy': 'energy',
    #          'gradient': 'gradient',
    #          'hessian': 'hessian',
    #          'properties': 'prop',
    #         }[driver]

    runtyp = {
        0: 'energy',
        1: 'gradient',
        2: 'hessian',
        #'properties': 'prop',
    }[dertype]

    ropts.require('GAMESS', 'contrl__runtyp', runtyp, accession=accession, verbose=verbose)

    if lowername == 'gms-gamess':
        pass

    elif lowername in ['gms-scf', 'gms-hf']:
        ropts.require('GAMESS', 'contrl__mplevl', 0, accession=accession, verbose=verbose)
        ropts.require('GAMESS', 'contrl__cityp', 'none', accession=accession, verbose=verbose)
        ropts.require('GAMESS', 'contrl__cctyp', 'none', accession=accession, verbose=verbose)

    elif lowername == 'gms-mp2':
        ropts.require('GAMESS', 'contrl__mplevl', 2, accession=accession, verbose=verbose)
        ropts.require('GAMESS', 'contrl__cityp', 'none', accession=accession, verbose=verbose)
        ropts.require('GAMESS', 'contrl__cctyp', 'none', accession=accession, verbose=verbose)

    elif lowername == 'gms-ccsd':
        ropts.require('GAMESS', 'contrl__mplevl', 0, accession=accession, verbose=verbose)
        ropts.require('GAMESS', 'contrl__cityp', 'none', accession=accession, verbose=verbose)
        ropts.require('GAMESS', 'contrl__cctyp', 'ccsd', accession=accession, verbose=verbose)

    elif lowername == 'gms-ccsd(t)':
        ropts.require('GAMESS', 'contrl__mplevl', 0, accession=accession, verbose=verbose)
        ropts.require('GAMESS', 'contrl__cityp', 'none', accession=accession, verbose=verbose)
        ropts.require('GAMESS', 'contrl__cctyp', 'ccsd(t)', accession=accession, verbose=verbose)

    elif lowername == 'gms-cr-ccl':
        ropts.require('GAMESS', 'contrl__mplevl', 0, accession=accession, verbose=verbose)
        ropts.require('GAMESS', 'contrl__cityp', 'none', accession=accession, verbose=verbose)
        ropts.require('GAMESS', 'contrl__cctyp', 'cr-ccl', accession=accession, verbose=verbose)

    elif lowername == 'gms-ccsd(tq)':
        ropts.require('GAMESS', 'contrl__mplevl', 0, accession=accession, verbose=verbose)
        ropts.require('GAMESS', 'contrl__cityp', 'none', accession=accession, verbose=verbose)
        ropts.require('GAMESS', 'contrl__cctyp', 'ccsd(tq)', accession=accession, verbose=verbose)

    elif lowername == 'gms-fci':
        ropts.require('GAMESS', 'contrl__mplevl', 0, accession=accession, verbose=verbose)
        ropts.require('GAMESS', 'contrl__cityp', 'aldet', accession=accession, verbose=verbose)
        ropts.require('GAMESS', 'contrl__cctyp', 'none', accession=accession, verbose=verbose)

        ropts.suggest('GAMESS', 'cidet__ncore', sysinfo["fc"]['ncore'], accession=accession, verbose=verbose)
        ropts.suggest('GAMESS', 'cidet__nact', sysinfo["fc"]['nact'], accession=accession, verbose=verbose)
        ropts.suggest('GAMESS', 'cidet__nels', sysinfo["fc"]['nels'], accession=accession, verbose=verbose)
        # TODO FC hack!!!

    # unused from Nuwan


#    elif lowername == 'gms-dft':
#        if dertype == 0:
#            options ['GAMESS']['GAMESS_CONTRL_RUNTYP']['value'] = 'energy'
#            options ['GAMESS']['GAMESS_CONTRL_DFTTYP']['value'] = 'b3lyp'
#        elif dertype == 1:
#            options ['GAMESS']['GAMESS_CONTRL_RUNTYP']['value'] = 'gradient'
#            options ['GAMESS']['GAMESS_CONTRL_DFTTYP']['value'] = 'b3lyp'
#        elif dertype == 2:
#            options ['GAMESS']['GAMESS_CONTRL_RUNTYP']['value'] = 'hessian'
#            options ['GAMESS']['GAMESS_CONTRL_DFTTYP']['value'] = 'b3lyp'
#
#    elif lowername == 'gms-eom-ccsd':
#        if dertype == 0:
#            options ['GAMESS']['GAMESS_CONTRL_RUNTYP']['value'] = 'energy'
#            options ['GAMESS']['GAMESS_CONTRL_CCTYP']['value']  = 'eom-ccsd'
#        elif dertype == 1:
#            options ['GAMESS']['GAMESS_CONTRL_RUNTYP']['value'] = 'gradient'
#            options ['GAMESS']['GAMESS_CONTRL_CCTYP']['value']  = 'eom-ccsd'
#        elif dertype == 2:
#            options ['GAMESS']['GAMESS_CONTRL_RUNTYP']['value'] = 'hessian'
#            options ['GAMESS']['GAMESS_CONTRL_CCTYP']['value']  = 'eom-ccsd'
#
#    elif lowername == 'gms-cis':
#        if dertype == 0:
#            options ['GAMESS']['GAMESS_CONTRL_RUNTYP']['value'] = 'energy'
#            options ['GAMESS']['GAMESS_CONTRL_CITYP']['value']  = 'cis'
#        elif dertype == 1:
#            options ['GAMESS']['GAMESS_CONTRL_RUNTYP']['value'] = 'gradient'
#            options ['GAMESS']['GAMESS_CONTRL_CITYP']['value']  = 'cis'
#        elif dertype == 2:
#            options ['GAMESS']['GAMESS_CONTRL_RUNTYP']['value'] = 'hessian'
#            options ['GAMESS']['GAMESS_CONTRL_CITYP']['value']  = 'cis'
#
#    elif lowername == 'gms-efp':
#        if dertype == 0:
#            options ['GAMESS']['GAMESS_CONTRL_RUNTYP']['value'] = 'energy'
#            options ['GAMESS']['GAMESS_CONTRL_COORD']['value']  = 'fragonly'
#        elif dertype == 1:
#            options ['GAMESS']['GAMESS_CONTRL_COORD']['value']  = 'gradient'
#            options ['GAMESS']['GAMESS_CONTRL_COORD']['value']  = 'fragonly'
#        elif dertype == 2:
#            options ['GAMESS']['GAMESS_CONTRL_RUNTYP']['value'] = 'hessian'
#            options ['GAMESS']['GAMESS_CONTRL_COORD']['value']  = 'fragonly'

    else:
        raise ValidationError("""Requested GAMESS computational methods %d is not available.""" % (lowername))


def muster_inherited_keywords(ropts: 'Keywords', sysinfo: Dict, verbose: int = 1) -> None:
    accession = uuid.uuid4()

    kwgs = {'accession': accession, 'verbose': verbose}
    do_translate = ropts.scroll['QCDB']['TRANSLATE_QCDB'].value

    # qcdb/memory [B] --> gamess/system__mwords [M QW]
    qopt = ropts.scroll['QCDB']['MEMORY']
    if do_translate or qopt.is_required():
        mem = int(qopt.value / 8e6)
        print('\n\nMEMORY', mem, '\n\n')
        ropts.suggest('GAMESS', 'system__mwords', mem, **kwgs)

    # qcdb/reference --> gamess/contrl__scftyp
    # TODO ref or scf__ref?
    qref = ropts.scroll['QCDB']['SCF__REFERENCE'].value
    if qref in ['RHF', 'UHF', 'ROHF']:
        #ref = {'RHF': 'RHF',
        #       'UHF': 'UHF',
        #       'ROHF': 'ROHF'}[ropts.scroll['QCDB']['REFERENCE'].value]
        ropts.suggest('GAMESS', 'contrl__scftyp', qref, **kwgs)

    # qcdb/scf__d_convergence --> gamess/scf__conv
    qopt = ropts.scroll['QCDB']['SCF__D_CONVERGENCE']
    if qopt.disputed():
        conv = conv_float2negexp(qopt.value)
        ropts.suggest('GAMESS', 'scf__conv', conv, **kwgs)

    # --> gamess/ccinp__iconv
    qopt = ropts.scroll['QCDB']['E_CONVERGENCE']
    if qopt.disputed():
        conv = conv_float2negexp(qopt.value)
        ropts.suggest('GAMESS', 'ccinp__iconv', conv, **kwgs)

    # qcdb/freeze_core --> gamess/
    qopt = ropts.scroll["QCDB"]["FREEZE_CORE"]
    if qopt.disputed():
        if qopt.value is True:
            ncore = sysinfo["fc"]["ncore"]
        elif qopt.value is False:
            ncore = 0
        ropts.suggest("GAMESS", "ccinp__ncore", ncore, accession=accession, verbose=verbose)
        ropts.suggest("GAMESS", "mp2__nacore", ncore, accession=accession, verbose=verbose)
