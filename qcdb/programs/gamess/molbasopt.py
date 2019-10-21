import uuid
from typing import Dict

import qcelemental as qcel

from ...molecule import Molecule


def muster_and_format_molecule_and_basis_for_gamess(molrec: Dict, ropts: 'Keywords', qbs: 'BasisSet',
                                                    verbose: int = 1) -> str:
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
