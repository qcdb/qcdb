#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2017 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

import collections
from typing import Any, Dict, List
from decimal import Decimal as Dm


def _difference(args):
    minuend, subtrahend = args
    return minuend - subtrahend


def _product(args):
    multiplicand, multiplier = args
    return multiplicand * multiplier


def _spin_component_scaling(args):
    os_scale, ss_scale, tot_corl, ss_corl = args
    return os_scale * (tot_corl - ss_corl) + ss_scale * ss_corl


def _spin_component_scaling_wsing(args):
    os_scale, ss_scale, tot_corl, ss_corl, s_corl = args
    return os_scale * (tot_corl - ss_corl - s_corl) + ss_scale * ss_corl + s_corl


def _dispersion_weighting_wsing(args):
    omega, hb_mtd, dd_mtd, s_corl = args
    return omega * hb_mtd + (1.0 - omega) * dd_mtd + s_corl


def omega(args):
    alpha, beta, ratio = args
    return 0.5 * (1.0 + np.tanh(alpha + beta * ratio))


def _linear(args):
    return sum(c * v for c, v in zip(args))


def _solve_in_turn(args: List, coeff: List) -> List[Dict[str, Any]]:
    """Expands a description of a linear equation of variables `args` and weights `coeff` into all rearrangements of the equation."""

    assert len(args) == len(coeff)
    pv0 = []

    for itgt in range(len(args)):
        non_target_args = args[:]
        non_target_args.pop(itgt)

        non_target_coeff = coeff[:]
        non_target_coeff.pop(itgt)
        solve_by = -1 // coeff[itgt]
        non_target_coeff = [solve_by * c for c in non_target_coeff]

        pv0.append({
            'form': args[itgt],
            'func': lambda vv, cc=non_target_coeff: sum(c * v for c, v in zip(vv, cc)),
            'args': non_target_args
        })

    return pv0


def wfn_psivars() -> List[Dict[str, Any]]:
    """Define QCVariable identity equations (e.g., method total = method correlation + HF)."""

    pv0 = []

    # MP2
    pv0.extend(_solve_in_turn(args=['MP2 TOTAL ENERGY', 'HF TOTAL ENERGY', 'MP2 CORRELATION ENERGY'], coeff=[-1, 1,
                                                                                                             1]))
    pv0.extend(
        _solve_in_turn(
            args=['MP2 DOUBLES ENERGY', 'MP2 SAME-SPIN CORRELATION ENERGY', 'MP2 OPPOSITE-SPIN CORRELATION ENERGY'],
            coeff=[-1, 1, 1]))
    pv0.extend(
        _solve_in_turn(args=['MP2 CORRELATION ENERGY', 'MP2 DOUBLES ENERGY', 'MP2 SINGLES ENERGY'], coeff=[-1, 1, 1]))

    pv0.extend(
        _solve_in_turn(args=['CUSTOM SCS-MP2 TOTAL ENERGY', 'CUSTOM SCS-MP2 CORRELATION ENERGY', 'HF TOTAL ENERGY'],
                       coeff=[-1, 1, 1]))

    # SCS-MP2
    pv0.append({
        'form': 'SCS-MP2 CORRELATION ENERGY',
        'func': _spin_component_scaling_wsing,
        'args': [Dm(6) / Dm(5), Dm(1) / Dm(3), 'MP2 CORRELATION ENERGY', 'MP2 SAME-SPIN CORRELATION ENERGY', 'MP2 SINGLES ENERGY']})  # yapf: disable
    pv0.append({
        'form': 'SCS-MP2 TOTAL ENERGY',
        'func': sum,
        'args': ['HF TOTAL ENERGY', 'SCS-MP2 CORRELATION ENERGY']
    })

    # SCS(N)-MP2
    pv0.append({
        'form': 'SCS(N)-MP2 CORRELATION ENERGY',
        'func': _spin_component_scaling_wsing,
        'args': [Dm(0), Dm(1.76), 'MP2 CORRELATION ENERGY', 'MP2 SAME-SPIN CORRELATION ENERGY', 'MP2 SINGLES ENERGY']})  # yapf: disable
    pv0.append({
        'form': 'SCS(N)-MP2 TOTAL ENERGY',
        'func': sum,
        'args': ['HF TOTAL ENERGY', 'SCS(N)-MP2 CORRELATION ENERGY']
    })

    # DW-MP2
    # only defined at the (IE) reaction level (like SAPT)
    #    dwmp2['DW-MP2 OMEGA'][mt] = \
    #        rxnm_contract_expand(df.xs('HF TOTAL ENERGY', level='psivar').xs(mtl[0], level='meta')) / \
    #        rxnm_contract_expand(df.xs('MP2 TOTAL ENERGY', level='psivar').xs(mtl[1], level='meta'))
    #    df_omega = omega([0.15276, 1.89952, df_omega])

    #pv0.append({
    #    'form': 'DW-MP2 CORRELATION ENERGY',
    #    'func': _dispersion_weighting_wsing,
    #    'args': ['DW_MP2 OMEGA', 'MP2 CORRELATION ENERGY', 'SCS-MP2 CORRELATION ENERGY', 'MP2 SINGLES ENERGY']})
    #pv0.append({
    #    'form': 'DW-MP2 TOTAL ENERGY',
    #    'func': sum,
    #    'args': ['HF TOTAL ENERGY', 'DW-MP2 CORRELATION ENERGY']})

    # MPN
    for mpn in range(3, 20):
        pv0.extend(
            _solve_in_turn(
                args=['MP{} TOTAL ENERGY'.format(mpn), 'HF TOTAL ENERGY', 'MP{} CORRELATION ENERGY'.format(mpn)],
                coeff=[-1, 1, 1]))

    # CCSD
    pv0.extend(
        _solve_in_turn(args=['CCSD TOTAL ENERGY', 'HF TOTAL ENERGY', 'CCSD CORRELATION ENERGY'], coeff=[-1, 1, 1]))
    pv0.extend(
        _solve_in_turn(
            args=['CCSD DOUBLES ENERGY', 'CCSD SAME-SPIN CORRELATION ENERGY', 'CCSD OPPOSITE-SPIN CORRELATION ENERGY'],
            coeff=[-1, 1, 1]))
    pv0.extend(
        _solve_in_turn(args=['CCSD CORRELATION ENERGY', 'CCSD DOUBLES ENERGY', 'CCSD SINGLES ENERGY'],
                       coeff=[-1, 1, 1]))

    # CCSD(T)
    pv0.extend(
        _solve_in_turn(args=['CCSD(T) CORRELATION ENERGY', 'CCSD CORRELATION ENERGY', '(T) CORRECTION ENERGY'],
                       coeff=[-1, 1, 1]))
    pv0.extend(
        _solve_in_turn(args=['CCSD(T) TOTAL ENERGY', 'HF TOTAL ENERGY', 'CCSD(T) CORRELATION ENERGY'],
                       coeff=[-1, 1, 1]))
    pv0.extend(
        _solve_in_turn(args=['CCSD(T) CORRELATION ENERGY', 'CCSD CORRELATION ENERGY', '(T) CORRECTION ENERGY'],
                       coeff=[-1, 1, 1]))  # duplicate of first so that all minimal combinations covered

    # CCSD[T]
    pv0.extend(
        _solve_in_turn(args=['CCSD[T] TOTAL ENERGY', 'HF TOTAL ENERGY', 'CCSD[T] CORRELATION ENERGY'],
                       coeff=[-1, 1, 1]))
    pv0.extend(
        _solve_in_turn(args=['CCSD[T] CORRELATION ENERGY', 'CCSD CORRELATION ENERGY', '[T] CORRECTION ENERGY'],
                       coeff=[-1, 1, 1]))

    # FCI
    pv0.extend(_solve_in_turn(args=['FCI TOTAL ENERGY', 'HF TOTAL ENERGY', 'FCI CORRELATION ENERGY'], coeff=[-1, 1,
                                                                                                             1]))

    # Generics
    pv0.extend(_solve_in_turn(args=['CC TOTAL ENERGY', 'HF TOTAL ENERGY', 'CC CORRELATION ENERGY'], coeff=[-1, 1, 1]))
    pv0.extend(_solve_in_turn(args=['CI TOTAL ENERGY', 'HF TOTAL ENERGY', 'CI CORRELATION ENERGY'], coeff=[-1, 1, 1]))

    # DFT

    #   fctl
    #       TODO want B97 here?
    for fctl in [
            'B3LYP', 'B3LYP5', 'WPBE', 'PBE', 'CAM-B3LYP', 'B97', 'WB97X', 'SVWN', 'PW91', 'BLYP', 'PW86PBE', 'FT97',
            'BOP', 'MPWPW', 'SOGGA11', 'BP86', 'B86BPBE', 'PW6B95', 'PBE0', 'SOGGA11-X', 'MN15'
    ]:
        pv0.extend(
            _solve_in_turn(args=['{} TOTAL ENERGY'.format(fctl), '{} FUNCTIONAL TOTAL ENERGY'.format(fctl)],
                           coeff=[-1, 1]))

    #   fctl + D
    for dash in ['-D2', '-D3', '-D3(BJ)', '-D3M', '-D3M(BJ)']:
        for fctl in ['B3LYP', 'B3LYP5', 'PBE', 'B97', 'BLYP', 'BP86', 'PBE0', 'WPBE']:
            pv0.extend(
                _solve_in_turn(args=[
                    '{}{} TOTAL ENERGY'.format(fctl, dash),
                    '{} FUNCTIONAL TOTAL ENERGY'.format('B97-D' if fctl == 'B97' else fctl),
                    '{}{} DISPERSION CORRECTION ENERGY'.format(fctl, dash)
                ],
                               coeff=[-1, 1, 1]))

    #   fctl + DH
    for fctl in ['B2PLYP', 'DSD-PBEP86', 'PBE0-2', 'B2GPPLYP', 'PTPSS', 'PWPB95', 'DSD-BLYP', 'PBE0-DH']:
        pv0.extend(
            _solve_in_turn(args=[
                '{} TOTAL ENERGY'.format(fctl), '{} FUNCTIONAL TOTAL ENERGY'.format(fctl),
                '{} DOUBLE-HYBRID CORRECTION ENERGY'.format(fctl)
            ],
                           coeff=[-1, 1, 1]))

    #   fctl + D + DH
    #   no psivar for fctl + dh, which would be the more restrictive def
    for dash in ['-D2', '-D3', '-D3(BJ)', '-D3M', '-D3M(BJ)']:
        for fctl in ['B2PLYP']:
            pv0.extend(
                _solve_in_turn(args=[
                    '{}{} TOTAL ENERGY'.format(fctl, dash), '{} TOTAL ENERGY'.format(fctl),
                    '{}{} DISPERSION CORRECTION ENERGY'.format(fctl, dash)
                ],
                               coeff=[-1, 1, 1]))

    #   misc.
    pv0.extend(
        _solve_in_turn(args=[
            'DLDF+D09 TOTAL ENERGY', 'DLDF+D09 FUNCTIONAL TOTAL ENERGY', 'DLDF-DAS2009 DISPERSION CORRECTION ENERGY'
        ],
                       coeff=[-1, 1, 1]))

    pv0.extend(
        _solve_in_turn(
            args=['WB97X-D TOTAL ENERGY', 'WB97X FUNCTIONAL TOTAL ENERGY', 'WB97-CHG DISPERSION CORRECTION ENERGY'],
            coeff=[-1, 1, 1]))

    # misc.
    pv0.extend(
        _solve_in_turn(args=['CURRENT ENERGY', 'CURRENT REFERENCE ENERGY', 'CURRENT CORRELATION ENERGY'],
                       coeff=[-1, 1, 1]))

    return pv0


def sapt_psivars():
    """Returns dictionary of QCVariable definitions for SAPT."""

    # yapf: disable
    pv1 = collections.OrderedDict()
    pv1['SAPT EXCHSCAL1'] = {'func': lambda x: 1.0 if x[0] < 1.0e-5 else x[0] / x[1],
                             'args': ['SAPT EXCH10 ENERGY', 'SAPT EXCH10(S^2) ENERGY']}  # special treatment in pandas
    pv1['SAPT EXCHSCAL3'] = {'func': lambda x: x[0] ** 3,
                             'args': ['SAPT EXCHSCAL1']}
    pv1['SAPT EXCHSCAL'] = {'func': lambda x: x[0] ** x[1],
                            'args': ['SAPT EXCHSCAL1', 'SAPT ALPHA']}
    pv1['SAPT HF(2) ALPHA=0.0 ENERGY'] = {'func': lambda x: x[0] - (x[1] + x[2] + x[3] + x[4]),
                                          'args': ['SAPT HF TOTAL ENERGY', 'SAPT ELST10,R ENERGY', 'SAPT EXCH10 ENERGY',
                                                   'SAPT IND20,R ENERGY', 'SAPT EXCH-IND20,R ENERGY']}

    pv1['SAPT HF(2) ENERGY'] = {'func': lambda x: x[1] + (1.0 - x[0]) * x[2],
                                'args': ['SAPT EXCHSCAL', 'SAPT HF(2) ALPHA=0.0 ENERGY', 'SAPT EXCH-IND20,R ENERGY']}
    pv1['SAPT HF(3) ENERGY'] = {'func': lambda x: x[1] - (x[2] + x[0] * x[3]),
                                'args': ['SAPT EXCHSCAL', 'SAPT HF(2) ENERGY', 'SAPT IND30,R ENERGY', 'SAPT EXCH-IND30,R ENERGY']}
    pv1['SAPT MP2(2) ENERGY'] = {'func': lambda x: x[1] - (x[2] + x[3] + x[4] + x[0] * (x[5] + x[6] + x[7] + x[8])),
                                 'args': ['SAPT EXCHSCAL', 'SAPT MP2 CORRELATION ENERGY', 'SAPT ELST12,R ENERGY',  # MP2 CORRELATION ENERGY renamed here from pandas since this is IE  # renamed again SA --> SAPT
                                          'SAPT IND22 ENERGY', 'SAPT DISP20 ENERGY', 'SAPT EXCH11(S^2) ENERGY',
                                          'SAPT EXCH12(S^2) ENERGY', 'SAPT EXCH-IND22 ENERGY', 'SAPT EXCH-DISP20 ENERGY']}
    pv1['SAPT MP2(3) ENERGY'] = {'func': lambda x: x[1] - (x[2] + x[0] * x[3]),
                                 'args': ['SAPT EXCHSCAL', 'SAPT MP2(2) ENERGY', 'SAPT IND-DISP30 ENERGY', 'SAPT EXCH-IND-DISP30 ENERGY']}
    pv1['SAPT MP4 DISP'] = {'func': lambda x: x[0] * x[1] + x[2] + x[3] + x[4] + x[5],
                            'args': ['SAPT EXCHSCAL', 'SAPT EXCH-DISP20 ENERGY', 'SAPT DISP20 ENERGY',
                                     'SAPT DISP21 ENERGY', 'SAPT DISP22(SDQ) ENERGY', 'SAPT EST.DISP22(T) ENERGY']}
    pv1['SAPT CCD DISP'] = {'func': lambda x: x[0] * x[1] + x[2] + x[3] + x[4],
                            'args': ['SAPT EXCHSCAL', 'SAPT EXCH-DISP20 ENERGY', 'SAPT DISP2(CCD) ENERGY',
                                     'SAPT DISP22(S)(CCD) ENERGY', 'SAPT EST.DISP22(T)(CCD) ENERGY']}
    pv1['SAPT0 ELST ENERGY'] = {'func': sum, 'args': ['SAPT ELST10,R ENERGY']}
    pv1['SAPT0 EXCH ENERGY'] = {'func': sum, 'args': ['SAPT EXCH10 ENERGY']}
    pv1['SAPT0 IND ENERGY'] = {'func': lambda x: x[1] + x[2] + x[0] * x[3],
                                'args': ['SAPT EXCHSCAL', 'SAPT HF(2) ENERGY', 'SAPT IND20,R ENERGY', 'SAPT EXCH-IND20,R ENERGY']}
    pv1['SAPT0 DISP ENERGY'] = {'func': lambda x: x[0] * x[1] + x[2],
                                'args': ['SAPT EXCHSCAL', 'SAPT EXCH-DISP20 ENERGY', 'SAPT DISP20 ENERGY']}
    pv1['SAPT0 TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT0 ELST ENERGY', 'SAPT0 EXCH ENERGY', 'SAPT0 IND ENERGY', 'SAPT0 DISP ENERGY']}
    pv1['SSAPT0 ELST ENERGY'] = {'func': sum, 'args': ['SAPT0 ELST ENERGY']}
    pv1['SSAPT0 EXCH ENERGY'] = {'func': sum, 'args': ['SAPT0 EXCH ENERGY']}
    pv1['SSAPT0 IND ENERGY'] = {'func': lambda x: x[1] + (x[0] - 1.0) * x[2],
                                 'args': ['SAPT EXCHSCAL3', 'SAPT0 IND ENERGY', 'SAPT EXCH-IND20,R ENERGY']}
    pv1['SSAPT0 DISP ENERGY'] = {'func': lambda x: x[0] * x[1] + x[2],
                                 'args': ['SAPT EXCHSCAL3', 'SAPT EXCH-DISP20 ENERGY', 'SAPT DISP20 ENERGY']}
    pv1['SSAPT0 TOTAL ENERGY'] = {'func': sum, 'args': ['SSAPT0 ELST ENERGY', 'SSAPT0 EXCH ENERGY', 'SSAPT0 IND ENERGY', 'SSAPT0 DISP ENERGY']}
    pv1['SCS-SAPT0 ELST ENERGY'] = {'func': sum, 'args': ['SAPT0 ELST ENERGY']}
    pv1['SCS-SAPT0 EXCH ENERGY'] = {'func': sum, 'args': ['SAPT0 EXCH ENERGY']}
    pv1['SCS-SAPT0 IND ENERGY'] = {'func': sum, 'args': ['SAPT0 IND ENERGY']}
    pv1['SCS-SAPT0 DISP ENERGY'] = {'func': lambda x: (x[0] - x[3]) * (x[1] + x[2]) + x[3] * (x[4] + x[5]),
                                    'args': [0.66, 'SAPT SAME-SPIN EXCH-DISP20 ENERGY', 'SAPT SAME-SPIN DISP20 ENERGY',
                                             1.2, 'SAPT EXCH-DISP20 ENERGY', 'SAPT DISP20 ENERGY']}  # note no xs for SCS disp
    pv1['SCS-SAPT0 TOTAL ENERGY'] = {'func': sum, 'args': ['SCS-SAPT0 ELST ENERGY', 'SCS-SAPT0 EXCH ENERGY', 'SCS-SAPT0 IND ENERGY', 'SCS-SAPT0 DISP ENERGY']}
    pv1['SAPT2 ELST ENERGY'] = {'func': sum, 'args': ['SAPT ELST10,R ENERGY', 'SAPT ELST12,R ENERGY']}
    pv1['SAPT2 EXCH ENERGY'] = {'func': lambda x: x[1] + x[0] * (x[2] + x[3]),
                                'args': ['SAPT EXCHSCAL', 'SAPT EXCH10 ENERGY', 'SAPT EXCH11(S^2) ENERGY', 'SAPT EXCH12(S^2) ENERGY']}
    pv1['SAPT2 IND ENERGY'] = {'func': lambda x: x[1] + x[2] + x[0] * x[3] + x[4] + x[0] * x[5],
                                'args': ['SAPT EXCHSCAL', 'SAPT HF(2) ENERGY', 'SAPT IND20,R ENERGY', 'SAPT EXCH-IND20,R ENERGY',
                                         'SAPT IND22 ENERGY', 'SAPT EXCH-IND22 ENERGY']}
    pv1['SAPT2 DISP ENERGY'] = {'func': lambda x: x[0] * x[1] + x[2],
                                'args': ['SAPT EXCHSCAL', 'SAPT EXCH-DISP20 ENERGY', 'SAPT DISP20 ENERGY']}
    pv1['SAPT2 TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT2 ELST ENERGY', 'SAPT2 EXCH ENERGY', 'SAPT2 IND ENERGY', 'SAPT2 DISP ENERGY']}
    pv1['SAPT2+ ELST ENERGY'] = {'func': sum, 'args': ['SAPT ELST10,R ENERGY', 'SAPT ELST12,R ENERGY']}
    pv1['SAPT2+ EXCH ENERGY'] = {'func': lambda x: x[1] + x[0] * (x[2] + x[3]),
                                 'args': ['SAPT EXCHSCAL', 'SAPT EXCH10 ENERGY', 'SAPT EXCH11(S^2) ENERGY', 'SAPT EXCH12(S^2) ENERGY']}
    pv1['SAPT2+ IND ENERGY'] = {'func': lambda x: x[1] + x[2] + x[0] * x[3] + x[4] + x[0] * x[5],
                                 'args': ['SAPT EXCHSCAL', 'SAPT HF(2) ENERGY', 'SAPT IND20,R ENERGY', 'SAPT EXCH-IND20,R ENERGY',
                                          'SAPT IND22 ENERGY', 'SAPT EXCH-IND22 ENERGY']}
    pv1['SAPT2+ DISP ENERGY'] = {'func': sum, 'args': ['SAPT MP4 DISP']}
    pv1['SAPT2+ TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT2+ ELST ENERGY', 'SAPT2+ EXCH ENERGY', 'SAPT2+ IND ENERGY', 'SAPT2+ DISP ENERGY']}
    pv1['SAPT2+(CCD) ELST ENERGY'] = {'func': sum, 'args': ['SAPT2+ ELST ENERGY']}
    pv1['SAPT2+(CCD) EXCH ENERGY'] = {'func': sum, 'args': ['SAPT2+ EXCH ENERGY']}
    pv1['SAPT2+(CCD) IND ENERGY'] = {'func': sum, 'args': ['SAPT2+ IND ENERGY']}
    pv1['SAPT2+(CCD) DISP ENERGY'] = {'func': sum, 'args': ['SAPT CCD DISP']}
    pv1['SAPT2+(CCD) TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT2+(CCD) ELST ENERGY', 'SAPT2+(CCD) EXCH ENERGY', 'SAPT2+(CCD) IND ENERGY', 'SAPT2+(CCD) DISP ENERGY']}
    pv1['SAPT2+DMP2 ELST ENERGY'] = {'func': sum, 'args': ['SAPT2+ ELST ENERGY']}
    pv1['SAPT2+DMP2 EXCH ENERGY'] = {'func': sum, 'args': ['SAPT2+ EXCH ENERGY']}
    pv1['SAPT2+DMP2 IND ENERGY'] = {'func': sum, 'args': ['SAPT2+ IND ENERGY', 'SAPT MP2(2) ENERGY']}
    pv1['SAPT2+DMP2 DISP ENERGY'] = {'func': sum, 'args': ['SAPT2+ DISP ENERGY']}
    pv1['SAPT2+DMP2 TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT2+DMP2 ELST ENERGY', 'SAPT2+DMP2 EXCH ENERGY', 'SAPT2+DMP2 IND ENERGY', 'SAPT2+DMP2 DISP ENERGY']}
    pv1['SAPT2+(CCD)DMP2 ELST ENERGY'] = {'func': sum, 'args': ['SAPT2+ ELST ENERGY']}
    pv1['SAPT2+(CCD)DMP2 EXCH ENERGY'] = {'func': sum, 'args': ['SAPT2+ EXCH ENERGY']}
    pv1['SAPT2+(CCD)DMP2 IND ENERGY'] = {'func': sum, 'args': ['SAPT2+DMP2 IND ENERGY']}
    pv1['SAPT2+(CCD)DMP2 DISP ENERGY'] = {'func': sum, 'args': ['SAPT2+(CCD) DISP ENERGY']}
    pv1['SAPT2+(CCD)DMP2 TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT2+(CCD)DMP2 ELST ENERGY', 'SAPT2+(CCD)DMP2 EXCH ENERGY', 'SAPT2+(CCD)DMP2 IND ENERGY', 'SAPT2+(CCD)DMP2 DISP ENERGY']}
    pv1['SAPT2+(3) ELST ENERGY'] = {'func': sum, 'args': ['SAPT ELST10,R ENERGY', 'SAPT ELST12,R ENERGY', 'SAPT ELST13,R ENERGY']}
    pv1['SAPT2+(3) EXCH ENERGY'] = {'func': lambda x: x[1] + x[0] * (x[2] + x[3]),
                                    'args': ['SAPT EXCHSCAL', 'SAPT EXCH10 ENERGY', 'SAPT EXCH11(S^2) ENERGY', 'SAPT EXCH12(S^2) ENERGY']}
    pv1['SAPT2+(3) IND ENERGY'] = {'func': lambda x: x[1] + x[2] + x[0] * x[3] + x[4] + x[0] * x[5],
                                    'args': ['SAPT EXCHSCAL', 'SAPT HF(2) ENERGY', 'SAPT IND20,R ENERGY', 'SAPT EXCH-IND20,R ENERGY',
                                             'SAPT IND22 ENERGY', 'SAPT EXCH-IND22 ENERGY']}
    pv1['SAPT2+(3) DISP ENERGY'] = {'func': sum, 'args': ['SAPT MP4 DISP', 'SAPT DISP30 ENERGY']}
    pv1['SAPT2+(3) TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT2+(3) ELST ENERGY', 'SAPT2+(3) EXCH ENERGY', 'SAPT2+(3) IND ENERGY', 'SAPT2+(3) DISP ENERGY']}
    pv1['SAPT2+(3)(CCD) ELST ENERGY'] = {'func': sum, 'args': ['SAPT2+(3) ELST ENERGY']}
    pv1['SAPT2+(3)(CCD) EXCH ENERGY'] = {'func': sum, 'args': ['SAPT2+(3) EXCH ENERGY']}
    pv1['SAPT2+(3)(CCD) IND ENERGY'] = {'func': sum, 'args': ['SAPT2+(3) IND ENERGY']}
    pv1['SAPT2+(3)(CCD) DISP ENERGY'] = {'func': sum, 'args': ['SAPT CCD DISP', 'SAPT DISP30 ENERGY']}
    pv1['SAPT2+(3)(CCD) TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT2+(3)(CCD) ELST ENERGY', 'SAPT2+(3)(CCD) EXCH ENERGY', 'SAPT2+(3)(CCD) IND ENERGY', 'SAPT2+(3)(CCD) DISP ENERGY']}
    pv1['SAPT2+(3)DMP2 ELST ENERGY'] = {'func': sum, 'args': ['SAPT2+(3) ELST ENERGY']}
    pv1['SAPT2+(3)DMP2 EXCH ENERGY'] = {'func': sum, 'args': ['SAPT2+(3) EXCH ENERGY']}
    pv1['SAPT2+(3)DMP2 IND ENERGY'] = {'func': sum, 'args': ['SAPT2+(3) IND ENERGY', 'SAPT MP2(2) ENERGY']}
    pv1['SAPT2+(3)DMP2 DISP ENERGY'] = {'func': sum, 'args': ['SAPT2+(3) DISP ENERGY']}
    pv1['SAPT2+(3)DMP2 TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT2+(3)DMP2 ELST ENERGY', 'SAPT2+(3)DMP2 EXCH ENERGY', 'SAPT2+(3)DMP2 IND ENERGY', 'SAPT2+(3)DMP2 DISP ENERGY']}
    pv1['SAPT2+(3)(CCD)DMP2 ELST ENERGY'] = {'func': sum, 'args': ['SAPT2+(3) ELST ENERGY']}
    pv1['SAPT2+(3)(CCD)DMP2 EXCH ENERGY'] = {'func': sum, 'args': ['SAPT2+(3) EXCH ENERGY']}
    pv1['SAPT2+(3)(CCD)DMP2 IND ENERGY'] = {'func': sum, 'args': ['SAPT2+(3)DMP2 IND ENERGY']}
    pv1['SAPT2+(3)(CCD)DMP2 DISP ENERGY'] = {'func': sum, 'args': ['SAPT2+(3)(CCD) DISP ENERGY']}
    pv1['SAPT2+(3)(CCD)DMP2 TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT2+(3)(CCD)DMP2 ELST ENERGY', 'SAPT2+(3)(CCD)DMP2 EXCH ENERGY', 'SAPT2+(3)(CCD)DMP2 IND ENERGY', 'SAPT2+(3)(CCD)DMP2 DISP ENERGY']}
    pv1['SAPT2+3 ELST ENERGY'] = {'func': sum, 'args': ['SAPT ELST10,R ENERGY', 'SAPT ELST12,R ENERGY', 'SAPT ELST13,R ENERGY']}
    pv1['SAPT2+3 EXCH ENERGY'] = {'func': lambda x: x[1] + x[0] * (x[2] + x[3]),
                                  'args': ['SAPT EXCHSCAL', 'SAPT EXCH10 ENERGY', 'SAPT EXCH11(S^2) ENERGY', 'SAPT EXCH12(S^2) ENERGY']}
    pv1['SAPT2+3 IND ENERGY'] = {'func': lambda x: x[1] + x[2] + x[0] * x[3] + x[4] + x[0] * x[5] + x[6] + x[0] * x[7],
                                  'args': ['SAPT EXCHSCAL', 'SAPT HF(3) ENERGY', 'SAPT IND20,R ENERGY', 'SAPT EXCH-IND20,R ENERGY',
                                           'SAPT IND22 ENERGY', 'SAPT EXCH-IND22 ENERGY', 'SAPT IND30,R ENERGY', 'SAPT EXCH-IND30,R ENERGY']}
    pv1['SAPT2+3 DISP ENERGY'] = {'func': lambda x: x[1] + x[2] + x[0] * x[3] + x[4] + x[0] * x[5],
                                  'args': ['SAPT EXCHSCAL', 'SAPT MP4 DISP', 'SAPT DISP30 ENERGY', 'SAPT EXCH-DISP30 ENERGY',
                                           'SAPT IND-DISP30 ENERGY', 'SAPT EXCH-IND-DISP30 ENERGY']}
    pv1['SAPT2+3 TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT2+3 ELST ENERGY', 'SAPT2+3 EXCH ENERGY', 'SAPT2+3 IND ENERGY', 'SAPT2+3 DISP ENERGY']}
    pv1['SAPT2+3(CCD) ELST ENERGY'] = {'func': sum, 'args': ['SAPT2+3 ELST ENERGY']}
    pv1['SAPT2+3(CCD) EXCH ENERGY'] = {'func': sum, 'args': ['SAPT2+3 EXCH ENERGY']}
    pv1['SAPT2+3(CCD) IND ENERGY'] = {'func': sum, 'args': ['SAPT2+3 IND ENERGY']}
    pv1['SAPT2+3(CCD) DISP ENERGY'] = {'func': lambda x: x[1] + x[2] + x[0] * x[3] + x[4] + x[0] * x[5],
                                       'args': ['SAPT EXCHSCAL', 'SAPT CCD DISP', 'SAPT DISP30 ENERGY', 'SAPT EXCH-DISP30 ENERGY',
                                                'SAPT IND-DISP30 ENERGY', 'SAPT EXCH-IND-DISP30 ENERGY']}
    pv1['SAPT2+3(CCD) TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT2+3(CCD) ELST ENERGY', 'SAPT2+3(CCD) EXCH ENERGY', 'SAPT2+3(CCD) IND ENERGY', 'SAPT2+3(CCD) DISP ENERGY']}
    pv1['SAPT2+3DMP2 ELST ENERGY'] = {'func': sum, 'args': ['SAPT2+3 ELST ENERGY']}
    pv1['SAPT2+3DMP2 EXCH ENERGY'] = {'func': sum, 'args': ['SAPT2+3 EXCH ENERGY']}
    pv1['SAPT2+3DMP2 IND ENERGY'] = {'func': sum, 'args': ['SAPT2+3 IND ENERGY', 'SAPT MP2(3) ENERGY']}
    pv1['SAPT2+3DMP2 DISP ENERGY'] = {'func': sum, 'args': ['SAPT2+3 DISP ENERGY']}
    pv1['SAPT2+3DMP2 TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT2+3DMP2 ELST ENERGY', 'SAPT2+3DMP2 EXCH ENERGY', 'SAPT2+3DMP2 IND ENERGY', 'SAPT2+3DMP2 DISP ENERGY']}
    pv1['SAPT2+3(CCD)DMP2 ELST ENERGY'] = {'func': sum, 'args': ['SAPT2+3 ELST ENERGY']}
    pv1['SAPT2+3(CCD)DMP2 EXCH ENERGY'] = {'func': sum, 'args': ['SAPT2+3 EXCH ENERGY']}
    pv1['SAPT2+3(CCD)DMP2 IND ENERGY'] = {'func': sum, 'args': ['SAPT2+3DMP2 IND ENERGY']}
    pv1['SAPT2+3(CCD)DMP2 DISP ENERGY'] = {'func': sum, 'args': ['SAPT2+3(CCD) DISP ENERGY']}
    pv1['SAPT2+3(CCD)DMP2 TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT2+3(CCD)DMP2 ELST ENERGY', 'SAPT2+3(CCD)DMP2 EXCH ENERGY', 'SAPT2+3(CCD)DMP2 IND ENERGY', 'SAPT2+3(CCD)DMP2 DISP ENERGY']}
    # yapf: enable

    return pv1
