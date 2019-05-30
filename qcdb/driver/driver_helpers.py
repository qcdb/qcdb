import os
import re

import numpy as np

from ..exceptions import *
from ..molecule import Molecule
from ..util import import_ignorecase
#from .driver import options
from . import pe


def _set_convergence_criterion(ptype, method_name, scf_Ec, pscf_Ec, scf_Dc, pscf_Dc, gen_Ec, verbose=1):
    r"""
    This function will set local SCF and global energy convergence criterion
    to the defaults listed at:
    http://www.psicode.org/psi4manual/master/scf.html#convergence-and-
    algorithm-defaults. SCF will be converged more tightly if a post-SCF
    method is select (pscf_Ec, and pscf_Dc) else the looser (scf_Ec, and
    scf_Dc convergence criterion will be used).

    Parameters
    ----------
    ptype : str
        Procedure type (energy, gradient, etc). Nearly always test on
        procedures['energy'] since that's guaranteed to exist for a method.
    method_name : str
        Name of the method.
    scf_Ec : int or float
        E convergence criterion for scf target method.
    pscf_Ec : int or float
        E convergence criterion for scf of post-scf target method
    scf_Dc : int or float
        D convergence criterion for scf target method
    pscf_Dc : int or float
        D convergence criterion for scf of post-scf target method
    gen_Ec : int or float
        E convergence criterion for post-scf target method

    """
    optstash = p4util.OptionsState(
        ['SCF', 'E_CONVERGENCE'],
        ['SCF', 'D_CONVERGENCE'],
        ['E_CONVERGENCE'])

    # Kind of want to move this out of here
    _method_exists(ptype, method_name)

    if verbose >= 2:
        print('      Setting convergence', end=' ')
    # Set method-dependent scf convergence criteria, check against energy routines
    if not core.has_option_changed('SCF', 'E_CONVERGENCE'):
        if procedures['energy'][method_name] == proc.run_scf:
            core.set_local_option('SCF', 'E_CONVERGENCE', scf_Ec)
            if verbose >= 2:
                print(scf_Ec, end=' ')
        else:
            core.set_local_option('SCF', 'E_CONVERGENCE', pscf_Ec)
            if verbose >= 2:
                print(pscf_Ec, end=' ')
    else:
        if verbose >= 2:
            print('CUSTOM', core.get_option('SCF', 'E_CONVERGENCE'), end=' ')

    if not core.has_option_changed('SCF', 'D_CONVERGENCE'):
        if procedures['energy'][method_name] == proc.run_scf:
            core.set_local_option('SCF', 'D_CONVERGENCE', scf_Dc)
            if verbose >= 2:
                print(scf_Dc, end=' ')
        else:
            core.set_local_option('SCF', 'D_CONVERGENCE', pscf_Dc)
            if verbose >= 2:
                print(pscf_Dc, end=' ')
    else:
        if verbose >= 2:
            print('CUSTOM', core.get_option('SCF', 'D_CONVERGENCE'), end=' ')

    # Set post-scf convergence criteria (global will cover all correlated modules)
    if not core.has_global_option_changed('E_CONVERGENCE'):
        if procedures['energy'][method_name] != proc.run_scf:
            core.set_global_option('E_CONVERGENCE', gen_Ec)
            if verbose >= 2:
                print(gen_Ec, end=' ')
    else:
        if procedures['energy'][method_name] != proc.run_scf:
            if verbose >= 2:
                print('CUSTOM', core.get_global_option('E_CONVERGENCE'), end=' ')

    if verbose >= 2:
        print('')
    return optstash




def _parse_arbitrary_order(name):
    r"""Function to parse name string into a method family like CI or MRCC and specific
    level information like 4 for CISDTQ or MRCCSDTQ.

    Separates `name` into (name, None) mostly, but into ('mp', 3) for `name='mp3'` and `('mrcc', dict)` for `name='mrccsd')

    """
    name = name.lower()
    mtdlvl_mobj = re.match(r"""\A(?P<method>[a-z]+)(?P<level>\d+)\Z""", name)

    # matches 'mrccsdt(q)'
    if name.startswith('mrcc'):

        # avoid undoing fn's good work when called twice
        if name == 'mrcc':
            return name, None

        # grabs 'sdt(q)'
        ccfullname = name[4:]

        # A negative order indicates perturbative method
        methods = {
            'sd'          : { 'method': 1, 'order':  2, 'fullname': 'CCSD'         },
            'sdt'         : { 'method': 1, 'order':  3, 'fullname': 'CCSDT'        },
            'sdtq'        : { 'method': 1, 'order':  4, 'fullname': 'CCSDTQ'       },
            'sdtqp'       : { 'method': 1, 'order':  5, 'fullname': 'CCSDTQP'      },
            'sdtqph'      : { 'method': 1, 'order':  6, 'fullname': 'CCSDTQPH'     },
            'sd(t)'       : { 'method': 3, 'order': -3, 'fullname': 'CCSD(T)'      },
            'sdt(q)'      : { 'method': 3, 'order': -4, 'fullname': 'CCSDT(Q)'     },
            'sdtq(p)'     : { 'method': 3, 'order': -5, 'fullname': 'CCSDTQ(P)'    },
            'sdtqp(h)'    : { 'method': 3, 'order': -6, 'fullname': 'CCSDTQP(H)'   },
            'sd(t)_l'     : { 'method': 4, 'order': -3, 'fullname': 'CCSD(T)_L'    },
            'sdt(q)_l'    : { 'method': 4, 'order': -4, 'fullname': 'CCSDT(Q)_L'   },
            'sdtq(p)_l'   : { 'method': 4, 'order': -5, 'fullname': 'CCSDTQ(P)_L'  },
            'sdtqp(h)_l'  : { 'method': 4, 'order': -6, 'fullname': 'CCSDTQP(H)_L' },
            'sdt-1a'      : { 'method': 5, 'order':  3, 'fullname': 'CCSDT-1a'     },
            'sdtq-1a'     : { 'method': 5, 'order':  4, 'fullname': 'CCSDTQ-1a'    },
            'sdtqp-1a'    : { 'method': 5, 'order':  5, 'fullname': 'CCSDTQP-1a'   },
            'sdtqph-1a'   : { 'method': 5, 'order':  6, 'fullname': 'CCSDTQPH-1a'  },
            'sdt-1b'      : { 'method': 6, 'order':  3, 'fullname': 'CCSDT-1b'     },
            'sdtq-1b'     : { 'method': 6, 'order':  4, 'fullname': 'CCSDTQ-1b'    },
            'sdtqp-1b'    : { 'method': 6, 'order':  5, 'fullname': 'CCSDTQP-1b'   },
            'sdtqph-1b'   : { 'method': 6, 'order':  6, 'fullname': 'CCSDTQPH-1b'  },
            '2'           : { 'method': 7, 'order':  2, 'fullname': 'CC2'          },
            '3'           : { 'method': 7, 'order':  3, 'fullname': 'CC3'          },
            '4'           : { 'method': 7, 'order':  4, 'fullname': 'CC4'          },
            '5'           : { 'method': 7, 'order':  5, 'fullname': 'CC5'          },
            '6'           : { 'method': 7, 'order':  6, 'fullname': 'CC6'          },
            'sdt-3'       : { 'method': 8, 'order':  3, 'fullname': 'CCSDT-3'      },
            'sdtq-3'      : { 'method': 8, 'order':  4, 'fullname': 'CCSDTQ-3'     },
            'sdtqp-3'     : { 'method': 8, 'order':  5, 'fullname': 'CCSDTQP-3'    },
            'sdtqph-3'    : { 'method': 8, 'order':  6, 'fullname': 'CCSDTQPH-3'   }
        }  # yapf: disable

        # looks for 'sdt(q)' in dictionary
        if ccfullname in methods:
            return 'mrcc', methods[ccfullname]
        else:
            raise ValidationError('Invalid MRCC method ({})'.format(name))

    elif mtdlvl_mobj:
        namestump = mtdlvl_mobj.group('method')
        namelevel = int(mtdlvl_mobj.group('level'))

        if namestump in ['mp', 'zapt', 'ci']:
            # let mp2, mp3, mp4 pass through to select functions
            if namestump == 'mp' and namelevel in [2, 3, 4]:
                return name, None
            # otherwise return method and order
            else:
                return namestump, namelevel
        else:
            return name, None
    else:
        return name, None


def set_molecule(molinit, name='default'):

    if molinit.startswith('db:'):
        db, rxn = molinit[3:].strip().split('-', 1)

        libraryPath = os.sep.join([pe.data_dir, 'databases'])
        dbPath = os.path.abspath('.') + \
                 ':' + ':'.join([os.path.abspath(x) for x in os.environ.get('PSIPATH', '').split(':')]) + \
                 ':' + libraryPath

        dbmod = import_ignorecase(db, lenv=dbPath.split(':'))
        if dbmod is None:
            raise ImportError('Python module loading problem for database ({}): {}'.
                format(db, dbPath))

        molecule = dbmod.GEOS[dbmod.dbse + '-' + rxn]
        # let KeyError on fail

        molecule.update_geometry()
        pe.active_molecule = molecule
        return molecule

    molecule = Molecule(molinit)
    pe.active_molecule = molecule
    return molecule


def activate(mol):
    pe.active_molecule = mol


def set_options(options_dict):
    """
    Sets Psi4 global options from an input dictionary.
    """
    optionre = re.compile(r'\A((?P<silo>(cfour|psi4|nwchem|gamess|dftd3|resp))_)?(?P<module>\w+__)?(?P<option>\w+)\Z', re.IGNORECASE)

    if len(pe.nu_options.scroll) == 0:
        #print('EMPTY OPT')
        pe.load_nu_options()

    for k, v, in options_dict.items():
        mobj = optionre.match(k.strip())
        try:
            v = v.strip()
        except AttributeError:
            pass

        if mobj:
            silo = mobj.group('silo').upper() if mobj.group('silo') else 'QCDB'
            module = mobj.group('module').upper() if mobj.group('module') else ''
            option = mobj.group('option').upper()

            print('SET_OPTIONS: [{}][{}] = {}'.format(silo, module + option, v)) 
            pe.nu_options.require(silo, module + option, v, accession=pe.nu_options.mark_of_the_user)
        else:
            raise ValidationError('Option not in {space}?_{module}?__{option} format: {}'.format(k))

    #print('after SET_OPT')
    #print(pe.nu_options)


def old_set_options(options_dict):
    """
    Sets Psi4 global options from an input dictionary.
    """
    optionre = re.compile(r'\A((?P<silo>(cfour|psi))_)?(?P<module>\w+__)?(?P<option>\w+)\Z', re.IGNORECASE)

    for k, v, in options_dict.items():
        mobj = optionre.match(k.strip())
        try:
            v = v.strip()
        except AttributeError:
            pass

        if mobj:
            silo = mobj.group('silo').upper() if mobj.group('silo') else 'GLOBALS'
            module = mobj.group('module').upper() if mobj.group('module') else ''
            option = mobj.group('option').upper()

            pe.active_options[silo][module + option]['value'] = v
            print('SET_OPTIONS: [{}][{}] = {}'.format(silo, module + option, v)) 
        else:
            raise ValidationError('Option not in {space}?_{module}?__{option} format: {}'.format(k))
            

def get_variable(key):
    ukey = key.upper()
    if ukey in pe.active_qcvars:
        return pe.active_qcvars[ukey].data
    else:
        # TODO this matches psi, as None confuses compare_values, but is it right soln?
        #return 0.
        raise ValidationError('No such var as {}'.format(key))


def get_active_molecule():
    return pe.active_molecule


def get_active_options():
    return pe.nu_options


def print_variables(qcvars=None):
    """Form a printable representation of qcvariables.

    Parameters
    ----------
    qcvars : dict of QCAspect, optional
        Group of QCAspect objects to print. If `None`, will use `qcdb.pe.active_qcvars`.

    Returns
    -------
    str
        Printable string representation of label, data, and unit in QCAspect-s.

    """
    text = []
    text.append('\n  Variable Map:')
    text.append('  ----------------------------------------------------------------------------')

    if qcvars is None:
        qcvars = pe.active_qcvars

    if len(qcvars) == 0:
        text.append('  (none)')
        return '\n'.join(text)

    largest_key = max(len(k) for k in qcvars) + 2  # for quotation marks
    for k, qca in sorted(qcvars.items()):
        if k != qca.label:
            raise ValidationError(f'Huh? {k} != {qca.label}')

        if isinstance(qca.data, np.ndarray):
            data = np.array_str(qca.data, max_line_width=120, precision=8, suppress_small=True)
            data = '\n'.join('        ' + ln for ln in data.splitlines())
            text.append("""  {:{keywidth}} => {:{width}} [{}]""".
                format('"' + k + '"', '', qca.units, keywidth=largest_key, width=20))
            text.append(data)
        else:
            text.append("""  {:{keywidth}} => {:{width}.{prec}f} [{}]""".
                format('"' + k + '"', qca.data, qca.units, keywidth=largest_key, width=20, prec=12))

    text.append('')
    return '\n'.join(text)


