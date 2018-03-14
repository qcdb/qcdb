import copy
import uuid
import collections

from ..exceptions import *
from . import parsers


#def read_options(options):
def load_qcdb_defaults(options):

    options.add('qcdb', RottenOption(
            keyword='translate_qcdb',
            default=True,
            validator=parsers.boolean,
            glossary='assert always translated. should qcdb or program defaults prevail?'))

    options.add('qcdb', RottenOption(  # true global
            keyword='memory',
            default='700 mb',
            validator=parsers.parse_memory,
            glossary='Total memory allocation in bytes.'))

    options.add('qcdb', RottenOption(  # true global
            keyword='basis',
            default='',
            validator=lambda x: x.upper(),
            glossary='Primary orbital basis set.'))

    #options.add('qcdb', RottenOption(
    #        keyword='e_convergence',
    #        default=1.e-6,
    #        validator=parsers.parse_convergence,
    #        glossary='Convergence criterion for energy.'))

    options.add('qcdb', RottenOption(  # derived shorthand global
            keyword='scf__e_convergence',
            default=1.e-6,
            validator=parsers.parse_convergence,
            glossary='Convergence criterion for SCF energy.'))

    options.add('qcdb', RottenOption(  # derived shorthand global
            keyword='scf__d_convergence',
            default=1.e-6,
            validator=parsers.parse_convergence,
            glossary='Convergence criterion for SCF density (psi: which is defined as the RMS value of the orbital gradient.'))

    options.add('qcdb', RottenOption(  # true global
            keyword='puream',
            default=True,
            validator=parsers.sphcart,
            glossary="""Do use pure angular momentum basis functions?
  If not explicitly set, the default comes from the basis set.
  **Cfour Interface:** Keyword translates into |cfour__cfour_spherical|."""))

    options.add('qcdb', RottenOption(
            keyword='reference',  # TODO don't want higher and local
            default='',
            validator=lambda x: x.upper(),  # TODO
            glossary="""Reference wavefunction type.
    **Cfour Interface:** Keyword translates into |cfour__cfour_reference|."""))
    #options.add_str("REFERENCE", "RHF", "RHF ROHF UHF CUHF RKS UKS")

    options.add('qcdb', RottenOption(
            keyword='scf__reference',
            default='',
            validator=lambda x: x.upper(),  # TODO
            glossary="""Reference wavefunction type.
    **Cfour Interface:** Keyword translates into |cfour__cfour_reference|."""))
    #options.add_str("REFERENCE", "RHF", "RHF ROHF UHF CUHF RKS UKS")

    options.add('qcdb', RottenOption(
            keyword='scf_type',  # TODO ditto, 2-leveled
            default='',
            validator=lambda x: x.upper(),  # TODO
            glossary="""What algorithm to use for the SCF computation. See Table :ref:`SCF
    Convergence & Algorithm <table:conv_scf>` for default algorithm for
    different calculation types."""))
    #options.add_str("SCF_TYPE", "PK", "DIRECT DF PK OUT_OF_CORE CD GTFOCK");

    options.add('qcdb', RottenOption(
            keyword='scf__scf_type',
            default='',
            validator=lambda x: x.upper(),  # TODO
            glossary="""What algorithm to use for the SCF computation. See Table :ref:`SCF
    Convergence & Algorithm <table:conv_scf>` for default algorithm for
    different calculation types."""))
    #options.add_str("SCF_TYPE", "PK", "DIRECT DF PK OUT_OF_CORE CD GTFOCK");

    options.add('qcdb', RottenOption(
            keyword='scf__maxiter',
            default=100,
            validator=parsers.positive_integer,
            glossary="""Maximum number of iterations.
    **Cfour Interface:** Keyword translates into |cfour__cfour_scf_maxcyc|."""))

    options.add('qcdb', RottenOption(
            keyword='scf__damping_percentage',
            default=0.0,
            validator=parsers.percentage,
            glossary="""The amount (percentage) of damping to apply to the early density updates.
        0 will result in a full update, 100 will completely stall the update. A
        value around 20 (which corresponds to 20\% of the previous iteration's
        density being mixed into the current density)
        could help to solve problems with oscillatory convergence."""))

    options.add('qcdb', RottenOption(
            keyword='mp2__mp2_type',
            default='',
            validator=lambda x: x.upper(),  # TODO
            glossary="""What algorithm to use for MP2 computation. See 
        :ref:`Cross-module Redundancies <table:managedmethods>` for details."""))
    #options.add_str("MP2_TYPE", "DF", "DF CONV CD");

    options.add('qcdb', RottenOption(
            keyword='mp2_ss_scale',
            default=1.0,
            validator=lambda x: float(x),
            glossary="""MP2 same-spin scaling value. Default produces canonical MP2, not canonical SCS-MP2."""))

    options.add('qcdb', RottenOption(
            keyword='mp2_os_scale',
            default=1.0,
            validator=lambda x: float(x),
            glossary="""MP2 opposite-spin scaling value. Default produces canonical MP2, not canonical SCS-MP2."""))

    options.add('qcdb', RottenOption(
            keyword='writer_file_label',
            default='',
            validator=lambda x: x,
            glossary="""Base filename for text files written by QCDB, such as the
  MOLDEN output file, the Hessian file, the internal coordinate file, etc.
"""))

    options.add('qcdb', RottenOption(
            keyword='geom_maxiter',
            default=20,
            validator=parsers.positive_integer,
            glossary="""Maximum number of geometry optimization steps."""))


    options.add('qcdb', RottenOption(
            keyword='g_convergence',
            default='qchem',
            validator=parsers.enum("QCHEM MOLPRO GAU GAU_LOOSE GAU_TIGHT INTERFRAG_TIGHT GAU_VERYTIGHT TURBOMOLE CFOUR NWCHEM_LOOSE"),
            glossary="""Set of optimization criteria. Specification of any MAX_*_G_CONVERGENCE
      or RMS_*_G_CONVERGENCE options will append to overwrite the criteria set here
      unless |optking__flexible_g_convergence| is also on.      See Table :ref:`Geometry Convergence <table:optkingconv>` for details."""))

    options.add('qcdb', RottenOption(
            keyword='max_force_g_convergence',
            default=3.e-4,
            validator=parsers.parse_convergence,
            glossary="""Convergence criterion for geometry optmization: maximum force
      (internal coordinates, atomic units)."""))

      #/*- Convergence criterion for geometry optmization: rms force
      #(internal coordinates, atomic units). -*/
      #options.add_double("RMS_FORCE_G_CONVERGENCE", 3.0e-4);
      #/*- Convergence criterion for geometry optmization: maximum energy change. -*/
      #options.add_double("MAX_ENERGY_G_CONVERGENCE", 1.0e-6);
      #/*- Convergence criterion for geometry optmization: maximum displacement
      #(internal coordinates, atomic units). -*/
      #options.add_double("MAX_DISP_G_CONVERGENCE", 1.2e-3);
      #/*- Convergence criterion for geometry optmization: rms displacement
      #(internal coordinates, atomic units). -*/
      #options.add_double("RMS_DISP_G_CONVERGENCE", 1.2e-3);
      #/*- Even if a user-defined threshold is set, allow for normal, flexible convergence criteria -*/
      #options.add_bool("FLEXIBLE_G_CONVERGENCE", false);

#    options.add('qcdb', RottenOption(
#            keyword='',
#            default=,
#            validator=,
#            glossary="""."""))

#    options.add('qcdb', RottenOption(
#            keyword='',
#            default=,
#            validator=,
#            glossary="""."""))

#    options.add('qcdb', RottenOption(
#            keyword='',
#            default=,
#            validator=,
#            glossary="""."""))

#    options.add('qcdb', RottenOption(
#            keyword='',
#            default=,
#            validator=,
#            glossary="""."""))

#    options.add('qcdb', RottenOption(
#            keyword='',
#            default=,
#            validator=,
#            glossary="""."""))



class RottenOptions(object):
    mark_of_the_user = '00000000'
    mark_of_the_default = 'ffffffff'

    def __init__(self):
        self.scroll = collections.defaultdict(dict)

    def __str__(self):
        text = []
        for pkg in self.scroll:
            text.append('  <<<  {}  >>>'.format(pkg))
            for opt, oopt in sorted(self.scroll[pkg].items()):
                text.append(str(oopt))

        return '\n'.join(text)

    def print_changed(self):
        text = []
        for pkg in self.scroll:
            text.append('  <<<  {}  >>>'.format(pkg))
            for opt, oopt in sorted(self.scroll[pkg].items()):
                #if not oopt.is_default():
                if len(oopt.history) > 1:
                    text.append(str(oopt))

        return '\n'.join(text)

    def add(self, package, opt):
        up = package.upper()
        if up in ['QCDB', 'PSI4', 'CFOUR', 'DFTD3', 'NWCHEM', 'GAMESS', 'RESP']:
            self.scroll[up][opt.keyword] = opt
        else:
            raise ValidationError('Domain not supported: {}'.format(package))
    
    def require(self, package, option, value, accession, verbose=1):
        self._set(True, package, option, value, accession, verbose)

    def suggest(self, package, option, value, accession, verbose=1):
        self._set(False, package, option, value, accession, verbose)

    def _set(self, imperative, package, option, value, accession, verbose):
        count = 0
        for ropt, oropt in self.scroll[package.upper()].items():
            #if ropt.endswith(option.upper()):
            if ropt == option.upper() or ropt.endswith('__' + option.upper()):  # psi wants
                overlap = len(option)
                if imperative:
                    oropt.require(value, overlap=overlap, accession=accession, verbose=verbose)
                else:
                    oropt.suggest(value, overlap=overlap, accession=accession, verbose=verbose)
                count += 1
        if count == 0:
            raise ValidationError('Option ({}) does not exist in domain ({}).'.format(option, package))

    def unwind_by_accession(self, accession):
        for pkg in self.scroll:
            for ropt, oropt in self.scroll[pkg].items():
               oropt.history = [entry for entry in oropt.history if entry[3] != accession]
            

class RottenOption(object):
    mark_of_the_user = '00000000'
    mark_of_the_default = 'ffffffff'

    def __init__(self, keyword, default, validator, glossary='', expert=False):
        self.keyword = keyword.upper()
        self.glossary = glossary
        self.validator = validator
        self.history = []  # list of quads (value, required, overlap, accession)
        self.suggest(default, accession=self.mark_of_the_default, verbose=0)
        self.has_changed = False
        self.expert = expert

    def __str__(self):
        text = []
        text.append('  {:50} {:>30} {} {}'.format(
                                        self.keyword + ':',
                                        str(self.value),
                                        '  ' if self.is_default() else '<>',
                                        '(' + str(self.history[0][0]) + ')'))
        if self.disputed():
            text.extend(['         ' + str(entry) for entry in self.history])
                #format(self.keyword, added[0], added[2] + 100 * int(added[1]), added[3]))
        return '\n'.join(text)

    #@property
    #def value(self):
    def _compute(self):
        """The all-important `self.value` is read-only and computed on-the-fly from `self.history`."""

        scores = [cand[2] + 100 * int(cand[1]) for cand in self.history]
        max_score = max(scores)

        # only catch user and driver reqd of highest relevance and most recent vintage
        user = None
        for score, candidate in zip(reversed(scores), reversed(self.history)):
            if score == max_score and candidate[3] == self.mark_of_the_user:
                user = candidate
                break

        driver = None
        for score, candidate in zip(reversed(scores), reversed(self.history)):
            if score == max_score and candidate[3] != self.mark_of_the_user:
                driver = candidate
                break

        if user is None and driver is None:
            raise OptionReconciliationError('No info')
        elif user is None and driver is not None:
            hist = driver
        elif user is not None and driver is None:
            hist = user
        elif user is not None and driver is not None:
            if user[0] == driver[0]:
                hist = user
            else:
                raise OptionReconciliationError(
                    'Conflicting option requirements btwn user ({}) and driver ({})'.
                    format(user[0], driver[0]))
        
        #self.score = max_score
        return hist[0], max_score, hist
        
    @property
    def value(self):
        val, score, hist = self._compute()
        return val

    def inherit(self, other, transform, suggests_too=True): #score_cutoff=100):
        assert isinstance(other, RottenOption)
        #oval, oscore, ohist = other._compute()
        _, _, ohist = other._compute()
        oval, oimperative, ooverlap, oaccession = ohist
        
        # if `other` has no required and mere suggestions not wanted, done
        #if not ohist[1] and not suggests_too:
        if not oimperative and not suggests_too:
            return

        #self.set(oimperative, transform(oval), ???, oaccession)

    def suggest(self, value, overlap=None, accession=None, verbose=1):
        self._set(False, value, overlap, accession, verbose)
        
    def require(self, value, overlap=None, accession=None, verbose=1):
        """

        Parameters
        ----------
        value
            Asserted value for option `self`. Will be checked against
            `self.validator` before function returns. May still be incompatible with
            other `require` calls of same priority, but that won't be checked until
            `self.value` is accessed.
        overlap : int, optional
            Specificity of assertion. If `self.keyword='CC_MAXITER'` and `value`
            is set for `MAXITER`, `overlap=7`, whereas if set for `CC_MAXITER`,
            `overlap=10`.
        accession
            Tag of who is setting this option.

        """
        self._set(True, value, overlap, accession, verbose)

    def _set(self, imperative, value, overlap, accession, verbose=1):
        if overlap is None:
            overlap = len(self.keyword)
        if accession is None:
            accession = uuid.uuid4()

        self.history.append((self._check(value), imperative, overlap, accession))

        if verbose >= 1:
            added = self.history[-1]
            # TODO add back somehow
            #print('Setting {} to {} priority {} accession {}'.
            #    format(self.keyword, added[0], added[2] + 100 * int(added[1]), added[3]))

    def _check(self, val):
        """Common function to check `val` against `self.validator` for setting, defaulting, etc."""

        try:
            nuval = self.validator(val)
        except Exception as err:
            raise OptionValidationError(
                'Option ({}) value ({}) does not pass.'.format(self.keyword, val)) from err
        return nuval

    def is_default(self):
        return self.value == self.history[0][0]

    def disputed(self):
        return len(self.history) > 1

    def is_required(self, score_cutoff=100):
        val, score, hist = self._compute()
        return score >= score_cutoff

