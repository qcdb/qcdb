import uuid
import collections

from ..exceptions import KeywordValidationError, KeywordReconciliationError, ValidationError


class Keywords(object):
    mark_of_the_user = '00000000'
    mark_of_the_default = 'ffffffff'

    def __init__(self):
        self.domains = ['QCDB', 'PSI4', 'CFOUR', 'DFTD3', 'NWCHEM', 'GAMESS', 'RESP']
        self.aliases = collections.defaultdict(dict)
        self.scroll = collections.defaultdict(dict)

    def __str__(self):
        text = []
        for pkg in self.scroll:
            text.append('  <<<  {}  >>>'.format(pkg))
            for opt, oopt in sorted(self.scroll[pkg].items()):
                text.append(str(oopt))

        return '\n'.join(text)

    def print_changed(self, history=True):
        text = []
        for pkg in self.scroll:
            text.append('  <<<  {}  >>>'.format(pkg))
            for opt, oopt in sorted(self.scroll[pkg].items()):
                #if not oopt.is_default():
                if len(oopt.history) > 1:
                    if history:
                        text.append(str(oopt))
                    else:
                        text.append(oopt.shortstr())

        return '\n'.join(text)

    def add(self, package, opt):
        up = package.upper()
        if up in self.domains:
            self.scroll[up][opt.keyword] = opt
        else:
            raise ValidationError(f'Domain not supported: {package}')

    def add_alias(self, package, opt):
        up = package.upper()
        if up in self.domains:
            if opt.alias in self.scroll[up]:
                raise ValidationError(f'Keyword alias must not share a name with keyword proper: {opt.alias}')
            if opt.target not in self.scroll[up]:
                raise ValidationError(
                    f'Keyword alias must point to existing keyword proper: {opt.alias} --/--> {opt.target}')
            self.aliases[up][opt.alias] = opt
        else:
            raise ValidationError(f'Domain not supported: {package}')

    def require(self, package, option, value, accession, verbose=1):
        self._set(True, package, option, value, accession, verbose)

    def suggest(self, package, option, value, accession, verbose=1):
        self._set(False, package, option, value, accession, verbose)

    def _set(self, imperative, package, option, value, accession, verbose):
        count = 0
        acount = 0
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
            for aopt, oaopt in self.aliases[package.upper()].items():
                if aopt == option.upper():
                    oropt = self.scroll[package.upper()][oaopt.target]
                    overlap = len(oropt.keyword)
                    if imperative:
                        oropt.require(value, overlap=overlap, accession=accession, verbose=verbose)
                    else:
                        oropt.suggest(value, overlap=overlap, accession=accession, verbose=verbose)
                    acount += 1

        if count == 0 and acount == 0:
            raise ValidationError('Keyword ({}) does not exist in domain ({}).'.format(option, package))

    def unwind_by_accession(self, accession):
        for pkg in self.scroll:
            for ropt, oropt in self.scroll[pkg].items():
                oropt.history = [entry for entry in oropt.history if entry[3] != accession]


class AliasKeyword(object):
    def __init__(self, alias, target):
        self.alias = alias.upper()
        self.target = target.upper()


class Keyword(object):
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
        text.append('  {:50} {:>30} {} {}'.format(self.keyword + ':', str(self.value),
                                                  '  ' if self.is_default() else '<>',
                                                  '(' + str(self.history[0][0]) + ')'))
        if self.disputed():
            text.extend(['         ' + str(entry) for entry in self.history])
        return '\n'.join(text)

    def shortstr(self):
        text = []
        text.append('  {:50} {:>30} {} {}'.format(self.keyword + ':', str(self.value),
                                                  '  ' if self.is_default() else '<>',
                                                  '(' + str(self.history[0][0]) + ')'))
        return '\n'.join(text)

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
            raise KeywordReconciliationError('No info')
        elif user is None and driver is not None:
            hist = driver
        elif user is not None and driver is None:
            hist = user
        elif user is not None and driver is not None:
            if user[0] == driver[0]:
                hist = user
            else:
                raise KeywordReconciliationError(
                    f'Conflicting option requirements btwn user ({user[0]}) and driver ({driver[0]}) for {self.keyword}'
                )

        #self.score = max_score
        return hist[0], max_score, hist

    @property
    def value(self):
        val, score, hist = self._compute()
        return val

    def inherit(self, other, transform, suggests_too=True):  #score_cutoff=100):
        assert isinstance(other, KeywordOption)
        _, _, ohist = other._compute()
        oval, oimperative, ooverlap, oaccession = ohist

        # if `other` has no required and mere suggestions not wanted, done
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

        if verbose >= 2:
            added = self.history[-1]
            print(f'Setting {self.keyword} to {added[0]} priority {added[2] + 100 * int(added[1])} accession {added[3]}')

    def _check(self, val):
        """Common function to check `val` against `self.validator` for setting, defaulting, etc."""

        try:
            nuval = self.validator(val)
        except Exception as err:
            raise KeywordValidationError('Keyword ({}) value ({}) does not pass.'.format(self.keyword, val)) from err
        return nuval

    def is_default(self):
        return self.value == self.history[0][0]

    def disputed(self):
        return len(self.history) > 1

    def is_required(self, score_cutoff=100):
        val, score, hist = self._compute()
        return score >= score_cutoff
