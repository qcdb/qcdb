import uuid
from collections import defaultdict
from typing import Any, Optional

from ..exceptions import KeywordReconciliationError, KeywordValidationError, ValidationError


class Keywords:
    mark_of_the_user = "00000000"
    mark_of_the_default = "ffffffff"

    def __init__(self):
        self.domains = ["QCDB", "PSI4", "CFOUR", "DFTD3", "NWCHEM", "GAMESS", "RESP"]
        self.aliases = defaultdict(dict)
        self.scroll = defaultdict(dict)

    def __str__(self) -> str:
        text = []
        for pkg in self.scroll:
            text.append(f"  <<<  {pkg}  >>>")
            for key, okey in sorted(self.scroll[pkg].items()):
                text.append(str(okey))

        return "\n".join(text)

    def print_changed(self, history: bool = True) -> str:
        """Compose a printout of default & current values for all packages & all keywords.
        If ``history=True``, includes intermediate values, too."""

        text = []
        for pkg in self.scroll:
            text.append(f"  <<<  {pkg}  >>>")
            for key, okey in sorted(self.scroll[pkg].items()):
                # if not okey.is_default():
                if len(okey.history) > 1:
                    if history:
                        text.append(str(okey))
                    else:
                        text.append(okey.shortstr())

        return "\n".join(text)

    def add(self, package: str, key: "Keyword") -> None:
        """Register single new keyword ``key`` of domain ``package`` into the keyword set."""

        pkg = package.upper()
        if pkg in self.domains:
            self.scroll[pkg][key.keyword] = key
        else:
            raise ValidationError(f"Domain not supported: {package}")

    def add_alias(self, package: str, key: "AliasKeyword") -> None:
        """Register single new alias keyword ``key`` of domain ``package`` into the keyword set."""
        pkg = package.upper()
        if pkg in self.domains:
            if key.alias in self.scroll[pkg]:
                raise ValidationError(f"Keyword alias must not share a name with keyword proper: {key.alias}")
            if key.target not in self.scroll[pkg]:
                raise ValidationError(
                    f"Keyword alias must point to existing keyword proper: {key.alias} --/--> {key.target}")
            self.aliases[pkg][key.alias] = key
        else:
            raise ValidationError(f"Domain not supported: {package}")

    def require(self, package: str, key: str, value: str, accession: str, verbose: int = 1) -> None:
        """Propose that Keyword with name ``key`` in domain ``package`` should have ``value``.
        Used by driver or harnesses to translate specifications into package options.
        Use this function to insist upon ``value`` such that contrary should raise an error.
        See ``sugguest`` for lighter touch.

        Parameters
        ----------
        package
            Domain (QC program) containing keyword for which to propose a value.
        key
            Keyword name for which to propose a value.
        value
            Proposed keyword value. Will be checked against the validator before function returns.
            May still be incompatible with other ``require`` calls of the same priority, but that
            won't be checked until `self.scroll["package"]["key"].value` is accessed.
        accession
            Tag of who (function, person, etc.) is setting this keyword.
        verbose
            If >1, print line upon being added to history.

        """
        self._set(True, package, key, value, accession, verbose)

    def suggest(self, package: str, key: str, value: str, accession: str, verbose: int = 1) -> None:
        """Propose that Keyword with name ``key`` in domain ``package`` should have ``value``.
        Used by driver or harnesses to translate specifications into package options.
        Use this function to urge ``value`` but contrary should not raise an error.
        See ``require`` for stronger touch.

        """
        self._set(False, package, key, value, accession, verbose)

    def _set(self, imperative: bool, package: str, key: str, value: Any, accession: str, verbose: int) -> None:
        """Backend to ``require()`` and ``suggest()`` to add to the history."""
        pkg = package.upper()
        ukey = key.upper()
        count = 0
        acount = 0
        for ropt, oropt in self.scroll[pkg].items():
            #if ropt.endswith(ukey):
            if ropt == ukey or ropt.endswith("__" + ukey):  # psi wants
                overlap = len(key)
                if imperative:
                    oropt.require(value, overlap=overlap, accession=accession, verbose=verbose)
                else:
                    oropt.suggest(value, overlap=overlap, accession=accession, verbose=verbose)
                count += 1
        if count == 0:
            for aopt, oaopt in self.aliases[pkg].items():
                if aopt == ukey:
                    oropt = self.scroll[pkg][oaopt.target]
                    overlap = len(oropt.keyword)
                    if imperative:
                        oropt.require(value, overlap=overlap, accession=accession, verbose=verbose)
                    else:
                        oropt.suggest(value, overlap=overlap, accession=accession, verbose=verbose)
                    acount += 1

        if count == 0 and acount == 0:
            raise ValidationError(f"Keyword ({key}) does not exist in domain ({package}).")

    def unwind_by_accession(self, accession):
        for pkg in self.scroll:
            for rkey, orkey in self.scroll[pkg].items():
                orkey.history = [entry for entry in orkey.history if entry[3] != accession]


class AliasKeyword:
    def __init__(self, alias, target):
        self.alias = alias.upper()
        self.target = target.upper()


class Keyword:
    mark_of_the_user = "00000000"
    mark_of_the_default = "ffffffff"

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

    def require(self, value: Any, overlap: Optional[int] = None, accession: Optional[str] = None, verbose: int = 1) -> None:
        """

        Parameters
        ----------
        value
            Asserted value for keyword `self`. Will be checked against
            `self.validator` before function returns. May still be incompatible with
            other `require` calls of same priority, but that won't be checked until
            `self.value` is accessed.
        overlap
            Specificity of assertion. If `self.keyword='CC_MAXITER'` and `value`
            is set for `MAXITER`, `overlap=7`, whereas if set for `CC_MAXITER`,
            `overlap=10`.
        accession
            Tag of who is setting this option.

        """
        self._set(True, value, overlap, accession, verbose)

    def _set(self, imperative: bool, value: Any, overlap: Optional[int], accession: Optional[str], verbose: int = 1) -> None:
        """Backend to ``suggest()`` and ``require()``. Fills in ``overlap`` and ``accession`` defaults, validates
        ``value``, and appends entry to ``history``. Prints if ``verbose``."""

        if overlap is None:
            overlap = len(self.keyword)
        if accession is None:
            accession = uuid.uuid4()

        self.history.append((self._check(value), imperative, overlap, accession))

        if verbose >= 2:
            added = self.history[-1]
            print(
                f"Setting {self.keyword} to {added[0]} priority {added[2] + 100 * int(added[1])} accession {added[3]}")

    def _check(self, val: Any) -> Any:
        """Common function to check `val` against `self.validator` for setting, defaulting, etc."""

        try:
            nuval = self.validator(val)
        except Exception as err:
            raise KeywordValidationError(f"Keyword ({self.keyword}) value ({val}) does not pass.") from err
        return nuval

    def is_default(self) -> bool:
        """Whether the present evaluated value matches the initial value (set at definition time, so default)"""
        return self.value == self.history[0][0]

    def disputed(self) -> bool:
        """Whether value candidates other than the initial value have been proposed."""
        return len(self.history) > 1

    def is_required(self, score_cutoff: int = 100) -> bool:
        """Whether the present evaluated value was suggested or required (approximately)."""
        val, score, hist = self._compute()
        return score >= score_cutoff
