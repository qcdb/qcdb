import pytest

import qcdb
from qcdb.keywords import AliasKeyword, Keyword, Keywords, parsers, register_kwds


def validator(v):
    if v > 5:
        return v
    else:
        raise qcdb.KeywordValidationError


def test_1():
    subject = Keyword(keyword='opt1', glossary='does stuff', validator=validator, default=6)
    assert subject.value == 6
    #assert subject.has_changed is False
    assert subject.is_default() is True
    assert subject.keyword == 'OPT1'

    #subject.value = 6
    subject.require(6)
    #assert subject.has_changed is True
    assert subject.is_default() is True

def test_2():
    with pytest.raises(qcdb.KeywordValidationError):
        subject = Keyword(keyword='opt2', glossary='does stuff', validator=validator, default=4)


def test_3():
    subject = Keyword(keyword='opt1', glossary='does stuff', validator=validator, default=6)
    #subject.value = 10
    subject.require(10)

    assert subject.value == 10
    #assert subject.has_changed is True
    assert subject.is_default() is False


def test_4_conv():
    subject = Keyword(keyword='opt1', glossary='does stuff', validator=parsers.parse_convergence, default=6)

    assert subject.value == 1.e-6
    #assert subject.has_changed is False
    assert subject.is_default() is True 


def test_5_conv():
    subject = Keyword(keyword='opt1', glossary='does stuff', validator=parsers.parse_convergence, default=6)

    with pytest.raises(qcdb.KeywordValidationError):
        #subject.value = 5.5
        subject.require(5.5)


def test_6_conv():
    subject = Keyword(keyword='opt1', glossary='does stuff', validator=parsers.parse_convergence, default=6)

    #subject.value = 1
    subject.require(1)
    assert subject.value == 1.e-1
    #assert subject.has_changed is True
    assert subject.is_default() is False

def test_7_conv():
    subject = Keyword(keyword='opt1', glossary='does stuff', validator=parsers.parse_convergence, default=6)

    with pytest.raises(qcdb.KeywordValidationError):
        #subject.value = -1
        subject.require(-1)

def test_8_conv():
    subject = Keyword(keyword='opt1', glossary='does stuff', validator=parsers.parse_convergence, default=6)

    with pytest.raises(qcdb.KeywordValidationError):
        #subject.value = -1.0
        subject.require(-1.0)


def test_9_mem():
    subject = Keyword(keyword='memory', default='700 mb', validator=parsers.parse_memory)

    assert subject.value == 700000000

    subject.require(800000000)
    assert subject.value == 800000000

    subject.require('.6 Gb')
    assert subject.value == 600000000

    with pytest.raises(qcdb.KeywordValidationError):
        subject.require('8 dimms')


def test_20():
    subjects = Keywords()
    subjects.add('qcdb', Keyword(keyword='memory', default='700 mb', validator=parsers.parse_memory))
    subjects.add('qcdb', Keyword(keyword='e_convergence', default=7, validator=parsers.parse_convergence))
    subjects.add('dftd3', Keyword(keyword='Opt1', default=4, validator=lambda x: x))
    subjects.add('dftd3', Keyword(keyword='Opt1', default='cat', validator=lambda x: x))

    print(subjects)
    #assert False


def test_21a():
    subjects = Keywords()

    with pytest.raises(qcdb.ValidationError):
        subjects.add('random', Keyword(keyword='memory', default='700 mb', validator=parsers.parse_memory))

def test_21b():
    subjects = Keywords()
    with pytest.raises(qcdb.ValidationError):
        subjects.require('qcdb', 'mmry', '4 gb', 1234)


def test_22a():
    subjects = Keywords()
    subjects.add('qcdb', Keyword(keyword='memory', default='700 mb', validator=parsers.parse_memory))
    subjects.require('qcdb', 'memory', 9000000000, 22342345)
    subjects.suggest('qcdb', 'memory', 4000000000, 12342345)
    subjects.require('qcdb', 'memory', '9 gb', '00000000')

    assert subjects.scroll['QCDB']['MEMORY'].value == 9000000000


def test_22b():
    subjects = Keywords()
    subjects.add('qcdb', Keyword(keyword='memory', default='700 mb', validator=parsers.parse_memory))
    subjects.require('qcdb', 'memory', 9000000000, 22342345)
    subjects.suggest('qcdb', 'memory', 4000000000, 12342345)
    subjects.require('qcdb', 'memory', '8 gb', '00000000')

    with pytest.raises(qcdb.KeywordReconciliationError):
        assert subjects.scroll['QCDB']['MEMORY'].value == 8000000000


def test_22c():
    subjects = Keywords()
    subjects.add('qcdb', Keyword(keyword='memory', default='700 mb', validator=parsers.parse_memory))
    subjects.require('qcdb', 'memory', 9000000000, 22342345)
    subjects.suggest('qcdb', 'memory', 4000000000, 12342345)
    subjects.require('qcdb', 'memory', '8 gb', 555)  # no user signal so trumps 2234

    assert subjects.scroll['QCDB']['MEMORY'].value == 8000000000
    assert subjects.scroll['QCDB']['MEMORY'].is_default() is False

    import json
    s = json.dumps(subjects, sort_keys=True, indent=2, default=lambda x: x.__dict__)
    print(s)  # all but validator


def test_22d():
    subjects = Keywords()
    subjects.add('qcdb', Keyword(keyword='memory', default='700 mb', validator=parsers.parse_memory))
    subjects.suggest('qcdb', 'memory', 4000000000, 12342345)

    assert subjects.scroll['QCDB']['MEMORY'].value == 4000000000
    assert subjects.scroll['QCDB']['MEMORY'].is_default() is False


def test_22e():
    subjects = Keywords()
    subjects.add('qcdb', Keyword(keyword='memory', default='700 mb', validator=parsers.parse_memory))

    assert subjects.scroll['QCDB']['MEMORY'].value == 700000000
    assert subjects.scroll['QCDB']['MEMORY'].is_default() is True


def test_23a():
    subjects = Keywords()
    subjects.add('qcdb', Keyword(keyword='scf_e_conv', default=5, validator=parsers.parse_convergence))
    subjects.suggest('qcdb', 'scf_e_conv', 1.e-6, 1234)

    assert subjects.scroll['QCDB']['SCF_E_CONV'].value == 1.e-6


@pytest.mark.xfail(True, reason='have not yet healed namespaced options', run=True, strict=True)
def test_23b():
    subjects = Keywords()
    subjects.add('qcdb', Keyword(keyword='scf_e_conv', default=5, validator=parsers.parse_convergence))
    subjects.suggest('qcdb', 'e_conv', 1.e-6, 1234)

    assert subjects.scroll['QCDB']['SCF_E_CONV'].value == 1.e-5


@pytest.mark.xfail(True, reason='have not yet healed namespaced options', run=True, strict=True)
def test_23c():
    subjects = Keywords()
    subjects.add('qcdb', Keyword(keyword='scf_e_conv', default=5, validator=parsers.parse_convergence))
    subjects.require('qcdb', 'e_conv', 1.e-6, 1234)

    assert subjects.scroll['QCDB']['SCF_E_CONV'].value == 1.e-6


@pytest.mark.xfail(True, reason='have not yet healed namespaced options', run=True, strict=True)
def test_23d():
    subjects = Keywords()
    subjects.add('qcdb', Keyword(keyword='scf_e_conv', default=5, validator=parsers.parse_convergence))
    subjects.require('qcdb', 'scf_e_conv', 7, 1234)
    subjects.require('qcdb', 'e_conv', 1.e-6, 1234)

    assert subjects.scroll['QCDB']['SCF_E_CONV'].value == 1.e-7


def test_24():
    subjects = Keywords()
    subjects.add('qcdb', Keyword(keyword='scf_e_conv', default=5, validator=parsers.parse_convergence))

    @register_kwds(subjects)
    def energy(count=1, **kwargs):
        if count > 3:
            assert subjects.scroll['QCDB']['SCF_E_CONV'].value == 1.e-4
            return
        else:
            count += 1
            #subjects.require('qcdb', 'scf_c_conv', count, accession=kwargs['accession'])
            subjects.require('qcdb', 'scf_e_conv', count, accession=kwargs['accession'])
            proc(count)

    @register_kwds(subjects)
    def proc(count, **kwargs):
        energy(count)

    assert subjects.scroll['QCDB']['SCF_E_CONV'].value == 1.e-5
    energy()
    assert subjects.scroll['QCDB']['SCF_E_CONV'].value == 1.e-5


@pytest.fixture
def alias_setup():
    subjects = Keywords()
    subjects.add('qcdb', Keyword(keyword='freeze__core', default=0, validator=parsers.nonnegative_integer))
    subjects.add_alias('qcdb', AliasKeyword(alias='freeze', target='freeze__core'))

    return subjects


def test_alias_a(alias_setup):
    alias_setup.require('qcdb', 'freeze__core', 4, Keywords.mark_of_the_user)

    assert alias_setup.scroll['QCDB']['FREEZE__CORE'].value == 4
    assert alias_setup.scroll['QCDB']['FREEZE__CORE'].is_default() is False


def test_alias_b(alias_setup):
    alias_setup.require('qcdb', 'freeze', 4, Keywords.mark_of_the_user)

    assert alias_setup.scroll['QCDB']['FREEZE__CORE'].value == 4
    assert alias_setup.scroll['QCDB']['FREEZE__CORE'].is_default() is False


def test_alias_c(alias_setup):

    with pytest.raises(qcdb.KeywordValidationError) as e:
        alias_setup.require('qcdb', 'freeze', -1, Keywords.mark_of_the_user)

    assert 'Keyword (FREEZE__CORE) value (-1) does not pass' in str(e.value)


def test_alias_d(alias_setup):

    with pytest.raises(qcdb.ValidationError) as e:
        alias_setup.add_alias('qcdb', AliasKeyword(alias='melt', target='melt__core'))

    assert 'Keyword alias must point to existing keyword proper' in str(e.value)


def test_alias_e(alias_setup):

    with pytest.raises(qcdb.ValidationError) as e:
        alias_setup.add_alias('qcdb', AliasKeyword(alias='freeze__core', target='melt__core'))

    assert 'Keyword alias must not share a name with keyword proper' in str(e.value)
