import copy

from . import parsers


def read_options(module):
    opts = {}

    if module == 'GLOBALS':
        opts['GLOBALS'] = {}
    
        opts['GLOBALS']['MEMORY'] = RottenOption(
            default='700 mb',
            validator=parsers.parse_memory,
            glossary='Total memory allocation in bytes.')

        opts['GLOBALS']['BASIS'] = RottenOption(
            default='',
            validator=lambda x: True,
            glossary='Primary orbital basis set.')

#        opts['GLOBALS'][''] = RottenOption(
#            default=,
#            validator=,
#            glossary='')
#
#        opts['GLOBALS'][''] = RottenOption(
#            default=,
#            validator=,
#            glossary='')
#
#        opts['GLOBALS'][''] = RottenOption(
#            default=,
#            validator=,
#            glossary='')

    return opts


class RottenOption(object):

    def __init__(self, default, validator, glossary, expert=False):
        self.glossary = glossary
        self.validator = validator
        self.value = default
        self.default = copy.deepcopy(self.value)
        self.has_changed = False
        self.expert = expert

    @property
    def value(self):
        return self._value
    
    @value.setter
    def value(self, val):
        if self.validator(val):
            self._value = val
            self.has_changed = True
        else:
            raise ValueError('Value {} does not pass validator {}'.format(val, self.validator))

    def is_default(self):
        return self.value == self.default


if __name__ == '__main__':

    validator = lambda x: x > 5

    asdf = RottenOption(glossary='does stuff', validator=validator, default=6)
    assert asdf.value == 6
    assert asdf.has_changed is False
    assert asdf.is_default() is True

    try:
        asdf2 = RottenOption(glossary='does stuff', validator=validator, default=4)
    except ValueError:
        assert True
    else:
        assert False

    asdf.value = 10
    assert asdf.value == 10
    assert asdf.has_changed is True
    assert asdf.is_default() is False

    try:
        asdf.value = -10
    except ValueError:
        assert True
    else:
        assert False

    #print(read_options('GLOBALS'))
