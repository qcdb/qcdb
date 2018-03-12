import collections

#QCAspect = collections.namedtuple('QCAspect', 'lbl unit data comment')


class QCAspect(collections.namedtuple('QCAspect', 'lbl units data comment doi glossary')):
    #"""Facilitates the application of the simple transformation operations
    #defined by namedtuple of arrays as recipe to the data structures
    #describing Cartesian molecular coordinates. Attaches functions to
    #transform the geometry, element list, gradient, etc. to the
    #AlignmentRecipe. When `mirror` attribute (defaults to False) active,
    #then molecular system can be substantively changed by procedure.

    #"""

    def __new__(cls, lbl, units, data, comment='', doi=None, glossary=''):
        return super(QCAspect, cls).__new__(cls, lbl, units, data, comment, doi, glossary)

    def __str__(self, label=''):
        width = 40
        text = []
        text.append('-' * width)
        text.append('{:^{width}}'.format('QCAspect ' + self.lbl, width=width))
        if label:
            text.append('{:^{width}}'.format(label))
        text.append('-' * width)
        text.append('Data:     {}'.format(self.data))
        text.append('Units:    [{}]'.format(self.units))
        text.append('doi:      {}'.format(self.doi))
        text.append('Comment:  {}'.format(self.comment))
        text.append('Glossary: {}'.format(self.glossary))
        #text.append('Rotation:')
        #text.append('{}'.format(self.rotation))
        text.append('-' * width)
        return ('\n'.join(text))


