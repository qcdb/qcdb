import collections

class QCAspect(collections.namedtuple('QCAspect', 'lbl units data comment doi glossary')):
    """Generic value plus metadata storage.

    Attributes
    ----------
    lbl : str
        Official label for `data`, often qcvar. May contain spaces.
    units : str
        ASCII, LaTeX-like representation of units, without square brackets.
    data : float or :py:class:`numpy.ndarray`
        Value for `lbl`.
    comment : str, optional
        Additional notes.
    doi : str, optional
        Literature citation or definition DOI link.
    glossary : str, optional
        Extended description or definition.

    """
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
        text.append('-' * width)
        return ('\n'.join(text))


