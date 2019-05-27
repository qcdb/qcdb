import collections

import numpy as np

class QCAspect(collections.namedtuple('QCAspect', 'label units data comment doi glossary')):
    """Facilitates the storage of quantum chemical results by labeling them with basic metadata.

    Attributes
    ----------
    label : str
        Official label for `data`, often qcvar. May contain spaces.
    units : str
        ASCII, LaTeX-like representation of units, without square brackets.
    data : float or :py:class:`numpy.ndarray`
        Value for `label`.
    comment : str, optional
        Additional notes.
    doi : str, optional
        Literature citation or definition DOI link.
    glossary : str, optional
        Extended description or definition.

    """
    def __new__(cls, label, units, data, comment='', doi=None, glossary=''):
        return super(QCAspect, cls).__new__(cls, label, units, data, comment, doi, glossary)

    def __str__(self, label=''):
        width = 40
        text = []
        text.append('-' * width)
        text.append('{:^{width}}'.format('QCAspect ' + self.label, width=width))
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


    def to_dict(self):
        dicary = dict(self._asdict())  # dict, not OrderedDict
        for d in ['doi', 'comment', 'glossary']:
            dicary.pop(d)
        if isinstance(self.data, (np.ndarray, np.number)):
            if self.data.dtype == np.complex:
                dicary['data'] = [dicary['data'].real.tolist(), dicary['data'].imag.tolist()]
            else:
                dicary['data'] = dicary['data'].tolist()
        elif isinstance(self.data, (complex, np.complex)):
            dicary['data'] = [self.data.real, self.data.imag]

        return dicary
