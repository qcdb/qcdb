import collections

class UgghOption(collections.namedtuple('UgghOption', 'value has_changed, clobber, superclobber')):

    def __new__(cls, value, has_changed, clobber=None, superclobber=None):
        return super(UgghOption, cls).__new__(cls, value, has_changed, clobber, superclobber)

    def __str__(self, label=''):
        width = 40
        text = []
        text.append('-' * width)
        text.append('{:^{width}}'.format('AlignmentMill', width=width))
        if label:
            text.append('{:^{width}}'.format(label))
        text.append('-' * width)
        text.append('Value:        {}'.format(self.value))
        text.append('has_changed:  {}'.format(self.has_changed))
        text.append('Clobber:      {}'.format(self.clobber))
        text.append('Superclobber: {}'.format(self.superclobber))
        text.append('-' * width)
        return ('\n'.join(text))

    def align_coordinates(self, geom, reverse=False):
        """suitable for geometry or displaced geometry"""

        algeom = np.copy(geom)
        if reverse:
            algeom = algeom.dot(self.rotation)
            algeom = algeom + self.shift
            if self.mirror:


