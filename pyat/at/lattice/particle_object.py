from .constants import e_mass, p_mass


class Particle(object):
    """Particle object

    This object defines the properties of the particles circulating in a ring

    Particle(name, **params)

    PARAMETERS
        name        Particle name. 'electron' and 'proton' are predefined. For
                    other particles, the mass and charge must be provided.

    KEYWORDS
        rest_energy Particle rest energy [ev]
        charge      Particle charge [elementary charge]
        *           Other keywords will be set as attributes of the particle
    """
    _known = dict(
        relativistic=dict(rest_energy=0.0, charge=-1.0),
        electron=dict(rest_energy=e_mass, charge=-1.0),
        positron=dict(rest_energy=e_mass, charge=1.0),
        proton=dict(rest_energy=p_mass, charge=1.0)
    )

    def __init__(self, name, **kwargs):
        if name in self._known:
            kwargs.update(self._known[name])
        self.name = name
        for key in ('rest_energy', 'charge'):
            if key not in kwargs:
                raise KeyError('"{}" is undefined'.format(key))
        for (key, val) in kwargs.items():
            setattr(self, key, val)

    def __repr__(self):
        if self.name in self._known:
            return "Particle('{0}')".format(self.name)
        else:
            attrs = vars(self).copy()
            name = attrs.pop('name')
            args = ', '.join('{0}={1!r}'.format(k, v) for k, v in attrs.items())
            return "Particle('{0}', {1})".format(name, args)

    def __str__(self):
        if self.name in self._known:
            return self.name
        else:
            return self.__repr__()
