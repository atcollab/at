"""Description of particles"""
from ..constants import e_mass, p_mass, mu_mass
import numpy
from warnings import warn
from typing import Optional, Dict


class Particle(object):
    """
    Particle object: it defines the properties of the particles circulating
    in a ring
    """
    _known = dict(
        relativistic=dict(rest_energy=0.0, charge=-1.0),
        electron=dict(rest_energy=e_mass, charge=-1.0),
        positron=dict(rest_energy=e_mass, charge=1.0),
        proton=dict(rest_energy=p_mass, charge=1.0),
        antiproton=dict(rest_energy=mu_mass, charge=-1.0),
        posmuon=dict(rest_energy=mu_mass, charge=1.0),
        negmuon=dict(rest_energy=mu_mass, charge=-1.0),
    )

    def __init__(self, name: Optional[str] = 'relativistic', **kwargs):
        """

        Parameters:
            name:   Particle name. 'electron', 'positron', 'proton', 'posmuon',
              'negmuon' are predefined. For other particles, the rest energy
              and charge must be provided as keywords.

        Keyword Arguments:
            rest_energy:    Particle rest energy [ev]
            charge:         Particle charge [elementary charge]
            *:              Other keywords will be set as attributes of the
                              particle
        """
        if name != 'relativistic':
            warn(UserWarning("AT tracking still assumes beta==1\n"
                             "Make sure your particle is ultra-relativistic"))
        if name in self._known:
            kwargs.update(self._known[name])
        self.name = name
        # Use a numpy scalar to allow division by zero
        self._rest_energy = numpy.array(kwargs.pop('rest_energy'), dtype=float)
        self._charge = kwargs.pop('charge')
        for (key, val) in kwargs.items():
            setattr(self, key, val)

    def to_dict(self) -> Dict:
        attrs = vars(self).copy()
        attrs['rest_energy'] = attrs.pop('_rest_energy')
        attrs['charge'] = attrs.pop('_charge')
        return attrs

    def __repr__(self):
        if self.name in self._known:
            return "Particle('{0}')".format(self.name)
        else:
            attrs = self.to_dict()
            name = attrs.pop('name')
            args = ', '.join('{0}={1!r}'.format(k, v)
                             for k, v in attrs.items())
            return "Particle('{0}', {1})".format(name, args)

    def __str__(self):
        if self.name in self._known:
            return self.name
        else:
            return self.__repr__()

    # Use properties so that they are read-only
    @property
    def rest_energy(self) -> numpy.ndarray:
        """Particle rest energy [eV]"""
        return self._rest_energy

    @property
    def charge(self) -> float:
        """Particle charge [elementary charge]"""
        return self._charge
