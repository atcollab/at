"""Description of particles"""
from ..constants import e_mass, p_mass
import numpy
from warnings import warn
from typing import Optional, Dict


class Particle(object):
    """
    Particle object: it defines the properties of the particles circulating
    in a ring and the beam properties
    """
    _known = dict(
        relativistic=dict(rest_energy=0.0, charge=-1.0),
        electron=dict(rest_energy=e_mass, charge=-1.0),
        positron=dict(rest_energy=e_mass, charge=1.0),
        proton=dict(rest_energy=p_mass, charge=1.0)
    )

    def __init__(self, name: Optional[str] = 'relativistic', **kwargs):
        """

        Parameters:
            name:   Particle name. 'electron', 'positron and 'proton' are
              predefined. For other particles, the rest energy and charge
              must be provided as keywords.

        Keyword Arguments:
            rest_energy:    Particle rest energy [ev]
            charge:         Particle charge [elementary charge]
            beam_current    Total beam current [A]
            fillpattern     boolean array to define filled buckets
                              or a scalar to define the number of equidistant
                              bunches (default=``numpy.array([True])``) .
            weights         vector of double to define bunch weights
                              (default=``numpy.ones(len(fillpattern))``)
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
        # Load parameters of the beam
        self.beam_current = kwargs.pop('beam_current', 0.0)
        self._fillpattern = kwargs.pop('fillpattern', numpy.ones(1))
        self._weights = kwargs.pop('weights',
                                   numpy.ones(len(self._fillpattern)))
        # load remaining keyword arguments
        for (key, val) in kwargs.items():
            setattr(self, key, val)

    def to_dict(self) -> Dict:
        attrs = vars(self).copy()
        attrs['rest_energy'] = attrs.pop('_rest_energy')
        attrs['charge'] = attrs.pop('_charge')
        attrs['nbunch'] = self.nbunch
        attrs = dict((k, v) for (k, v) in attrs.items()
                     if not k.startswith('_'))
        return attrs

    def __repr__(self):
        if self.name in self._known:
            return ("Particle('{0}', beam_current={1}, nbunch={2})"
                    .format(self.name, self.beam_current, self.nbunch))
        else:
            attrs = self.to_dict()
            name = attrs.pop('name')
            args = ', '.join('{0}={1!r}'.format(k, v)
                             for k, v in attrs.items)
            return "Particle('{0}', {1})".format(name, args)

    # Use properties so that they are read-only
    @property
    def rest_energy(self) -> numpy.ndarray:
        """Particle rest energy [eV]"""
        return self._rest_energy

    @property
    def charge(self) -> float:
        """Particle charge [elementary charge]"""
        return self._charge
        
    @property
    def fillpattern(self):
        return numpy.array(self._fillpattern, dtype=bool)   

    def set_fillpattern(self, bunches=1, harmonic_number=None):
        """Changes the fill pattern

        Keyword Arguments:
            bunches:    Sets the fillpattern.
                        Can be either a vector containing zeros, ones
                        or booleans or a scalar. Used only if the
                        harmonic number is provided. A scalar will
                        generate equidistant bunches while a vector
                        allows to generate arbitrary bunch patterns
                        Changing the length of ``fillpattern`` will reset
                        all weights to 1
            harmonic_number:   number of buckets for the ring, required to
                               generate multi-bunch beam
        .. Note::

           The fillpattern can be set directly byt the ``Lattice`` object
           using ``ring.set_fillpattern(bunches)``, in which the ``ring``
           harmonic number is automatically used.
           Assigning ``Particle`` object with
           ``len(fillpattern)!=ring.harmonic_number`` to ``ring`` will
           reset its fillpattern to the default=numpy.array([True])
           which corresponds to a single bunch
        """
        if (harmonic_number is None):
            warn(UserWarning("Harmonic number not provided, "
                             "nbunch and fillpattern set to 1"))
            fp = numpy.ones(1)
        elif numpy.isscalar(bunches):
            if bunches == 1:
                fp = numpy.ones(1)
            else:
                bs = int(harmonic_number/bunches)
                fp = numpy.zeros((harmonic_number,))
                fp[0::bs] = 1
        else:
            assert len(bunches) == harmonic_number, \
                'Fill pattern has to be of shape ({0},)' \
                .format(harmonic_number)
            assert numpy.all((bunches == 0) | (bunches == 1)), \
                'Fill pattern can only contain 0 or 1'
            fp = numpy.array(bunches)
        self._fillpattern = numpy.array(fp)    
        if len(self.weights) != len(self._fillpattern):
            self.weights = numpy.ones(len(self._fillpattern))

    @property
    def weights(self):
        return self._weights

    @weights.setter
    def weights(self, weights):
        assert len(weights) == len(self._fillpattern), \
            'Weights and fill pattern must have the same length'
        self._weights = numpy.array(weights)

    @property
    def nbunch(self):
        return numpy.sum(self._fillpattern)

    @property
    def bunch_currents(self):
        fw = self.weights[self.fillpattern]
        return numpy.squeeze(self.beam_current*fw/numpy.sum(fw))
