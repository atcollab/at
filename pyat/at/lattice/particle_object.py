"""Description of particles"""
from ..constants import e_mass, p_mass
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
        self.beam_current = kwargs.pop('beam_current',0.0)
        self._harmn = kwargs.pop('harmonic_number', None)
        self._fillpattern = kwargs.pop('fillpattern', numpy.ones(1))
        self._weights = kwargs.pop('weights', numpy.ones(1))
        #load remaining keyword arguments
        for (key, val) in kwargs.items():
            setattr(self, key, val)

    def to_dict(self) -> Dict:
        attrs = vars(self).copy()
        attrs['rest_energy'] = attrs.pop('_rest_energy')
        attrs['charge'] = attrs.pop('_charge')
        attrs['nbunch'] = self.nbunch
        attrs['bunch_currents'] = self.bunch_currents
        return attrs

    def __repr__(self):
        attrs = dict((k, v) for (k, v) in self.to_dict().items()
                     if not k.startswith('_'))
        return '{0}({1})'.format(self.__class__.__name__, attrs)

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
    def harmonic_number(self):
        return self._harmn
        
    @harmonic_number.setter    
    def harmonic_number(self, harmn):
        self._harmn  = harmn   
      
    @property    
    def fillpattern(self):
        return self._fillpattern

    @fillpattern.setter
    def fillpattern(self, bunches):
        if (self._harmn is None):
            warn(UserWarning("Harmonic number not set in Beam(), "
                             "nbunch and fillpattern set to 1"))    
            fp = numpy.ones(1)        
        elif numpy.isscalar(bunches):
            if bunches == 1:
                fp = numpy.ones(1)
            else:
                bs = int(self._harmn/bunches)
                fp = numpy.zeros((self._harmn,))
                fp[0::bs] = 1
        else:
            assert len(bunches) == self._harmn, \
                'Fill pattern has to be of shape ({0},)'.format(self._harmn)
            assert numpy.all((bunches==0) | (bunches==1)), \
                'Fill pattern can only contain 0 or 1'                       
            fp = numpy.array(bunches)
        self._fillpattern = fp
        
    @property    
    def weights(self):
        return self._weights

    @weights.setter
    def weights(self, weights):
        assert len(weights)==len(self._fillpattern), \
            'Weights and fill pattern must have the same length'                   
        self._weights = numpy.array(weights)
        
    @property
    def nbunch(self):
        return numpy.sum(self._fillpattern)
        
    @property
    def bunch_currents(self):
        return numpy.squeeze(self.current*self._weights/numpy.sum(self._weights))
