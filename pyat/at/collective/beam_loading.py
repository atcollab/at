import numpy
from enum import IntEnum
from ..lattice import Lattice
from at.lattice import Element, RFCavity, Collective
from at.physics import get_timelag_fromU0
from at.constants import clight


class BLMode(IntEnum):
    WAKE = 1
    PHASOR = 2


def add_beamloading(ring: Lattice, index: int, *args, **kwargs):
    """Function to add beam loading to a cavity element, the cavity
    element is changed to a beam loading element that combines the energy
    kick from both the cavity and the resonator

    Parameters:
        ring:           Lattice object
        index:          Index of the cavity element
        qfactor:        Q factor
        rshunt:         Shunt impedance, [:math:`\Omega`] for longitudinal,
            [:math:`\Omega/m`] for transverse

    Keyword Arguments:
        Nslice (int):       Number of slices per bunch. Default: 101
        Nturns (int):       Number of turn for the wake field. Default: 1
        ZCuts:              Limits for fixed slicing, default is adaptive
        NormFact (Tuple[float,...]):    Normalization for the 3 planes,
            to account for beta function at the observation point for
            example. Default: (1,1,1)
        mode (BLMode):  method for beam loading calculation BLMode.PHASOR
            (default) uses the phasor method, BLMode.WAKE uses the wake
            function. For high Q resonator, the phasor method should be
            used
    Returns:
        bl_elem (Element): beam loading element
    """
    c = ring[index]
    assert isinstance(c, RFCavity), \
        'Beam loading can only be assigned to a cavity element'
    bl_elem = BeamLoadingElement(ring, c, *args, **kwargs)
    ring[index] = bl_elem
    return bl_elem


def remove_beamloading(ring, index):
    """Function to remove beam loading from a cavity element, the beam
    loading
    element is changed to a beam loading element that combines the energy
    kick from both the cavity and the resonator

    Parameters:
        ring:           Lattice object
        index:          Index of the cavity element

    Returns:
        cav_elem (Element): cavity element
    """
    c = ring[index]
    assert isinstance(c, BeamLoadingElement), \
        'Cannot remove beam loading: not a beam loading element'
    family_name = c.FamName.replace('_BL', '')
    length = c.Length
    voltage = c.Voltage
    energy = c.Energy
    frequency = c.Frequency
    timelag = getattr(c, 'TimeLag', 0.0)
    phaselag = getattr(c, 'PhaseLag', 0.0)
    harm = numpy.round(frequency*ring.circumference/clight)
    cav_elem = RFCavity(family_name, length, voltage, frequency,
                        harm, energy, TimeLag=timelag,
                        PhaseLag=phaselag)
    ring[index] = cav_elem
    return cav_elem


class BeamLoadingElement(Collective, Element):
    """Class to generate a beamloading element, inherits from Element
       additional argument are ring, cavity, qfactor, rshunt
    """
    _BUILD_ATTRIBUTES = Element._BUILD_ATTRIBUTES
    default_pass = {False: 'IdentityPass', True: 'RFCavityBeamLoadingPass'}
    _conversions = dict(Element._conversions, _nslice=int, _nturns=int,
                        NumParticles=float, Circumference=float,
                        NormFact=float, WakeFact=float)

    def __init__(self, ring, cavelem, qfactor, rshunt, mode=BLMode.PHASOR,
                 **kwargs):
        """
        Parameters:
            ring:           Lattice object
            index:          Index of the cavity element
            qfactor:        Q factor
            rshunt:         Shunt impedance, [:math:`\Omega`] for longitudinal,
                [:math:`\Omega/m`] for transverse

            Keyword Arguments:
            Nslice (int):       Number of slices per bunch. Default: 101
            Nturns (int):       Number of turn for the wake field. Default: 1
            ZCuts:              Limits for fixed slicing, default is adaptive
            NormFact (Tuple[float,...]):    Normalization for the 3 planes,
                to account for beta function at the observation point for
                example. Default: (1,1,1)
            mode (BLMode):  method for beam loading calculation BLMode.PHASOR
                (default) uses the phasor method, BLMode.WAKE uses the wake
                function. For high Q resonator, the phasor method should be
                used
        Returns:
            bl_elem (Element): beam loading element
        """
        kwargs.setdefault('PassMethod', self.default_pass[True])
        assert isinstance(mode, BLMode), \
            'Beam loading mode has to be an instance of BLMode'
        zcuts = kwargs.pop('ZCuts', None)
        phil = kwargs.pop('phil', 0)
        family_name = cavelem.FamName + '_BL'
        self.Length = cavelem.Length
        self.Voltage = cavelem.Voltage
        self.Energy = cavelem.Energy
        self.Frequency = cavelem.Frequency
        self.TimeLag = getattr(cavelem, 'TimeLag', 0.0)
        self.PhaseLag = getattr(cavelem, 'PhaseLag', 0.0)
        self.Rshunt = rshunt
        self.Qfactor = qfactor
        self.NormFact = kwargs.pop('NormFact', 1.0)
        self.PhaseGain = kwargs.pop('PhaseGain', 1.0)
        self.VoltGain = kwargs.pop('VoltGain', 1.0)
        self._mode = int(mode)
        self._beta = ring.beta
        self._wakefact = - ring.circumference/(clight *
                                               ring.energy*ring.beta**3)
        self._nslice = kwargs.pop('Nslice', 101)
        self._nturns = kwargs.pop('Nturns', 1)
        self._turnhistory = None    # Defined here to avoid warning
        self._vbunch = None
        _, ts = get_timelag_fromU0(ring)
        self._phis = 2*numpy.pi*self.Frequency*(ts+self.TimeLag)/clight
        self._vbeam_phasor = numpy.zeros(2)
        self._vbeam = numpy.zeros(2)
        self._vgen = numpy.zeros(2)
        self._vcav = numpy.array([self.Voltage,
                                  numpy.pi/2-self._phis-phil])
        if zcuts is not None:
            self.ZCuts = zcuts
        self.clear_history(ring=ring)
        super(BeamLoadingElement, self).__init__(family_name, **kwargs)

    def clear_history(self, ring=None):
        if ring is None:
            current = 0.0
            nbunch = 1
        else:
            current = ring.beam_current
            nbunch = ring.nbunch
        tl = self._nturns*self._nslice*nbunch
        self._turnhistory = numpy.zeros((tl, 4), order='F')
        self._vbunch = numpy.zeros((nbunch, 2), order='F')
        self._init_bl_params(current)

    def _init_bl_params(self, current):
        theta = -self._vcav[1]+numpy.pi/2
        vb = 2*current*self.Rshunt
        a = self.Voltage*numpy.cos(theta-self._phis)
        b = self.Voltage*numpy.sin(theta-self._phis)-vb*numpy.cos(theta)
        psi = numpy.arcsin(b/numpy.sqrt(a**2+b**2))
        vgen = self.Voltage*numpy.cos(psi) + \
            vb*numpy.cos(psi)*numpy.sin(self._phis)
        self._vgen = numpy.array([vgen, psi])
        self._vbeam = numpy.array([2*current*self.Rshunt*numpy.cos(psi),
                                   numpy.pi-psi])

    @property
    def ResFrequency(self):
        return self.Frequency/(1-numpy.tan(self.Vgen[1])/(2*self.Qfactor))

    @property
    def Vbeam(self):
        return self._vbeam
        
    @property
    def Vbunch(self):
        return self._vbunch

    @property
    def Vcav(self):
        return self._vcav

    @property
    def Vgen(self):
        return self._vgen

    @Vgen.setter
    def Vgen(self, value):
        self._vgen = value

    def __repr__(self):
        """Simplified __repr__ to avoid errors due to arguments
        not defined as attributes
        """
        attrs = dict((k, v) for (k, v) in self.items()
                     if not k.startswith('_'))
        return '{0}({1})'.format(self.__class__.__name__, attrs)
