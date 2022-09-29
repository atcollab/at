import numpy
from enum import IntEnum
from at.lattice import Element, RFCavity, Collective
from at.physics import get_timelag_fromU0
from at.constants import qe, clight


def get_qs_beamloading(qs0, vbeam, volt, phis, psi):
    '''
    This function computes the analytical synchrotron tun shift
    '''
    xi = numpy.sin(psi)
    qs = qs0*(numpy.sqrt(volt*numpy.cos(phis)+vbeam*numpy.cos(xi) *
              numpy.sin(xi))/numpy.sqrt(volt*numpy.cos(phis)))
    return qs


def get_params_beamloading(frf, current, volt, qfactor,
                           rshunt, phil, phis, vb=None):
    '''
    This function computes some relevant quantities for beam loading
    Reference: Wilson, P B. Beam loading in high-energy storage rings.
    United States: N. p., 1974. Web. doi:10.2172/7166384.
    INPUTS
    frf     RF frequency
    current Total beam current
    volt    total RF voltage
    qfactor cavity Q factor
    rshunt  cavity shunt impedance
    phil    cavity loaded angle
    phis    synchronous phase
    OUTPUTS
    vgen   generator voltage
    fres   resonator frequency
    psi    cavity phase offset (tuning angle)
    vcav   cavity voltage
    '''
    theta = phis+phil
    if vb is None:
        vb = 2*current*rshunt
    a = volt*numpy.cos(phil)
    b = volt*numpy.sin(phil)-vb*numpy.cos(theta)
    psi = numpy.arcsin(b/numpy.sqrt(a**2+b**2))
    vgen = volt*numpy.cos(psi)+vb*numpy.cos(psi)*numpy.sin(phis)
    fres = frf/(1-numpy.tan(psi)/(2*qfactor))
    vcav = (vgen*numpy.sin(theta-psi) -
            vb*numpy.cos(psi)**2)/numpy.sin(phis)
    return vgen, fres, psi, vcav


def add_beamloading(ring, index, *args, **kwargs):
    c = ring[index]
    assert isinstance(c, RFCavity), \
        'Beam loading can only be assigned to a cavity element'
    bl_elem = BeamLoadingElement(ring, c, *args, **kwargs)
    ring[index] = bl_elem
    return bl_elem


def remove_beamloading(ring, index):
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

    _BUILD_ATTRIBUTES = Element._BUILD_ATTRIBUTES
    default_pass = {False: 'IdentityPass', True: 'RFCavityBeamLoadingPass'}
    _conversions = dict(Element._conversions, _nslice=int, _nturns=int,
                        NumParticles=float, Circumference=float,
                        NormFact=float, WakeFact=float)

    def __init__(self, ring, cavelem, qfactor, rshunt, **kwargs):
        kwargs.setdefault('PassMethod', self.default_pass[True])
        zcuts = kwargs.pop('ZCuts', None)
        family_name = cavelem.FamName + '_BL'
        self.Length = cavelem.Length
        self.Voltage = cavelem.Voltage
        self.Energy = cavelem.Energy
        self.Frequency = cavelem.Frequency
        self.TimeLag = getattr(cavelem, 'TimeLag', 0.0)
        self.PhaseLag = getattr(cavelem, 'PhaseLag', 0.0)
        self.Rshunt = rshunt
        self.Qfactor = qfactor
        self._beta = ring.beta
        self._wakefact = - ring.circumference/(clight *
                                               ring.energy*ring.beta**3)
        self._nslice = kwargs.pop('Nslice', 101)
        self._nturns = kwargs.pop('Nturns', 1)
        self.Phil = kwargs.pop('Phil', 0)
        self._turnhistory = None    # Defined here to avoid warning
        self._vbunch = None
        self._bl_params = numpy.array([self.Frequency, 0.0,
                                       self.Voltage, self.Voltage, 0.0])
        _, ts = get_timelag_fromU0(ring)
        self._phis = 2*numpy.pi*self.Frequency*(ts+self.TimeLag)/clight
        self.clear_history(ring=ring)
        self.NormFact = kwargs.pop('NormFact', 1.0)
        self.PhaseGain = kwargs.pop('PhaseGain', 1.0)
        self.VoltGain = kwargs.pop('VoltGain', 1.0)
        if zcuts is not None:
            self.ZCuts = zcuts
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
        self._vbunch = numpy.zeros((nbunch, ), order='F')
        self._init_bl_params(current)

    def _init_bl_params(self, current):
        vgen, fres, psi, vcav = get_params_beamloading(self.Frequency,
                                                       current,
                                                       self.Voltage,
                                                       self.Qfactor,
                                                       self.Rshunt,
                                                       self.Phil,
                                                       self._phis)
        self._bl_params = numpy.array([fres, 2*current*self.Rshunt,
                                       vcav, vgen, psi])

    @property
    def ResFrequency(self):
        return self._bl_params[0]

    @ResFrequency.setter
    def ResFrequency(self, freq):
        self._bl_params[0] = freq

    @property
    def Vbeam(self):
        return self._bl_params[1]
        
    @property
    def Vbunch(self):
        return self._vbunch

    @property
    def Vcav(self):
        return self._bl_params[2]

    @property
    def Vgen(self):
        return self._bl_params[3]

    @Vgen.setter
    def Vgen(self, volt):
        self._bl_params[3] = volt

    @property
    def Psi(self):
        return self._bl_params[4]

    @Psi.setter
    def Psi(self, psi):
        self._bl_params[4] = psi

    def __repr__(self):
        """Simplified __repr__ to avoid errors due to arguments
        not defined as attributes
        """
        attrs = dict((k, v) for (k, v) in self.items()
                     if not k.startswith('_'))
        return '{0}({1})'.format(self.__class__.__name__, attrs)
