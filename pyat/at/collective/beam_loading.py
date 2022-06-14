import numpy
from at.lattice import Element, RFCavity
from at.physics import ELossMethod
from at.constants import qe, clight


def get_anal_qs(qs0, vbeam, volt, u0, psi):
    '''
    This function computes the analytical synchrotron tun shift
    '''
    phis = numpy.arcsin(u0/volt)
    xi = numpy.sin(psi)
    qs = qs0*(numpy.sqrt(volt*numpy.cos(phis)+vbeam*numpy.cos(xi) * \
        numpy.sin(xi))/numpy.sqrt(volt*numpy.cos(phis)))
    return qs


def get_anal_values_phasor(frf, current, volt, qfactor, rshunt, phil, u0, vb=None):
    '''
    This function computes some relevant quantities for beam loading
    Reference: Wilson, P B. Beam loading in high-energy storage rings.
    United States: N. p., 1974. Web. doi:10.2172/7166384.
    INPUTS
    frf     RF frequency
    volt    total RF voltage
    qfactor cavity Q factor
    rshunt  cavity shunt impedance
    phil    cavity loaded angle
    u0      energy loss per turn
    OUTPUTS
    vgen   generator voltage
    fres   resonator frequency
    psi    cavity phase offset (tuning angle)
    '''
    phis = numpy.arcsin(u0/volt)
    theta = phis+phil
    if vb is None:
        vb = 2*current*rshunt
    a = volt*numpy.cos(phis)*numpy.cos(theta)+u0*numpy.sin(theta)
    b = volt*numpy.cos(phis)*numpy.sin(theta)-(vb+u0)*numpy.cos(theta)
    x = a*b/numpy.sqrt(a**2*(a**2+b**2))
    psi = numpy.arcsin(x)
    vgen = volt*numpy.cos(psi)+vb*numpy.cos(psi)*numpy.sin(phis)
    dff=-numpy.tan(psi)/(2*qfactor);
    fres=frf/(dff+1)
    vcav = (vgen*numpy.sin(theta-psi) - vb*numpy.cos(psi)**2)/numpy.sin(theta)
    return vgen, fres, psi, vcav
  

def add_beamloading(ring, index, *args, **kwargs):
    c = ring[index]
    assert isinstance(c, RFCavity), \
        'Beam loading can only be assigned to a cavity element'
    bl_elem = BeamLoadingElement(ring, c, *args, **kwargs)
    ring[index] = bl_elem
    return bl_elem  
  
    
class BeamLoadingElement(Element):

    REQUIRED_ATTRIBUTES = Element.REQUIRED_ATTRIBUTES

    _conversions = dict(Element._conversions, _nslice=int, _nturns=int,
                        NumParticles=float, Circumference=float,
                        NormFact=float, WakeFact=float)
                        
    def __init__(self, ring, cavelem, qfactor, rshunt, **kwargs):
        kwargs.setdefault('PassMethod', 'RFCavityBeamLoadingPass')
        zcuts = kwargs.pop('ZCuts', None)
        family_name = cavelem.FamName + '_BL'
        self.Length = cavelem.Length
        self.Voltage = cavelem.Voltage
        self.Energy = cavelem.Energy
        self.Frequency = cavelem.Frequency
        self.Timelag = getattr(cavelem, 'TimeLag', 0.0)
        self.Phaselag = getattr(cavelem, 'PhaseLag', None)
        self.Rshunt = rshunt
        self.Qfactor = qfactor
        self._beta = ring.beta
        self._charge2current = clight*ring.beta*qe/ring.circumference
        self._wakefact = -qe/(ring.energy*ring.beta**2)
        self.NumParticles = kwargs.pop('NumParticles', 0.0)
        self._nslice = kwargs.pop('Nslice', 101)
        self._nturns = kwargs.pop('Nturns', 1)
        self.Phil = kwargs.pop('Phil', 0)
        self._turnhistory = None    # Defined here to avoid warning
        self._bl_params = None
        self._u0 = ring.get_energy_loss(method=ELossMethod.TRACKING)
        self.clear_history()
        self.init_bl_params()
        self.NormFact = kwargs.pop('NormFact', 1.0)
        if zcuts is not None:
            self.ZCuts = zcuts
        super(BeamLoadingElement, self).__init__(family_name, **kwargs)
        
    @property
    def Current(self):
        return self.NumParticles*self._charge2current

    @Current.setter
    def Current(self, current):
        self.NumParticles = current/self._charge2current
        self.init_bl_params()
        
    def clear_history(self):
        self._turnhistory = numpy.zeros((self._nturns*self._nslice, 4),
                                        order='F')
                                        
    def init_bl_params(self):
        vgen, fres, psi, vcav = get_anal_values_phasor(self.Frequency, self.Current,
                                                       self.Voltage, self.Qfactor,
                                                       self.Rshunt, self.Phil, self._u0)
        self._bl_params = numpy.array([fres, 2*self.Current*self.Rshunt, vcav, vgen, psi])

    @property  
    def ResFrequency(self):
        return self._bl_params[0]
             
    @property  
    def Vbeam(self):
        return self._bl_params[1] 
        
    @property
    def Vcav(self):
        return self._bl_params[2] 
        
    @property
    def Vgen(self):
        return self._bl_params[3] 
        
    @property
    def Psi(self):
        return self._bl_params[4]   
