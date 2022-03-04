import numpy
from at.collective.wake_elements import LongResonatorElement
from at.collective.wake_object import Wake, WakeType, WakeComponent
from at.physics import ELossMethod
from at.lattice import AtError, get_cells, checktype, RFCavity


def _select_cav(ring, cavpts):
    """Select the cavities"""
    if cavpts is None:
        try:
            cavpts = ring.cavpts
        except AttributeError:
            cavpts = get_cells(ring, checktype(RFCavity))
    return cavpts


def get_anal_values_phasor(frf, current, volt, qfactor, rshunt, phil, u0):
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
    theta = numpy.pi/2-phis-phil
    vb = 2*current*rshunt
    a = volt*numpy.cos(phis)*numpy.sin(theta)+u0*numpy.cos(theta)
    b = volt*numpy.cos(phis)*numpy.cos(theta)-(vb+u0)*numpy.sin(theta)
    x = a*b/numpy.sqrt(a**2*(a**2+b**2))
    psi = numpy.arcsin(x)
    vgr=(u0+vb*(1-x**2))/(numpy.cos(psi)*numpy.cos(psi+theta));
    vgen = rshunt*numpy.cos(psi)*(volt/rshunt+2*current*numpy.sin(phis))/numpy.cos(phil)
    dff=-numpy.tan(psi)/(2*qfactor);
    fres=frf/(dff+1)
    psi = numpy.arctan(2*qfactor*(fres-frf)/fres)
    return vgen,fres,psi


class BeamLoadingElement(LongResonatorElement):
    """Class to generate a longitudinal resonator, inherits from WakeElement
       additonal argument are frequency, qfactor, rshunt
    """
    def __init__(self, family_name, ring, srange, qfactor, rshunt,
                 cavpts=None, **kwargs):
        self._cavpts = _select_cav(ring, cavpts)
        self._ring = ring
        self._cavs = ring.select(cavpts)
        self._v0 = ring.get_rf_voltage(cavpts=cavpts)
        self._frf0 = ring.get_rf_frequency(cavpts=cavpts)
        self._u0 = ring.get_energy_loss(method=ELossMethod.TRACKING)
        self._phil = kwargs.pop('PhiL',0.0)
        self._beta = ring.beta
        super(BeamLoadingElement, self).__init__(family_name, self._ring, srange, 
                                                 self._frf0, qfactor, rshunt,
                                                 **kwargs)

    @property
    def Current(self):
        return self.NumParticles*self._charge2current

    @Current.setter
    def Current(self, current):
        self.NumParticles = current/self._charge2current
        self.update_gen_res(current)        

    def update_gen_res(self, current):
        vgen, fres, psi = get_anal_values_phasor(self._frf0, current, self._v0,
                                                 self._qfactor, self._rshunt,
                                                 self._phil, self._u0)                                       
        self.ResFrequency = fres    
        for c in self._ring.select(self._cavpts):
            c.PhaseLag = -psi
            c.Voltage = vgen          
