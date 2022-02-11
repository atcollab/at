import numpy
from at.collective.wake_elements import LongResonatorElement
from at.collective.wake_object import Wake, WakeType, WakeComponent
from at.physics import ELossMethod
from at.lattice.cavity_access import _select_cav
from at.lattice import AtError


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
        if not ring.radiation:
            raise AtError('BeamLoading needs an active cavity')
        self.ring = ring
        self.cavpts = _select_cav(self.ring, cavpts)
        self.v0 = self.ring.get_rf_voltage(cavpts=self.cavpts)
        self.frf = self.ring.get_rf_frequency(cavpts=self.cavpts)
        self.u0 = self.ring.get_energy_loss(method=ELossMethod.TRACKING)
        self.srange = srange
        self.phil = kwargs.pop('PhiL',0.0)
        self.beta = ring.beta
        super(BeamLoadingElement, self).__init__(family_name, self.ring, self.srange, 
                                                 self.frf, qfactor, rshunt, **kwargs)

    # noinspection PyPep8Naming
    @property
    def Current(self):
        return self.Intensity*self.int2curr

    # noinspection PyPep8Naming
    @Current.setter
    def Current(self, current):
        self.Intensity = current/self.int2curr
        self.update_gen_res(current)        

    def update_gen_res(self, current):
        vgen, fres, psi = get_anal_values_phasor(self.frf, current, self.v0,
                                                 self.QFactor, self.RShunt,
                                                 self.phil, self.u0)
        wake = Wake(self.srange)
        wake.add(WakeType.RESONATOR, WakeComponent.Z,
                 fres, self.QFactor, self.RShunt, self.beta)
        self.WakeZ = wake.Z
        vfact = vgen/v0
        for c in self.ring.select(self.cavpts):
            c.update({'PhaseLag': -psi})
            c.update({'Voltage': c.Voltage*vfact})           
