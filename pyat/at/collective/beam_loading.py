import numpy
from .wake_elements import LongResonatorElement
from .wake_object import Wake, WakeType, WakeComponent
from at.physics import ELossMethod
from at.lattice import AtError, get_cells, checktype, RFCavity
from .wake_functions import long_resonator_wf
from ..lattice.constants import qe


def _select_cav(ring, cavpts):
    """Select the cavities"""
    if cavpts is None:
        cavpts = get_cells(ring, checktype(RFCavity))
    return cavpts


def _check_unique(vals, attr):
    vu = numpy.unique(vals)
    if len(vu) > 1:
        raise AtError('Not all cavities {0} are identical: '
                      'please define one beamloading element '
                      'per cavity type'.format(attr))


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


class BeamLoadingElement(LongResonatorElement):
    """Class to generate a longitudinal resonator, inherits from WakeElement
       additonal argument are frequency, qfactor, rshunt
    """
    def __init__(self, family_name, ring, srange, qfactor, rshunt,
                 cavpts=None, **kwargs):
        self._cavpts = _select_cav(ring, cavpts)
        self._ncavs = ring.refcount(self._cavpts)
        self._ring = ring
        self._cavs = ring.select(self._cavpts)
        va = ring.get_rf_voltage(cavpts=self._cavpts, array=True)
        _check_unique(va, 'Voltage')
        frf = ring.get_rf_frequency(cavpts=self._cavpts, array=True)
        _check_unique(frf, 'Frequency')
        self._v0 = numpy.sum(va)
        self._frf0 = numpy.unique(frf)[0]
        self._u0 = ring.get_energy_loss(method=ELossMethod.TRACKING)
        self.Phil = kwargs.pop('PhiL', 0.0)
        self._beta = ring.beta
        self._psi = 0.0
        self._vgen = numpy.sum(va)
        self._vcav = numpy.sum(va)
        super(BeamLoadingElement, self).__init__(family_name, self._ring,
                                                 srange, self._frf0,
                                                 qfactor, rshunt, **kwargs)

    @property
    def Current(self):
        return self.NumParticles*self._charge2current

    @Current.setter
    def Current(self, current):
        self.NumParticles = current/self._charge2current
        self.update_gen_res(current)
     
    @property  
    def Vbeam(self):
        thist = self._turnhistory[-1,2] - self._turnhistory[:,2]
        swz = long_resonator_wf(thist, self.ResFrequency, self.Qfactor,
                                self.Rshunt, self._beta)
        swz *= self._turnhistory[:,3]
        swz = numpy.sum(swz)*self.NumParticles*qe 
        return numpy.absolute(swz)/numpy.cos(self._psi)**2  
        
    @property
    def Vcav(self):
        return self._vcav 
        
    @property
    def Vgen(self):
        return self._vgen
        
    @property
    def Psi(self):
        return self._psi   

    def update_gen_res(self, current, vb=None):
        vgen, fres, psi, vcav = get_anal_values_phasor(self._frf0, current, self._v0,
                                                       self._qfactor, self._rshunt,
                                                       self.Phil, self._u0, vb=vb)
        self.ResFrequency = fres
        self._psi = psi
        self._vgen = vgen
        self._vcav = vcav
        for c in self._ring.select(self._cavpts):
            c.PhaseLag = -psi
            c.Voltage = vgen/self._ncavs
