import numpy
from .wake_elements import LongResonatorElement
from .wake_object import Wake, WakeType, WakeComponent
from at.physics import ELossMethod
from at.lattice import AtError, get_cells, checktype, RFCavity


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


def get_anal_qs(qs0, current, volt, rshunt, phil, u0):
    '''
    This function computes the analytical synchrtron tun shift
    '''
    phis = numpy.arcsin(u0/volt)
    theta = numpy.pi/2-phis-phil
    vb = 2*current*rshunt
    a = volt*numpy.cos(phis)*numpy.sin(theta)+u0*numpy.cos(theta)
    b = volt*numpy.cos(phis)*numpy.cos(theta)-(vb+u0)*numpy.sin(theta)
    xi = numpy.arcsin(a*b/numpy.sqrt(a**2*(a**2+b**2)))
    qs = qs0*(numpy.sqrt(volt*numpy.cos(phis)-vb*numpy.cos(xi) *
              (numpy.cos(xi+phis) * numpy.sin(phis)-numpy.sin(xi+phis) *
              numpy.cos(phis)))/numpy.sqrt(volt*numpy.cos(phis)))
    return qs


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
    vgen = rshunt*numpy.cos(psi)*(volt/rshunt+2*current *
                                  numpy.sin(phis))/numpy.cos(phil)
    dff = -numpy.tan(psi)/(2*qfactor)
    fres = frf/(dff+1)
    return vgen, fres, psi


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

    def update_gen_res(self, current):
        vgen, fres, psi = get_anal_values_phasor(self._frf0, current, self._v0,
                                                 self._qfactor, self._rshunt,
                                                 self.Phil, self._u0)
        self.ResFrequency = fres
        for c in self._ring.select(self._cavpts):
            c.PhaseLag = -psi
            c.Voltage = vgen/self._ncavs
