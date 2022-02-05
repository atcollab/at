import at
import numpy
# noinspection PyProtectedMember
from at.lattice.elements import Element, _array
from at.lattice.constants import clight, qe, e_mass
from at.collective.wake_object import Wake, WakeComponent, WakeType


# noinspection PyPep8Naming
class WakeElement(Element):
    """Class to generate an AT wake element using the
    passmethod WakeFieldPass
    args:  family name, ring, wake object
    kwargs: Intensity  (default=0) bunch intensity
            Passmethod (default=WakeFieldPass)
            Nslice     (default=101) number of slices
                       per bunch
            Nturns     (default=1) number of turn for
                       the wake field
            ZCuts      (default=None)limits for fixed
                       slicing, default is adaptative
            NormFact   (default=[1,1,1]) normalization
                       for the 3 planes, to account for
                       beta function at the observation
                       point for example
    """
    REQUIRED_ATTRIBUTES = Element.REQUIRED_ATTRIBUTES

    _conversions = dict(Element._conversions, _nslice=int, _nturns=int,
                        _nelem=int, Intensity=float, Circumference=float,
                        NormFact=lambda v: _array(v, (3,)),
                        WakeFact=float,
                        _wakeDX=lambda v: _array(v),
                        _wakeDY=lambda v: _array(v),
                        _wakeQX=lambda v: _array(v),
                        _wakeQY=lambda v: _array(v),
                        _wakeZ=lambda v: _array(v),
                        _wakeT=lambda v: _array(v))

    def __init__(self, family_name, ring, wake, **kwargs):
        kwargs.setdefault('PassMethod', 'WakeFieldPass')
        self.Intensity = kwargs.pop('Intensity', 0.0)
        self._nslice = kwargs.pop('Nslice', 101)
        self._nturns = kwargs.pop('Nturns', 1)
        self.clear_history()
        self.NormFact = kwargs.pop('NormFact', numpy.ones(3, order='F'))
        self.Circumference = ring.circumference
        self.Wakefact = self.get_wakefact(ring)
        self.int2curr = self.get_int2curr(ring)
        self._wakeT = wake.get_srange()
        self._nelem = len(self._wakeT)
        zcuts = kwargs.pop('ZCuts', None)
        if zcuts is not None:
            self.ZCuts = zcuts
        if wake.Z is not None:
            self._wakeZ = wake.Z
        if wake.DX is not None:
            self._wakeDX = wake.DX
        if wake.DY is not None:
            self._wakeDY = wake.DY
        if wake.QX is not None:
            self._wakeQX = wake.QX
        if wake.QY is not None:
            self._wakeQY = wake.QY
        super(WakeElement, self).__init__(family_name, **kwargs)

    def get_wakefact(self, ring):
        betrel = numpy.sqrt(1.0-(e_mass/ring.energy)**2)
        return -qe/(ring.energy*betrel**2)

    def get_int2curr(self, ring):
        betrel = numpy.sqrt(1.0-(e_mass/ring.energy)**2)
        return clight*betrel*qe/ring.circumference

    def clear_history(self):
        self._turnhistory = numpy.zeros((self._nturns*self._nslice, 4),
                                       order='F')

    def set_normfactxy(self, ring):
        l0, _, _ = at.get_optics(ring)
        self.NormFact[0] = 1/l0['beta'][0]
        self.NormFact[1] = 1/l0['beta'][1]

    @property
    def Nslice(self):
        return self._nslice

    @Nslice.setter
    def Nslice(self, nslice):
        self._nslice = nslice
        self.clear_history()

    @property
    def Nturns(self):
        return self._nslice

    @Nturns.setter
    def Nturns(self, nslice):
        self._nslice = nslice
        self.clear_history()

    @property
    def Current(self):
        return self.Intensity*self.int2curr

    @Current.setter
    def Current(self, current):
        self.Intensity = current/self.int2curr

    def __repr__(self):
        """Simplified __repr__ to avoid errors due to arguments
        not defined as attributes
        """
        attrs = dict(self.items())
        return '{0}({1})'.format(self.__class__.__name__, attrs)


class LongResonatorElement(WakeElement):
    """Class to generate a longitudinal resonator, inherits from WakeElement
       additonal argument are frequency, qfactor, rshunt
    """
    def __init__(self, family_name, ring, srange, frequency, qfactor,
                 rshunt, **kwargs):
        self.Frequency = frequency
        self.QFactor = qfactor
        self.RShunt = rshunt
        beta = numpy.sqrt(1.0-(e_mass/ring.energy)**2)
        wake = Wake(srange)
        wake.add(WakeType.RESONATOR, WakeComponent.Z,
                 frequency, qfactor, rshunt, beta)
        super(LongResonatorElement, self).__init__(family_name, ring=ring,
                                                   wake=wake, **kwargs)


class TransResonatorElement(WakeElement):
    """Class to generate a transverse resonator, inherits from WakeElement
       additonal argument are frequency, qfactor, rshunt
    """
    def __init__(self, family_name, ring, srange, wakecomp, frequency, qfactor,
                 rshunt, **kwargs):
        self.Frequency = frequency
        self.QFactor = qfactor
        self.RShunt = rshunt
        beta = numpy.sqrt(1.0-(e_mass/ring.energy)**2)
        wake = Wake(srange)
        wake.add(WakeType.RESONATOR, wakecomp,
                 frequency, qfactor, rshunt, beta)
        super(TransResonatorElement, self).__init__(family_name, ring=ring,
                                                    wake=wake, **kwargs)


class ResWallElement(WakeElement):
    """Class to generate a resistive wall element, inherits from WakeElement
       additonal argument are yokoya_factor, length, pipe radius, conductivity
    """
    def __init__(self, family_name, ring, srange, yokoya_factor, length,
                 rvac, conduc, wakecomp, **kwargs):
        self.yokoya_factor = yokoya_factor
        self.length = length
        self.rvac = rvac
        self.conduc = conduc
        beta = numpy.sqrt(1.0-(e_mass/ring.energy)**2)
        wake = Wake(srange)
        wake.add(WakeType.RESWALL, wakecomp,
                 length, rvac, conduc, beta, yokoya_factor)
        super(ResWallElement, self).__init__(family_name, ring=ring,
                                             wake=wake, **kwargs)
