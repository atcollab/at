import numpy
# noinspection PyProtectedMember
from ..lattice.elements import Element, _array
from ..lattice.constants import clight, qe
from .wake_object import Wake, WakeComponent


# noinspection PyPep8Naming
class WakeElement(Element):
    """Class to generate an AT wake element using the passmethod WakeFieldPass
    args:  family name, ring, wake object
    kwargs: PassMethod=WakeFieldPass
            Current=0   Bunch current [A]
            Nslice=101  Number of slices per bunch
            Nturns=1    Number of turn for the wake field
            ZCuts=None  Limits for fixed slicing, default is adaptive
            NormFact    (default=[1,1,1]) normalization for the 3 planes,
                        to account for beta function at the observation
                        point for example
    """
    REQUIRED_ATTRIBUTES = Element.REQUIRED_ATTRIBUTES

    _conversions = dict(Element._conversions, _nslice=int, _nturns=int,
                        _nelem=int, NumParticles=float, Circumference=float,
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
        zcuts = kwargs.pop('ZCuts', None)
        betrel = ring.beta
        self._charge2current = clight*betrel*qe/ring.circumference
        self._wakefact = -qe/(ring.energy*betrel**2)
        self.NumParticles = kwargs.pop('NumParticles', 0.0)
        self._nslice = kwargs.pop('Nslice', 101)
        self._nturns = kwargs.pop('Nturns', 1)
        self._turnhistory = None    # Defined here to avoid warning
        self.clear_history()
        self.NormFact = kwargs.pop('NormFact', numpy.ones(3, order='F'))
        self._build(wake)
        if zcuts is not None:
            self.ZCuts = zcuts
        super(WakeElement, self).__init__(family_name, **kwargs)

    def _build(self, wake):
        self._wakeT = wake.srange
        self._nelem = len(self._wakeT)
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

    def rebuild_wake(self, wake):
        self._build(wake)

    def clear_history(self):
        self._turnhistory = numpy.zeros((self._nturns*self._nslice, 4),
                                        order='F')

    def set_normfactxy(self, ring):
        l0, _, _ = ring.get_optics(ring)
        self.NormFact[0] = 1/l0['beta'][0]
        self.NormFact[1] = 1/l0['beta'][1]

    @property
    def WakeT(self):
        return self._wakeT

    @property
    def WakeZ(self):
        return getattr(self, '_wakeZ', None)

    @property
    def WakeDX(self):
        return getattr(self, '_wakeDX', None)

    @property
    def WakeDY(self):
        return getattr(self, '_wakeDY', None)

    @property
    def WakeQX(self):
        return getattr(self, '_wakeQX', None)

    @property
    def WakeQY(self):
        return getattr(self, '_wakeQY', None)

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
        return self.NumParticles*self._charge2current

    @Current.setter
    def Current(self, current):
        self.NumParticles = current/self._charge2current

    def __repr__(self):
        """Simplified __repr__ to avoid errors due to arguments
        not defined as attributes
        """
        attrs = dict((k, v) for (k, v) in self.items()
                     if not k.startswith('_'))
        return '{0}({1})'.format(self.__class__.__name__, attrs)


class ResonatorElement(WakeElement):
    """Class to generate a resonator, inherits from WakeElement
       additonal argument are frequency, qfactor, rshunt
    """
    def __init__(self, family_name, ring, srange, wakecomp, frequency, qfactor,
                 rshunt, yokoya_factor=1, **kwargs):
        self._resfrequency = frequency
        self._qfactor = qfactor
        self._rshunt = rshunt
        self._yokoya = yokoya_factor
        self._wakecomponent = wakecomp
        self._beta = ring.beta
        wake = Wake.resonator(srange, wakecomp, frequency, qfactor, rshunt,
                              ring.beta, yokoya_factor=yokoya_factor)
        super(ResonatorElement, self).__init__(family_name, ring, wake,
                                               **kwargs)

    def rebuild_wake(self):
        wake = Wake.resonator(self.WakeT, self._wakecomponent,
                              self._resfrequency, self._qfactor,
                              self._rshunt, self._beta, self._yokoya)
        self._build(wake)

    @property
    def ResFrequency(self):
        return self._resfrequency

    @ResFrequency.setter
    def ResFrequency(self, frequency):
        self._resfrequency = frequency
        self.rebuild_wake()

    @property
    def Qfactor(self):
        return self._qfactor

    @Qfactor.setter
    def Qfactor(self, qfactor):
        self._qfactor = qfactor
        self.rebuild_wake()

    @property
    def Rshunt(self):
        return self._rshunt

    @Rshunt.setter
    def Rshunt(self, rshunt):
        self._rshunt = rshunt
        self.rebuild_wake()

    @property
    def Yokoya(self):
        return self._yokoya

    @Yokoya.setter
    def Yokoya(self, yokoya):
        self._yokoya = yokoya
        self.rebuild_wake()


class LongResonatorElement(ResonatorElement):
    """Class to generate a longitudinal resonator, inherits from WakeElement
       additional argument are frequency, qfactor, rshunt
    """
    def __init__(self, family_name, ring, srange, frequency, qfactor, rshunt,
                 **kwargs):
        super(LongResonatorElement, self).__init__(family_name, ring, srange,
                                                   WakeComponent.Z, frequency,
                                                   qfactor, rshunt, **kwargs)


class ResWallElement(WakeElement):
    """Class to generate a resistive wall element, inherits from WakeElement
       additional argument are yokoya_factor, length, pipe radius, conductivity
    """
    def __init__(self, family_name, ring, srange, wakecomp, rwlength,
                 rvac, conduc, yokoya_factor=1, **kwargs):
        self._wakecomponent = wakecomp
        self._rwlength = rwlength
        self._rvac = rvac
        self._conductivity = conduc
        self._yokoya = yokoya_factor
        self._beta = ring.beta
        wake = Wake.resistive_wall(srange, wakecomp, rwlength, rvac, conduc,
                                   ring.beta, yokoya_factor=yokoya_factor)
        super(ResWallElement, self).__init__(family_name, ring, wake, **kwargs)

    def rebuild_wake(self):
        wake = Wake.resistive_wall(self.WakeT, self._wakecomp, self._rwlength,
                                   self._rvac, self._conduc, self._beta,
                                   yokoya_factor=self._yokoya)
        self._build(wake)

    @property
    def RWLength(self):
        return self._rwlength

    @RWLength.setter
    def RWLength(self, length):
        self._rwlength = length
        self.rebuild_wake()

    @property
    def Conductivity(self):
        return self._conductivity

    @Conductivity.setter
    def Conductivity(self, conduct):
        self._conductivity = conduct
        self.rebuild_wake()

    @property
    def Rvac(self):
        return self._rvac

    @Rvac.setter
    def Rvac(self, rvac):
        self._rvac = rvac
        self.rebuild_wake()

    @property
    def Yokoya(self):
        return self._yokoya

    @Yokoya.setter
    def Yokoya(self, yokoya):
        self._yokoya = yokoya
        self.rebuild_wake()
