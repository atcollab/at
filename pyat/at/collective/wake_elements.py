import numpy
# noinspection PyProtectedMember
from ..lattice.elements import Element, _array
from ..lattice.constants import clight, qe
from .wake_object import resonator, longresonator, reswall


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
        self.Current = kwargs.pop('Current', 0.0)
        self._nslice = kwargs.pop('Nslice', 101)
        self._nturns = kwargs.pop('Nturns', 1)
        self._turnhistory = None    # Defined here to avoid warning
        self.clear_history()
        self.NormFact = kwargs.pop('NormFact', numpy.ones(3, order='F'))
        self._wakeT = wake.srange
        self._nelem = len(self._wakeT)
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
        
    def rebuild_wake(wake):
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

    def clear_history(self):
        self._turnhistory = numpy.zeros((self._nturns*self._nslice, 4),
                                        order='F')

    def set_normfactxy(self, ring):
        l0, _, _ = ring.get_optics(ring)
        self.NormFact[0] = 1/l0['beta'][0]
        self.NormFact[1] = 1/l0['beta'][1]
        
    def get_opt_attr(self, attrname):
        if hasattr(self, attrname):
            return getattr(self, attrname)
        else:
            return None
        
    @property
    def WakeT(self):
        return self._wakeT
        
    @property
    def WakeZ(self):
        return self.get_opt_attr('_wakeZ')
        
    @property
    def WakeDX(self):
        return self.get_opt_attr('_wakeDX')

    @property
    def WakeDY(self):
        return self.get_opt_attr('_wakeDY')
        
    @property
    def WakeQX(self):
        return self.get_opt_attr('_wakeQX')

    @property
    def WakeQY(self):
        return self.get_opt_attr('_wakeQY')

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
        self._wakecomponent = wakecomp
        self._frequency = frequency
        self._qfactor = qfactor
        self._rshunt = rshunt
        self._yokoya_factor = yokoya_factor
        self._beta = beta
        wake = resonator(srange, wakecomp, frequency, qfactor, rshunt,
                         ring.beta, yokoya_factor=yokoya_factor)
        super(ResonatorElement, self).__init__(family_name, ring, wake,
                                               **kwargs)    
                                               
    def rebuild_wake(self):          
        wake = resonator(self.WakeT, self._wakecomponent, self.ResFrequency,
                         self.Qfactor, self.Rshunt, self._beta, self.Yokoya)    
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
                                               
    @property
    def Qfactor(self):
        return self._qfactor
            
    @Qfactor.setter
    def Qfactor(self, qfactor):
        self._qfactor = qfactor
        self.rebuild_wake()
            
                                              


class LongResonatorElement(WakeElement):
    """Class to generate a longitudinal resonator, inherits from WakeElement
       additional argument are frequency, qfactor, rshunt
    """
    def __init__(self, family_name, ring, srange, frequency, qfactor, rshunt,
                 **kwargs):
        self.WParams{'Frequency': frequency,
                     'Qfactor': qfactor,
                     'Rshunt': rshunt
                     'Beta': ring.beta}
        wake = longresonator(srange, frequency, qfactor, rshunt, ring.beta)
        super(LongResonatorElement, self).__init__(family_name, ring, wake,
                                                   **kwargs)


class ResWallElement(WakeElement):
    """Class to generate a resistive wall element, inherits from WakeElement
       additional argument are yokoya_factor, length, pipe radius, conductivity
    """
    def __init__(self, family_name, ring, srange, wakecomp, rwlength,
                 rvac, conduc, yokoya_factor=1, **kwargs):
        self.WParams{'WakeComponent': wakecomp,
                     'Length': rwlength,
                     'Rvac': rvac,
                     'Conductivity': conduct,
                     'Yokoya': yokoya_factor,
                     'Beta': ring.beta}
        wake = reswall(srange, wakecomp, rwlength, rvac, conduc, ring.beta,
                       yokoya_factor=yokoya_factor)
        super(ResWallElement, self).__init__(family_name, ring, wake, **kwargs)
