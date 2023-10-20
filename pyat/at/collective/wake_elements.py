import numpy
from ..lattice import Lattice
# noinspection PyProtectedMember
from ..lattice.elements import Element, Collective, _array
from ..constants import clight, qe
from .wake_object import Wake, WakeComponent


# noinspection PyPep8Naming
class WakeElement(Collective, Element):
    """Class to generate an AT wake element using the passmethod WakeFieldPass
    """
    _BUILD_ATTRIBUTES = Element._BUILD_ATTRIBUTES
    default_pass = {False: 'IdentityPass', True: 'WakeFieldPass'}
    _conversions = dict(Element._conversions, _nslice=int, _nturns=int,
                        _nelem=int, _wakeFact=float,
                        NormFact=lambda v: _array(v, (3,)),
                        ZCuts=lambda v: _array(v),
                        _wakeDX=lambda v: _array(v),
                        _wakeDY=lambda v: _array(v),
                        _wakeQX=lambda v: _array(v),
                        _wakeQY=lambda v: _array(v),
                        _wakeZ=lambda v: _array(v),
                        _wakeT=lambda v: _array(v))

    def __init__(self, family_name: str, ring: Lattice, wake: Wake, **kwargs):
        """
        Parameters:
            family name:    Element name
            ring:           Lattice in which the element will be inserted
            wake:           :py:class:`.Wake` object

        Keyword Arguments:
            Nslice (int):       Number of slices per bunch. Default: 101
            Nturns (int):       Number of turn for the wake field. Default: 1
            ZCuts:              Limits for fixed slicing, default is adaptive
            NormFact (Tuple[float,...]):    Normalization for the 3 planes,
              to account for beta function at the observation point for
              example. Default: (1,1,1)
"""
        kwargs.setdefault('PassMethod', self.default_pass[True])
        zcuts = kwargs.pop('ZCuts', None)
        self._wakefact = - ring.circumference/(clight *
                                               ring.energy*ring.beta**3)
        self._nslice = kwargs.pop('Nslice', 101)
        self._nturns = kwargs.pop('Nturns', 1)
        self._turnhistory = None    # Defined here to avoid warning
        self.clear_history(ring=ring)
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

    def clear_history(self, ring=None):
        if ring is not None:
            self._nbunch = ring.nbunch
        tl = self._nturns*self._nslice*self._nbunch
        self._turnhistory = numpy.zeros((tl, 4), order='F')

    def set_normfactxy(self, ring):
        l0, _, _ = ring.get_optics()
        self.NormFact[0] = 1/l0['beta'][0]
        self.NormFact[1] = 1/l0['beta'][1]

    @property
    def WakeT(self):
        return self._wakeT

    @property
    def WakeZ(self):
        """Longitudinal component"""
        return getattr(self, '_wakeZ', None)

    @property
    def WakeDX(self):
        """Dipole X component"""
        return getattr(self, '_wakeDX', None)

    @property
    def WakeDY(self):
        """Dipole Y component"""
        return getattr(self, '_wakeDY', None)

    @property
    def WakeQX(self):
        """Quadrupole X component"""
        return getattr(self, '_wakeQX', None)

    @property
    def WakeQY(self):
        """Quadrupole Y component"""
        return getattr(self, '_wakeQY', None)

    @property
    def Nslice(self):
        """Number of slices per bunch"""
        return self._nslice

    @Nslice.setter
    def Nslice(self, nslice):
        self._nslice = nslice
        self.clear_history()

    @property
    def Nturns(self):
        """Number of turn for the wake field"""
        return self._nturns

    @Nturns.setter
    def Nturns(self, nturns):
        self._nturns = nturns
        self.clear_history()

    @property
    def TurnHistory(self):
        """Turn histroy of the slices center of mass"""
        return self._turnhistory

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
    def __init__(self, family_name: str, ring: Lattice, srange,
                 wakecomp: WakeComponent,
                 frequency: float, qfactor: float, rshunt: float,
                 yokoya_factor=1, **kwargs):
        r"""
        Parameters:
            family name:    Element name
            ring:           Lattice in which the element will be inserted
            srange:         Vector of s position where to sample the wake
            wakecomp:       Wake component
            frequency:      Resonator frequency [Hz]
            qfactor:        Q factor
            rshunt:         Shunt impedance, [:math:`\Omega`] for longitudinal,
              [:math:`\Omega/m`] for transverse
            yokoya_factor:  Yokoya factor

        Keyword Arguments:
            Nslice (int):       Number of slices per bunch. Default: 101
            Nturns (int):       Number of turn for the wake field. Default: 1
            ZCuts:              Limits for fixed slicing, default is adaptive
            NormFact (Tuple[float,...]):    Normalization for the 3 planes,
              to account for beta function at the observation point for
              example. Default: (1,1,1)
"""
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
        """Resonator frequency [Hz]"""
        return self._resfrequency

    @ResFrequency.setter
    def ResFrequency(self, frequency):
        self._resfrequency = frequency
        self.rebuild_wake()

    @property
    def Qfactor(self):
        """Resonator Q factor"""
        return self._qfactor

    @Qfactor.setter
    def Qfactor(self, qfactor):
        self._qfactor = qfactor
        self.rebuild_wake()

    @property
    def Rshunt(self):
        r"""Resonator shunt impedance, [:math:`\Omega`] for longitudinal,
        [:math:`\Omega/m`] for transverse"""
        return self._rshunt

    @Rshunt.setter
    def Rshunt(self, rshunt):
        self._rshunt = rshunt
        self.rebuild_wake()

    @property
    def Yokoya(self):
        """Resonator Yokoya factor"""
        return self._yokoya

    @Yokoya.setter
    def Yokoya(self, yokoya):
        self._yokoya = yokoya
        self.rebuild_wake()


class LongResonatorElement(ResonatorElement):
    """Class to generate a longitudinal resonator, inherits from WakeElement
       additional argument are frequency, qfactor, rshunt
    """
    def __init__(self, family_name: str, ring: Lattice, srange,
                 frequency: float, qfactor: float, rshunt: float,
                 **kwargs):
        r"""
        Parameters:
            family name:    Element name
            ring:           Lattice in which the element will be inserted
            srange:         Vector of s position where to sample the wake
            frequency:      Resonator frequency [Hz]
            qfactor:        Q factor
            rshunt:         Shunt impedance, [:math:`\Omega`] for longitudinal,
              [:math:`\Omega/m`] for transverse

        Keyword Arguments:
            yokoya_factor:  Yokoya factor
            Nslice (int):       Number of slices per bunch. Default: 101
            Nturns (int):       Number of turn for the wake field. Default: 1
            ZCuts:              Limits for fixed slicing, default is adaptive
            NormFact (Tuple[float,...]):    Normalization for the 3 planes,
              to account for beta function at the observation point for
              example. Default: (1,1,1)
"""
        super(LongResonatorElement, self).__init__(family_name, ring, srange,
                                                   WakeComponent.Z, frequency,
                                                   qfactor, rshunt, **kwargs)


class ResWallElement(WakeElement):
    """Class to generate a resistive wall element, inherits from WakeElement
       additional argument are yokoya_factor, length, pipe radius, conductivity
    """
    def __init__(self, family_name: str, ring: Lattice, srange,
                 wakecomp: WakeComponent, rwlength: float,
                 rvac: float, conduc: float, yokoya_factor=1, **kwargs):
        """
        Parameters:
            family name:    Element name
            ring:           Lattice in which the element will be inserted
            srange:         Vector of s position where to sample the wake
            wakecomp:       Wake component
            rwlength:
            rvac:
            conduc:
            yokoya_factor:  Yokoya factor

        Keyword Arguments:
            Nslice (int):       Number of slices per bunch. Default: 101
            Nturns (int):       Number of turn for the wake field. Default: 1
            ZCuts:              Limits for fixed slicing, default is adaptive
            NormFact (Tuple[float,...]):    Normalization for the 3 planes,
              to account for beta function at the observation point for
              example. Default: (1,1,1)
        """
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
        wake = Wake.resistive_wall(self.WakeT, self._wakecomponent,
                                   self._rwlength, self._rvac,
                                   self._conductivity, self._beta,
                                   yokoya_factor=self._yokoya)
        self._build(wake)

    @property
    def RWLength(self):
        """Length of the resistive wall"""
        return self._rwlength

    @RWLength.setter
    def RWLength(self, length):
        self._rwlength = length
        self.rebuild_wake()

    @property
    def Conductivity(self):
        """Conductivity of the beam pipe [S/m]"""
        return self._conductivity

    @Conductivity.setter
    def Conductivity(self, conduct):
        self._conductivity = conduct
        self.rebuild_wake()

    @property
    def Rvac(self):
        """Radius of the beam pipe [m]"""
        return self._rvac

    @Rvac.setter
    def Rvac(self, rvac):
        self._rvac = rvac
        self.rebuild_wake()

    @property
    def Yokoya(self):
        """Yokoya factor for the reistive wall"""
        return self._yokoya

    @Yokoya.setter
    def Yokoya(self, yokoya):
        self._yokoya = yokoya
        self.rebuild_wake()
