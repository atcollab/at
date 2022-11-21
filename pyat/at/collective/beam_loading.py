import numpy
from enum import IntEnum
from ..lattice import Lattice
from at.lattice import RFCavity, Collective
from at.lattice.elements import _array
from at.lattice.utils import Refpts, uint32_refpts, make_copy
from at.physics import get_timelag_fromU0
from at.constants import clight
from typing import Sequence, Optional, Union


class BLMode(IntEnum):
    WAKE = 1
    PHASOR = 2


def add_beamloading(ring: Lattice, qfactor: Union[float, Sequence[float]],
                    rshunt: Union[float, Sequence[float]],
                    cavpts: Refpts = None, copy: Optional[bool] = False,
                    **kwargs):
    r"""Function to add beam loading to a cavity element, the cavity
    element is changed to a beam loading element that combines the energy
    kick from both the cavity and the resonator

    Parameters:
        ring:           Lattice object
        qfactor:        Q factor. Scalar or array of float are accepted
        rshunt:         Shunt impedance, [:math:`\Omega`] for longitudinal
                        Scalar or array of float are accepted

    Keyword Arguments:
        cavpts (Refpts):    refpts of the cavity. If None (default) apply
                            to all cavity in the ring
        Nslice (int):       Number of slices per bunch. Default: 101
        Nturns (int):       Number of turn for the wake field. Default: 1
        ZCuts:              Limits for fixed slicing, default is adaptive
        NormFact (float):   Normalization factor
        mode (BLMode):  method for beam loading calculation BLMode.PHASOR
            (default) uses the phasor method, BLMode.WAKE uses the wake
            function. For high Q resonator, the phasor method should be
            used
        copy:       If True, returns a shallow copy of ring with new
                    beam loading elements. Otherwise, modify ring in-place
    """
    @make_copy(copy)
    def apply(ring, cavpts, newelems):
        for ref, elem in zip(cavpts, newelems):
            ring[ref] = elem

    if cavpts is None:
        cavpts = ring.get_refpts(RFCavity)
    else:
        cavpts = uint32_refpts(cavpts, len(ring))
    qfactor = numpy.broadcast_to(qfactor, (len(cavpts), ))
    rshunt = numpy.broadcast_to(rshunt, (len(cavpts), ))
    new_elems = []
    for ref, qf, rs in zip(cavpts, qfactor, rshunt):
        cav = ring[ref]
        assert isinstance(cav, RFCavity), \
            'Beam loading can only be assigned to a cavity element'
        new_elems.append(BeamLoadingElement.build_from_cav(cav, ring, qf,
                                                           rs, **kwargs))
    return apply(ring, cavpts, new_elems)


def remove_beamloading(ring, cavpts: Refpts = None,
                       copy: Optional[bool] = False):
    """Function to remove beam loading from a cavity element, the beam
    loading element is changed to a beam loading element that
    combines the energy kick from both the cavity and the resonator

    Parameters:
        ring:           Lattice object

    Keyword Arguments:
        cavpts (Refpts):    refpts of the beam loading Elements.
                            If None (default) apply to all elements
        copy:       If True, returns a shallow copy of ring with new
                      cavity elements. Otherwise, modify ring in-place
    """
    @make_copy(copy)
    def apply(ring, cavpts, newelems):
        for ref, elem in zip(cavpts, newelems):
            ring[ref] = elem

    if cavpts is None:
        cavpts = ring.get_refpts(BeamLoadingElement)
    else:
        cavpts = uint32_refpts(cavpts, len(ring))
    new_elems = []
    for ref in cavpts:
        bl = ring[ref]
        assert isinstance(bl, BeamLoadingElement), \
            'Cannot remove beam loading: not a beam loading element'
        family_name = bl.FamName.replace('_BL', '')
        harm = numpy.round(bl.Frequency*ring.circumference/clight)
        new_elems.append(RFCavity(family_name, bl.Length, bl.Voltage,
                                  bl.Frequency, harm, bl.Energy,
                                  TimeLag=getattr(bl, 'TimeLag', 0.0),
                                  PhaseLag=getattr(bl, 'PhaseLag', 0.0)))
    return apply(ring, cavpts, new_elems)


class BeamLoadingElement(RFCavity, Collective):
    """Class to generate a beamloading element, inherits from Element
       additional argument are ring, cavity, qfactor, rshunt
    """
    default_pass = {False: 'DriftPass', True: 'BeamLoadingCavityPass'}
    _conversions = dict(RFCavity._conversions,
                        Rshunt=float, Qfactor=float, NormFact=float,
                        PhaseGain=float, VoltGain=float, _mode=int,
                        _beta=float, _wakefact=float, _nslice=int,
                        ZCuts=lambda v: _array(v),
                        _nturns=int, _phis=float,
                        _turnhistory=lambda v: _array(v),
                        _vbunch=lambda v: _array(v),
                        _vbeam_phasor=lambda v: _array(v, shape=(2,)),
                        _vbeam=lambda v: _array(v, shape=(2,)),
                        _vcav=lambda v: _array(v, shape=(2,)),
                        _vgen=lambda v: _array(v, shape=(2,)))

    def __init__(self, family_name: str, length: float, voltage: float,
                 frequency: float, ring: Lattice, qfactor: float,
                 rshunt: float, mode: Optional[BLMode] = BLMode.PHASOR,
                 **kwargs):
        r"""
        Parameters:
            ring:            Lattice object
            length:          Length of the cavity
            voltage:         Cavity voltage [V]
            frequency:       Cavity frequency [Hz]
            qfactor:         Q factor
            rshunt:          Shunt impedance, [:math:`\Omega`]

        Keyword Arguments:
            Nslice (int):       Number of slices per bunch. Default: 101
            Nturns (int):       Number of turn for the wake field. Default: 1
            ZCuts:              Limits for fixed slicing, default is adaptive
            NormFact (float):   Normalization factor
            mode (BLMode):  method for beam loading calculation BLMode.PHASOR
                (default) uses the phasor method, BLMode.WAKE uses the wake
                function. For high Q resonator, the phasor method should be
                used
        Returns:
            bl_elem (Element): beam loading element
        """
        kwargs.setdefault('PassMethod', self.default_pass[True])
        assert isinstance(mode, BLMode), \
            'Beam loading mode has to be an instance of BLMode'
        zcuts = kwargs.pop('ZCuts', None)
        phil = kwargs.pop('phil', 0)
        energy = ring.energy
        harmonic_number = numpy.round(frequency*ring.circumference/clight)
        self.Rshunt = rshunt
        self.Qfactor = qfactor
        self.NormFact = kwargs.pop('NormFact', 1.0)
        self.PhaseGain = kwargs.pop('PhaseGain', 1.0)
        self.VoltGain = kwargs.pop('VoltGain', 1.0)
        self._mode = int(mode)
        self._beta = ring.beta
        self._wakefact = - ring.circumference/(clight *
                                               ring.energy*ring.beta**3)
        self._nslice = kwargs.pop('Nslice', 101)
        self._nturns = kwargs.pop('Nturns', 1)
        self._turnhistory = None    # Defined here to avoid warning
        self._vbunch = None
        if zcuts is not None:
            self.ZCuts = zcuts
        super(BeamLoadingElement, self).__init__(family_name, length,
                                                 voltage, frequency,
                                                 harmonic_number,
                                                 energy, **kwargs)
        _, ts = get_timelag_fromU0(ring)
        self._phis = 2*numpy.pi*self.Frequency*(ts+self.TimeLag)/clight
        self._vbeam_phasor = numpy.zeros(2)
        self._vbeam = numpy.zeros(2)
        self._vgen = numpy.zeros(2)
        self._vcav = numpy.array([self.Voltage,
                                  numpy.pi/2-self._phis-phil])
        self.clear_history(ring=ring)

    def is_compatible(self, other):
        return False

    def clear_history(self, ring=None):
        if ring is None:
            current = 0.0
            nbunch = 1
        else:
            current = ring.beam_current
            nbunch = ring.nbunch
        tl = self._nturns*self._nslice*nbunch
        self._turnhistory = numpy.zeros((tl, 4), order='F')
        self._vbunch = numpy.zeros((nbunch, 2), order='F')
        self._init_bl_params(current)

    def _init_bl_params(self, current):
        theta = -self._vcav[1]+numpy.pi/2
        vb = 2*current*self.Rshunt
        a = self.Voltage*numpy.cos(theta-self._phis)
        b = self.Voltage*numpy.sin(theta-self._phis)-vb*numpy.cos(theta)
        psi = numpy.arcsin(b/numpy.sqrt(a**2+b**2))
        vgen = self.Voltage*numpy.cos(psi) + \
            vb*numpy.cos(psi)*numpy.sin(self._phis)
        self._vgen = numpy.array([vgen, psi])
        self._vbeam = numpy.array([2*current*self.Rshunt*numpy.cos(psi),
                                   numpy.pi-psi])

    @property
    def ResFrequency(self):
        """Resonator frequency"""
        return self.Frequency/(1-numpy.tan(self.Vgen[1])/(2*self.Qfactor))

    @property
    def Vbeam(self):
        """Beam phasor (amplitude, phase)"""
        return self._vbeam

    @property
    def Vbunch(self):
        """Bunch phasor (amplitude, phase)"""
        return self._vbunch

    @property
    def Vcav(self):
        """Cavity phasor (amplitude, phase)"""
        return self._vcav

    @property
    def Vgen(self):
        """Generator phasor (amplitude, phase)"""
        return self._vgen

    @Vgen.setter
    def Vgen(self, value):
        self._vgen = value

    @staticmethod
    def build_from_cav(cav: RFCavity, ring: Sequence,
                       qfactor: float, rshunt: float,
                       mode: Optional[BLMode] = BLMode.PHASOR,
                       **kwargs):
        r"""Function to build the BeamLoadingElement from a cavity
        the FamName, Length, Voltage, Frequency and HarmNumber are
        taken from the cavity element

        Parameters:
            ring:            Lattice object
            qfactor:         Q factor
            rshunt:          Shunt impedance, [:math:`\Omega`]

        Keyword Arguments:
            Nslice (int):       Number of slices per bunch. Default: 101
            Nturns (int):       Number of turn for the wake field. Default: 1
            ZCuts:              Limits for fixed slicing, default is adaptive
            NormFact (float):   Normalization factor
            mode (BLMode):  method for beam loading calculation BLMode.PHASOR
                (default) uses the phasor method, BLMode.WAKE uses the wake
                function. For high Q resonator, the phasor method should be
                used
        Returns:
            bl_elem (Element): beam loading element
        """
        _CAV_ATTRIBUTES = ['Length', 'Voltage', 'Frequency']
        _EXCL_ATTRIBUTES = ['PassMethod', 'Energy', 'HarmNumber']
        cav_attrs = dict(cav.items())
        family_name = cav_attrs.pop('FamName') + '_BL'
        [cav_attrs.pop(k) for k in _EXCL_ATTRIBUTES]
        cav_args = [cav_attrs.pop(k, getattr(cav, k)) for k in
                    _CAV_ATTRIBUTES]
        return BeamLoadingElement(family_name, *cav_args, ring,
                                  qfactor, rshunt, mode=mode,
                                  **cav_attrs, **kwargs)

    def __repr__(self):
        """Simplified __repr__ to avoid errors due to arguments
        not defined as attributes
        """
        attrs = dict((k, v) for (k, v) in self.items()
                     if not k.startswith('_'))
        return '{0}({1})'.format(self.__class__.__name__, attrs)
