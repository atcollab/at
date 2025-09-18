from __future__ import annotations

import warnings
from enum import IntEnum
from typing import Sequence, Optional, Union

import numpy

from ..constants import clight
from ..lattice import Lattice, AtWarning, AtError, RFCavity, Collective
from ..lattice import Refpts, uint32_refpts, make_copy
from ..lattice.elements.conversions import _array
from ..physics import get_timelag_fromU0


class BLMode(IntEnum):
    WAKE = 1
    PHASOR = 2


class CavityMode(IntEnum):
    """
    CavityMode.ACTIVE is an active feedback loop, acting on both
    detuning and generator voltages to compensate the beam induced
    voltage.

    CavityMode.PASSIVE takes a fixed cavity frequency and assumes
    zero generator voltage.

    CavityMode.PASSIVE_SETVOLTAGE is a feedback loop for a passive harmonic
    cavity. The frequency of the cavity is varied in order to maintain
    a voltage setpoint.
    """

    ACTIVE = 1
    PASSIVE = 2
    PASSIVE_SETVOLTAGE = 3


class FeedbackMode(IntEnum):
    ONETURN = 1
    WINDOW = 2


def add_beamloading(
    ring: Lattice,
    qfactor: Union[float, Sequence[float]],
    rshunt: Union[float, Sequence[float]],
    cavpts: Refpts = None,
    copy: Optional[bool] = False,
    **kwargs,
):
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
        blmode (BLMode):  method for beam loading calculation BLMode.PHASOR
            (default) uses the phasor method, BLMode.WAKE uses the wake
            function. For high Q resonator, the phasor method should be
            used
        copy:       If True, returns a shallow copy of ring with new
                    beam loading elements. Otherwise, modify ring in-place
        cavitymode (CavityMode):     Define PASSIVE or ACTIVE cavity
        buffersize (int):  Size of the history buffer for vbeam, vgen, vbunch
                (default 0)
    """

    @make_copy(copy)
    def apply(ring, cavpts, newelems):
        for ref, elem in zip(cavpts, newelems):
            ring[ref] = elem

    if cavpts is None:
        cavpts = ring.get_refpts(RFCavity)
    else:
        cavpts = uint32_refpts(cavpts, len(ring))
    qfactor = numpy.broadcast_to(qfactor, (len(cavpts),))
    rshunt = numpy.broadcast_to(rshunt, (len(cavpts),))
    new_elems = []
    for ref, qf, rs in zip(cavpts, qfactor, rshunt):
        cav = ring[ref]
        if not isinstance(cav, RFCavity):
            raise TypeError(
                "Beam loading can only be assigned" + "to an RFCavity element"
            )
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
        if not isinstance(bl, BeamLoadingElement):
            raise TypeError("Cannot remove beam loading: " +
                            "not a BeamLoadingElement")
        family_name = bl.FamName.replace("_BL", "")
        harm = numpy.round(bl.Frequency * ring.circumference / clight)
        new_elems.append(
            RFCavity(
                family_name,
                bl.Length,
                bl.Voltage,
                bl.Frequency,
                harm,
                bl.Energy,
                TimeLag=getattr(bl, "TimeLag", 0.0),
                PhaseLag=getattr(bl, "PhaseLag", 0.0),
            )
        )
    return apply(ring, cavpts, new_elems)


class BeamLoadingElement(RFCavity, Collective):
    """Class to generate a beamloading element, inherits from Element
    additional argument are ring, cavity, qfactor, rshunt
    """

    default_pass = {False: "DriftPass", True: "BeamLoadingCavityPass"}
    _conversions = dict(
        RFCavity._conversions,
        Rshunt=float,
        Qfactor=float,
        NormFact=float,
        PhaseGain=float,
        VoltGain=float,
        _blmode=int,
        _beta=float,
        _wakefact=float,
        _nslice=int,
        ZCuts=lambda v: _array(v),
        _cavitymode=int,
        _nturns=int,
        _phis=float,
        _buffersize=int,
        _vbeam_phasor=lambda v: _array(v, shape=(2,)),
        _vbeam=lambda v: _array(v, shape=(2,)),
        _vcav=lambda v: _array(v, shape=(2,)),
        _vgen=lambda v: _array(v, shape=(3,)),
    )

    def __init__(
        self,
        family_name: str,
        length: float,
        voltage: float,
        frequency: float,
        ring: Lattice,
        qfactor: float,
        rshunt: float,
        detune: Optional[float] = 0.0,
        blmode: Optional[BLMode] = BLMode.PHASOR,
        cavitymode: Optional[CavityMode] = CavityMode.ACTIVE,
        fbmode: Optional[FeedbackMode] = FeedbackMode.ONETURN,
        buffersize: Optional[int] = 0,
        **kwargs,
    ):
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
            detune [Hz] (float):     Define how much to detune the cavity from
                resonance in unints of Hz
            blmode (BLMode):    method for beam loading calculation
                BLMode.PHASOR (default) uses the phasor method, BLMode.WAKE
                uses the wake function. For high Q resonator, the phasor
                method should be used.
            cavitymode (CavityMode):  Is cavity ACTIVE (default), PASSIVE or
                PASSIVE_SETVOLTAGE (Passive with a voltage feedback).
                For PASSIVE_SETVOLTAGE, the voltage setpoint is specified with
                passive_voltage
            passive_voltage [V] (float): Voltage setpoint with the passive
                cavity with feedback.
            PhaseGain (float):  Used for cavity feedbacks. States the gain on
                the phase correction factor to be applied.
            VoltGain (float):  Used for cavity feedbacks. States the gain on
                the phase correction factor to be applied.
            buffersize (int):  Size of the history buffer for vbeam, vgen,
                vbunch (default 0)
            feedback_angle_offset:      Fixed detuning from optimal tuning
                angle [rad]. For a negative slope of the RF voltage at the
                synchronous position, the optimum detuning is negative.
                Applying a positive feedback_angle_offset will therefore
                reduce the detuning. The reverse is true for positive RF
                slope.
            ts (float):        The timelag of the synchronous particle in the
                full RF system [m]. If not specified, it will be calculated
                using get_timelag_fromU0. Defines the expected position of the
                beam to be used for the beam loading setpoints.
            fbmode (FeedbackMode): States the mode for the parameter
                calculation to be used for the loop. ONETURN (default) takes
                only the current turn, compared to WINDOW which takes a
                sliding window.
            windowlength (int): for WINDOW feedback mode, states the length
                [turns] for the sliding window. Must be smaller than
                buffersize.
            system_harmonic (float): Used to compute the nominal rf frequency
                for the given system. e.g. third of fourth harmonic of
                rf_frequency. If None, then will be computed to the nearest
                integer multiple of rf_frequency.
        Returns:
            bl_elem (Element): beam loading element
        """
        kwargs.setdefault("PassMethod", self.default_pass[True])
        if not isinstance(blmode, BLMode):
            raise TypeError("blmode mode has to be an " +
                            "instance of BLMode")
        if not isinstance(cavitymode, CavityMode):
            raise TypeError("cavitymode has to be an " +
                            "instance of CavityMode")

        zcuts = kwargs.pop("ZCuts", None)
        ts = kwargs.pop("ts", None)
        self.system_harmonic = kwargs.pop(
            "system_harmonic", int(numpy.round(frequency / ring.rf_frequency))
        )
        self.detune = detune

        check_frequency = numpy.abs(frequency -
                                    self.system_harmonic * ring.rf_frequency)

        if check_frequency > 1.0:  # 1 Hz is the limit for the float check

            error_string = (
                "Cavity must be an integer of rf_frequency, otherwise" +
                "the phi_s computation will be wrong. Please use the detune" +
                "argument when adding beamloading to a cavity that is an" +
                "integer harmonic."
            )
            raise AtError(error_string)
            
        self.circumference = ring.circumference
        self.bunch_spos = ring.bunch_spos
        energy = ring.energy
        harmonic_number = self.system_harmonic * ring.harmonic_number
        self.feedback_angle_offset = kwargs.pop("feedback_angle_offset", 0)
        self.Rshunt = rshunt
        self.Qfactor = qfactor
        self.NormFact = kwargs.pop("NormFact", 1.0)
        self.PhaseGain = kwargs.pop("PhaseGain", 1.0)
        self.VoltGain = kwargs.pop("VoltGain", 1.0)
        self._blmode = int(blmode)
        self._cavitymode = int(cavitymode)

        if self._cavitymode == 1:
            if not isinstance(fbmode, FeedbackMode):
                raise TypeError(
                    "For an active cavity, fbmode has to be defined and an " +
                    "instance of FeedbackMode"
                )
            self._fbmode = int(fbmode)
        else:
            self._fbmode = 0

        if self.detune == 0 and self._cavitymode == 3:
            raise AtError(
                "Cannot start passive cavity feedback from zero detuning." +
                "You must decide at the beginning which polarity you want." +
                "This problem arises because the loop does not know which " +
                "polarity to take!"
            )

        self._passive_vset = kwargs.pop("passive_voltage", 0.0)
        self._beta = ring.beta
        self._wakefact = (-ring.circumference /
                          (clight * ring.energy * ring.beta**3))
        self._nslice = kwargs.pop("Nslice", 101)
        self._nturns = kwargs.pop("Nturns", 1)
        self._nbunch = ring.nbunch
        self._turnhistory = None  # Defined here to avoid warning
        self._vbunch = None
        self._buffersize = buffersize
        self._windowlength = kwargs.pop("windowlength", 0)

        if self._windowlength > self._buffersize:
            raise ValueError("The windowlength must be smaller" +
                             "than the buffersize")
        self._vgen_buffer = numpy.zeros(1)
        self._vbeam_buffer = numpy.zeros(1)
        self._vbunch_buffer = numpy.zeros(1)
        if zcuts is not None:
            self.ZCuts = zcuts
        super(BeamLoadingElement, self).__init__(
            family_name, length, voltage, frequency, harmonic_number,
            energy, **kwargs)
        if ts is None:
            _, ts = get_timelag_fromU0(ring)
        self._ts = ts
        self._phis = (2 * numpy.pi * self.Frequency *
                      (self._ts + self.TimeLag) / clight)
        self._vbeam_phasor = numpy.zeros(2)
        self._vbeam = numpy.zeros(2)
        self._vgen = numpy.zeros(3)

        cavity_voltage = self.Voltage
        if self._cavitymode == 3:
            cavity_voltage = self._passive_vset

        self._vcav = numpy.array([cavity_voltage, numpy.pi / 2 - self._phis])
        self.clear_history(ring=ring)

    def is_compatible(self, other):
        return False

    def clear_history(self, ring=None):
        if ring is not None:
            self._nbunch = ring.nbunch
            current = ring.beam_current
            nbunch = ring.nbunch
            self._vbunch = numpy.zeros((nbunch, 2), order="F")
            self._init_bl_params(current)
        tl = self._nturns * self._nslice * self._nbunch
        self._turnhistory = numpy.zeros((tl, 4), order="F")
        if self._buffersize > 0:
            self._vgen_buffer = numpy.zeros((3, self._buffersize), order="F")
            self._vbeam_buffer = numpy.zeros((2, self._buffersize), order="F")
            self._vbunch_buffer = numpy.zeros(
                (self._nbunch, 2, self._buffersize), order="F"
            )

    def _init_bl_params(self, current):
        if (self._cavitymode == 1) and (current > 0.0):
            theta = -self._vcav[1] + numpy.pi / 2 #just self._phis
            vb = 2 * current * self.Rshunt
            a = self.Voltage * numpy.cos(theta - self._phis)
            b = (self.Voltage * numpy.sin(theta - self._phis) -
                 vb * numpy.cos(theta))
            psi = numpy.arcsin(b / numpy.sqrt(a**2 + b**2))
            if numpy.isnan(psi):
                psi = 0.0
                warning_string = (
                    "Unusual cavity configuration found." +
                    "Setting initial psi to 0 to avoid NaNs"
                )
                warnings.warn(AtWarning(warning_string))
            psi += self.feedback_angle_offset
            vgen = (self.Voltage * numpy.cos(psi) +
                    vb * numpy.cos(psi) * numpy.sin(self._phis))

        elif self._cavitymode == 2 or self._cavitymode == 3:
            vgen = 0
            psi = numpy.arctan(
                2 * self.Qfactor *
                (1 - self.Frequency / (self.Frequency + self.detune))
            )
        else:
            vgen = self.Voltage
            psi = 0

        self._vbeam = numpy.array(
            [2 * current * self.Rshunt * numpy.cos(psi), numpy.pi + psi]
        )

        self._vgen = numpy.array([vgen, psi - self._phis, psi])

        omr = self.ResFrequency * 2 * numpy.pi
        vbp_amp = 2 * current * self.Rshunt * numpy.cos(psi)
        vbp_phase = numpy.pi + psi
        vbp_complex = vbp_amp*numpy.cos(vbp_phase) + 1j*vbp_amp*numpy.sin(vbp_phase)
        
        dt = (self.circumference - self.bunch_spos[0])/clight    
        vbp_complex *= numpy.exp((1j*omr-omr/(2*self.Qfactor))*dt)
    
    
        self._vbeam_phasor = numpy.array(
            [numpy.abs(vbp_complex), numpy.angle(vbp_complex)]
        )
                


    @property
    def Buffersize(self):
        return self._buffersize

    @Buffersize.setter
    def Buffersize(self, value):
        self._buffersize = value
        self.clear_history()

    @property
    def Vgen_buffer(self):
        """Stored generator voltage data"""
        return self._vgen_buffer

    @property
    def Vbeam_buffer(self):
        """Stored beam induced voltage data"""
        return self._vbeam_buffer

    @property
    def Vbunch_buffer(self):
        """Stored bunch induced voltage data"""
        return numpy.moveaxis(self._vbunch_buffer, 0, -1)

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
        """Turn history of the slices center of mass"""
        return self._turnhistory

    @property
    def ResFrequency(self):
        """Resonator frequency"""
        return (self.Frequency /
                (1 - numpy.tan(self.Vgen[2]) / (2 * self.Qfactor)))

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
    def build_from_cav(
        cav: RFCavity,
        ring: Sequence,
        qfactor: float,
        rshunt: float,
        blmode: Optional[BLMode] = BLMode.PHASOR,
        cavitymode: Optional[CavityMode] = CavityMode.ACTIVE,
        buffersize: Optional[int] = 0,
        **kwargs,
    ):
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
            blmode (BLMode):  method for beam loading calculation BLMode.PHASOR
                (default) uses the phasor method, BLMode.WAKE uses the wake
                function. For high Q resonator, the phasor method should be
                used
            cavitymode (CavityMode):  type of beam loaded cavity ACTIVE
                (default) for a cavity with active compensation, or
                PASSIVE to only include the beam induced voltage
            buffersize (int):  Size of the history buffer for vbeam, vgen,
                vbunch (default 0)

        Returns:
            bl_elem (Element): beam loading element
        """
        _CAV_ATTRIBUTES = ["Length", "Voltage", "Frequency"]
        _EXCL_ATTRIBUTES = ["PassMethod", "Energy", "HarmNumber"]
        cav_attrs = dict(cav.items())
        family_name = cav_attrs.pop("FamName") + "_BL"
        [cav_attrs.pop(k) for k in _EXCL_ATTRIBUTES]
        cav_args = [cav_attrs.pop(k, getattr(cav, k)) for k in _CAV_ATTRIBUTES]
        if cavitymode == CavityMode.PASSIVE:
            if cav_args[1] != 0.0:
                warnings.warn(AtWarning("Setting Cavity Voltage to 0"))
            cav_args[1] = 0.0
        return BeamLoadingElement(
            family_name,
            *cav_args,
            ring,
            qfactor,
            rshunt,
            blmode=blmode,
            cavitymode=cavitymode,
            buffersize=buffersize,
            **cav_attrs,
            **kwargs,
        )

    def __repr__(self):
        """Simplified __repr__ to avoid errors due to arguments
        not defined as attributes
        """
        att = dict((k, v) for (k, v) in self.items() if not k.startswith("_"))
        return "{0}({1})".format(self.__class__.__name__, att)
