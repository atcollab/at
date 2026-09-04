from __future__ import annotations

import warnings
from enum import IntEnum
from collections.abc import Sequence

import numpy as np

from ..constants import clight
from ..lattice import Lattice, AtWarning, AtError, RFCavity, Collective
from ..lattice import Refpts, uint32_refpts, make_copy
from ..lattice.elements.conversions import _array
from ..physics import get_timelag_fromU0


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
    """
    FeedbackMode.PROP means the feedback is only using the most recent
    turn, and the delta is applied each turn multiplied by VoltGain or
    PhaseGain.

    FeedbackMode.PROP_INTEGRAL means
    
    """

    PROP = 1
    PROP_INTEGRAL = 2


def add_beamloading(
    ring: Lattice,
    qfactor: float | Sequence[float],
    rshunt: float | Sequence[float],
    cavitybeta: float | Sequence[float],
    cavpts: Refpts = None,
    copy: bool | None = False,
    **kwargs,
):
    r"""Function to add beam loading to a cavity element, the cavity
    element is changed to a beam loading element that combines the energy
    kick from both the cavity and the resonator.

    Parameters:
        ring:           Lattice object
        qfactor:        Unloaded Q. Scalar or array of float are accepted
        rshunt:         Unloaded Shunt impedance, [:math:`\Omega`] for longitudinal
                        Scalar or array of float are accepted
        cavitybeta:     Cavity coupling factor. 

    Keyword Arguments:
        cavpts (Refpts):    refpts of the cavity. If None (default) apply
                            to all cavity in the ring
        Nslice (int):       Number of slices per bunch. Default: 101
        Nturns (int):       Number of turn for the wake field. Default: 1
        ZCuts:              Limits for fixed slicing, default is adaptive
        NormFact (float):   Normalization factor
        copy:       If True, returns a shallow copy of ring with new
                    beam loading elements. Otherwise, modify ring in-place
        cavitymode (CavityMode):     Define PASSIVE or ACTIVE cavity
        buffersize (int):  Size of the history buffer for vbeam, vgen, vbunch
                (default 0)
    """

    @make_copy(copy)
    def apply(ring, cavpts, newelems):
        for ref, elem in zip(cavpts, newelems, strict=True):
            ring[ref] = elem

    if cavpts is None:
        cavpts = ring.get_refpts(RFCavity)
    else:
        cavpts = uint32_refpts(cavpts, len(ring))
    qfactor = np.broadcast_to(qfactor, (len(cavpts),))
    rshunt = np.broadcast_to(rshunt, (len(cavpts),))
    cavitybeta = np.broadcast_to(cavitybeta, (len(cavpts),))
    new_elems = []
    for ref, qf, rs, be in zip(cavpts, qfactor, rshunt, cavitybeta, strict=True):
        cav = ring[ref]
        if not isinstance(cav, RFCavity):
            raise TypeError(
                "Beam loading can only be assigned" + "to an RFCavity element"
            )
        new_elems.append(BeamLoadingElement.build_from_cav(cav, ring, qf, rs, be, **kwargs))
    return apply(ring, cavpts, new_elems)


def remove_beamloading(ring, cavpts: Refpts = None, copy: bool | None = False):
    """Function to remove beam loading from a cavity element, the beam
    loading element is changed to a beam loading element that
    combines the energy kick from both the cavity and the resonator.

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
        for ref, elem in zip(cavpts, newelems, strict=True):
            ring[ref] = elem

    if cavpts is None:
        cavpts = ring.get_refpts(BeamLoadingElement)
    else:
        cavpts = uint32_refpts(cavpts, len(ring))
    new_elems = []
    for ref in cavpts:
        bl = ring[ref]
        if not isinstance(bl, BeamLoadingElement):
            raise TypeError("Cannot remove beam loading: " + "not a BeamLoadingElement")
        family_name = bl.FamName.replace("_BL", "")
        harm = np.round(bl.Frequency * ring.circumference / clight)
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
    additional argument are ring, cavity, qfactor, rshunt.
    """

    default_pass = {False: "DriftPass", True: "BeamLoadingCavityPass"}
    _conversions = dict(
        RFCavity._conversions,
        Rshunt=float,
        Qfactor=float,
        CavityBeta=float,
        NormFact=float,
        TunerGain=float,
        Gain=lambda v: _array(v, shape=(2,)),
        IIR_cutoff=float,
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
        delay=int,
        every=int,
        samplenum=int,
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
        cavitybeta: float,
        detune: float | None = 0.0,
        cavitymode: CavityMode | None = CavityMode.ACTIVE,
        fbmode: FeedbackMode | None = FeedbackMode.PROP,
        IIR_cutoff: float = 0.0,
        **kwargs,
    ):
        r"""
        Parameters:
            ring:            Lattice object
            length:          Length of the cavity
            voltage:         Cavity voltage [V]
            frequency:       Cavity frequency [Hz]
            qfactor:         Unloaded Q factor
            rshunt:          Unloaded Shunt impedance, [:math:`\Omega`].
            cavitybeta:      Cavity coupling factor
            
        Keyword Arguments:
            Nslice (int):       Number of slices per bunch. Default: 101
            Nturns (int):       Number of turn for the wake field. Default: 1
            ZCuts:              Limits for fixed slicing, default is adaptive
            NormFact (float):   Normalization factor
            detune [Hz] (float):     Define how much to detune the cavity from
                resonance in unints of Hz
            cavitymode (CavityMode):  Is cavity ACTIVE (default), PASSIVE or
                PASSIVE_SETVOLTAGE (Passive with a voltage feedback).
                For PASSIVE_SETVOLTAGE, the voltage setpoint is specified with
                passive_voltage
            passive_voltage [V] (float): Voltage setpoint with the passive
                cavity with feedback.
            Gain ([float, float]):  Used for cavity feedbacks. If FBMode is PROP
                then Gain[0] is the amplitude gain, and Gain[1] is the phase gain.
                If FBMode is PROP_INTEGRAL, then Gain[0] is Prop gain and Gain[1]
                is Integral gain. 
            buffersize (int):  Size of the history buffer for vbeam, vgen,
                vbunch (default 0)
            TunerOffset:      Fixed detuning from optimal tuning
                angle [rad]. For a negative slope of the RF voltage at the
                synchronous position, the optimum detuning is negative.
                Applying a positive TunerOffset will therefore
                reduce the detuning. The reverse is true for positive RF
                slope.
            TunerGain (float): Used for detuning of the cavities. States the gain
                of the correction factor to be applied.
                
            TunerAveragingPeriod (int): default 1
            
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
            IIR_cutoff: cutoff frequency of the IIR filter [Hz]. If 0,
                a cutoff frequency of infinity is assumed. 
            Delay: Loop delay. If FBMode is PROP, then the units are in turns,
                if FBMode is PROP_INTEGRAL, then the units are in buckets.
            Every: Every what
            samplenum: Sample whhatt?
        Returns:
            bl_elem (Element): beam loading element
        """
        kwargs.setdefault("PassMethod", self.default_pass[True])

        # Initialise ring parameters
        self.circumference = ring.circumference
        self.bunch_spos = ring.bunch_spos
        energy = ring.energy
        self.system_harmonic = kwargs.pop(
            "system_harmonic", int(np.round(frequency / ring.rf_frequency))
        )
        harmonic_number = self.system_harmonic * ring.harmonic_number #cavity harmonic number
        self.ring_harmonic_number = ring.harmonic_number #ring harmonic number (nbuckets)        
        self._nbunch = ring.nbunch
        self._beta = ring.beta #particle velocity (NOT cavity beta)
        self._cavitymode = int(cavitymode)        
        
        zcuts = kwargs.pop("ZCuts", None) #fixed or adaptive slicing
        if zcuts is not None:
            self.ZCuts = zcuts
        

        # Initialise resonator parameters
        self.CavityBeta = cavitybeta
        self.Rshunt = rshunt/(1+self.CavityBeta)
        self.Rshunt_unloaded = rshunt
        self.Qfactor = qfactor/(1+self.CavityBeta)
        self.Qfactor_unloaded = qfactor

        # Initialise wake computation parameters
        self._wakefact = -ring.circumference / (clight * ring.energy * ring.beta**3)
        self._nslice = kwargs.pop("Nslice", 101)
        self._nturns = kwargs.pop("Nturns", 1)
        self._turnhistory = None  # Defined here to avoid warning
        self._vbunch = None
        self.NormFact = kwargs.pop("NormFact", 1.0)


        # Initialise tuner parameters
        self.TunerGain = kwargs.pop("TunerGain", 0.01)
        self.TunerOffset = kwargs.pop("TunerOffset", 0)
        self.TunerAveragingPeriod = kwargs.pop("TunerAveragingPeriod", 1)
        if self.TunerAveragingPeriod <1:
            raise AtError('Tuner Averaging Period must be >=1)')
        self._TunerParams = np.array([0.0, 0.0]) #TunerCount and TunerDiff
            
        # Initialise common regulator parameters
        self.Gain = kwargs.pop("Gain", [1e-3,1e-3])
        self.delay = kwargs.pop("delay", 1)
        if self.delay <= 0:
            raise AtError('Attribute delay must be >= 1')  


        # Initialise FBMode=PROP buffers
        self.VoltDelay = np.ones(self.delay)
        self.PhaseDelay = np.ones(self.delay)

        # Initialise FBMode=PROP_INTEGRAL buffers                  
        self.every = kwargs.pop("every", 1)
        self.samplenum = kwargs.pop("samplenum", 1)
        self.samplelist_length = int(np.ceil(ring.harmonic_number/self.every))
        self.recordsize = int(np.ceil(self.delay / self.every))

        self.cutoff = kwargs.pop("IIRcutoff",0.0)
        self.FF = kwargs.pop("FF", 1) #bool
        self._IIRcoef = np.zeros(1)
        self._IIRout = np.zeros(2)
        self._FFconst = np.zeros(2)        
        self._I_record = np.zeros(2)
        self.OpenLoop = int(kwargs.pop("OpenLoop", 0))

        # Initlise CavityMode=Passive parameters
        self.detune = detune
        
        # Initialise CavityMode=Passive_SetVoltage
        self._passive_vset = kwargs.pop("passive_voltage", 0.0)  
              
              
        ####################################
        ### Next we perform all the checks #
        ####################################
        if not isinstance(cavitymode, CavityMode):
            error_string = ("cavitymode has to be an "
                            "instance of CavityMode")
            raise TypeError(error_string)

        if self._cavitymode == 1:
            if not isinstance(fbmode, FeedbackMode):
                err_string = (
                    "For an active cavity, fbmode has to be defined as an "
                    "instance of FeedbackMode"
                )
                raise TypeError(err_string)
            self._fbmode = int(fbmode)
        else:
            self._fbmode = 0


        # Here we make checks relating to the passive cavity with FB
        if self.detune == 0 and self._cavitymode == 3:
            err_string = (
                "Cannot start passive cavity feedback from zero detuning."
                "You must decide at the beginning which polarity you want."
                "This problem arises because the loop does not know which "
                "polarity to take!"
            )
            raise AtError(err_string)

        # Check that the cavity frequency provided is system_harmonic * ring.rf_frequency
        check_frequency = np.abs(frequency - self.system_harmonic * ring.rf_frequency)
        if check_frequency > 1.0:  # 1 Hz is the limit for the float check
            error_string = (
                "Cavity frequency must be system_harmonic*rf_frequency, otherwise"
                "the phi_s computation will be wrong. Please use the detune"
                "argument but keep the resonant frequency on resonance."
            )
            raise AtError(error_string)


        # REMOVE WINDOWLENGTH 
        # buffer size and windowlength verification
        self._windowlength = kwargs.pop("windowlength", 0)
        self._buffersize = kwargs.pop("buffersize", 0) #is it still needed?
        if self._windowlength > self._buffersize:
            err_string = "The windowlength must be smaller than the buffersize"
            raise ValueError(err_string)


        # Initlise the buffers before super. Redefined later.
        self._vgen_buffer = np.zeros(1)
        self._vbeam_buffer = np.zeros(1)
        self._vbunch_buffer = np.zeros(1)
        


        super().__init__(
            family_name, length, voltage, frequency, harmonic_number, energy, **kwargs
        )
        
        # ts says where the beam will be. Is either provided or computed.
        ts = kwargs.pop("ts", None)
        if ts is None:
            _, ts = get_timelag_fromU0(ring)
        self._ts = ts

        # this one has to come after the super in order to properly initialise
        cavity_voltage = self.Voltage
        if self._cavitymode == 3:
            cavity_voltage = self._passive_vset

        # phis must be defined in such a way that the beam loading compensation
        # works for any value of TimeLag. 
        self._phis = 2 * np.pi * self.Frequency * (-self._ts - self.TimeLag) / clight

        # The below is needed because atan2 returns phases between -pi and pi
        # this prevents a setpoint being provided that is impossible
        # to achieve
        while self._phis < -np.pi:
            self._phis += 2 * np.pi
        while self._phis > np.pi:
            self._phis -= 2 * np.pi

        # Initlise some parameters needed
        self._vbeam_phasor = np.zeros(2)
        self._vbeam = np.zeros(2)
        self._vgen = np.zeros(3)

        # Here we define the cavity setpoints, finally
        self._vcav = np.array([cavity_voltage, self._phis])
        self.VoltDelay *= self._vcav[0] # PROP 
        self.PhaseDelay *= self._vcav[1] # PROP
        self.clear_history(ring=ring)

        # these big chaps are all for the memory assignment in C
        self._Ig2Vg_vec = np.zeros(ring.harmonic_number*2)
        self._Ig2Vg_tmp = np.zeros(ring.harmonic_number*2)
        self._ig_phasor = np.zeros(ring.harmonic_number*2)
        self._ig_phasor_record = np.zeros(ring.harmonic_number*2)
        self._dot_output = np.zeros(ring.harmonic_number*2)
        self._generator_phasor_record = np.zeros(ring.harmonic_number*2)
        self._beam_phasor_record = np.zeros(ring.harmonic_number*2)                
        self._cavity_phasor_record = np.zeros(ring.harmonic_number*2)        
        
        self._Ig2Vg_mat = np.zeros(ring.harmonic_number**2 * 2)
        self._vc_previous = np.zeros(self.samplenum*2)
        self._diff_record = np.zeros(self.recordsize*2)
        self._samplelist = np.zeros(self.samplelist_length)

        self._vc_list = np.zeros((ring.harmonic_number + self.samplenum)*2)        

        
    def is_compatible(self, other):
        return False

    def clear_history(self, ring=None):
        if ring is not None:
            current = ring.beam_current
            self._vbunch = np.zeros((self.ring_harmonic_number, 2), order="F")
            self._init_bl_params(current)

        tl = self._nturns * self._nslice * self._nbunch
        self._turnhistory = np.zeros((tl, 4), order="F")


        if self._buffersize > 0:
            self._vgen_buffer = np.zeros((3, self._buffersize), order="F")
            self._vbeam_buffer = np.zeros((2, self._buffersize), order="F")
            self._vbunch_buffer = np.zeros(
                (self.ring_harmonic_number, 2, self._buffersize), order="F"
            )
            
    def _set_optimum_detuning(self, current):
        psi = np.arctan(- 2 * current * self.Rshunt / self.Voltage * np.cos(self._phis))
        return psi
       
    def _compute_generator_parameters(self, current, psi):
        """
         Thanks to MBTRACK2 for the formula. Taken from 
            [1] Wilson, P. B. (1994). Fundamental-mode rf design in e+ e− storage ring
                factories. In Frontiers of Particle Beams: Factories with e+ e-Rings
                (pp. 293-311). Springer, Berlin, Heidelberg.
        """
        # Generator power [W] - Eq. (4.1.2) [1] corrected with factor
        # (1+beta)**2 instead of (1+beta**2)
        
        Pg = self.Voltage**2 * (1 + self.CavityBeta)**2 / (
            2 * self.Rshunt_unloaded * 4 * self.CavityBeta * np.cos(psi)**2) * (
                (-np.sin(self._phis) + 2 * current * self.Rshunt_unloaded /
                 (self.Voltage * (1 + self.CavityBeta)) * np.cos(psi)**2)**2 +
                (np.cos(self._phis) + 2 * current * self.Rshunt_unloaded /
                 (self.Voltage *
                  (1 + self.CavityBeta)) * np.cos(psi) * np.sin(psi))**2)

                  
        # Generator voltage at resonance [V] - Eq. (3.2.2) [1]
        Vgr = 2 * self.CavityBeta**(1 / 2) / (1 + self.CavityBeta) * (
            2 * self.Rshunt_unloaded * Pg)**(1 / 2)

        # Generator voltage [V]
        Vg = Vgr * np.cos(psi)
            
        vbr = 2 * current * self.Rshunt   
                 
        theta_g = np.arctan2(
            (self.Voltage * np.cos(self._phis) +
             vbr * np.cos(psi) * np.sin(psi)),
            (-self.Voltage * np.sin(self._phis) +
             vbr * np.cos(psi)**2)) - np.pi/2
             
        return np.array([Vg, theta_g])
            
            
    def _init_bl_params(self, current):
        if (self._cavitymode == 1) and (current > 0.0):
            psi = self._set_optimum_detuning(current)
            if np.isnan(psi):
                psi = 0.0
                warning_string = (
                    "Unusual cavity configuration found."
                    "Setting initial psi to 0 to avoid NaNs"
                )
                warnings.warn(AtWarning(warning_string), stacklevel=2)
            psi += self.TunerOffset
            vgen, theta_g = self._compute_generator_parameters(current, psi)
            
        elif self._cavitymode in {2, 3}:
            vgen = 0
            psi = np.arctan(
                2 * self.Qfactor * (1 - self.Frequency / (self.Frequency + self.detune))
            )
        else:
            vgen = self.Voltage
            psi = 0

        self._vbeam = np.array([2 * current * self.Rshunt * np.cos(psi), np.pi + psi])

        self._vgen = np.array([vgen, theta_g, psi])

        vbp_amp = 2 * current * self.Rshunt * np.cos(psi)
        vbp_phase = np.pi + psi
        vbp_complex = vbp_amp * np.cos(vbp_phase) + 1j * vbp_amp * np.sin(vbp_phase)

        self._vbeam_phasor = np.array([np.abs(vbp_complex), np.angle(vbp_complex)])

    @property
    def Buffersize(self):
        return self._buffersize

    @Buffersize.setter
    def Buffersize(self, value):
        self._buffersize = value
        self.clear_history()

    @property
    def Vgen_buffer(self):
        """Stored generator voltage data."""
        return self._vgen_buffer

    @property
    def Vbeam_buffer(self):
        """Stored beam induced voltage data."""
        return self._vbeam_buffer

    @property
    def Vbunch_buffer(self):
        """Stored bunch induced voltage data."""
        return np.moveaxis(self._vbunch_buffer, 0, -1)

    @property
    def Nslice(self):
        """Number of slices per bunch."""
        return self._nslice

    @Nslice.setter
    def Nslice(self, nslice):
        self._nslice = nslice
        self.clear_history()

    @property
    def Nturns(self):
        """Number of turn for the wake field."""
        return self._nturns

    @Nturns.setter
    def Nturns(self, nturns):
        self._nturns = nturns
        self.clear_history()

    @property
    def TurnHistory(self):
        """Turn history of the slices center of mass."""
        return self._turnhistory

    @property
    def ResFrequency(self):
        """Resonator frequency."""
        delta = (
            self.Frequency * np.tan(self.Vgen[2]) / self.Qfactor
        ) ** 2 + 4 * self.Frequency**2
        freqres = (
            self.Frequency * np.tan(self.Vgen[2]) / self.Qfactor + np.sqrt(delta)
        ) / 2

        return freqres

    @property
    def Vbeam(self):
        """Beam phasor (amplitude, phase)."""
        return self._vbeam

    @property
    def Vbunch(self):
        """Bunch phasor (amplitude, phase)."""
        return self._vbunch

    @property
    def Vcav(self):
        """Cavity phasor (amplitude, phase)."""
        return self._vcav

    @property
    def Vgen(self):
        """Generator phasor (amplitude, phase)."""
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
        cavitybeta: float,
        cavitymode: CavityMode | None = CavityMode.ACTIVE,
        buffersize: int | None = 0,
        **kwargs,
    ):
        r"""Function to build the BeamLoadingElement from a cavity
        the FamName, Length, Voltage, Frequency and HarmNumber are
        taken from the cavity element.

        Parameters:
            ring:            Lattice object
            qfactor:         Q factor
            rshunt:          Shunt impedance, [:math:`\Omega`]

        Keyword Arguments:
            Nslice (int):       Number of slices per bunch. Default: 101
            Nturns (int):       Number of turn for the wake field. Default: 1
            ZCuts:              Limits for fixed slicing, default is adaptive
            NormFact (float):   Normalization factor
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
                warnings.warn(AtWarning("Setting Cavity Voltage to 0"), stacklevel=2)
            cav_args[1] = 0.0
        return BeamLoadingElement(
            family_name,
            *cav_args,
            ring,
            qfactor,
            rshunt,
            cavitybeta,
            cavitymode=cavitymode,
            buffersize=buffersize,
            **cav_attrs,
            **kwargs,
        )

    def __repr__(self):
        """Simplified __repr__ to avoid errors due to arguments
        not defined as attributes.
        """
        att = {k: v for (k, v) in self.items() if not k.startswith("_")}
        return f"{self.__class__.__name__}({att})"
