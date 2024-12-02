from __future__ import annotations

from datetime import datetime
from enum import IntEnum

import numpy as np

from .elements import Element, _array
from .utils import AtError


class ACMode(IntEnum):
    """Class to define the excitation types."""

    SINE = 0
    WHITENOISE = 1
    ARBITRARY = 2


class VariableMultipole(Element):
    """Class to generate an AT variable thin multipole element."""

    _BUILD_ATTRIBUTES = Element._BUILD_ATTRIBUTES
    _conversions = dict(
        Element._conversions,
        Mode=int,
        AmplitudeA=_array,
        AmplitudeB=_array,
        FrequencyA=float,
        FrequencyB=float,
        PhaseA=float,
        PhaseB=float,
        Seed=int,
        NSamplesA=int,
        NSamplesB=int,
        FuncA=_array,
        FuncB=_array,
        Ramps=_array,
        Periodic=bool,
    )

    def __init__(self, family_name: str, **kwargs: dict[str, any]):
        # noinspection PyUnresolvedReferences
        r"""
        Create VariableMultipole.

        This function creates a thin multipole of order and type defined by the
        amplitude components; the polynoms PolynomA and PolynomB are calculated
        on every turn depending on the chosen mode.

        Keep in mind that this element varies on every turn, therefore, any ring
        containing a variable element may change after tracking n turns.

        There are three different modes implemented: SINE, WHITENOISE and ARBITRARY.

        The SINE mode requires amplitude, frequency and phase of at least one of the
        two polynoms A or B. The j-th component of the polynom on the n-th turn
        is given by:
          amplitude_j*sin[ 2\pi*frequency*(nth_turn*T0 + c\tau_k) + phase],
        where T0 is the revolution period of the ideal ring, and c\tau_k is the delay
        of the particle i.e. the sixth coordinate.
        The following is an example of the SINE mode of an skew quad:
            eleskew = at.VariableMultipole('VAR_SKEW',
                AmplitudeA=[0,skewa2],FrequencyA=freqA,PhaseA=phaseA)
        The WHITENOISE mode requires the amplitude.
        THe ARBITRARY mode requires the Amplitude


        Parameters:
            family_name(str):    Element name
            mode(ACMode): defines the evaluation grid. Default ACMode.SINE
              * :py:attr:`.ACMode.SINE`: sine function
              * :py:attr:`.ACMode.WHITENOISE`: gaussian white noise
              * :py:attr:`.ACMode.ARBITRARY`: user defined turn-by-turn kick list
        Keyword Arguments:
            AmplitudeA(list,float): Amplitude of the excitation for PolynomA.
                Default None
            AmplitudeB(list,float): Amplitude of the excitation for PolynomB.
                Default None
            FrequencyA(float): Frequency of the sine excitation for PolynomA
            FrequencyB(float): Frequency of the sine excitation for PolynomB
            PhaseA(float): Phase of the sine excitation for PolynomA. Default 0
            PhaseB(float): Phase of the sine excitation for PolynomB. Default 0
            MaxOrder(int): Order of the multipole for scalar amplitude. Default 0
                It is overwritten with the length of Amplitude(A,B)
            Seed(int): Seed of the random number generator for white
                       noise excitation. Default datetime.now()
            FuncA(list): User defined tbt kick list for PolynomA
            FuncB(list): User defined tbt kick list for PolynomB
            Periodic(bool): If True (default) the user defined kick is repeated
            Ramps(list): Vector (t0, t1, t2, t3) in turn number to define
                        the ramping of the excitation

              * ``t<t0``: excitation amplitude is zero
              * ``t0<t<t1``: excitation amplitude is linearly ramped up
              * ``t1<t<t2``: excitation amplitude is constant
              * ``t2<t<t3``: excitation amplitude is linearly ramped down
              * ``t3<t``: excitation amplitude is zero
        Raises:
            AtError if none of AmplitudeA or AmplitudeB is passed.
            AtError if ramp is not vector of length 4 when using Ramps
            AtError when Frequency is not defined if using Mode ``ACMode.SINE``
            AtError when Funct is not defined if using Mode ``ACMode.ARBITRARY``


        Examples:

            >>> acmpole = at.VariableMultipole('ACMPOLE', AmplitudeB=amp, FrequencyB=frequency)
            >>> acmpole = at.VariableMultipole('ACMPOLE', AmplitudeB=amp, mode=at.ACMode.WHITENOISE)
            >>> acmpole = at.VariableMultipole('ACMPOLE', AmplitudeB=amp, FuncB=fun, mode=at.ACMode.ARBITRARY)

        .. note::

            * If no parameters are given it will be initialized as IdentityPass.
            * For all excitation modes at least one amplitudes (A or B) has
              to be provided.
              The default excitation is ``ACMode.SINE``
            * For ``mode=ACMode.SINE`` the ``Frequency(A,B)`` corresponding to the
              ``Amplitude(A,B)`` has to be provided.
            * For ``mode=ACMode.ARBITRARY`` the ``Func(A,B)`` corresponding to the
              ``Amplitude(A,B)`` has to be provided.
        """
        if len(kwargs) > 0:
            self.FamName = family_name
            if "AmplitudeA" not in kwargs and "AmplitudeB" not in kwargs:
                raise AtError("Please provide at least one amplitude for A or B")
            # start setting up Amplitudes and modes
            # fist modes are called differently
            modepyatinput = kwargs.pop("mode", ACMode.SINE)
            modefromdict = kwargs.pop("Mode", None)
            if modefromdict is not None:
                mode = modefromdict  # matlab class
            else:
                mode = modepyatinput
            self.Mode = int(mode)
            # MaxOrder is later overwritten by lenth of amplitudes
            self.MaxOrder = kwargs.get("MaxOrder", 0)
            amplitudea = None
            amplitudeb = None
            if "AmplitudeA" in kwargs:
                amplitudea = kwargs.pop("AmplitudeA", None)
                amplitudea = self._set_params(amplitudea, mode, "A", **kwargs)
            if "AmplitudeB" in kwargs:
                amplitudeb = kwargs.pop("AmplitudeB", None)
                amplitudeb = self._set_params(amplitudeb, mode, "B", **kwargs)
            if amplitudea is not None:
                self.AmplitudeA = amplitudea
            if amplitudeb is not None:
                self.AmplitudeB = amplitudeb
            # end setting up Amplitudes and modes
            kwargs.setdefault("PassMethod", "VariableThinMPolePass")
            # this overwrites MaxOrder
            self._setmaxorder(amplitudea, amplitudeb)
            if mode == ACMode.WHITENOISE and "Seed" not in kwargs:
                kwargs.update({"Seed": datetime.now().timestamp()})
            self.Periodic = kwargs.pop("Periodic", True)
            self.PolynomA = kwargs.pop("PolynomA", np.zeros(self.MaxOrder + 1))
            self.PolynomB = kwargs.pop("PolynomB", np.zeros(self.MaxOrder + 1))
            # check ramps
            ramps = kwargs.pop("Ramps", None)
            if ramps is not None:
                if len(ramps) != 4:
                    raise AtError("Ramps has to be a vector with 4 elements")
                self.Ramps = ramps
        # fill in super class
        super().__init__(family_name, **kwargs)

    def _setmaxorder(self, ampa: np.ndarray, ampb: np.ndarray):
        mxa, mxb = 0, 0
        if ampa is not None:
            mxa = np.max(np.append(np.nonzero(ampa), 0))
        if ampb is not None:
            mxb = np.max(np.append(np.nonzero(ampb), 0))
        self.MaxOrder = max(mxa, mxb)
        if ampa is not None:
            delta = self.MaxOrder - len(ampa)
            if delta > 0:
                ampa = np.pad(ampa, (0, delta))
            self.AmplitudeA = ampa
        if ampb is not None:
            delta = self.MaxOrder + 1 - len(ampb)
            if delta > 0:
                ampb = np.pad(ampb, (0, delta))
            self.AmplitudeB = ampb

    def _set_params(
        self, amplitude: int or str, mode, a_b: str, **kwargs: dict[str, any]
    ):
        if amplitude is not None:
            if np.isscalar(amplitude):
                amp = np.zeros(self.MaxOrder)
                amplitude = np.append(amp, amplitude)
            if mode == ACMode.SINE:
                self._set_sine(a_b, **kwargs)
            if mode == ACMode.ARBITRARY:
                self._set_arb(a_b, **kwargs)
        return amplitude

    def _set_sine(self, a_b: str, **kwargs: dict[str, any]):
        frequency = kwargs.pop("Frequency" + a_b, None)
        phase = kwargs.pop("Phase" + a_b, 0)
        if frequency is None:
            raise AtError("Please provide a value for Frequency" + a_b)
        setattr(self, "Frequency" + a_b, frequency)
        setattr(self, "Phase" + a_b, phase)

    def _set_arb(self, a_b: str, **kwargs: dict[str, any]):
        func = kwargs.pop("Func" + a_b, None)
        nsamp = len(func)
        if func is None:
            raise AtError("Please provide a value for Func" + a_b)
        setattr(self, "Func" + a_b, func)
        setattr(self, "NSamples" + a_b, nsamp)
