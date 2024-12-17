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


def _array(value, shape=(-1,), dtype=np.float64):
    # Ensure proper ordering(F) and alignment(A) for "C" access in integrators
    return np.require(value, dtype=dtype, requirements=["F", "A"]).reshape(
        shape, order="F"
    )


class VariableMultipole(Element):
    """Class to generate an AT variable thin multipole element."""

    _BUILD_ATTRIBUTES = Element._BUILD_ATTRIBUTES + ["Mode"]
    _conversions = dict(
        Element._conversions,
        amplitudea=_array,
        amplitudeb=_array,
        frequencya=float,
        frequencyb=float,
        phasea=float,
        phaseb=float,
        seed=int,
        seedinitstate=int,
        nsamplesa=int,
        nsamplesb=int,
        funca=_array,
        funcb=_array,
        ramps=_array,
        periodic=bool,
    )

    def __init__(self, family_name: str, mode: int, **kwargs: dict[str, any]):
        r"""
        Create VariableMultipole.

        This function creates a thin multipole any order (1 to k) and type (Normal
        or Skew) defined by the amplitude A or B components; the polynoms PolynomA
        and PolynomB are calculated on every turn depending on the chosen mode. All
        modes could be ramped.

        Keep in mind that this element varies on every turn, therefore, any ring
        containing a variable element may change after tracking n turns.

        There are three different modes implemented: SINE, WHITENOISE and ARBITRARY.
        By default all are periodic.

        The SINE mode requires amplitude, frequency and phase of at least one of the
        two polynoms A or B. The j-th component of the kth order polynom on the
        n-th turn is given by:
          amplitude_j*sin[ 2\pi*frequency*(nth_turn*T0 + \tau_p) + phase],
        where T0 is the revolution period of the ideal ring, and \tau_p is the delay
        of the pth particle i.e. the sixth coordinate over the speed of light. Also,
        note that the position of the element on the ring has no effect, the phase
        should be used to add any delay due to the position along s.
        The following is an example of the SINE mode of an skew quad:
            eleskew = at.VariableMultipole('VAR_SKEW',ACMode.SINE,
                AmplitudeA=[0,skewa2],FrequencyA=freqA,PhaseA=phaseA)

        The WHITENOISE mode requires the amplitude.
        THe ARBITRARY mode requires the amplitude


        Parameters:
            family_name(str):    Element name
            mode(ACMode): defines the evaluation grid. Default ACMode.SINE
              * :py:attr:`.ACMode.SINE`: sine function
              * :py:attr:`.ACMode.WHITENOISE`: gaussian white noise
              * :py:attr:`.ACMode.ARBITRARY`: user defined turn-by-turn kick list
        Keyword Arguments:
            amplitudeA(list,float): Amplitude of the excitation for PolynomA.
                Default None
            amplitudeB(list,float): Amplitude of the excitation for PolynomB.
                Default None
            frequencyA(float): Frequency of the sine excitation for PolynomA
            frequencyB(float): Frequency of the sine excitation for PolynomB
            phaseA(float): Phase of the sine excitation for PolynomA. Default 0 rad
            phaseB(float): Phase of the sine excitation for PolynomB. Default 0 rad
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

            >>> acmpole = at.VariableMultipole('ACMPOLE', at.ACMode.SINE, amplitudeb=amp, frequencyb=frequency)
            >>> acmpole = at.VariableMultipole('ACMPOLE', at.ACMode.WHITENOISE, amplitudeb=amp)
            >>> acmpole = at.VariableMultipole('ACMPOLE', at.ACMode.ARBITRARY, amplitudeb=amp, funcb=fun)

            * For ``mode=ACMode.ARBITRARY`` the ``Func(A,B)`` corresponding to the
              ``Amplitude(A,B)`` has to be provided.
        """
        kwargs["Mode"] = kwargs.get("Mode", mode)
        kwargs.setdefault("PassMethod", "VariableThinMPolePass")
        if len(kwargs) > 2:
            amp_aux = self._check_amplitudes(**kwargs)
            for k, v in amp_aux.items():
                if v is not None:
                    kwargs["Amplitude" + k] = self._set_amplitude(v)
                    self._check_mode(mode, k, **kwargs)
            maxorder = self._getmaxorder(amp_aux["A"], amp_aux["B"])
            kwargs["MaxOrder"] = kwargs.get("MaxOrder", maxorder)
            for key in amp_aux.keys():
                kwargs["Polynom" + key] = kwargs.get(
                    "Polynom" + key, np.zeros(maxorder + 1)
                )
            ramps = self._check_ramp(**kwargs)
            if ramps is not None:
                kwargs["Ramps"] = ramps
        super().__init__(family_name, **kwargs)

    def _check_amplitudes(self, **kwargs: dict[str, any]):
        amp_aux = {"A": None, "B": None}
        all_amplitudes_are_none = True
        for key in amp_aux.keys():
            if "Amplitude" + key in kwargs and "amplitude" + key in kwargs:
                raise AtError("Duplicated amplitude " + key + "parameters.")
            lower_case_kwargs = {k.lower(): v for k, v in kwargs.items()}
            amp_aux[key] = lower_case_kwargs.pop("amplitude" + key.lower(), None)
            if amp_aux[key] is not None:
                all_amplitudes_are_none = False
        if all_amplitudes_are_none:
            raise AtError("Please provide at least one amplitude for A or B")
        return amp_aux

    def _set_amplitude(self, amplitude: float or _array or None):
        if np.isscalar(amplitude):
            amplitude = [amplitude]
        amplitude = np.asarray(amplitude)
        return amplitude

    def _getmaxorder(self, ampa: np.ndarray, ampb: np.ndarray):
        mxa, mxb = 0, 0
        if ampa is not None:
            mxa = np.max(np.append(np.nonzero(ampa), 0))
        if ampb is not None:
            mxb = np.max(np.append(np.nonzero(ampb), 0))
        return max(mxa, mxb)

    def _check_mode(self, mode, a_b: str, **kwargs: dict[str, any]):
        if mode == ACMode.WHITENOISE:
            if "seed" not in kwargs:
                kwargs["Seed"] = kwargs.get("Seed", datetime.now().timestamp())
            if "seedinitstate" not in kwargs:
                kwargs["Seedinitstate"] = kwargs.get("Seedinitstate", 1)
            if 'Mean' not in kwargs:
                kwargs['Mean'] = kwargs.get("Mean", 0)
            if 'Std' not in kwargs:
                kwargs['Std'] = kwargs.get("Std", 1)
        if mode == ACMode.SINE:
            self._check_sine(a_b, **kwargs)
        if mode == ACMode.ARBITRARY:
            self._set_arb(a_b, **kwargs)
            kwargs["Periodic"] = kwargs.get("Periodic", True)

    def _check_sine(self, a_b: str, **kwargs: dict[str, any]):
        frequency = kwargs.pop("Frequency" + a_b, None)
        if frequency is None:
            raise AtError("Please provide a value for Frequency" + a_b)
        kwargs["Phase" + a_b] = kwargs.get("Phase" + a_b, 0)

    def _set_arb(self, a_b: str, **kwargs: dict[str, any]):
        func = kwargs.pop("Func" + a_b, None)
        if func is None:
            raise AtError("Please provide a value for Func" + a_b)
        nsamp = len(func)
        setattr(self, "Func" + a_b, func)
        setattr(self, "NSamples" + a_b, nsamp)

    def _check_ramp(self, **kwargs: dict[str, any]):
        ramps = kwargs.get("Ramps", None)
        if ramps is not None:
            if len(ramps) != 4:
                raise AtError("Ramps has to be a vector with 4 elements")
            ramps = np.asarray(ramps)
        return ramps
