from __future__ import annotations

from datetime import datetime
from enum import IntEnum

import numpy as np

from .elements import Element, _array
from .utils import AtError

__all__ = ["ACMode","VariableMultipole"]

class ACMode(IntEnum):
    """Class to define the excitation types."""

    SINE = 0
    WHITENOISE = 1
    ARBITRARY = 2


class VariableMultipole(Element):
    """Class to generate an AT variable thin multipole element."""

    _BUILD_ATTRIBUTES = Element._BUILD_ATTRIBUTES + ["Mode"]
    _conversions = dict(
        Element._conversions,
        AmplitudeA=_array,
        AmplitudeB=_array,
        FrequencyA=float,
        FrequencyB=float,
        PhaseA=float,
        PhaseB=float,
        SinlimitA=float,
        SinlimitB=float,
        NsamplesA=int,
        NsamplesB=int,
        FuncA=_array,
        FuncAderiv1=_array,
        FuncAderiv2=_array,
        FuncAderiv3=_array,
        FuncAderiv4=_array,
        FuncB=_array,
        FuncBderiv1=_array,
        FuncBderiv2=_array,
        FuncBderiv3=_array,
        FuncBderiv4=_array,
        Ramps=_array,
        Periodic=bool,
    )

    def __init__(self, family_name: str, mode: int, **kwargs: dict[str, any]):
        r"""
        Create VariableMultipole.

        This function creates a thin multipole any order (1 to k) and type (Normal
        or Skew) defined by the amplitude A or B components; the polynoms PolynomA
        and PolynomB are calculated on every turn depending on the chosen mode, and
        for some modes on the particle time delay. All modes could be ramped.

        Keep in mind that this element varies on every turn, therefore, any ring
        containing a variable element may change after tracking n turns.

        There are three different modes implemented:
            SINE = 0, WHITENOISE = 1 and ARBITRARY = 2.

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
        The values of the sin function could be limited to be above a defined
        threshold using `Sinlimit`. For example, you could create a half-sin
        by setting `Sinlimit` to zero.

        The WHITENOISE mode requires the amplitude of either A or B. For example
            elenoise = at.VariableMultipole('MYNOISE',ACMode.WHITENOISE,
                AmplitudeA=[noiseA1])
        creates a gaussian vertical noise of amplitude noiseA1. The gaussian
        distribution is generated with zero-mean and one standard deviation from
        a pseudo-random stream pcg32. The pcg32 seed is fixed by the tracking
        function, therefore using the same stream on all trackings (sequencial or
        parallel). See https://github.com/atcollab/at/discussions/879 for more
        details on the pseudo random stream.  Additionally, the user could set
        the mean and std values by setting MeanA, MeanB, StdA, StdB.

        The ARBITRARY mode requires the definition of a custom turn-by-turn function.
        The user defines the value of the function and its Taylor expansion with
        respect to \tau up to fourth order.
            value = f(turn) + f'(turn)*tau + 0.5*f''(turn)*tau**2
                    + 1/6*f'''(turn)*tau**3 + 1/24*f''''(turn)*tau**4
        where f is an array of values, f',f'',f''',f'''', are arrays containing
        the derivatives wrt \tau, and \tau is the time delay of the particle, i.e.
        the the sixth coordinate divided by the speed of light.
        For example, the following is a positive vertical kick in the first turn,
        negative on the second turn, and zero on the third turn.
            funca = [1 -1 0];
            elesinglekick = at.VariableMultipole('CUSTOMFUNC','ARBITRARY',
            AmplitudeA=1e-4,FuncA=funca,Periodic=True)
        by default the array is assumed periodic. If Periodic is set to False
        any turn exceeding the defined length of the function is assumed to have
        no effect on the beam.


        Parameters:
            family_name(str):  Element name
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
            PhaseA(float): Phase of the sine excitation for PolynomA. Default 0 rad
            PhaseB(float): Phase of the sine excitation for PolynomB. Default 0 rad
            SinlimitA(float): Default -1. Values of the sin function will be above
                SinlimitA or zero.
            SinlimitB(float): Default -1. Values of the sin function will be above
                SinlimitB or zero.
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
                    kwargs = self._set_mode(mode, k, **kwargs)
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

    def _set_mode(self, mode, a_b: str, **kwargs: dict[str, any]):
        if mode == ACMode.WHITENOISE:
            kwargs = self._set_white_noise(a_b, **kwargs)
        if mode == ACMode.SINE:
            kwargs = self._set_sine(a_b, **kwargs)
        if mode == ACMode.ARBITRARY:
            kwargs = self._set_arb(a_b, **kwargs)
        return kwargs

    def _set_white_noise(self, a_b: str, **kwargs: dict[str, any]):
        if "Mean" + a_b not in kwargs:
            kwargs["Mean" + a_b] = kwargs.get("Mean" + a_b, 0)
        if "Std" not in kwargs:
            kwargs["Std" + a_b] = kwargs.get("Std" + a_b, 1)
        return kwargs

    def _set_sine(self, a_b: str, **kwargs: dict[str, any]):
        frequency = kwargs.get("Frequency" + a_b, None)
        if frequency is None:
            raise AtError("Please provide a value for Frequency" + a_b)
        kwargs["Phase" + a_b] = kwargs.get("Phase" + a_b, 0)
        kwargs["Sinlimit" + a_b] = kwargs.get("Sinlimit" + a_b, -1)
        return kwargs

    def _set_arb(self, a_b: str, **kwargs: dict[str, any]):
        func = kwargs.get("Func" + a_b, None)
        if func is None:
            raise AtError("Please provide a value for Func" + a_b)
        nsamp = len(func)
        kwargs["Func" + a_b + "deriv1"] = kwargs.get(
            "Func" + a_b + "deriv1", np.zeros(nsamp)
        )
        kwargs["Func" + a_b + "deriv2"] = kwargs.get(
            "Func" + a_b + "deriv2", np.zeros(nsamp)
        )
        kwargs["Func" + a_b + "deriv3"] = kwargs.get(
            "Func" + a_b + "deriv3", np.zeros(nsamp)
        )
        kwargs["Func" + a_b + "deriv4"] = kwargs.get(
            "Func" + a_b + "deriv4", np.zeros(nsamp)
        )
        kwargs["Func" + a_b + "TimeDelay"] = kwargs.get(
            "Func" + a_b + "TimeDelay", 0
        )
        kwargs["NSamples" + a_b] = nsamp
        kwargs["Periodic"] = kwargs.get("Periodic", True)
        return kwargs

    def _check_ramp(self, **kwargs: dict[str, any]):
        ramps = kwargs.get("Ramps", None)
        if ramps is not None:
            if len(ramps) != 4:
                raise AtError("Ramps has to be a vector with 4 elements")
            ramps = np.asarray(ramps)
        return ramps
