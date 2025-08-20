"""Time-dependent Multipole"""

from __future__ import annotations

from datetime import datetime
from enum import IntEnum

import numpy as np

from .conversions import _array
from .element_object import Element
from ..exceptions import AtError


class ACMode(IntEnum):
    """Class to define the excitation types"""

    SINE = 0
    WHITENOISE = 1
    ARBITRARY = 2

def _anyarray(value):
    # Ensure proper ordering(F) and alignment(A) for "C" access in integrators
    return np.require(value, dtype=np.float64, requirements=["F", "A"])


class VariableThinMultipole(Element):
    """Turn by turn variable thin multipole."""

    _BUILD_ATTRIBUTES = Element._BUILD_ATTRIBUTES + ["Mode"]
    _conversions = dict(
        Element._conversions,
        AmplitudeA=_array,
        AmplitudeB=_array,
        FrequencyA=float,
        FrequencyB=float,
        PhaseA=float,
        PhaseB=float,
        SinAabove=float,
        SinBabove=float,
        NsamplesA=int,
        NsamplesB=int,
        KtaylorA=int,
        KtaylorB=int,
        FuncA=_anyarray,
        FuncB=_anyarray,
        Ramps=_array,
        Periodic=bool,
    )

    def __init__(self, family_name: str, mode: int or ACMode, **kwargs):
        r"""VariableThinMultipole initialization.

        Default pass method: ``VariableThinMPolePass``.

        This Class creates a thin multipole of any order (dipole kick, quadrupole,
        sextupole, etc.) and type (Normal or Skew) defined by AmplitudeA and/or
        AmplitudeB components; the polynoms PolynomA and PolynomB are calculated
        on every turn depending on the chosen mode, and for some modes on the
        particle time delay. All modes could be ramped.

        Keep in mind that as this element varies on every turn, and at the end of
        the tracking PolynomA and PolynomB are set to zero.

        Passing arrays of zeros as amplitude will initialize the MaxOrder to
        zero, and the polynom to a single zero.

        There are three different modes that could be set:
            SINE = 0, WHITENOISE = 1 and ARBITRARY = 2. See ACMode.
        For example, use at.ACMode.SINE or 0 to create an sin function element.

        The **SINE** mode requires amplitude and frequency for A and/or B.
        The value of the jth component of the polynom (A or B) at the nth turn
        is given by
            Amplitude[j]*sin[TWOPI*frequency*(n*T0 + \tau_p) + phase],
        where T0 is the revolution period of the ideal ring, and \tau_p is the delay
        of the pth particle i.e. the sixth coordinate over the speed of light. Also,
        note that the position of the element on the ring has no effect, the phase
        could be used to add any delay due to the position along s. The following is
        an example of the SINE mode of an skew quad
            eleskew = at.VariableThinMultipole('VAR_SKEW',at.ACMode.SINE,
                AmplitudeA=[0,skewa2],FrequencyA=freqA,PhaseA=phaseA)
        The values of the sin function could be limited to be above a defined
        threshold using ``Sin[AB]above``. For example, you could create a positive
        half-sin by setting ``Sin[AB]above`` to zero. You could also create a
        negative half-sin function by setting the amplitude to -1.

        The **WHITENOISE** mode requires the amplitude of A and/or B. The gaussian
        distribution is generated with zero-mean and one standard deviation from
        a pseudo-random stream pcg32. The pcg32 seed is fixed by the tracking
        function, therefore using the same stream on all trackings (sequencial or
        parallel). See https://github.com/atcollab/at/discussions/879 for more
        details on the pseudo random stream. For example
            elenoise = at.VariableThinMultipole('MYNOISE',at.ACMode.WHITENOISE,
                AmplitudeA=[noiseA1])
        creates a vertical kick as gaussian noise of amplitude noiseA1.

        The **ARBITRARY** mode requires the definition of a custom discrete function
        to be sampled at every turn. The function and its Taylor expansion with
        respect to \tau up to any given order is
            value = f(turn) + f'(turn)*tau + 0.5*f''(turn)*tau**2
                    + 1/6*f'''(turn)*tau**3 + 1/24*f''''(turn)*tau**4 ...
        f is an array of values, f',f'',f''',f'''', are arrays containing
        the derivatives wrt \tau, and \tau is the time delay of the particle, i.e.
        the the sixth coordinate divided by the speed of light. Therefore, the
        function func is a 2D-array with columns corresponding to the f and its
        derivatives and rows to turns. For example, a positive value on the first turn
        f(1)=1 with positive derivative f'(1)=0.1 followed by a negative value in the
        second turn f(2)=-1 with negative derivative f'(2)-0.2 would be
        Func=[[1,-1],[0.1,-0.2]]. Use FuncA or FuncB, and AmplitudeA or AmplitudeB
        accordingly.
        tau could be offset using ``FuncATimeDelay`` or ``FuncBTimeDelay``.
            tau = tau - Func[AB]TimeDelay
        The function `value` is then **multiplied by Amplitude A and/or B**.
        For example, the following is a positive vertical kick in the first turn,
        negative on the second turn, and zero on the third turn.
            elesinglekick = at.VariableThinMultipole('CUSTOMFUNC',at.ACMODE.ARBITRARY,
                AmplitudeA=1e-4,FuncA=[1,-1,0],Periodic=True)
        By default the array is assumed non periodic, the function has no effect
        on the particle in turns exceeding the function definition. If
        ``Periodic`` is set to True, the sequence is repeated.

        One could use the method inspect_polynom_values to check the polynom values
        used in every turn.


        Parameters:
            family_name(str):  Element name
            mode(at.ACMode): defines the mode. Default ACMode.SINE:

              * :py:attr:`.ACMode.SINE`: sine function
              * :py:attr:`.ACMode.WHITENOISE`: gaussian white noise
              * :py:attr:`.ACMode.ARBITRARY`: user defined turn-by-turn kick list
        Keyword Arguments:
            AmplitudeA(list,float): Amplitude of the excitation for PolynomA.
                Default None
            AmplitudeB(list,float): Amplitude of the excitation for PolynomB.
                Default None
            FrequencyA(float): Frequency of the PolynomA sine excitation. Unit Hz
            FrequencyB(float): Frequency of the PolynomB sine excitation. Unit Hz
            PhaseA(float): Phase of the sine excitation for PolynomA. Default 0 rad
            PhaseB(float): Phase of the sine excitation for PolynomB. Default 0 rad
            SinAabove(float): Limit the sin function to be above. Default -1.
            SinBabove(float): Limit the sin function to be above. Default -1.
            FuncA(Sequence[float]):   User defined tbt list for PolynomA
            FuncB(Sequence[float]):   User defined tbt list for PolynomB
            FuncATimeDelay(float): generate a time offset on the function FUNCA.
                It only has an effect if any of the derivatives is not zero.
            FuncBTimeDelay(float): generate a time offset on the function FUNCB.
                It only has an effect if any of the derivatives is not zero.
            Periodic(bool): If True the user definition is repeated. Default False.
            Ramps(list): Four values (t0,t1,t2,t3) defining the ramping steps:

              * ``t<=t0``: excitation amplitude is zero.
              * ``t0<t<=t1``: excitation amplitude is linearly ramped up.
              * ``t1<t<=t2``: excitation amplitude is constant.
              * ``t2<t<=t3``: excitation amplitude is linearly ramped down.
              * ``t3<t``: excitation amplitude is zero.
        Raises:
            AtError: if none of ``AmplitudeA`` or ``AmplitudeB`` is passed.
            AtError: if ramp is not vector of length 4 when using ``Ramps``.
            AtError: when frequency is not defined if using Mode ``ACMode.SINE``.
            AtError: when the custom function is not defined if using mode
                ``ACMode.ARBITRARY``.


        Examples:

            >>> acmpole = at.VariableThinMultipole('ACSKEW', at.ACMode.SINE, AmplitudeA=[0,amp], FrequencyA=freq)
            >>> acmpole = at.VariableThinMultipole('ACHKICK', at.ACMode.WHITENOISE, AmplitudeB=amp)
            >>> acmpole = at.VariableThinMultipole('ACMPOLE', at.ACMode.ARBITRARY, AmplitudeB=[0,0,amp], FuncB=fun)

        """

        def _check_amplitudes(**kwargs) -> dict[str, Any]:
            amp_aux = {"A": None, "B": None}
            all_amplitudes_are_none = True
            for key in amp_aux:
                if "amplitude" + key in kwargs:
                    raise AtWarning(
                        "amplitude" + key + " should be Amplitude" + key + "."
                    )
                lower_case_kwargs = {k.lower(): v for k, v in kwargs.items()}
                amp_aux[key] = lower_case_kwargs.pop("amplitude" + key.lower(), None)
                if amp_aux[key] is not None:
                    all_amplitudes_are_none = False
            if all_amplitudes_are_none:
                raise AtError("Please provide at least one amplitude for A or B")
            return amp_aux

        def _getmaxorder(ampa: np.ndarray, ampb: np.ndarray) -> int:
            mxa, mxb = 0, 0
            if ampa is not None:
                mxa = np.max(np.append(np.nonzero(ampa), 0))
            if ampb is not None:
                mxb = np.max(np.append(np.nonzero(ampb), 0))
            return max(mxa, mxb)

        def _set_thismode(mode: int, a_b: str, **kwargs) -> dict[str, Any]:
            if mode == ACMode.SINE:
                kwargs = _set_sine(a_b, **kwargs)
            if mode == ACMode.ARBITRARY:
                kwargs = _set_arb(a_b, **kwargs)
            return kwargs

        def _set_sine(a_b: str, **kwargs) -> dict[str, Any]:
            frequency = kwargs.get("Frequency" + a_b, None)
            if frequency is None:
                raise AtError("Please provide a value for Frequency" + a_b)
            kwargs.setdefault("Phase" + a_b, 0)
            kwargs.setdefault("Sin" + a_b + "above", -1)
            return kwargs

        def _set_arb(a_b: str, **kwargs) -> dict[str, Any]:
            func = kwargs.get("Func" + a_b, None)
            if func is None:
                raise AtError("Please provide a value for Func" + a_b)
            if len(np.squeeze(func.shape)) == 1:
                nsamp = len(func)
                ktaylor = 1
            else:
                ktaylor, nsamples = func.shape
            kwargs.setdefault("Func" + a_b + "TimeDelay", 0)
            kwargs["NSamples" + a_b] = nsamples
            kwargs["Ktaylor" + a_b] = ktaylor
            kwargs.setdefault("Periodic", False)
            return kwargs

        def _check_ramp(**kwargs) -> _array:
            ramps = kwargs.get("Ramps", None)
            if ramps is not None:
                if len(ramps) != 4:
                    raise AtError("Ramps has to be a vector with 4 elements")
                ramps = np.asarray(ramps)
            return ramps

        # init begins
        kwargs.setdefault("Mode", mode)
        kwargs.setdefault("PassMethod", "VariableThinMPolePass")
        if len(kwargs) > 2:
            amp_aux = _check_amplitudes(**kwargs)
            for key, value in amp_aux.items():
                if value is not None:
                    kwargs["Amplitude" + key] = value
                    kwargs = _set_thismode(mode, key, **kwargs)
            maxorder = _getmaxorder(amp_aux["A"], amp_aux["B"])
            kwargs["MaxOrder"] = kwargs.get("MaxOrder", maxorder)
            for key in amp_aux:
                kwargs["Polynom" + key] = kwargs.get(
                    "Polynom" + key, np.zeros(maxorder + 1)
                )
            ramps = _check_ramp(**kwargs)
            if ramps is not None:
                kwargs["Ramps"] = ramps
        super().__init__(family_name, **kwargs)

    def inspect_polynom_values(self, **kwargs) -> dict[str, list]:
        """
        Get the polynom values per turn.

        Translations (T1,T2) and Rotations (R1,R2) in the element are ignored.

        Keyword arguments:
            turns(int): Default 1. Number of turns to calculate.
            T0(float): revolution time in seconds. Use only in SINE mode.
            tparticle(float): Default 0. Time of the particle in seconds.

        Returns:
            Dictionary with a list of PolynomA and PolynomB per turn.
        """
        turns = kwargs.setdefault("turns", 1)
        mode = self.Mode
        timeoffset = 0
        if mode == 0:
            # revolution time
            trevol = float(kwargs["T0"])
            tparticle = float(kwargs.setdefault("tparticle", 0))
            timeoffset = trevol + tparticle
        elif mode == 2:
            # particle time
            timeoffset = float(kwargs.setdefault("tparticle", 0))
        ramps = getattr(self, "Ramps", 0)
        periodic = getattr(self, "Periodic", False)
        maxorder = self.MaxOrder

        pola = np.full(maxorder + 1, np.nan)
        polb = np.full(maxorder + 1, np.nan)

        listpola = []
        listpolb = []

        for turn in range(turns):
            for order in range(maxorder + 1):
                if hasattr(self, "AmplitudeA"):
                    pola[order] = self._get_pol(
                        "A", ramps, mode, timeoffset * turn, turn, order, periodic
                    )
                if hasattr(self, "AmplitudeB"):
                    polb[order] = self._get_pol(
                        "B", ramps, mode, timeoffset * turn, turn, order, periodic
                    )
            listpola.append(np.copy(pola))
            listpolb.append(np.copy(polb))
        return {"PolynomA": listpola, "PolynomB": listpolb}

    def _get_amp(self, amp: float, ramps: _array, _time: float) -> float:
        """get_amp returns the input value `amp` when ramps is False.

        If ramps is True, it returns a value linearly interpolated
        accoding to the ramping turn.

        Parameters:
            amp: amplitude component.
            ramps: array containing the turns that define the ramp
            _time: turn

        Returns:
            amp if no ramp.
            amp multiplied by the ramp state.
        """
        ampt = amp
        if ramps != 0:
            if _time <= ramps[0]:
                ampt = 0.0
            elif _time <= ramps[1]:
                ampt = amp * (_time - ramps[0]) / (ramps[1] - ramps[0])
            elif _time <= ramps[2]:
                ampt = amp
            elif _time <= ramps[3]:
                ampt = amp - amp * (_time - ramps[2]) / (ramps[3] - ramps[2])
            else:
                ampt = 0.0
        return ampt

    def _get_pol(
        self,
        a_b: str,
        ramps: _array,
        mode: int,
        _time: float,
        turn: int,
        order: int,
        periodic: bool,
    ) -> float:
        """
        Return the polynom component of a given order.

        Parameters:
            a_b: either 'A' or 'B' indicating the polynom.
            ramps: array containing the ramp definition.
            mode: value to specify the type of variable element.
            _time: time for this mode
            turn: turn to check
            order: order of the polynom
            periodic: whether the sequence is periodic or not.

        Returns:
            the amplitude for the polynom component
        """
        allamp = getattr(self, "Amplitude" + a_b)
        amp = allamp[order]
        ampout = 0
        # check if amp is zero
        if amp == 0:
            return ampout

        # get the ramp value
        ampout = self._get_amp(amp, ramps, turn)

        if mode == 0:
            # sin mode parameters
            whole_sin_above = getattr(self, "Sin" + a_b + "above")
            freq = getattr(self, "Frequency" + a_b)
            phase = getattr(self, "Phase" + a_b)
            sinval = np.sin(2 * np.pi * freq * _time + phase)
            if sinval >= whole_sin_above:
                ampout = ampout * sinval
            else:
                ampout = 0
        elif mode == 1:
            ampout = np.nan
        elif mode == 2:
            nsamples = getattr(self, "NSamples" + a_b)
            if periodic or turn < nsamples:
                func = getattr(self, "Func" + a_b)
                functdelay = float(getattr(self, "Func" + a_b + "TimeDelay"))
                turnidx = np.mod(turn, nsamples)

                _time = _time - functdelay
                _time_2 = _time * _time
            else:
                ampout = 0.0
        else:
            ampout = 0.0
        return ampout

