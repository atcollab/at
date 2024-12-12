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

    _BUILD_ATTRIBUTES = Element._BUILD_ATTRIBUTES + ['Mode']
    _conversions = dict(
        Element._conversions,
        amplitudea=_array,
        amplitudeb=_array,
        frequencya=float,
        frequencyb=float,
        phasea=float,
        phaseb=float,
        seed=int,
        nsamplesa=int,
        nsamplesb=int,
        funca=_array,
        funcb=_array,
        ramps=_array,
        periodic=bool,
    )

    def __init__(self, family_name: str, mode:int, **kwargs: dict[str, any]):
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
        of the kth particle i.e. the sixth coordinate over the speed of light. Also,
        note that the position of the element on the ring has no effect, the phase
        should be used to add any delay due to the position along s.
        The following is an example of the SINE mode of an skew quad:
            eleskew = at.VariableMultipole('VAR_SKEW',
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
        print(kwargs)
        print(vars(self))
        if len(kwargs) > 0:
            self.FamName = family_name
            self.Mode = int(mode)
            # start setting up Amplitudes
            # amplitude names need to be different in pyat and matlab in order to avoid
            # problems between args and kwargs in python when loading a .mat file
            # as a dictionary.
            # I chose to use lower case in the pyat args, and A and B when reading
            # from kwargs.
            amp_aux = {'A':None, 'B':None}
            all_amplitudes_are_none = True
            for key in amp_aux.keys():
                if 'Amplitude'+key in kwargs and 'amplitude'+key in kwargs:
                    raise AtError('Duplicated amplitude '+key+'parameters.')
                    some_amplitude_exists = 1
                lower_case_kwargs = {k.lower(): v for k, v in kwargs.items()}
                amp_aux[key] = lower_case_kwargs.pop('amplitude'+key.lower(), None)
                if amp_aux[key] is not None:
                    all_amplitudes_are_none = False
            if all_amplitudes_are_none:
                raise AtError("Please provide at least one amplitude for A or B")
            for k,v in amp_aux.items():
                amp_aux[k] = self._set_amplitude(v)
                if amp_aux[k] is not None:
                    setattr(self, 'Amplitude'+k, amp_aux[k])
                    self._set_params(mode, k, **kwargs)
            kwargs.setdefault("PassMethod", "VariableThinMPolePass")
            self._setmaxorder(amp_aux['A'], amp_aux['B'])
            if mode == ACMode.WHITENOISE and "seed" not in kwargs:
                kwargs.update({"Seed": datetime.now().timestamp()})
            self.Periodic = kwargs.pop("Periodic", True)
            for key in amp_aux.keys():
                if amp_aux[key] is not None:
                    setattr(self, 'Polynom'+key, amp_aux[key])
                else:
                    setattr(self, 'Polynom'+key, np.zeros(self.MaxOrder+1))
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

    def _set_amplitude(self, amplitude: float or _array or None):
        if amplitude is not None:
            if np.isscalar(amplitude):
                amplitude = [amplitude]
            amplitude = np.asarray(amplitude)
        return amplitude

    def _set_params( self, mode, a_b: str, **kwargs: dict[str, any]):
        if mode == ACMode.SINE:
            self._set_sine(a_b, **kwargs)
        if mode == ACMode.ARBITRARY:
            self._set_arb(a_b, **kwargs)

    def _set_sine(self, a_b: str, **kwargs: dict[str, any]):
        frequency = kwargs.pop("Frequency" + a_b, None)
        if frequency is None:
            raise AtError("Please provide a value for Frequency" + a_b)
        phase = kwargs.pop("Phase" + a_b, 0)
        setattr(self, "Frequency" + a_b, frequency)
        setattr(self, "Phase" + a_b, phase)

    def _set_arb(self, a_b: str, **kwargs: dict[str, any]):
        func = kwargs.pop("Func" + a_b, None)
        if func is None:
            raise AtError("Please provide a value for Func" + a_b)
        nsamp = len(func)
        setattr(self, "Func" + a_b, func)
        setattr(self, "NSamples" + a_b, nsamp)
