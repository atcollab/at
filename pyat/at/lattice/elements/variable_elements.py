"""Time-dependent thin multipole."""

from __future__ import annotations

from datetime import datetime
from enum import IntEnum

import numpy as np

from ..exceptions import AtError
from .conversions import _anyarray, _array
from .element_object import Element


class ACMode(IntEnum):
    """Class to define the excitation types."""

    SINE = 0
    WHITENOISE = 1
    ARBITRARY = 2


class VariableThinMultipole(Element):
    """Class to generate an AT variable thin multipole element."""

    _BUILD_ATTRIBUTES = Element._BUILD_ATTRIBUTES
    _conversions = dict(
        Element._conversions,
        Mode=int,
        ModeName=str,
        AmplitudeA=_array,
        AmplitudeB=_array,
        FrequencyA=float,
        FrequencyB=float,
        PhaseA=float,
        PhaseB=float,
        Sinmin=float,
        Sinmax=float,
        Seed=int,
        NSamplesA=int,
        NSamplesB=int,
        FuncA=_anyarray,
        FuncB=_anyarray,
        Ramps=_array,
        Periodic=bool,
    )

    def __init__(
        self, family_name, mode=ACMode.SINE, AmplitudeA=None, AmplitudeB=None, **kwargs
    ):
        # noinspection PyUnresolvedReferences,SpellCheckingInspection
        r"""Create a variable thin multipole.

        Parameters:
            family_name(str):    Element name
            mode(ACMode):  one of the following at.ACMode. Default at.ACMode.SINE.

              * :py:attr:`at.ACMode.SINE`: sine function
              * :py:attr:`at.ACMode.WHITENOISE`: gaussian white noise
              * :py:attr:`at.ACMode.ARBITRARY`: user defined turn-by-turn kick list

        Keyword Arguments:
            AmplitudeA(list,float): Amplitude of the excitation for PolynomA.
              Default None
            AmplitudeB(list,float): Amplitude of the excitation for PolynomB.
              Default None

            FrequencyA(float): Frequency of the sine excitation for PolynomA
            FrequencyB(float): Frequency of the sine excitation for PolynomB
            PhaseA(float): Phase of the sine excitation for PolynomA. Default 0
            PhaseB(float): Phase of the sine excitation for PolynomB. Default 0
            Sinmin(float): Sine function min limit. Default -1.1
            Sinmax(float): Sine function max limit. Default +1.1
            MaxOrder(int): Order of the multipole for scalar amplitude. Default 0
            Seed(int): Seed of the random number generator for white
                       noise excitation. Default datetime.now()
            FuncA(ndarray): User defined tbt kick list for PolynomA
            FuncB(ndarray): User defined tbt kick list for PolynomB
            FuncDelay(float): Value to substract from the particles 6th coordinate
            Periodic(bool): If True (default) the user defined kick is repeated
            Ramps(list): Vector (t0, t1, t2, t3) in turn number to define the ramping
                         of the excitation

              * ``t<t0``: excitation amplitude is zero
              * ``t0<t<t1``: excitation amplitude is linearly ramped up
              * ``t1<t<t2``: excitation amplitude is constant
              * ``t2<t<t3``: excitation amplitude is linearly ramped down
              * ``t3<t``: excitation amplitude is zero

        Examples:

            >>> acmpole = at.VariableThinMultipole(
            ...     "ACMPOLE", at.ACMode.SINE, AmplitudeB=amp, FrequencyB=frequency
            ... )
            >>> pos_halfsine = at.VariableThinMultipole(
            ...     "PHALFSINE", at.ACMode.SINE, AmplitudeB=amp, FrequencyB=frequency,
            ...     Sinmin=0)
            >>> sine_saturated, at.VariableThinMultipole(
            ...     "SATSINE", at.ACMode.SINE, AmplitudeB=amp, FrequencyB=frequency,
            ...     Sinmax=0.9)
            >>> acmpole = at.VariableThinMultipole(
            ...     "ACMPOLE", at.ACMode.WHITENOISE, AmplitudeB=amp, ... )
            >>> acmpole = at.VariableThinMultipole(
            ...     "ACMPOLE", at.ACMode.ARBITRARY, AmplitudeB=amp, FuncB=fun, ... )
            >>> customf = at.VariableThinMultiple(
            ...     "ACFUNC", at.ACMoode.ARBITRARY, AmplitudA=amp, ...
            ...     FuncA=np.array([[1,.1,0.01],[0.2,0.02,0.002]]), ...
            ...     FuncDelay=0.2)

        .. note::

            * At least AmplitudeA or AmplitudeB has to be provided.
            * For ``mode=at.ACMode.SINE`` the ``Frequency(A,B)`` corresponding to the
              ``Amplitude(A,B)`` has to be provided
            * For ``mode=at.ACMode.ARBITRARY`` the ``Func(A,B)`` corresponding to the
              ``Amplitude(A,B)`` has to be provided
            * Func(A,B) could be an array of size (m,n) with n coefficients in the first
              row for the function over n turns, and other m-1 rows with higher order
              derivatives with respect to ctau
                i.e. on the kth turn the ith component of the Polynom(A/B) seen by a
                particle with 6th coordinate ctau is calculated as,
                   amp(A/B)[i] * func(0,k) + (ctau-delay) * func(1,k) + ...
                                        + (ctau-delay)**m * func(m,k)
        """

        def _default_amplitudes(ampa, ampb):
            if ampa is None and ampb is None:
                ampb = np.array([0])
            if np.isscalar(ampa):
                ampa = np.array([ampa])
            if np.isscalar(ampb):
                ampb = np.array([ampb])
            return ampa, ampb

        def _getmaxorder(ampa, ampb):
            mxa, mxb = 0, 0
            if ampa is not None:
                mxa = np.max(np.append(np.nonzero(ampa), 0))
            if ampb is not None:
                mxb = np.max(np.append(np.nonzero(ampb), 0))
            return max(mxa, mxb)

        self.Mode = kwargs.get("Mode", mode.value)
        self.ModeName = kwargs.get("ModeName", mode.name)
        kwargs.setdefault("PassMethod", "VariableThinMPolePass")
        AmplitudeA, AmplitudeB = _default_amplitudes(AmplitudeA, AmplitudeB)
        # MaxOrder is set finally by the user if given
        max_order_ampab = _getmaxorder(AmplitudeA, AmplitudeB)
        self.MaxOrder = kwargs.get("MaxOrder", max_order_ampab)
        # after the definition of MaxOrder we can create Amplitudes
        self._set_amplitudes(AmplitudeA, AmplitudeB)
        self._set_params(AmplitudeB, "B", **kwargs)
        self._set_params(AmplitudeA, "A", **kwargs)
        if self.Mode == ACMode.WHITENOISE:
            self.Seed = kwargs.pop("Seed", datetime.now().timestamp())
        self.PolynomA = kwargs.get("PolynomA", np.zeros(self.MaxOrder + 1))
        self.PolynomB = kwargs.get("PolynomB", np.zeros(self.MaxOrder + 1))
        ramps = kwargs.pop("Ramps", None)
        if ramps is not None:
            assert len(ramps) == 4, "Ramps has to be a vector with 4 elements"
            self.Ramps = ramps
        super().__init__(family_name, **kwargs)

    def _set_amplitudes(self, ampa, ampb):
        if ampa is not None:
            delta = self.MaxOrder + 1 - len(ampa)
            if delta > 0:
                ampa = np.pad(ampa, (0, delta))
            self.AmplitudeA = ampa
        if ampb is not None:
            delta = self.MaxOrder + 1 - len(ampb)
            if delta > 0:
                ampb = np.pad(ampb, (0, delta))
            self.AmplitudeB = ampb

    def _set_params(self, amplitude, ab, **kwargs):
        if amplitude is not None:
            if self.Mode == ACMode.SINE:
                self._set_sine(ab, **kwargs)
            if self.Mode == ACMode.ARBITRARY:
                self._set_arb(ab, **kwargs)

    def _set_sine(self, ab, **kwargs):
        frequency = kwargs.pop("Frequency" + ab, 0)
        phase = kwargs.pop("Phase" + ab, 0)
        sinmin = kwargs.pop("Sinmin", -1.1)
        sinmax = kwargs.pop("Sinmax", 1.1)
        setattr(self, "Frequency" + ab, frequency)
        setattr(self, "Phase" + ab, phase)
        self.Sinmin = sinmin
        self.Sinmax = sinmax

    def _set_arb(self, ab, **kwargs):
        func = kwargs.pop("Func" + ab, None)
        if func is None:
            raise AtError("Please provide a value for Func" + ab)
        if np.any(np.isnan(func)):
            raise AtError("Function" + ab + " contains nan.")
        if not np.all(np.isreal(func)):
            raise AtError("Function" + ab + " contains non real values.")
        func = np.squeeze(func)
        if func.ndim == 1:
            nsamples = len(func)
            rows = 1
        else:
            rows, nsamples = func.shape
        setattr(self, "Func" + ab, func)
        setattr(self, "NSamples" + ab, nsamples)
        setattr(self, "Dorder" + ab, rows - 1)
        self.FuncDelay = kwargs.setdefault("FuncDelay", 0)
        self.Periodic = kwargs.setdefault("Periodic", False)
