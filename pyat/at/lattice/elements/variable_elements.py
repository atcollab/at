"""Time-dependent thin multipole."""

from __future__ import annotations

from enum import IntEnum
from warnings import warn

import numpy as np

from .conversions import _anyarray, _array
from .element_object import Element


class ACMode(IntEnum):
    """Class to define the excitation types."""

    SINE = 0
    WHITENOISE = 1
    ARBITRARY = 2
    INTERPOLATION_TABLE = 3


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
        NSamplesA=int,
        NSamplesB=int,
        FuncA=_anyarray,
        FuncB=_anyarray,
        FinterpolateA=_array,
        FinterpolateB=_array,
        TinterpolateA=_array,
        TinterpolateB=_array,
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
              * :py:attr:`at.ACMode.INTERPOLATION_TABLE`: linear interpolation from user curve

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
            FuncA(list): User defined tbt kick list for PolynomA
            FuncB(list): User defined tbt kick list for PolynomB
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
            >>> fvst = at.VariableThinMultipole(
            ...     "FvsT", at.ACMode.INTERPOLATION_TABLE, AmplitudeA=amp, FuncB=func, ...)

        .. note::

            * At least AmplitudeA or AmplitudeB has to be provided.
            * For ``mode=at.ACMode.SINE`` the ``Frequency(A,B)`` corresponding to the
              ``Amplitude(A,B)`` has to be provided
            * For ``mode=at.ACMode.ARBITRARY`` the ``Func(A,B)`` corresponding to the
              ``Amplitude(A,B)`` has to be provided
            * For ``mode=at.ACMode.INTERPOLATION_TABLE`` the ``Fnc(A,B)`` corresponding to the
              ``Amplitude(A,B)`` needs to be of shape (2, n) with n >= 2. The first row
              is time in seconds, and the second row is the function value.
        """

        def _default_amplitudes(ampa, ampb):
            if ampa is None and ampb is None:
                ampb = np.array([0])
            if np.ndim(ampa) == 0 and ampa is not None:
                ampa = np.array([float(ampa)])
            if np.ndim(ampb) == 0 and ampb is not None:
                ampb = np.array([float(ampb)])
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
        AmplitudeA = kwargs.get("AmplitudeA", AmplitudeA)
        AmplitudeB = kwargs.get("AmplitudeB", AmplitudeB)
        self.ModeName = kwargs.get("ModeName", mode.name)
        kwargs.setdefault("PassMethod", "VariableThinMPolePass")
        AmplitudeA, AmplitudeB = _default_amplitudes(AmplitudeA, AmplitudeB)
        # MaxOrder is set finally by the user if given
        max_order_ampab = _getmaxorder(AmplitudeA, AmplitudeB)
        self.MaxOrder = kwargs.get("MaxOrder", max_order_ampab)
        # after the definition of MaxOrder we can create Amplitudes
        self._set_amplitudes(AmplitudeA, AmplitudeB)
        self.Periodic = kwargs.pop("Periodic", True)
        self._set_params(AmplitudeB, "B", **kwargs)
        self._set_params(AmplitudeA, "A", **kwargs)
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
            if self.Mode == ACMode.INTERPOLATION_TABLE:
                self._set_interpolate(ab, **kwargs)

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
        assert func is not None, "Please provide a value for Func" + ab
        nsamp = len(func)
        setattr(self, "Func" + ab, func)
        setattr(self, "NSamples" + ab, nsamp)

    def _set_interpolate(self, ab, **kwargs):
        interpolate = kwargs.get("Func" + ab)
        assert np.ndim(interpolate) == 2, "Func" + ab + " should be of 2 dimensions."
        _, nsamp = np.shape(interpolate)
        assert nsamp >= 2, "Func" + ab + " requires at least two points to interpolate."
        assert ~np.any(np.isnan(interpolate)), "Function has nan values."
        assert ~np.any(np.isinf(interpolate)), "Function has inf values."
        tsort = interpolate[0, :]
        assert len(tsort) == len(np.unique(tsort)), "Time array has repeated elements."
        idxsort = np.argsort(interpolate[0, :])
        fsort = interpolate[1, :]
        if ~np.all(np.diff(idxsort) == 1):
            warn(UserWarning("Time is not sorted. It will be rearanged."), stacklevel=2)
            tsort = interpolate[0, idxsort]
            fsort = interpolate[1, idxsort]
        assert (tsort[-1] - tsort[0]) > 0, "Zero time cannot be interpolated"
        setattr(self, "Tinterpolate" + ab, tsort)
        setattr(self, "Finterpolate" + ab, fsort)
        setattr(self, "NSamples" + ab, nsamp)
