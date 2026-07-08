"""Time-dependent thin multipole."""

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


class VariableThinMultipole(Element):
    """Class to generate an AT variable thin multipole element"""

    _BUILD_ATTRIBUTES = Element._BUILD_ATTRIBUTES + ["Mode"]
    _conversions = dict(
        Element._conversions,
        Mode=str,
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

    def __init__(
        self, family_name, mode, AmplitudeA=None, AmplitudeB=None, **kwargs
    ):
        # noinspection PyUnresolvedReferences,SpellCheckingInspection
        r"""
        Parameters:
            family_name(str):    Element name
            mode:  one of the following:

              * :py:attr:`SINE`: sine function
              * :py:attr:`WHITENOISE`: gaussian white noise
              * :py:attr:`ARBITRARY`: user defined turn-by-turn kick list

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
            Seed(int): Seed of the random number generator for white
                       noise excitation. Default datetime.now()
            FuncA(list): User defined tbt kick list for PolynomA
            FuncB(list): User defined tbt kick list for PolynomB
            Periodic(bool): If True (default) the user defined kick is repeated
            Ramps(list): Vector (t0, t1, t2, t3) in turn number to define the ramping
                         of the excitation

              * ``t<t0``: excitation amplitude is zero
              * ``t0<t<t1``: exciation amplitude is linearly ramped up
              * ``t1<t<t2``: exciation amplitude is constant
              * ``t2<t<t3``: exciation amplitude is linearly ramped down
              * ``t3<t``: exciation amplitude is zero

        Examples:

            >>> acmpole = at.VariableThinMultipole(
            ...     "ACMPOLE", "SINE", AmplitudeB=amp, FrequencyB=frequency
            ... )
            >>> acmpole = at.VariableThinMultipole(
            ...     "ACMPOLE", "WHITENOISE", AmplitudeB=amp, ... )
            >>> acmpole = at.VariableThinMultipole(
            ...     "ACMPOLE", "ARBITRARY", AmplitudeB=amp, FuncB=fun, ... )

        .. note::

            * At least AmplitudeA or AmplitudeB has to be provided.
            * For ``mode="SINE"`` the ``Frequency(A,B)`` corresponding to the
              ``Amplitude(A,B)`` has to be provided
            * For ``mode="ARBITRARY"`` the ``Func(A,B)`` corresponding to the
              ``Amplitude(A,B)`` has to be provided
        """
        kwargs.setdefault("PassMethod", "VariableThinMPolePass")
        self.MaxOrder = kwargs.pop("MaxOrder", 0)
        self.Periodic = kwargs.pop("Periodic", True)
        self.Mode = mode
        if AmplitudeA is None and AmplitudeB is None:
            raise AtError("Please provide at least one amplitude for A or B")
        AmplitudeB = self._set_params(AmplitudeB, mode, "B", **kwargs)
        AmplitudeA = self._set_params(AmplitudeA, mode, "A", **kwargs)
        self._setmaxorder(AmplitudeA, AmplitudeB)
        if mode == "WHITENOISE":
            self.Seed = kwargs.pop("Seed", datetime.now().timestamp())
        self.PolynomA = np.zeros(self.MaxOrder + 1)
        self.PolynomB = np.zeros(self.MaxOrder + 1)
        ramps = kwargs.pop("Ramps", None)
        if ramps is not None:
            assert len(ramps) == 4, "Ramps has to be a vector with 4 elements"
            self.Ramps = ramps
        super().__init__(family_name, **kwargs)

    def _setmaxorder(self, ampa, ampb):
        mxa, mxb = 0, 0
        if ampa is not None:
            mxa = np.max(np.nonzero(ampa))
        if ampb is not None:
            mxb = np.max(np.nonzero(ampb))
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

    def _set_params(self, amplitude, mode, ab, **kwargs):
        if amplitude is not None:
            if np.isscalar(amplitude):
                amp = np.zeros(self.MaxOrder)
                amplitude = np.append(amp, amplitude)
            if mode == "SINE":
                self._set_sine(ab, **kwargs)
            if mode == "ARBITRARY":
                self._set_arb(ab, **kwargs)
        return amplitude

    def _set_sine(self, ab, **kwargs):
        frequency = kwargs.pop("Frequency" + ab, None)
        phase = kwargs.pop("Phase" + ab, 0)
        assert frequency is not None, "Please provide a value for Frequency" + ab
        setattr(self, "Frequency" + ab, frequency)
        setattr(self, "Phase" + ab, phase)

    def _set_arb(self, ab, **kwargs):
        func = kwargs.pop("Func" + ab, None)
        nsamp = len(func)
        assert func is not None, "Please provide a value for Func" + ab
        setattr(self, "Func" + ab, func)
        setattr(self, "NSamples" + ab, nsamp)
