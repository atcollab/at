"""Basic :py:class:`.Element` classes."""

from __future__ import annotations

__all__ = [
    "M66",
    "Aperture",
    "BeamMoments",
    "Collimator",
    "Drift",
    "EnergyLoss",
    "LongElement",
    "LongtAperture",
    "Marker",
    "Monitor",
    "QuantumDiffusion",
    "RFCavity",
    "SimpleQuantDiff",
    "SimpleRadiation",
    "SliceMoments",
]

import warnings
from collections.abc import Iterable

import numpy as np

from ..exceptions import AtWarning
from .conversions import _array, _array66, _int
from .abstract_elements import LongtMotion, Radiative, _DictLongtMotion
from .element_object import Element

# AtWarning from this module should always be issued (not only on the first occurrence)
warnings.filterwarnings("always", category=AtWarning, module=__name__)


class LongElement(Element):
    """Base class for long elements."""

    _BUILD_ATTRIBUTES = [*Element._BUILD_ATTRIBUTES, "Length"]

    def __init__(self, family_name: str, length: float, *args, **kwargs):
        """
        Args:
            family_name:    Name of the element
            length:         Element length [m].

        Other arguments and keywords are given to the base class
        """
        kwargs.setdefault("Length", length)
        # Ancestor may be either Element or ThinMultipole
        # noinspection PyArgumentList
        super().__init__(family_name, *args, **kwargs)

    def _part(self, fr, sumfr):
        pp = self.copy()
        pp.Length = fr * self.get_parameter("Length")
        if hasattr(self, "KickAngle"):
            pp.KickAngle = fr / sumfr * self.KickAngle
        return pp

    def divide(self, frac) -> list[Element]:
        def popattr(element, attr):
            val = element.get_parameter(attr)  # get the parameter itself
            delattr(element, attr)
            return attr, val

        frac = np.asarray(frac, dtype=float)
        el = self.copy()
        # Remove entrance and exit attributes
        attrs = el.keys()
        fin = dict(popattr(el, key) for key in attrs if key in self._entrance_fields)
        fout = dict(popattr(el, key) for key in attrs if key in self._exit_fields)
        # Split element
        element_list = [el._part(f, np.sum(frac)) for f in frac]
        # Restore entrance and exit attributes
        for key, value in fin.items():
            setattr(element_list[0], key, value)
        for key, value in fout.items():
            setattr(element_list[-1], key, value)
        return element_list

    def is_compatible(self, other) -> bool:
        def compatible_field(fieldname):
            f1 = getattr(self, fieldname, None)
            f2 = getattr(other, fieldname, None)
            if f1 is None and f2 is None:  # no such field
                return True
            elif f1 is None or f2 is None:  # only one
                return False
            else:  # both
                return np.all(f1 == f2)

        if not (type(other) is type(self) and self.PassMethod == other.PassMethod):
            return False
        return all(compatible_field(fname) for fname in ("RApertures", "EApertures"))

    def merge(self, other) -> None:
        super().merge(other)
        self.Length += other.Length


class Marker(Element):
    """Marker element."""


class Monitor(Element):
    """Monitor element."""


class BeamMoments(Element):
    """Element to compute bunches mean and std."""

    def __init__(self, family_name: str, **kwargs):
        """
        Args:
            family_name:    Name of the element.

        Default PassMethod: ``BeamMomentsPass``
        """
        kwargs.setdefault("PassMethod", "BeamMomentsPass")
        self._stds = np.zeros((6, 1, 1), order="F")
        self._means = np.zeros((6, 1, 1), order="F")
        super().__init__(family_name, **kwargs)

    def set_buffers(self, nturns, nbunch):
        self._stds = np.zeros((6, nbunch, nturns), order="F")
        self._means = np.zeros((6, nbunch, nturns), order="F")

    @property
    def stds(self):
        """Beam 6d standard deviation."""
        return self._stds

    @property
    def means(self):
        """Beam 6d centre of mass."""
        return self._means


class SliceMoments(Element):
    """Element computing the mean and std of slices."""

    _BUILD_ATTRIBUTES = [*Element._BUILD_ATTRIBUTES, "nslice"]
    _conversions = dict(Element._conversions, nslice=_int)

    def __init__(self, family_name: str, nslice: int, **kwargs):
        """
        Args:
            family_name:    Name of the element
            nslice:         Number of slices.

        Keyword arguments:
            startturn:      Start turn of the acquisition (Default 0)
            endturn:        End turn of the acquisition (Default 1)

        Default PassMethod: ``SliceMomentsPass``
        """
        kwargs.setdefault("PassMethod", "SliceMomentsPass")
        self._startturn = kwargs.pop("startturn", 0)
        self._endturn = kwargs.pop("endturn", 1)
        super().__init__(family_name, nslice=nslice, **kwargs)
        self._nbunch = 1
        self.startturn = self._startturn
        self.endturn = self._endturn
        self._dturns = self.endturn - self.startturn
        self._stds = np.zeros((3, nslice, self._dturns), order="F")
        self._means = np.zeros((3, nslice, self._dturns), order="F")
        self._spos = np.zeros((nslice, self._dturns), order="F")
        self._weights = np.zeros((nslice, self._dturns), order="F")
        self.set_buffers(self._endturn, 1)

    def set_buffers(self, nturns, nbunch):
        self.endturn = min(self.endturn, nturns)
        self._dturns = self.endturn - self.startturn
        self._nbunch = nbunch
        self._stds = np.zeros((3, nbunch * self.nslice, self._dturns), order="F")
        self._means = np.zeros((3, nbunch * self.nslice, self._dturns), order="F")
        self._spos = np.zeros((nbunch * self.nslice, self._dturns), order="F")
        self._weights = np.zeros((nbunch * self.nslice, self._dturns), order="F")

    @property
    def stds(self):
        """Slices x,y,dp standard deviation."""
        return self._stds.reshape((3, self._nbunch, self.nslice, self._dturns))

    @property
    def means(self):
        """Slices x,y,dp centre of mass."""
        return self._means.reshape((3, self._nbunch, self.nslice, self._dturns))

    @property
    def spos(self):
        """Slices s position."""
        return self._spos.reshape((self._nbunch, self.nslice, self._dturns))

    @property
    def weights(self):
        """Slices weights in mA if beam current >0,
        otherwise fraction of total number of
        particles in the bunch.
        """
        return self._weights.reshape((self._nbunch, self.nslice, self._dturns))

    @property
    def startturn(self):
        """Start turn of the acquisition."""
        return self._startturn

    @startturn.setter
    def startturn(self, value):
        if value < 0:
            msg = "start-turn must be greater or equal to 0"
            raise ValueError(msg)
        if value >= self._endturn:
            msg = "start-turn must be smaller than endturn"
            raise ValueError(msg)
        self._startturn = value

    @property
    def endturn(self):
        """End turn of the acquisition."""
        return self._endturn

    @endturn.setter
    def endturn(self, value):
        if value <= 0:
            msg = "end-turn must be greater than 0"
            raise ValueError(msg)
        if value <= self._startturn:
            msg = "end-turn must be greater than startturn"
            raise ValueError(msg)
        self._endturn = value


class Aperture(Element):
    """Transverse aperture element."""

    _BUILD_ATTRIBUTES = [*Element._BUILD_ATTRIBUTES, "Limits"]
    _conversions = dict(Element._conversions, Limits=lambda v: _array(v, (4,)))

    def __init__(self, family_name, limits, **kwargs):
        """
        Args:
            family_name:    Name of the element
            limits:         (4,) array of physical aperture:
              [xmin, xmax, ymin, ymax]
        Default PassMethod: ``AperturePass``.
        """
        kwargs.setdefault("PassMethod", "AperturePass")
        super().__init__(family_name, Limits=limits, **kwargs)


class LongtAperture(Element):
    """Longitudinal aperture element."""

    _BUILD_ATTRIBUTES = [*Element._BUILD_ATTRIBUTES, "Limits"]
    _conversions = dict(Element._conversions, Limits=lambda v: _array(v, (4,)))

    def __init__(self, family_name, limits, **kwargs):
        """
        Args:
            family_name:    Name of the element
            limits:         (4,) array of physical aperture:
              [dpmin, dpmax, ctmin, ctmax]
        Default PassMethod: ``LongtAperturePass``.
        """
        kwargs.setdefault("PassMethod", "LongtAperturePass")
        super().__init__(family_name, Limits=limits, **kwargs)


class Drift(LongElement):
    """Drift space element."""

    def __init__(self, family_name: str, length: float, **kwargs):
        """
        Args:
            family_name:    Name of the element
            length:         Element length [m].

        Default PassMethod: ``DriftPass``
        """
        kwargs.setdefault("PassMethod", "DriftPass")
        super().__init__(family_name, length, **kwargs)

    def insert(
        self, insert_list: Iterable[tuple[float, Element | None]]
    ) -> list[Element]:
        # noinspection PyUnresolvedReferences
        """insert elements inside a drift.

        Arguments:
            insert_list: iterable, each item of insert_list is itself an
              iterable with 2 objects:

              1. the location where the centre of the element
                 will be inserted, given as a fraction of the Drift length.
              2. an element to be inserted at that location. If :py:obj:`None`,
                 the drift will be divided but no element will be inserted.

        Returns:
             elem_list: a list of elements.

        Drifts with negative lengths may be generated if necessary.

        Examples:

            >>> Drift("dr", 2.0).insert(((0.25, None), (0.75, None)))
            [Drift('dr', 0.5), Drift('dr', 1.0), Drift('dr', 0.5)]

            >>> Drift("dr", 2.0).insert(((0.0, Marker("m1")), (0.5, Marker("m2"))))
            [Marker('m1'), Drift('dr', 1.0), Marker('m2'), Drift('dr', 1.0)]

            >>> Drift("dr", 2.0).insert(((0.5, Quadrupole("qp", 0.4, 0.0)),))
            [Drift('dr', 0.8), Quadrupole('qp', 0.4), Drift('dr', 0.8)]
        """
        frac, elements = zip(*insert_list, strict=True)
        lg = [0.0 if el is None else el.Length for el in elements]
        fr = np.asarray(frac, dtype=float)
        lg = 0.5 * np.asarray(lg, dtype=float) / self.Length
        drfrac = np.hstack((fr - lg, 1.0)) - np.hstack((0.0, fr + lg))
        long_elems = drfrac != 0.0
        drifts = np.ndarray((len(drfrac),), dtype="O")
        drifts[long_elems] = self.divide(drfrac[long_elems])
        nline = len(drifts) + len(elements)
        line = [None] * nline  # type: list[Optional[Element]]
        line[::2] = drifts
        line[1::2] = elements
        return [el for el in line if el is not None]


class Collimator(Drift):
    """Collimator element."""

    _BUILD_ATTRIBUTES = [*LongElement._BUILD_ATTRIBUTES, "RApertures"]

    def __init__(self, family_name: str, length: float, limits, **kwargs):
        """
        Args:
            family_name:    Name of the element
            length:         Element length [m]
            limits:         (4,) array of physical aperture:
              [xmin, xmax, zmin, zmax] [m].

        Default PassMethod: ``DriftPass``
        """
        super().__init__(family_name, length, RApertures=limits, **kwargs)


class RFCavity(LongtMotion, LongElement):
    """RF cavity element."""

    _BUILD_ATTRIBUTES = [
        *LongElement._BUILD_ATTRIBUTES,
        "Voltage",
        "Frequency",
        "HarmNumber",
        "Energy",
    ]
    default_pass = {False: "DriftPass", True: "RFCavityPass"}
    _conversions = dict(
        LongElement._conversions,
        Voltage=float,
        Frequency=float,
        HarmNumber=int,
        TimeLag=float,
    )

    def __init__(
        self,
        family_name: str,
        length: float,
        voltage: float,
        frequency: float,
        harmonic_number: int,
        energy: float,
        **kwargs,
    ):
        """
        Args:
            family_name:    Name of the element
            length:         Element length [m]
            voltage:        RF voltage [V]
            frequency:      RF frequency [Hz]
            harmonic_number:
            energy:         ring energy [eV].

        Keyword Arguments:
            TimeLag=0:      Cavity time lag

        Default PassMethod: ``RFCavityPass``
        """
        kwargs.setdefault("TimeLag", 0.0)
        kwargs.setdefault("PassMethod", self.default_pass[True])
        super().__init__(
            family_name,
            length,
            Voltage=voltage,
            Frequency=frequency,
            HarmNumber=harmonic_number,
            Energy=energy,
            **kwargs,
        )

    def _part(self, fr, sumfr):
        pp = super()._part(fr, sumfr)
        pp.Voltage = fr * self.Voltage
        return pp

    def is_compatible(self, other) -> bool:
        return (
            super().is_compatible(other)
            and self.Frequency == other.Frequency
            and self.TimeLag == other.TimeLag
        )

    def merge(self, other) -> None:
        super().merge(other)
        self.Voltage += other.Voltage

    def _get_longt_motion(self):
        return self.PassMethod.endswith("CavityPass")

    # noinspection PyShadowingNames
    def set_longt_motion(self, enable, new_pass=None, **kwargs):
        if new_pass == "auto":
            new_pass = (
                self.default_pass[True]
                if enable
                else ("IdentityPass" if self.Length == 0 else "DriftPass")
            )
        return super().set_longt_motion(enable, new_pass=new_pass, **kwargs)


class M66(Element):
    """Linear (6, 6) transfer matrix."""

    _BUILD_ATTRIBUTES = [*Element._BUILD_ATTRIBUTES, "M66"]
    _conversions = dict(Element._conversions, M66=_array66)

    def __init__(self, family_name: str, m66=None, **kwargs):
        """
        Args:
            family_name:    Name of the element
            m66:            Transfer matrix. Default: Identity matrix.

        Default PassMethod: ``Matrix66Pass``
        """
        if m66 is None:
            m66 = np.identity(6)
        kwargs.setdefault("PassMethod", "Matrix66Pass")
        kwargs.setdefault("M66", m66)
        super().__init__(family_name, **kwargs)


class SimpleQuantDiff(_DictLongtMotion, Element):
    """
    Linear tracking element for a simplified quantum diffusion,
    radiation damping and energy loss.

    Note: The damping times are needed to compute the correct
    kick for the emittance. Radiation damping is NOT applied.
    """

    _BUILD_ATTRIBUTES = Element._BUILD_ATTRIBUTES
    default_pass = {False: "IdentityPass", True: "SimpleQuantDiffPass"}

    def __init__(
        self,
        family_name: str,
        betax: float = 1.0,
        betay: float = 1.0,
        emitx: float = 0.0,
        emity: float = 0.0,
        espread: float = 0.0,
        taux: float = 0.0,
        tauy: float = 0.0,
        tauz: float = 0.0,
        **kwargs,
    ):
        """
        Args:
            family_name:    Name of the element
            betax:         Horizontal beta function at element [m]
            betay:         Vertical beta function at element [m]
            emitx:         Horizontal equilibrium emittance [m.rad]
            emity:         Vertical equilibrium emittance [m.rad]
            espread:       Equilibrium energy spread
            taux:          Horizontal damping time [turns]
            tauy:          Vertical damping time [turns]
            tauz:          Longitudinal damping time [turns].

        Default PassMethod: ``SimpleQuantDiffPass``
        """
        kwargs.setdefault("PassMethod", self.default_pass[True])

        assert taux >= 0.0, "taux must be greater than or equal to 0"
        self.taux = taux

        assert tauy >= 0.0, "tauy must be greater than or equal to 0"
        self.tauy = tauy

        assert tauz >= 0.0, "tauz must be greater than or equal to 0"
        self.tauz = tauz

        assert emitx >= 0.0, "emitx must be greater than or equal to 0"
        self.emitx = emitx
        if emitx > 0.0:
            assert taux > 0.0, "if emitx is given, taux must be non zero"

        assert emity >= 0.0, "emity must be greater than or equal to 0"
        self.emity = emity
        if emity > 0.0:
            assert tauy > 0.0, "if emity is given, tauy must be non zero"

        assert espread >= 0.0, "espread must be greater than or equal to 0"
        self.espread = espread
        if espread > 0.0:
            assert tauz > 0.0, "if espread is given, tauz must be non zero"

        self.betax = betax
        self.betay = betay
        super().__init__(family_name, **kwargs)


class SimpleRadiation(_DictLongtMotion, Radiative, Element):
    """Simple radiation damping and energy loss."""

    _BUILD_ATTRIBUTES = Element._BUILD_ATTRIBUTES
    _conversions = dict(
        Element._conversions, U0=float, damp_mat_diag=lambda v: _array(v, shape=(6,))
    )

    default_pass = {False: "IdentityPass", True: "SimpleRadiationRadPass"}

    def __init__(
        self,
        family_name: str,
        taux: float = 0.0,
        tauy: float = 0.0,
        tauz: float = 0.0,
        U0: float = 0.0,
        **kwargs,
    ):
        """
        Args:
            family_name:    Name of the element
            taux:          Horizontal damping time [turns]
            tauy:          Vertical damping time [turns]
            tauz:          Longitudinal damping time [turns]
            U0:            Energy loss per turn [eV].

        Default PassMethod: ``SimpleRadiationRadPass``
        """
        assert taux >= 0.0, "taux must be greater than or equal to 0"
        dampx = 1 if taux == 0.0 else np.exp(-1 / taux)

        assert tauy >= 0.0, "tauy must be greater than or equal to 0"
        dampy = 1 if tauy == 0.0 else np.exp(-1 / tauy)

        assert tauz >= 0.0, "tauz must be greater than or equal to 0"
        dampz = 1 if tauz == 0.0 else np.exp(-1 / tauz)

        kwargs.setdefault("PassMethod", self.default_pass[True])
        kwargs.setdefault("U0", U0)
        kwargs.setdefault(
            "damp_mat_diag", np.array([dampx, dampx, dampy, dampy, dampz, dampz])
        )

        super().__init__(family_name, **kwargs)


class QuantumDiffusion(_DictLongtMotion, Element):
    _BUILD_ATTRIBUTES = [*Element._BUILD_ATTRIBUTES, "Lmatp"]
    default_pass = {False: "IdentityPass", True: "QuantDiffPass"}
    _conversions = dict(Element._conversions, Lmatp=_array66)

    def __init__(self, family_name: str, lmatp: np.ndarray, **kwargs):
        """Quantum diffusion element.

        Args:
            family_name:    Name of the element
            lmatp      :    Diffusion matrix for generation (see
              :py:func:`.gen_quantdiff_elem`)

        Default PassMethod: ``QuantDiffPass``
        """
        kwargs.setdefault("PassMethod", self.default_pass[True])
        super().__init__(family_name, Lmatp=lmatp, **kwargs)


class EnergyLoss(_DictLongtMotion, Element):
    _BUILD_ATTRIBUTES = [*Element._BUILD_ATTRIBUTES, "EnergyLoss"]
    _conversions = dict(Element._conversions, EnergyLoss=float)
    default_pass = {False: "IdentityPass", True: "EnergyLossRadPass"}

    def __init__(self, family_name: str, energy_loss: float, **kwargs):
        """Energy loss element.

        The :py:class:`EnergyLoss` element is taken into account in
        :py:func:`.radiation_parameters`: it adds damping by contributing to the
        :math:`I_2` integral, thus reducing the equilibrium emittance. But it does not
        generate any diffusion. This makes sense only if the losses summarised in
        the element occur in non-dispersive locations.

        It is a single thin, straight non-focusing radiative element that does not
        contribute to the diffusion. It's typical usage is to model the energy loss
        and contribution to the damping times from a thin wiggler located in a
        non dispersive region. More complex cases with focusing and / or  diffusion
        are not correctly handled by this element.

        Args:
            family_name:    Name of the element
            energy_loss:    Energy loss [eV]

        """
        kwargs.setdefault("PassMethod", self.default_pass[False])
        super().__init__(family_name, EnergyLoss=energy_loss, **kwargs)


Radiative.register(EnergyLoss)
