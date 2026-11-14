"""Definition of  the Xsuite-like element classes."""

from __future__ import annotations

__all__ = ["XsElement"]

import warnings
from math import pi
from collections.abc import Callable
from typing import Any, ClassVar, Protocol

import numpy as np

from ..lattice import AtWarning, ReferencePoint, constants as cst
from ..lattice import elements as elt
from .madx import poly_from_mad, poly_to_mad

_INDEX_TO_EDGE_MODEL = {
    -1: "suppressed",
    0: "linear",
    1: "full",
    2: "dipole-only",
}
_EDGE_MODEL_TO_INDEX = {k: v for v, k in _INDEX_TO_EDGE_MODEL.items()}


class _XsFactory(Protocol):
    """Base class for Xsuite element factories."""

    @classmethod
    def from_at(cls, atelem: elt.Element, **kwargs) -> XsElement: ...


class XsElement(dict, _XsFactory):
    """Base class for Xsuite elements."""

    # Class attributes
    _atClass: ClassVar[type[elt.Element]] = elt.Marker
    _xsuite2at_attr: ClassVar[dict[str, str]] = {"length": "Length"}
    _at2xsuite_model: ClassVar[dict[str | None, str]] = {}
    _xsuite2at_model: ClassVar[dict[str, str | None]] = {}
    _at_integrator: ClassVar[str | None] = None
    _allow_rot_shift: ClassVar[bool] = False

    # Instance attributes
    name: str  #: Element name
    origin: str  #: Xsuite class of the source element

    def __init__(self, name: str, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setdefault("length", 0.0)
        self.name = name
        self.origin = kwargs.get("__class__", "")

    def _warn(self, msg: str):
        """Issue a warning with the element name."""
        cname = type(self).__name__
        warnings.warn(AtWarning(f"Element {cname}({self.name}): {msg}"), stacklevel=3)

    @property
    def _integration_length(self):
        return self["length"]

    def _set_xs_model(self, atelem) -> None:
        """Set the model for the Xsuite element."""
        xs_model = self._at2xsuite_model.get(atelem.PassMethod, None)
        if xs_model is None:
            xs_model = self._at2xsuite_model.get(None, None)
        if xs_model is not None:
            self["model"] = xs_model

    def _set_xs_integrator(self, atelem) -> None:
        """Set the integrator for the Xsuite element."""
        if (integrator := self.get("_at_integrator", None)) is not None:
            self["integrator"] = integrator

    def _set_xs_transforms(self, atelem) -> None:
        """Generate Xsuite element displacements."""
        tr = {
            "shift_x": getattr(atelem, "_dx", 0.0),
            "shift_y": getattr(atelem, "_dy", 0.0),
            "shift_s": getattr(atelem, "_dz", 0.0),
            "rot_s_rad_no_frame": getattr(atelem, "_tilt", 0.0),
            "rot_x_rad": -getattr(atelem, "_pitch", 0.0),
            "rot_y_rad": getattr(atelem, "_yaw", 0.0),
        }
        misalign = {k: v for k, v in tr.items() if v != 0.0}
        if misalign and not self._allow_rot_shift:
            msg = "Misalignment and rotations are not avalaible in Xsuite, ignored"
            self._warn(msg)
            return
        refpoint = getattr(atelem, "_reference", ReferencePoint.CENTRE)
        rots_frame = getattr(atelem, "_tilt_frame", 0.0)

        if rots_frame != 0.0:
            if misalign:
                self._warn("alignment errors available only for H bends, ignored")
            misalign = {"rot_s_rad": rots_frame}
        elif misalign:
            if refpoint == ReferencePoint.ENTRANCE:
                misalign["rot_shift_anchor"] = 0.0
            else:
                misalign["rot_shift_anchor"] = self._integration_length / 2.0

        self.update(misalign)

    def _set_at_passmethod(self) -> dict[str, str]:
        """Get the AT passmethod."""
        passmethod = self._xsuite2at_model.get(self.get("model"), None)
        if passmethod is not None:
            passmethod = self._xsuite2at_model.get(None, None)
        return {} if passmethod is None else {"PassMethod": passmethod}

    def _set_at_transforms(self) -> dict:
        """Generate AT element displacements."""
        tr = {
            "dx": self.get("shift_x", 0.0),
            "dy": self.get("shift_y", 0.0),
            "dz": self.get("shift_s", 0.0),
            "pitch": -self.get("rot_x_rad", 0.0),
            "yaw": self.get("rot_y_rad", 0.0),
            "tilt": self.get("rot_s_rad_no_frame", 0.0),
        }
        transforms = {k: v for k, v in tr.items() if v != 0.0}
        ismisaligned = len(transforms) > 0

        rot_s_frame = self.get("rot_s_rad", 0.0)
        if rot_s_frame != 0.0:
            if ismisaligned:
                self._warn("alignment errors available only for H bends, ignored")
            transforms = {"tilt_frame": rot_s_frame}
        elif ismisaligned:
            length = self._integration_length
            rot_shift_anchor = self.get("rot_shift_anchor", 0.0)
            anchor = 0.5 if length == 0.0 else rot_shift_anchor / length
            if anchor == 0.0:
                reference = ReferencePoint.ENTRANCE
            elif abs(anchor - 0.5) < 1.0e-6:
                reference = ReferencePoint.CENTRE
            else:
                msg = (
                    f"Anchor point for rotation different from 0 or {length / 2.0},\n"
                    "setting reference to AT default: ReferencePoint.CENTRE"
                )
                self._warn(msg)
                reference = ReferencePoint.CENTRE
            transforms["reference"] = reference.value
        return transforms

    def _set_at_params(self, **kwargs) -> dict:
        """Translate Xsuite parameters to AT parameters.

        Args:
            beta0:          beta of the reference particle
        """
        atparams = {
            atk: v
            for k, v in self.items()
            if (atk := self._xsuite2at_attr.get(k, None)) is not None
        }
        atparams["FamName"] = self.name
        if self.origin:
            atparams["origin"] = self.origin
        atparams.update(self._set_at_passmethod())
        atparams.update(self._set_at_transforms())
        return atparams

    def _class_to_at(self, atparams: dict[str, Any]) -> type[elt.Element]:
        """Get the AT element class."""
        return self._atClass

    def to_dict(self) -> dict[str, Any]:
        """Build a dictionary representation of the XsElement.

        Returns:
            elem_dict:  dictionary representation of the :py:class:`XsElement`
        """
        return self

    def to_at(self, beta0: float = 1.0) -> elt.Element:
        """Generate the AT element.

        Args:
            beta0:          beta of the reference particle

        Returns:
            elem:   new :py:class:`.Element` object
        """
        elem_kwargs = self._set_at_params(beta0=beta0)
        cls = self._class_to_at(elem_kwargs)
        return cls.init_from_dict(**elem_kwargs)

    # noinspection PyUnusedLocal,unused-parameter
    @classmethod
    def from_dict(
        cls,
        xsparams: dict[str, Any],
        name: str = "?",
        do_warn: Callable[[str], bool] = lambda clname: True,
    ) -> XsElement:
        """Build a XsElement from its dictionary representation.

        Args:
            xsparams:   dictionary representation of the XsElement
            name:       Optional name of the element
            do_warn:    Function checking if a warning should be issued.

        Returns:
            xselement:  new :py:class:`XsElement` object
        """
        return cls(name=name, **xsparams)

    @classmethod
    def from_at(
        cls,
        atelem: elt.Element,
        *,
        match_model: bool = True,
        **kwargs,
    ) -> XsElement:
        """Build a XsElement element from an AT element.

        Args:
            atelem:         AT :py:class:`.Element`
            match_model:    If :py:obj:`True`, set the Xsuite model and integrator
              matching at best the AT PassMethod. Otherwise, use Xsuite defaults,
            do_warn:        Function checking if a warning should be issued.
            beta0:          beta of the reference particle
        Returns:
            xselement:  new :py:class:`XsElement` object
        """
        # Conversion of similar attributes
        xsparams = {
            kxs: v
            for kxs, kat in cls._xsuite2at_attr.items()
            if (v := getattr(atelem, kat, None)) is not None
        }
        # Set the Xsuite class
        xsparams["__class__"] = cls.__name__
        elem = cls(name=atelem.FamName, **xsparams)
        if match_model:
            elem._set_xs_model(atelem)
            elem._set_xs_integrator(atelem)
        elem._set_xs_transforms(atelem)
        return elem

    @staticmethod
    def static_from_dict(
        elem_dict: dict[str, Any],
        name: str = "?",
        do_warn: Callable[[str], bool] = lambda clname: True,
    ) -> XsElement:
        """Build a XsElement from its dictionary representation.

        Args:
            elem_dict:  dictionary representation of the XsElement
            name:       Optional name of the element
            do_warn:    Function checking if a warning should be issued.

        Returns:
            xselement:  new :py:class:`XsElement` object
        """
        xsclass = _xsclass.get(elem_dict.get("__class__", "Unknown"), NotInAT)
        return xsclass.from_dict(elem_dict, name=name, do_warn=do_warn)

    @staticmethod
    def static_from_at(
        atelem: elt.Element,
        *,
        beta0=1.0,
        match_model: bool = True,
        do_warn: Callable[[str], bool] = lambda x: True,
    ) -> XsElement:
        """Build a XsElement from an AT element.

        Args:
            atelem:         AT :py:class:`.Element`
            beta0:          beta of the reference particle
            match_model:    If :py:obj:`True`, set the Xsuite model and integrator
              matching at best the AT PassMethod. If :py:obj:`False`, the Xsuite default
              model and integrator are be used.
            do_warn:        Function to checking if a warning should be issued.

        Returns:
            xselement:  new :py:class:`XsElement` object
        """
        xsclass = _at2xsclass.get(type(atelem), NotInXsuite)
        return xsclass.from_at(
            atelem, beta0=beta0, match_model=match_model, do_warn=do_warn
        )


class Marker(XsElement):
    """Xsuite Marker element."""

    # Class variables
    _atClass: ClassVar[type[elt.Element]] = elt.Marker


class LimitRect(Marker):
    """Xsuite LimitRect element."""

    # Class variables
    _allow_rot_shift = True

    def _set_at_params(self, **kwargs) -> dict[str, Any]:
        atparams = super()._set_at_params(**kwargs)
        atparams["RApertures"] = [
            self.get("min_x", -np.inf),
            self.get("max_x", np.inf),
            self.get("min_y", -np.inf),
            self.get("max_y", np.inf),
        ]
        return atparams

    @classmethod
    def from_at(cls, atelem: elt.Element, **kwargs) -> XsElement:
        elem = super().from_at(atelem, **kwargs)
        ap = atelem.RApertures
        elem["min_x"] = ap[0]
        elem["max_x"] = ap[1]
        elem["min_y"] = ap[2]
        elem["max_y"] = ap[3]
        return elem


class LimitEllipse(Marker):
    """Xsuite LimitEllipse element."""

    # Class variables
    _allow_rot_shift = True

    def _set_at_params(self, **kwargs) -> dict[str, Any]:
        atparams = super()._set_at_params(**kwargs)
        atparams["EApertures"] = [self.get("a", np.inf), self.get("b", np.inf)]
        return atparams

    @classmethod
    def from_at(cls, atelem: elt.Element, **kwargs) -> XsElement:
        elem = super().from_at(atelem, **kwargs)
        ap = atelem.EApertures
        elem["a"] = ap[0]
        elem["b"] = ap[1]
        return elem


class Drift(XsElement):
    """Xsuite Drift element."""

    # Class variables
    _atClass = elt.Drift
    _at2xsuite_model = {"ExactDriftPass": "exact", None: "expanded"}
    _xsuite2at_model = {v: k for k, v in _at2xsuite_model.items()}


class Multipole(XsElement):
    """Base class for Xsuite magnetic elements."""

    # Class variables
    _atClass: ClassVar[type[elt.Element]] = elt.Multipole
    _at2xsuite_model = {
        "ExactMultipoleRadPass": "drift-kick-drift-exact",
        "ExactMultipolePass": "drift-kick-drift-exact",
        None: "drift-kick-drift-expanded",  # Default Xsuite model
        "StrMPoleSymplectic6Pass": None,  # Default AT passmethod
    }
    _xsuite2at_model = {v: k for k, v in _at2xsuite_model.items()}
    _at_integrator = "yoshida-6"
    _xsuite2at_attr = XsElement._xsuite2at_attr | {
        "order": "MaxOrder",
        "num_multipole_kicks": "NumIntKicks",
    }
    _allow_rot_shift = True

    def _set_at_poly(self) -> dict[str, Any]:
        """Generate the AT field expansion."""

        def xspoly(kmain: list[str], kerr: str) -> tuple[int, np.ndarray]:
            # Get main components
            pmain = np.array([self.get(k, 0.0) for k in kmain])
            # Get errors
            perr = np.array(self.get(kerr, np.zeros(4)))
            if length != 0.0:
                perr /= length
            # Pad to the same length
            lpoly = max(4, len(perr))
            poly = np.pad(pmain, (0, lpoly - 4)) + np.pad(perr, (0, lpoly - len(perr)))
            # Use the AT default order if possible
            porder = atorder if np.all(poly[atorder + 1 :] == 0.0) else xsorder
            # Convert to AT
            return porder, np.fromiter(poly_from_mad(poly), dtype=float, count=lpoly)

        length = self._integration_length
        xsorder: int = self.get("order", 0)
        atorder: int = getattr(self._atClass, "DefaultOrder", xsorder)
        aorder, polya = xspoly(["k0s", "k1s", "k2s", "k3s"], "ksl")
        border, polyb = xspoly(["k0", "k1", "k2", "k3"], "knl")
        maxorder = max(aorder, border)
        atparams = {
            "MaxOrder": maxorder,
            "PolynomB": polyb[: maxorder + 1],
            "PolynomA": polya[: maxorder + 1],
        }
        if (taper := self.get("delta_taper")) is not None:
            atparams["FieldScaling"] = 1.0 + taper
        return atparams

    def _set_at_fringe(self) -> dict[str, Any]:
        """generate the AT fringe field description."""
        atparams = {}
        if self.get("edge_entry_active", False):
            atparams["FringeQuadEntrance"] = 1
        if self.get("edge_exit_active", False):
            atparams["FringeQuadExit"] = 1
        return atparams

    def _set_xs_poly(self, atelem, match_model: bool = True) -> None:
        """Generate the AT field expansion."""

        def extract(poly, ord):
            try:
                v = poly[ord]
            except IndexError:
                v = 0.0
            else:
                poly[ord] = 0.0
            return v

        pata = getattr(atelem, "PolynomA", np.zeros(4))
        pola = np.fromiter(poly_to_mad(pata), dtype=float, count=pata.size)
        patb = getattr(atelem, "PolynomB", np.zeros(4))
        polb = np.fromiter(poly_to_mad(patb), dtype=float, count=patb.size)
        korder = getattr(self, "_mag_order", None)
        if korder == 0:  # dipole
            self["k1"] = extract(polb, 1)
            self["k2"] = extract(polb, 2)
        elif korder is not None:  # quadrupole, sextupole, octupole
            self["k" + str(korder)] = extract(polb, korder)
            self["k" + str(korder) + "s"] = extract(pola, korder)
        if (length := self._integration_length) > 0.0:
            polb *= length
            pola *= length
        if np.any(pola) or np.any(polb):
            self["knl"] = list(polb)
            self["ksl"] = list(pola)
        if (scaling := getattr(atelem, "FieldScaling", None)) is not None:
            self["delta_taper"] = scaling - 1.0
        self["_isthick"] = self["length"] != 0.0

    def _set_xs_fringe(self, atelem) -> None:
        self["edge_entry_active"] = bool(getattr(atelem, "FringeQuadEntrance", 0))
        self["edge_exit_active"] = bool(getattr(atelem, "FringeQuadExit", 0))

    def _class_to_at(self, atparams: dict[str, Any]) -> type[elt.Element]:
        if atparams.get("Length", 0.0) == 0.0:
            return elt.ThinMultipole
        else:
            return super()._class_to_at(atparams)

    def _set_at_params(self, **kwargs) -> dict[str, Any]:
        atparams = super()._set_at_params(**kwargs)
        atparams.update(self._set_at_poly())
        atparams.update(self._set_at_fringe())
        return atparams

    @classmethod
    def from_at(
        cls, atelem: elt.Element, *, match_model: bool = True, **kwargs
    ) -> XsElement:
        elem = super().from_at(atelem, match_model=match_model, **kwargs)
        elem._set_xs_poly(atelem, match_model=match_model)
        elem._set_xs_fringe(atelem)
        return elem


class Quadrupole(Multipole):
    """Xsuite Quadrupole element."""

    # Class variables
    _atClass = elt.Quadrupole
    _mag_order: ClassVar[int] = 1


class Sextupole(Multipole):
    """Xsuite Sextupole element."""

    # Class variables
    _atClass = elt.Sextupole
    _mag_order: ClassVar[int] = 2


class Octupole(Multipole):
    """Xsuite Octupole element."""

    # Class variables
    _atClass = elt.Octupole
    _mag_order: ClassVar[int] = 3


class Bend(Multipole):
    """Xsuite Bend element."""

    # Class variables
    _atClass = elt.Bend
    _at2xsuite_model = {
        "ExactSectorBendRadPass": "bend-kick-bend",
        "ExactSectorBendPass": "bend-kick-bend",
        None: "drift-kick-drift-expanded",  # Default Xsuite model
        "BndMPoleSymplectic6Pass": None,  # Default AT passmethod
    }
    _xsuite2at_model = {v: k for k, v in _at2xsuite_model.items()}
    _xsuite2at_attr = Multipole._xsuite2at_attr | {
        "angle": "BendingAngle",
        "edge_entry_fint": "FringeInt1",
        "edge_exit_fint": "FringeInt2",
        "edge_entry_angle": "EntranceAngle",
        "edge_exit_angle": "ExitAngle",
    }
    _mag_order = 0
    _default_bend_fringe: ClassVar[dict[bool, int]] = {
        True: 4,  # exact methods
        False: 1,  # non-exact methods
    }
    _at2xsuite_edge: ClassVar[dict[int, str]] = {
        0: "suppressed",
        1: "linear",
        2: "linear",
        3: "linear",
        4: "dipole-only",
    }
    _xsuite2at_edge: ClassVar[dict[str, tuple]] = {
        "suppressed": (0, None),
        "linear": (1, None),
        "dipole-only": (None, None),
        "full": (None, 1),
    }

    def _set_at_fringe(self) -> dict[str, Any]:
        def edge_model(xskey, bendkey, quadkey):
            [bend, quad] = self._xsuite2at_edge[self.get(xskey, "linear")]
            if bend is not None:
                atparams[bendkey] = bend
            if quad is not None:
                atparams[quadkey] = quad

        atparams = {}
        entry_hgap = self.get("edge_entry_hgap")
        exit_hgap = self.get("edge_exit_hgap")
        if entry_hgap != exit_hgap:
            self._warn("Entry and Exit gaps for dipole are different, keep entry")
        if entry_hgap is not None:
            atparams["FullGap"] = 2.0 * entry_hgap

        edge_model("edge_entry_model", "FringeBendEntrance", "FringeQuadEntrance")
        edge_model("edge_exit_model", "FringeBendExit", "FringeQuadExit")
        return atparams

    def _set_xs_fringe(self, atelem):
        def edge_model(bendkey, quadkey):
            at_quad_fringe = getattr(atelem, quadkey, 0)
            if at_quad_fringe > 0:
                return "full"
            else:
                return self._at2xsuite_edge[getattr(atelem, bendkey, def_bend_fringe)]

        if (gap := getattr(atelem, "FullGap", None)) is not None:
            self["edge_entry_hgap"] = 0.5 * gap
            self["edge_exit_hgap"] = 0.5 * gap
        exact = atelem.PassMethod.startswith("Exact")
        def_bend_fringe = self._default_bend_fringe[exact]
        self["edge_entry_model"] = edge_model(
            "FringeBendEntrance", "FringeQuadEntrance"
        )
        self["edge_exit_model"] = edge_model("FringeBendExit", "FringeQuadExit")

    def _set_at_params(self, **kwargs) -> dict[str, Any]:
        atparams = {"EntranceAngle": 0.0, "ExitAngle": 0.0}
        atparams.update(super()._set_at_params(**kwargs))
        if not self.get("k0_from_h", "k0" not in self):
            atparams["PolynomB"][0] -= self["angle"] / self["length"]
        return atparams

    @classmethod
    def from_at(
        cls, atelem: elt.Element, *, match_model: bool = True, **kwargs
    ) -> XsElement:
        elem = super().from_at(atelem, match_model=match_model, **kwargs)
        elem["k0_from_h"] = True
        return elem


class RBend(Bend):
    """Xsuite RBend element."""

    _at2xsuite_model = {
        None: "drift-kick-drift-exact",  # Default Xsuite model
        "ExactRectangularBendPass": None,  # Default AT passmethod,
    }
    _xsuite2at_model = {v: k for k, v in _at2xsuite_model.items()}

    @property
    def _integration_length(self):
        length = self["length"]
        hangle = 0.5 * self["angle"]
        return length * np.sinc(hangle / pi)

    def _set_at_passmethod(self) -> dict[str, str]:
        """Get the AT passmethod."""
        straight = self.get("rbend_model", "adaptive") == "straight-body"
        if straight:
            passmethod = self._xsuite2at_model.get(self.get("model"), None)
        else:
            passmethod = Bend._xsuite2at_model.get(self.get("model"), None)
        return {} if passmethod is None else {"PassMethod": passmethod}

    def _set_at_params(self, **kwargs) -> dict[str, Any]:
        atparams = super()._set_at_params(**kwargs)
        straight = self.get("rbend_model", "adaptive") == "straight-body"
        exact = atparams.get("PassMethod", "").startswith("Exact")
        hangle = 0.5 * self["angle"]
        length = self.get("length_straight", 0.0) / np.sinc(hangle / pi)
        if (
            straight
            and not exact
            and (
                np.any(atparams["PolynomA"] != 0.0)
                or np.any(atparams["PolynomB"] != 0.0)
            )
        ):
            msg = (
                "Multipole are present in the rectangular bend: "
                "PassMethod is forced to ExactRectangularBendPass"
            )
            self._warn(msg)
            atparams["PassMethod"] = "ExactRectangularBendPass"
        if self.get("rbend_angle_diff", 0.0) != 0.0:
            self._warn('"rbend_angle_diff" is not supported in AT, ignored.')
        atparams["Length"] = length
        atparams["EntranceAngle"] = self.get("edge_entry_angle", 0.0) + hangle
        atparams["ExitAngle"] = self.get("edge_exit_angle", 0.0) + hangle
        return atparams

    def to_at(self, **kwargs) -> elt.Element:
        atelem = super().to_at(**kwargs)
        atelem.rbendtune()
        return atelem

    @classmethod
    def from_at(
        cls, atelem: elt.Element, *, match_model: bool = True, **kwargs
    ) -> XsElement:

        def X0ref_to_rbend_shift(hangle, length, X0ref):
            return -X0ref - 0.125 * length * hangle * np.sinc(0.5 * hangle / pi) ** 2

        elem = super().from_at(atelem, match_model=match_model, **kwargs)
        hangle = 0.5 * elem["angle"]
        length = elem.get("length")
        X0ref = getattr(atelem, "X0ref", 0.0)
        elem["rbend_model"] = "straight-body"
        elem["edge_entry_angle"] -= hangle
        elem["edge_exit_angle"] -= hangle
        elem["length_straight"] = elem._integration_length
        elem["rbend_shift"] = X0ref_to_rbend_shift(hangle, length, X0ref)
        return elem


class Cavity(XsElement):
    # Class variables
    _atClass = elt.RFCavity
    _at_integrator = "yoshida-6"
    _xsuite2at_attr = XsElement._xsuite2at_attr | {
        "voltage": "Voltage",
        "frequency": "Frequency",
    }
    _allow_rot_shift: ClassVar[bool] = True

    def _set_xs_lag(self, atelem, beta0):
        freq = atelem.Frequency
        tl = getattr(atelem, "TimeLag", 0.0)
        phase = 2 * pi * (0.5 - tl * freq / beta0 / cst.clight)
        self["phase"] = phase

    def _set_at_params(self, **kwargs) -> dict[str, Any]:
        atparams = {"Voltage": 0.0, "Frequency": np.nan, "Energy": 0.0}
        atparams.update(super()._set_at_params(**kwargs))
        length = atparams.get("Length", 0.0)
        atparams["PassMethod"] = "IdentityPass" if length == 0.0 else "DriftPass"
        atparams["HarmNumber"] = self.get("harmonic", 0)
        if (phaselag := self.get("phase")) is not None:
            atparams["_phaselag"] = pi - phaselag
        return atparams

    @classmethod
    def from_at(cls, atelem: elt.Element, *, beta0: float = 1.0, **kwargs) -> XsElement:
        elem = super().from_at(atelem, beta0=beta0, **kwargs)
        elem._set_xs_lag(atelem, beta0)
        return elem


class CrabCavity(XsElement):
    # Class variables
    _atClass = elt.CrabCavity
    _at_integrator = "yoshida-6"
    _xsuite2at_attr = XsElement._xsuite2at_attr | {
        "frequency": "Frequency",
    }
    _allow_rot_shift: ClassVar[bool] = True

    def _set_xs_lag(self, atelem, beta0):
        freq = atelem.Frequency
        tl = getattr(atelem, "TimeLag", 0.0)
        phase = 2 * pi * tl * freq / beta0 / cst.clight
        self["phase"] = phase

    def _set_at_params(self, **kwargs) -> dict[str, Any]:
        atparams = {"Frequency": 0.0}
        atparams.update(super()._set_at_params(**kwargs))
        length = atparams.get("Length", 0.0)
        atparams["Voltages"] = [self.get("crab_voltage", 0.0), 0.0]
        atparams["PassMethod"] = "IdentityPass" if length == 0.0 else "DriftPass"
        atparams["HarmNumber"] = self.get("harmonic", 0)
        if (phaselag := self.get("phase")) is not None:
            atparams["_phaselag"] = phaselag
        return atparams

    @classmethod
    def from_at(cls, atelem: elt.Element, *, beta0: float = 1.0, **kwargs) -> XsElement:
        elem = super().from_at(atelem, beta0=beta0, **kwargs)
        elem._set_xs_lag(atelem, beta0)
        elem["crab_voltage"] = atelem.Voltages[0]
        if atelem.Voltages[1] != 0.0:
            msg = "Vertical voltage not available in Xsuite, ignored"
            elem._warn(msg)
        return elem


class FirstOrderTaylorMap(XsElement):
    r"""Class for Xsuite elements with first-order Taylor map.

    The coordinate systems used in transfer matrices are different between AT and
    Xsuite: :math:`(x, p_x, y, p_y, \delta, \beta c t)` in AT and
    :math:`(x, p_x, y, p_y, \tau, p_{\tau})` in Xsuite. This class makes the conversion,
    assuming :math:`p_{\tau} = \beta_0 \delta` and :math:`\tau = -\beta c t / \beta_0`.
    """

    _atClass = elt.M66

    def _set_at_params(self, beta0: float = 1.0, **kwargs) -> dict[str, Any]:
        atparams = super()._set_at_params(beta0=beta0, **kwargs)
        if (m1 := self.get("m1", None)) is not None:
            m66 = np.asfortranarray(m1)
            # Swap columns 4 and 5
            tmp = m66[:, 5].copy()
            m66[:, 5] = -m66[:, 4] / beta0
            m66[:, 4] = tmp * beta0
            # Swap lines 4 and 5
            tmp = m66[5, :].copy()
            m66[5, :] = -m66[4, :] * beta0
            m66[4, :] = tmp / beta0
            atparams["M66"] = m66
        if (m0 := self.get("m0", None)) is not None:
            m6 = np.asfortranarray(m0)
            # Swap items 4 and 5
            tmp = m6[5]
            m6[5] = -m6[4] * beta0
            m6[4] = tmp / beta0
            atparams["M6"] = m6
        return atparams

    @classmethod
    def from_at(cls, atelem: elt.Element, *, beta0: float = 1.0, **kwargs) -> XsElement:
        elem = super().from_at(atelem, beta0=beta0, **kwargs)
        m1 = atelem.M66.copy(order="C")
        # Swap columns 4 and 5
        tmp = m1[:, 5].copy()
        m1[:, 5] = m1[:, 4] / beta0
        m1[:, 4] = -tmp * beta0
        # Swap lines 4 and 5
        tmp = m1[5, :].copy()
        m1[5, :] = m1[4, :] * beta0
        m1[4, :] = -tmp / beta0
        elem["m1"] = m1
        if (m6 := getattr(atelem, "M6", None)) is not None:
            # Swap items 4 and 5
            m0 = m6.copy(order="C")
            tmp = m0[5]
            m0[5] = m0[4] * beta0
            m0[4] = -tmp / beta0
            elem["m0"] = m0
        return elem


class NotInAT:
    """Class for Xsuite elements without AT equivalent."""

    @classmethod
    def from_dict(
        cls,
        xsparams: dict[str, Any],
        name: str = "?",
        do_warn: Callable[[str], bool] = lambda clname: True,
    ) -> XsElement:
        origin = xsparams.get("__class__", "Unknown")
        xsclass = Marker if xsparams.get("length", 0.0) == 0.0 else Drift
        if do_warn(origin):
            repl = xsclass.__name__
            msg = f"Element {name!r}: unknown class {origin} replaced by {repl}"
            warnings.warn(AtWarning(msg), stacklevel=5)
        return xsclass.from_dict(xsparams, name=name, do_warn=do_warn)


class NotInXsuite:
    """Class for AT elements without Xsuite equivalent."""

    @classmethod
    def from_at(
        cls,
        atelem: elt.Element,
        *,
        do_warn: Callable[[str], bool] = lambda x: True,
        **kwargs,
    ) -> XsElement:
        atclass = atelem.__class__.__name__
        xsclass = Marker if atelem.Length == 0.0 else Drift
        if do_warn(atclass):
            name = atelem.FamName
            repl = xsclass.__name__
            msg = f"Element {name!r}: unknown class {atclass} replaced by {repl}"
            warnings.warn(AtWarning(msg), stacklevel=5)
        return xsclass.from_at(atelem, do_warn=do_warn, **kwargs)


class Dipole:
    """Class for handling AT dipoles."""

    @classmethod
    def from_at(cls, atelem: elt.Element, **kwargs) -> XsElement:
        if atelem.PassMethod in RBend._at2xsuite_model:
            return RBend.from_at(atelem, **kwargs)
        else:
            return Bend.from_at(atelem, **kwargs)


class Corrector:
    """Class for handling AT correctors."""

    _tr_keys: set[str] = {"dx", "dy", "dz", "tilt", "yaw", "pitch", "tilt_frame"}

    # noinspection PyPep8Naming
    @classmethod
    def from_at(cls, atelem: elt.Element, **kwargs) -> XsElement:
        at = {k: v for k in cls._tr_keys if (v := getattr(atelem, k, 0.0)) != 0.0}
        at["reference"] = atelem.ReferencePoint
        if (length := atelem.Length) == 0.0:
            atel = elt.ThinMultipole(atelem.FamName, np.zeros(1), np.zeros(1))
        else:
            atel = elt.Multipole(atelem.FamName, length, np.zeros(1), np.zeros(1))
        atel.KickAngle = atelem.KickAngle
        atel.transform(**at)
        return Multipole.from_at(atel, **kwargs)


class CheckAperture:
    @classmethod
    def from_at(cls, atelem: elt.Element, **kwargs) -> XsElement:
        if hasattr(atelem, "RApertures"):
            return LimitRect.from_at(atelem, **kwargs)
        elif hasattr(atelem, "EApertures"):
            return LimitEllipse.from_at(atelem, **kwargs)
        else:
            return Marker.from_at(atelem, **kwargs)


_xsclass: dict[str, type[XsElement]] = {
    "Marker": Marker,
    "Drift": Drift,
    "Cavity": Cavity,
    "Quadrupole": Quadrupole,
    "Sextupole": Sextupole,
    "Octupole": Octupole,
    "Multipole": Multipole,
    "Bend": Bend,
    "RBend": RBend,
    "LimitRect": LimitRect,
    "LimitEllipse": LimitEllipse,
    "FirstOrderTaylorMap": FirstOrderTaylorMap,
    "CrabCavity": CrabCavity,
}


_at2xsclass: dict[type[elt.Element], type[_XsFactory]] = {
    elt.Marker: CheckAperture,
    elt.Monitor: Marker,
    elt.Drift: Drift,
    elt.RFCavity: Cavity,
    elt.Quadrupole: Quadrupole,
    elt.Sextupole: Sextupole,
    elt.Octupole: Octupole,
    elt.Multipole: Multipole,
    elt.ThinMultipole: Multipole,
    elt.Dipole: Dipole,
    elt.Corrector: Corrector,
    elt.M66: FirstOrderTaylorMap,
    elt.CrabCavity: CrabCavity,
}
