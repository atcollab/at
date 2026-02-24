from __future__ import annotations

__all__ = [
    "XsElement",
    "XsLine",
    "lattice_from_line",
    "line_from_lattice",
    "load_xsuite",
    "save_xsuite",
]

import json
import warnings
from collections.abc import Generator, Iterable
from math import pi, sqrt
from pathlib import Path
from typing import Any, ClassVar, Protocol
import contextlib

import numpy as np

from ..lattice import AtWarning, ReferencePoint, constants as cst
from ..lattice import Lattice, Particle, elements as elt
from .madx import p_to_at

try:
    from xtrack import Line
except ImportError:
    # noinspection PyUnusedLocal
    def line_from_lattice(ring: Lattice, match_model: bool = False) -> Line:
        """Create a Xsuite :py:class:`.Line` from an AT :py:class:`.Lattice`.

        Args:
            ring:           AT lattice
            match_model:    if :py:obj:`True`, set the Xsuite model and integrator
              matching at best the AT PassMethod. Otherwise, use Xsuite defaults.

        Returns:
            line:           new :py:class:`.Line` object.
        """
        msg = "xtrack is not installed."
        raise ImportError(msg)
else:

    def line_from_lattice(ring: Lattice, match_model: bool = False) -> Line:
        """Create a Xsuite :py:class:`.Line` from an AT lattice.

        Args:
            ring:           AT lattice
            match_model:    if :py:obj:`True`, set the Xsuite model and integrator
              matching at best the AT PassMethod. Otherwise, use Xsuite defaults.

        Returns:
            line:           new :py:class:`.Line` object.
        """
        return Line.from_dict(XsLine.from_at(ring, match_model=match_model).to_dict())


# Default particle: electrons at 1 GeV
_default_energy = 1.0e9
_default_mass = cst.e_mass
_default_charge = -1
_default_gamma = _default_energy / _default_mass
_default_particle: dict[str, float | list[float]] = {
    "mass0": _default_mass,
    "gamma0": [_default_gamma],
    "beta0": [sqrt(1.0 - 1.0 / _default_gamma / _default_gamma)],
    "q0": _default_charge,
}

_INDEX_TO_EDGE_MODEL = {
    -1: "suppressed",
    0: "linear",
    1: "full",
    2: "dipole-only",
}
_EDGE_MODEL_TO_INDEX = {k: v for v, k in _INDEX_TO_EDGE_MODEL.items()}


class _XtLine(Protocol):
    """Defines an object having a to_dict method."""

    @classmethod
    def to_dict(cls) -> dict[str, Any]: ...


class _AtEncoder(json.JSONEncoder):
    """JSON encoder for specific AT types."""

    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, Particle):
            return obj.to_dict()
        elif isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        else:
            return super().default(obj)


class XsElement(dict):
    """Base class for Xsuite elements."""

    # Class attributes
    _atClass: ClassVar[type[elt.Element]] = elt.Marker
    _xsuite2at_attr: ClassVar[dict[str, str]] = {"length": "Length"}
    _at2xsuite_attr: ClassVar[dict[str, str]] = {
        v: k for k, v in _xsuite2at_attr.items()
    }
    _at2xsuite_model: ClassVar[dict[str | None, str]] = {}
    _xsuite2at_model: ClassVar[dict[str, str | None]] = {}
    _at_integrator: ClassVar[str | None] = None
    _index_to_model: ClassVar[dict[int, str]] = {
        0: "adaptive",
        1: "full",
        2: "bend-kick-bend",
        3: "rot-kick-rot",
        4: "mat-kick-mat",
        5: "drift-kick-drift-exact",
        6: "drift-kick-drift-expanded",
    }

    # Instance attributes
    name: str  #: Element name
    origin: str  #: Xsuite class of the source element

    def __init__(self, name: str, *args, origin: str = "", **kwargs):
        super().__init__(*args, **kwargs)
        self.setdefault("length", 0.0)
        self.name = name
        self.origin = origin

    def _get_model(self) -> str:
        try:
            # String model
            model = self.get("model")
        except KeyError:
            # Numeric model
            n_model = self.get("_model", -1)
            model = self._index_to_model.get(n_model, "unknown")
        return model

    def _params_to_at(self, **atparams) -> dict:
        """Translate Xsuite parameters to AT parameters."""
        for k, v in self.items():
            atk = self._xsuite2at_attr.get(k, None)
            if atk is not None:
                atparams[atk] = v
        atparams["FamName"] = self.name
        if self.origin:
            atparams["origin"] = self.origin
        atparams.update(self._model_to_at())
        return atparams

    def _class_to_at(self, atparams: dict[str, Any]) -> type[elt.Element]:
        """Get the AT element class."""
        return self._atClass

    def _model_to_at(self) -> dict[str, str]:
        """Get the AT passmethod."""
        passmethod = self._xsuite2at_model.get(self._get_model(), None)
        return {} if passmethod is None else {"PassMethod": passmethod}

    def to_dict(self) -> dict[str, Any]:
        """Build a dictionary representation of the XsElement.

        Returns:
            elem_dict:  dictionary representation of the :py:class:`XsElement`
        """
        return self

    def to_at(self) -> elt.Element:
        """Generate the AT element.

        Returns:
            elem:   new :py:class:`.Element` object
        """
        atparams = self._params_to_at()
        cls = self._class_to_at(atparams)
        return cls.from_file(atparams)

    @classmethod
    def from_dict(
        cls, xsparams: dict[str, Any], name: str = "?", origin: str = ""
    ) -> XsElement:
        """Build a XsElement from its dictionary representation.

        Args:
            xsparams:   dictionary representation of the XsElement
            name:       Optional name of the element
            origin:     Xsuite class of the source element

        Returns:
            xselement:  new :py:class:`XsElement` object
        """
        return cls(name, origin=origin, **xsparams)

    @classmethod
    def from_at(cls, FamName: str, match_model: bool = False, **atparams) -> XsElement:
        """Build a XsElement element from an AT element.

        Args:
            FamName:        Name of the element
            match_model:    If :py:obj:`True`, set the Xsuite model and integrator
              matching at best the AT PassMethod. Otherwise, use Xsuite defaults,
            atparams:       dictionary representation of the AT element.

        Returns:
            xselement:  new :py:class:`XsElement` object
        """
        # Conversion of similar attributes
        xsparams = {
            xsk: v
            for k, v in atparams.items()
            if (xsk := cls._at2xsuite_attr.get(k, None)) is not None
        }
        # Handling of model and integrator
        xs_model = cls._at2xsuite_model.get(atparams["PassMethod"], None)
        if match_model:
            if xs_model is None:
                xs_model = cls._at2xsuite_model.get(None, None)
            if (integrator := cls._at_integrator) is not None:
                xsparams["integrator"] = integrator
        else:
            xsparams.pop("num_multipole_kicks", None)
        if xs_model is not None:
            xsparams["model"] = xs_model
        # Set the Xsuite class
        xsparams["__class__"] = cls.__name__
        return cls(FamName, **xsparams)

    @staticmethod
    def static_from_dict(elem_dict: dict[str, Any], name: str = "?") -> XsElement:
        """Build a XsElement from its dictionary representation.

        Args:
            elem_dict:  dictionary representation of the XsElement
            name:       Optional name of the element

        Returns:
            xselement:  new :py:class:`XsElement` object
        """
        class_name = elem_dict.get("__class__", "Unknown")
        xsclass: type[XsElement] = _xsclass.get(class_name, Unknown)
        return xsclass.from_dict(elem_dict, name=name, origin=class_name)

    @staticmethod
    def static_from_at(atelem: elt.Element, match_model: bool = False) -> XsElement:
        """Build a XsElement from an AT element.

        Args:
            atelem:         AT :py:class:`.Element`
            match_model:    If :py:obj:`True`, set the Xsuite model and integrator
              matching at best the AT PassMethod. By default, the Xsuite default model
              and integrator will be used.

        Returns:
            xselement:  new :py:class:`XsElement` object
        """
        xsclass = _at2xsclass.get(type(atelem), XsElement)
        return xsclass.from_at(match_model=match_model, **atelem.to_file())


class Marker(XsElement):
    """Xsuite Marker element."""

    # Class variables
    _atClass: ClassVar[type[elt.Element]] = elt.Marker


class Drift(XsElement):
    """Xsuite Drift element."""

    # Class variables
    _atClass = elt.Drift
    _at2xsuite_model = {"ExactDriftPass": "exact", None: "expanded"}
    _xsuite2at_model = {v: k for k, v in _at2xsuite_model.items()}
    _index_to_model = {0: "adaptive", 1: "expanded", 2: "exact"}


class Multipole(XsElement):
    """Base class for Xsuite magnetic elements."""

    # Class variables
    _atClass = elt.Multipole
    _at2xsuite_model = {
        "ExactMultipoleRadPass": "drift-kick-drift-exact",
        "ExactMultipolePass": "drift-kick-drift-exact",
        None: "drift-kick-drift-expanded",
    }
    _xsuite2at_model = {v: k for k, v in _at2xsuite_model.items()}
    _at_integrator = "yoshida4"
    _xsuite2at_attr = XsElement._xsuite2at_attr | {
        "order": "MaxOrder",
        "num_multipole_kicks": "NumIntSteps",
    }
    _at2xsuite_attr = {v: k for k, v in _xsuite2at_attr.items()}

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
                msg = (
                    f"{self.name}:alignment errors available only for H bends, ignored"
                )
                warnings.warn(AtWarning(msg), stacklevel=2)
            transforms = {"tilt_frame": rot_s_frame}
        elif ismisaligned:
            anchor = self.pop("rot_shift_anchor", 0.0)
            if anchor == 0.0:
                reference = ReferencePoint.ENTRANCE
            elif anchor == 0.5:
                reference = ReferencePoint.CENTRE
            else:
                msg = (
                    "Anchor point for rotation different from 0 or 1, setting"
                    "to AT default: ReferencePoint.CENTRE"
                )
                warnings.warn(AtWarning(msg), stacklevel=2)
                reference = ReferencePoint.CENTRE
            transforms["Reference"] = reference
        return transforms

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
            return porder, p_to_at(list(poly))

        length = self.get("length", 0.0)
        xsorder = self["order"]
        atorder = getattr(self._atClass, "DefaultOrder", xsorder)
        aorder, polya = xspoly(["k0s", "k1s", "k2s", "k3s"], "ksl")
        border, polyb = xspoly(["k0", "k1", "k2", "k3"], "knl")
        aborder = max(aorder, border)
        return {
            "MaxOrder": aborder,
            "PolynomB": polyb[: aborder + 1],
            "PolynomA": polya[: aborder + 1],
        }

    def _set_at_fringe(self):
        """generate the AT fringe field description."""
        atparams = {}
        if self.get("edge_entry_active", 0) == 1:
            atparams["FringeQuadEntrance"] = 1
        if self.get("edge_exit_active", 0) == 1:
            atparams["FringeQuadExit"] = 1
        return atparams

    def _set_xs_transforms(self, atparams: dict) -> None:
        """Generate Xsuite element displacements."""
        tr = {
            "shift_x": atparams.pop("dx", 0.0),
            "shift_y": atparams.pop("dy", 0.0),
            "shift_s": atparams.pop("dz", 0.0),
            "rot_s_rad_no_frame": atparams.pop("tilt", 0.0),
            "rot_x_rad": -atparams.pop("pitch", 0.0),
            "rot_y_rad": atparams.pop("yaw", 0.0),
        }
        misalign = {k: v for k, v in tr.items() if v != 0.0}
        refpoint = atparams.pop("reference", ReferencePoint.ENTRANCE)
        rots_frame = atparams.pop("tilt_frame", 0.0)
        ismisaligned = len(misalign) > 0

        if rots_frame != 0.0:
            if ismisaligned:
                msg = (
                    f"{self.name}: alignment errors available only for H bends, ignored"
                )
                warnings.warn(AtWarning(msg), stacklevel=2)
            misalign = {"rot_s_rad": rots_frame}
        elif ismisaligned:
            if refpoint == ReferencePoint.ENTRANCE:
                misalign["rot_shift_anchor"] = 0.0
            else:
                misalign["rot_shift_anchor"] = atparams["Length"] / 2.0

        self.update(misalign)

    def _set_xs_poly(self, atparams: dict) -> None:
        """Generate the AT field expansion."""
        pola = self.p_to_xsuite(atparams.get("PolynomA", np.zeros(4)))
        polb = self.p_to_xsuite(atparams.get("PolynomB", np.zeros(4)))
        length = atparams.get("Length")
        korder = getattr(self, "_mag_order", None)
        if korder is not None:
            self["k" + str(korder)] = polb[korder]
            self["k" + str(korder) + "s"] = pola[korder]
            pola[korder] = 0.0
            polb[korder] = 0.0
        if length > 0.0:
            polb *= length
            pola *= length
        if np.any(pola) or np.any(polb):
            self["knl"] = list(polb)
            self["ksl"] = list(pola)
        if korder == 1:
            self["edge_entry_active"] = atparams.get("FringeQuadEntrance", 0)
            self["edge_exit_active"] = atparams.get("FringeQuadExit", 0)
        elif atparams.get("BendingAngle", 0.0) == 0.0:
            self["edge_exit_active"] = 0
            self["edge_entry_active"] = 0
        self["_isthick"] = length != 0.0

    def _set_xs_fringe(self, atparams: dict):
        """generate the Xsuite fringe field description."""

    @staticmethod
    def p_to_xsuite(a: Iterable[float]) -> np.ndarray:
        """Convert polynomials from XSUITE to AT."""

        def to_xsuite(x: Iterable[float], factor: float = 1.0) -> Generator[float]:
            f = 1.0
            for n, vx in enumerate(x):
                yield factor * float(vx * f)
                f *= n + 1

        return np.fromiter(to_xsuite(a), dtype=float)

    def _class_to_at(self, atparams: dict[str, Any]) -> type[elt.Element]:
        if atparams.get("Length", 0.0) == 0.0:
            return elt.ThinMultipole
        else:
            return super()._class_to_at(atparams)

    def _params_to_at(self, **atparams) -> dict[str, Any]:
        atparams = super()._params_to_at(**atparams)
        atparams.update(self._set_at_poly())
        atparams.update(self._set_at_fringe())
        atparams.update(self._set_at_transforms())
        return atparams

    @classmethod
    def from_at(cls, **atparams):
        elem = super().from_at(**atparams)
        elem._set_xs_poly(atparams)
        elem._set_xs_fringe(atparams)
        elem._set_xs_transforms(atparams)
        return elem


class Quadrupole(Multipole):
    """Xsuite Quadrupole element."""

    # Class variables
    _atClass = elt.Quadrupole
    _mag_order: ClassVar[int] = 1

    def _set_xs_fringe(self, atparams: dict):
        self["edge_entry_active"] = atparams.get("FringeQuadEntrance", 0)
        self["edge_exit_active"] = atparams.get("FringeQuadExit", 0)


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
        None: "drift-kick-drift-expanded",
    }
    _xsuite2at_model = {v: k for k, v in _at2xsuite_model.items()}
    _xsuite2at_attr = Multipole._xsuite2at_attr | {
        "angle": "BendingAngle",
        "edge_entry_fint": "FringeInt1",
        "edge_exit_fint": "FringeInt2",
        "edge_entry_angle": "EntranceAngle",
        "edge_exit_angle": "ExitAngle",
    }
    _at2xsuite_attr = {v: k for k, v in _xsuite2at_attr.items()}

    def _set_at_fringe(self) -> dict[str, Any]:

        def fringe(xs_inout, at_inout):
            prms = {}
            active = self.get("".join(("edge_", xs_inout, "_active")), 1)
            try:
                fringe_model = self.get("".join(("edge_", xs_inout, "_model")))
                fringe_model = _EDGE_MODEL_TO_INDEX[fringe_model]
            except KeyError:
                fringe_model = self.get("".join(("_edge_", xs_inout, "_model")), 0)
            if fringe_model == -1 or active == 0:
                prms["".join(("FringeBend", at_inout))] = 0
            elif fringe_model == 1:
                prms["".join(("FringeQuad", at_inout))] = 1
            return prms

        atparams = {}
        entry_hgap = self.get("edge_entry_hgap")
        exit_hgap = self.get("edge_exit_hgap")
        if entry_hgap != exit_hgap:
            msg = "Entry and Exit gaps for dipole are different, use entry"
            warnings.warn(AtWarning(msg), stacklevel=2)
        if entry_hgap is not None:
            atparams["FullGap"] = entry_hgap

        atparams.update(fringe("entry", "Entrance"))
        atparams.update(fringe("exit", "Exit"))
        return atparams

    def _set_xs_fringe(self, atparams: dict):

        def fringe(atbend: str, atquad: str) -> str:
            if atparams.get(atbend, 1) == 0:
                return "suppressed"
            elif atparams.get(atquad, 0) == 0:
                return "linear"
            else:
                return "full"

        if (gap := atparams.get("FullGap")) is not None:
            self["edge_entry_gap"] = gap
            self["edge_exit_gap"] = gap
        self["edge_entry_model"] = fringe("FringeBendEntrance", "FringeQuadEntrance")
        self["edge_exit_model"] = fringe("FringeBendExit", "FringeQuadExit")

    def _params_to_at(self, **atparams) -> dict[str, Any]:
        atparams = super()._params_to_at(EntranceAngle=0.0, ExitAngle=0.0, **atparams)
        return atparams

    @classmethod
    def from_at(cls, **atparams):
        elem = super().from_at(**atparams)
        elem["k0_from_h"] = True
        return elem


class RBend(Bend):
    """Xsuite RBend element."""

    _at2xsuite_model = {
        "ExactRectangularBendRadPass": "drift-kick-drift-exact",
        "ExactRectangularBendPass": "drift-kick-drift-exact",
    }
    _xsuite2at_model = {v: k for k, v in _at2xsuite_model.items()}

    def _model_to_at(self) -> dict[str, str]:
        """Get the AT passmethod."""
        straight = self.get("rbend_model", "adaptive") == "straight-body"
        if straight:
            passmethod = self._xsuite2at_model.get(self._get_model(), None)
        else:
            passmethod = Bend._xsuite2at_model.get(self._get_model(), None)
        return {} if passmethod is None else {"PassMethod": passmethod}

    def _params_to_at(self, **atparams) -> dict[str, Any]:
        atparams = super()._params_to_at(**atparams)
        straight = self.get("rbend_model", "adaptive") == "straight-body"
        exact = atparams.get("PassMethod", "").startswith("Exact")
        if (
            straight
            and not exact
            and not (
                np.all(atparams["PolynomA"] == 0.0)
                and np.all(atparams["PolynomB"] == 0.0)
            )
        ):
            msg = (
                "Multipole are present in the rectangular bend: "
                "PassMethod is forced to ExactRectangularBendPass"
            )
            warnings.warn(AtWarning(msg), stacklevel=3)
            atparams["PassMethod"] = "ExactRectangularBendPass"
        hangle = abs(0.5 * self.get("angle"))
        atparams["Length"] = self.get("length_straight", 0.0) / np.sinc(hangle / np.pi)
        atparams["EntranceAngle"] = self.get("edge_entry_angle", 0.0) + hangle
        atparams["ExitAngle"] = self.get("edge_exit_angle", 0.0) + hangle
        return atparams

    @classmethod
    def from_at(cls, **atparams):
        elem = super().from_at(**atparams)
        elem["rbend_model"] = "straight-body"
        hangle = 0.5 * elem["angle"]
        elem["edge_entry_angle"] -= hangle
        elem["edge_exit_angle"] -= hangle
        elem["length_straight"] = elem.pop("length") * np.sinc(hangle / np.pi)
        return elem


class Cavity(XsElement):
    # Class variables
    _atClass = elt.RFCavity
    _at_integrator = "yoshida4"
    _xsuite2at_attr = XsElement._xsuite2at_attr | {
        "voltage": "Voltage",
        "frequency": "Frequency",
        "harmonic": "HarmNumber",
    }
    _at2xsuite_attr = {v: k for k, v in _xsuite2at_attr.items()}

    def _set_xs_lag(self, atparams: dict):
        om = 2.0 * np.pi * atparams.get("Frequency")
        tl = atparams.get("TimeLag", 0.0)
        pl = tl * om / cst.clight * 360.0 + 180.0
        self["lag"] = pl

    def _params_to_at(self, **atparams) -> dict[str, Any]:
        atparams = super()._params_to_at(
            Voltage=0.0, Frequency=np.nan, HarmNumber=0, Energy=0.0, **atparams
        )
        length = atparams.get("Length", 0.0)
        atparams["PassMethod"] = "IdentityPass" if length == 0.0 else "DriftPass"
        if (phaselag := self.get("lag")) is not None:
            atparams["_phaselag"] = phaselag
        return atparams

    @classmethod
    def from_at(cls, HarmNumber=0, **atparams):
        # Discard HarmNumber
        elem = super().from_at(**atparams)
        elem._set_xs_lag(atparams)
        return elem


class Unknown:
    """Class for Xsuite elements without AT equivalent."""
    @classmethod
    def from_dict(cls, xsparams: dict[str, Any], name: str = "?"):
        origin = xsparams["__class__"]
        newcls = Marker if xsparams.get("length", 0.0) == 0.0 else Drift
        msg = f"Element {name!r}: unknown class {origin} replaced by {newcls.__name__}"
        warnings.warn(AtWarning(msg), stacklevel=2)
        return newcls.from_dict(xsparams, name=name, origin=origin)


class Dipole:
    """Class for handling AT dipoles."""
    @classmethod
    def from_at(cls, **atparams) -> XsElement:
        if atparams["PassMethod"] in {
            "ExactRectangularBendPass",
            "ExactRectangularBendRadPass",
        }:
            return RBend.from_at(**atparams)
        else:
            return Bend.from_at(**atparams)


_xsclass = {
    "Marker": Marker,
    "Drift": Drift,
    "Cavity": Cavity,
    "Quadrupole": Quadrupole,
    "Sextupole": Sextupole,
    "Octupole": Octupole,
    "Multipole": Multipole,
    "Bend": Bend,
    "RBend": RBend,
    "Unknown": Unknown,
}


_at2xsclass: dict[type[elt.Element], type[XsElement]] = {
    elt.Marker: Marker,
    elt.Monitor: Marker,
    elt.Drift: Drift,
    elt.RFCavity: Cavity,
    elt.Quadrupole: Quadrupole,
    elt.Sextupole: Sextupole,
    elt.Octupole: Octupole,
    elt.Multipole: Multipole,
    elt.Dipole: Dipole,
}


class XsLine:
    """Line object simulation a Xsuite Line."""

    # Class attributes
    _line_keys: ClassVar[set[str]] = {
        "name",
        "element_names",
        "particle_ref",
        "in_file",
    }

    # Instance attributes
    name: str = "test"
    element_names: list[str] = []
    particle_ref: dict[str, Any] = _default_particle
    in_file: str = ""

    def __init__(self, elements, **root):
        self.elements = {
            k: XsElement.static_from_dict(el, name=k) for k, el in elements.items()
        }
        for k in self._line_keys:
            if (value := root.get(k)) is not None:
                setattr(self, k, value)

    def to_dict(self) -> dict[str, Any]:
        """Create a Xsuite-compatible dictionary representation of the XsLine.

        Returns:
            dict[str, Any]: Dictionary representation of the XsLine.
        """
        return dict(
            ((k, v) for k in self._line_keys if (v := getattr(self, k)) is not None),
            elements={nm: self.elements[nm].to_dict() for nm in self.element_names},
            __class__="Line",
        )

    def to_at(self, **kwargs) -> Lattice:
        """Create an AT lattice from the XsLine.

        Keyword Args:
            name (str):         Name of the lattice. Default: Xsuite line name.
            particle(Particle): Circulating particle. Default: Xsuite reference
              particle
            energy (float):     Energy of the lattice [eV]. Default: from Xsuite
              reference particle
            periodicity(int):   Number of periods. Default: 1
            *:                  All other keywords will be set as Lattice attributes

        Returns:
            lattice:            New :py:class:`.Lattice`.
        """

        def json_generator(params: dict[str, Any]):

            element_names = self.element_names
            particle = self.particle_ref
            mass0 = float(particle["mass0"])
            gamma0 = float(particle["gamma0"][0])
            beta0 = float(particle["beta0"][0])
            energy = gamma0 * mass0

            beam_params = {
                "particle": Particle("", rest_energy=mass0, charge=particle["q0"]),
                "energy": energy,  # [eV]
                "periodicity": 1,
                "name": self.name,
            }

            for k, v in beam_params.items():
                params.setdefault(k, v)

            cavities = []
            cell_length = 0.0

            for nm in element_names:
                elem = self.elements[nm].to_at()
                if isinstance(elem, elt.RFCavity):
                    cavities.append(elem)
                cell_length += getattr(elem, "Length", 0.0)
                yield elem

            rev = beta0 * cst.clight / cell_length

            # Set the frequency of cavities in which it is not specified
            for cav in cavities:
                if np.isnan(cav.Frequency):
                    cav.Frequency = rev * cav.HarmNumber
                elif cav.HarmNumber == 0:
                    cav.HarmNumber = round(cav.Frequency / rev)
                lag = getattr(cav, "_phaselag", 0.0)
                with contextlib.suppress(AttributeError):
                    delattr(cav, "_phaselag")
                om = 2.0 * pi * cav.Frequency
                pl = lag - 180.0
                tl = pl / om * cst.clight / 360.0
                cav.Timelag = tl

        if in_file := getattr(self, "in_file", ""):
            kwargs["in_file"] = in_file
        return Lattice(iterator=json_generator, **kwargs)

    # "lattice" is an alias for "to_at"
    lattice = to_at

    def to_json(
        self, filename: str | Path | None = None, compact: bool = False
    ) -> None:
        """Save the XsLine into a JSON file.

        Parameters:
            filename:   Name of the JSON file. Default: outputs on
              :py:obj:`sys.stdout`
            compact:    If :py:obj:`False` (default), the JSON file is pretty-printed
              with line feeds and indentation. Otherwise, the output is a single line.
        """
        indent = None if compact else 2
        root = self.to_dict()
        if filename is None:
            print(json.dumps(root, cls=_AtEncoder, indent=indent))
        else:
            filename = Path(filename)
            with filename.open("w") as jsonfile:
                json.dump(root, jsonfile, cls=_AtEncoder, indent=indent)

    @classmethod
    def from_dict(
        cls, root: dict[str, Any], use: str | None = None, **kwargs
    ) -> XsLine:
        """Create an XsLine from a dictionary.

        Args:
            root:       Dictionary representation of the line or environment,
            use:        For environment, name of the sequence or line containing the
              desired lattice. Default: ``ring``.

        Returns:
            xsline:     New :py:class:`XsLine`
        """
        defclass = "Environment" if "lines" in root else "Line"
        if root.get("__class__", defclass) == "Line":
            return cls(**root)
        else:
            if use is None:
                lkeys = list(root["lines"].keys())
                use = lkeys[0] if len(lkeys) == 1 else "ring"
            line = root["lines"][use]
            return cls(name=use, elements=root["elements"], **line, **kwargs)

    @classmethod
    def from_json(
        cls, filename: str | Path, use: str | None = None, **kwargs
    ) -> XsLine:
        """Create a XsLine from a JSON file.

        Args:
            filename:   Name or path of the JSON file.
            use:        For environment, name of the sequence or line containing the
              desired lattice. Default: ``ring``

        Returns:
            xsline:     New :py:class:`XsLine`
        """
        filename = Path(filename)
        with filename.open() as f:
            return cls.from_dict(
                json.load(f), in_file=str(filename.resolve()), use=use, **kwargs
            )

    @classmethod
    def from_xsuite(cls, line, use: str | None = None) -> XsLine:
        """Create a XsLine from a Xsuite :py:class:`.Line` or :py:class:`.Environment`.

        Args:
            line:       Xsuite Line or Environment
            use:        For environment, name of the sequence or line containing the
              desired lattice. Default: ``ring``

        Returns:
            xsline:     New :py:class:`XsLine`
        """
        return cls.from_dict(line.to_dict(), use=use)

    @classmethod
    def from_at(cls, ring: Lattice, match_model: bool = False, **kwargs) -> XsLine:
        """Create a XsLine from an AT :py:class:`.Lattice`.

        Args:
            ring:           AT lattice
            match_model:    If :py:obj:`True`, set the Xsuite model and integrator
              matching at best the AT PassMethod. By default, the Xsuite default model
              and integrator will be used.

        Returns:
            xsline:     New :py:class:`XsLine`
        """

        def refpart(rng):
            prt = rng.particle
            if prt.name == "relativistic":
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=UserWarning)
                    prt = Particle("electron")
            gamma0 = rng.energy / prt.rest_energy
            beta0 = sqrt(1.0 - 1.0 / gamma0 / gamma0)
            return {
                "mass0": prt.rest_energy,
                "q0": prt.charge,
                "gamma0": [gamma0],
                "beta0": [beta0],
            }

        def name_gen(name: str, nmax: int = 10000) -> Generator[str, None, None]:
            """Generate a list of tentative names."""
            yield name
            for i in range(1, nmax):
                yield ".".join((name, str(i)))

        def check(name: str, elem: XsElement) -> tuple[bool, bool]:
            """Check if the element is already in the elements dictionary."""
            if name in elements:
                return elem == elements[name], False
            return True, True

        def store_elem(elem: elt.Element):
            """Store the element in the dictionary with a unique name."""
            xselem = XsElement.static_from_at(elem, match_model=match_model)
            for nm in name_gen(elem.FamName):
                valid, new = check(nm, xselem)
                if valid:
                    if new:
                        elements[nm] = xselem
                    return nm
            msg = f"Cannot store {elem.FamName}"
            raise NameError(msg)

        elements: dict[str, XsElement] = {}
        line = [store_elem(el) for el in ring]
        return cls(
            elements=elements,
            element_names=line,
            particle_ref=refpart(ring),
            name=ring.name,
            **kwargs,
        )


def load_xsuite(filename: str | Path, use: str | None = None, **kwargs) -> Lattice:
    """Create a :py:class:`.Lattice`  from a Xsuite JSON file.

    Parameters:
        filename:           Name or path of the JSON file
        use:                Name of the Xsuite line containing the desired lattice.
          Default: ``ring``

    Keyword Args:
        name (str):         Name of the lattice. Default: Xsuite line name.
        particle(Particle): Circulating particle. Default: Xsuite reference
          particle
        energy (float):     Energy of the lattice [eV]. Default: from Xsuite
          reference particle
        periodicity(int):   Number of periods. Default: 1
        *:                  All other keywords will be set as Lattice attributes

    Returns:
        lattice (Lattice):          New :py:class:`.Lattice` object
    """
    return XsLine.from_json(filename, use=use).to_at(**kwargs)


def save_xsuite(
    lattice: Lattice,
    filename: str | Path | None = None,
    match_model: bool = False,
    compact: bool = False,
) -> None:
    """Save a :py:class:`.Lattice` as a Xsuite JSON file.

    Parameters:
        lattice:        Lattice description
        filename:       Name of the JSON file. Default: outputs on
          :py:obj:`sys.stdout`
        match_model:    If :py:obj:`True`, set the Xsuite model matching at best
          the AT PassMethod. By default, the Xsuite default model will be used.
        compact:        If :py:obj:`False` (default), the JSON file is pretty-printed
          with line feeds and indentation. Otherwise, the output is a single line.
    """
    XsLine.from_at(lattice, match_model=match_model).to_json(filename, compact=compact)


def lattice_from_line(line: _XtLine, use: str | None = None, **kwargs) -> Lattice:
    """Create a :py:class:`.Lattice` object from a Xsuite :py:class:`.Line` or
    :py:class:`.Environment`.

    Args:
        line: Xsuite line or environment
        use (str):  Name of the line containing the desired lattice.
          Default: ``ring``.

    Keyword Args:
        name (str):         Name of the lattice. Default: Xsuite line name.
        particle(Particle): Circulating particle. Default: Xsuite reference
          particle
        energy (float):     Energy of the lattice [eV]. Default: from Xsuite
          reference particle
        periodicity(int):   Number of periods. Default: 1
        *:                  All other keywords will be set as Lattice attributes

    Returns:
        lattice:    New :py:class:`.Lattice` object.
    """
    return XsLine.from_xsuite(line, use=use).to_at(**kwargs)
