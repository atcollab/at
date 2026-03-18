r"""Using `Xsuite`_ .json files with PyAT.
==========================================

PyAT can read and write lattice descriptions in Xsuite .json files. If both PyAT and
Xsuite are installed, in-memory direct conversion is also possible in both directions.

However, because of intrinsic differences between PyAT and Xsuite, some
incompatibilities must be taken into account.

1. Differences between PyAT and Xsuite
-------------------------------------------

Magnet models
^^^^^^^^^^^^^

Xsuite separates the definition of the magnet modelling method (``model``) and the
definition of the integrator (``integrator``). It offers a large choice for both.
PyAT combines both in the ``PassMethod`` associated with the element. The conversion in
both directions is described below.

Magnet transformations
^^^^^^^^^^^^^^^^^^^^^^

PyAT uses a coherent definition of rotations: the positive sign corresponds to a
cork-screw rotation, or in other words an anti-clockwise rotation when facing the axis.
Xsuite has a special definition for the x-axis rotation, opposite to the standard. As a
consequence, the signs of the ``pitch`` in AT and ``rot-x-rad`` in Xsuite are
opposite. The conversion takes care of that.

Rectangular bending magnets
^^^^^^^^^^^^^^^^^^^^^^^^^^^
A Xsuite :py:class:`.RBend` magnet corresponds to a standard :py:class:`.Dipole` if its
``rbend_model`` is set to ``curved-body``. When converted to AT, its ``PassMethod``
will then be set to ``BndMPoleSymplectic4Pass`` or ``ExactSectorBendPass`` depending
on its ``model``.

It corresponds to a real rectangular magnet if its ``rbend_model`` is set to
``straight-body``. When converted to AT, its ``PassMethod`` will then be set to
``ExactRectangularBendPass``, the only one for rectangular bends.

2. Reading Xsuite files
-----------------------

PyAT can read both :py:class:`.Environment` and :py:class:`.Line` file descriptions.
For :py:class:`.Environment` files, a ``use`` keyword selects the desired
:py:class:`.Line`.

Xsuite elements converted to PyAT have an additional ``origin`` string attribute giving
the class of the original Xsuite element. Xsuite element without equivalent in PyAT are
converted to :py:class:`.Marker` or :py:class:`.Drift` depending on their length.

Model and integrator
^^^^^^^^^^^^^^^^^^^^

Xsuite models equivalent to PyAT "exact" methods are converted to the equivalent
``PassMethod``:

+------------------------------+--------------------------+----------------------------+
| Xsuite element               | Xsuite model             | AT PassMethod              |
+==============================+==========================+============================+
| Drift                        | "exact"                  | "ExactDriftPass"           |
+------------------------------+--------------------------+----------------------------+
| Straight magnets             | "drift-kick-drift-exact" | "ExactMultipolePass"       |
+------------------------------+--------------------------+----------------------------+
| Bend,                        | "bend-kick-bend"         | "ExactSectorBendPass"      |
| Rbend in "curved-body"       |                          |                            |
+------------------------------+--------------------------+----------------------------+
| Rbend in "straight-body"     | "drift-kick-drift-exact" | "ExactRectangularBendPass" |
+------------------------------+--------------------------+----------------------------+

All other models are converted to the default element ``PassMethod``.

Longitudinal motion
^^^^^^^^^^^^^^^^^^^
The generated AT lattice has longitudinal motion turned off: RF cavities are inactive
and synchrotron radiation is turned off.

3. Writing Xsuite files
-----------------------

Rectangular bending magnets
^^^^^^^^^^^^^^^^^^^^^^^^^^^

An AT :py:class:`.Dipole` will be converted to a :py:class:`.RBend` with ``rbend_model``
set to ``straight-body`` if its ``PassMethod`` is ``ExactRectangularBendPass``. In all
other cases, it will be converted to a :py:class:`.Bend`.

Model and integrator
^^^^^^^^^^^^^^^^^^^^
PyAT exact passmethods are converted to the equivalent Xsuite ``model``. Other methods
are converting according to a ``match_model`` keyword:

- ``match_model`` is False (the default): ``model``, ``integrator``,
  ``num_multipole_kicks`` are set to their default value. This way, the default
  behaviour in AT is turned into the default behaviour in Xsuite,
- ``match_model`` is True: a model matching at best the AT model is selected, the
  integrator is set to ``yoshida4`` and ``num_multipole_kicks`` is set equal to
  ``NumIntSteps``.

This is summarised in this table:

.. table::
   :align: center

   +-----------------+------------------------+--------------------------------------------+
   |                 |                        |Xsuite model, match_model=                  |
   |AT element       |AT PassMethod           +----------------+---------------------------+
   |                 |                        | False          | True                      |
   +=================+========================+================+===========================+
   |                 |ExactDriftPass          |"exact"                                     |
   |Drift            +------------------------+----------------+---------------------------+
   |                 |*default*               |"adaptive"      |"expanded"                 |
   +-----------------+------------------------+----------------+---------------------------+
   |                 |ExactMultipolePass      |"drift-kick-drift-exact"                    |
   |Straight magnet  +------------------------+----------------+---------------------------+
   |                 |*default*               |"adaptive"      |"drift-kick-drift-expanded"|
   +-----------------+------------------------+----------------+---------------------------+
   |                 |ExactSectorBendPass     |"bend-kick-bend"|                           |
   |                 +------------------------+----------------+---------------------------+
   |Dipole           |ExactRectangularBendPass|"drift-kick-drift-exact"                    |
   |                 +------------------------+----------------+---------------------------+
   |                 |*default*               |"adaptive"      |"rot-kick-rot"             |
   +-----------------+------------------------+----------------+---------------------------+

Longitudinal motion
^^^^^^^^^^^^^^^^^^^^
Whatever the 6D status of the AT lattice, the Xsuite lattice will have synchrotron
radiation turned off and active RF cavities.

Ignored attributes
^^^^^^^^^^^^^^^^^^

The following AT attributes are ignored in writing Xsuite .json files:

- ``RApertures``
- ``EApertures``
- ``X0ref``. This parameter defines the distance between the magnetic axis of a
  rectangular bend and the reference trajectory. It only matters for rectangular bends
  with focusing or higher order multipoles. Xsuite has a similar attribute
  ``rbend_shift`` but there is no analytical conversion possible between both.
- ``RefDZ``. This defines the correction of path lengthening for rectangular bends.


.. _xsuite: https://xsuite.readthedocs.io
"""  # fmt: skip # noqa: E501

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
from math import sqrt
from pathlib import Path
from collections.abc import Callable
from typing import Any, ClassVar
import contextlib

import numpy as np

from ..lattice import AtWarning, ReferencePoint, constants as cst
from ..lattice import Lattice, Particle, elements as elt
from .madx import poly_from_mad, poly_to_mad

try:
    from xtrack import Line
except ImportError:

    class Line:
        """Dummy class replacing xtrack.Line."""

        # noinspection PyUnusedLocal
        @classmethod
        def from_dict(cls, root: dict[str, Any]) -> Line:
            msg = "xtrack is not installed."
            raise ImportError(msg)

        # noinspection PyMethodMayBeStatic
        def to_dict(self) -> dict[str, Any]:
            msg = "xtrack is not installed."
            raise ImportError(msg)


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
    _at2xsuite_model: ClassVar[dict[str | None, str]] = {}
    _xsuite2at_model: ClassVar[dict[str, str | None]] = {}
    _at_integrator: ClassVar[str | None] = None

    # Instance attributes
    name: str  #: Element name
    origin: str  #: Xsuite class of the source element

    def __init__(self, name: str, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setdefault("length", 0.0)
        self.name = name
        self.origin = kwargs.get("__class__", "")

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
        passmethod = self._xsuite2at_model.get(self.get("model"), None)
        return {} if passmethod is None else {"PassMethod": passmethod}

    def to_dict(self) -> dict[str, Any]:
        """Build a dictionary representation of the XsElement.

        Returns:
            elem_dict:  dictionary representation of the :py:class:`XsElement`
        """
        return self

    def to_at(self) -> list[elt.Element]:
        """Generate the AT element.

        Returns:
            elem:   new :py:class:`.Element` object
        """
        atparams = self._params_to_at()
        cls = self._class_to_at(atparams)
        return [cls.from_file(atparams)]

    # noinspection PyUnusedLocal
    @classmethod
    def from_dict(
        cls,
        xsparams: dict[str, Any],
        name: str = "?",
        warn: Callable[[str], bool] = lambda clname: True,
    ) -> XsElement:
        """Build a XsElement from its dictionary representation.

        Args:
            xsparams:   dictionary representation of the XsElement
            name:       Optional name of the element
            warn:       warning trigger function

        Returns:
            xselement:  new :py:class:`XsElement` object
        """
        return cls(name=name, **xsparams)

    @classmethod
    def from_at(cls, match_model: bool = False, **atparams) -> XsElement:
        """Build a XsElement element from an AT element.

        Args:
            match_model:    If :py:obj:`True`, set the Xsuite model and integrator
              matching at best the AT PassMethod. Otherwise, use Xsuite defaults,
            atparams:       dictionary representation of the AT element.

        Returns:
            xselement:  new :py:class:`XsElement` object
        """
        # Conversion of similar attributes
        xsparams = {
            kxs: v
            for kxs, kat in cls._xsuite2at_attr.items()
            if (v := atparams.get(kat)) is not None
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
        return cls(name=atparams["FamName"], **xsparams)

    @staticmethod
    def static_from_dict(
        elem_dict: dict[str, Any],
        name: str = "?",
        warn: Callable[[str], bool] = lambda clname: True,
    ) -> XsElement:
        """Build a XsElement from its dictionary representation.

        Args:
            elem_dict:  dictionary representation of the XsElement
            name:       Optional name of the element
            warn:       warning trigger function

        Returns:
            xselement:  new :py:class:`XsElement` object
        """
        xsclass = _xsclass.get(elem_dict.get("__class__", "Unknown"), NotInAT)
        return xsclass.from_dict(elem_dict, name=name, warn=warn)

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
        xsclass = _at2xsclass.get(type(atelem), NotInXsuite)
        return xsclass.from_at(match_model=match_model, **atelem.to_file())

    @property
    def _l_anchor(self) -> float:
        return self.get("length", 0.0)


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


class Multipole(XsElement):
    """Base class for Xsuite magnetic elements."""

    # Class variables
    _atClass: ClassVar[type[elt.Element]] = elt.Multipole
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
            length = self._l_anchor
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
                warnings.warn(AtWarning(msg), stacklevel=2)
                reference = ReferencePoint.CENTRE
            transforms["reference"] = reference.value
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
            return porder, np.fromiter(poly_from_mad(poly), dtype=float, count=lpoly)

        length = self.get("length", 0.0)
        xsorder = self.get("order", 0)
        atorder = getattr(self._atClass, "DefaultOrder", xsorder)
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
        return {}

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
        refpoint = atparams.pop("reference", ReferencePoint.CENTRE)
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
        pata = atparams.get("PolynomA", np.zeros(4))
        pola = np.fromiter(poly_to_mad(pata), dtype=float, count=pata.size)
        patb = atparams.get("PolynomB", np.zeros(4))
        polb = np.fromiter(poly_to_mad(patb), dtype=float, count=patb.size)
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
        if (scaling := atparams.get("FieldScaling")) is not None:
            self["delta_taper"] = scaling - 1.0
        self["_isthick"] = length != 0.0

    def _set_xs_fringe(self, atparams: dict) -> None:
        """generate the Xsuite fringe field description."""

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

    def _set_at_fringe(self):
        """generate the AT fringe field description."""
        atparams = {}
        if self.get("edge_entry_active", 0):
            atparams["FringeQuadEntrance"] = 1
        if self.get("edge_exit_active", 0):
            atparams["FringeQuadExit"] = 1
        return atparams

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
    _mag_order = 1
    _edge_to_xs: ClassVar[dict[bool, dict[bool, str]]] = {
        True: {True: "full", False: "dipole-only"},
        False: {True: "linear", False: "linear"},
    }

    def _set_at_fringe(self) -> dict[str, Any]:
        atparams = {}
        entry_hgap = self.get("edge_entry_hgap")
        exit_hgap = self.get("edge_exit_hgap")
        if entry_hgap != exit_hgap:
            msg = "Entry and Exit gaps for dipole are different, use entry"
            warnings.warn(AtWarning(msg), stacklevel=2)
        if entry_hgap is not None:
            atparams["FullGap"] = entry_hgap

        if self.get("edge_entry_model", "linear") in ["linear", "full"]:
            atparams["FringeQuadEntrance"] = 1
        if self.get("edge_exit_model", "linear") in ["linear", "full"]:
            atparams["FringeQuadExit"] = 1
        return atparams

    def _set_xs_fringe(self, atparams: dict):
        if (gap := atparams.get("FullGap")) is not None:
            self["edge_entry_gap"] = gap
            self["edge_exit_gap"] = gap
        exact = atparams.get("PassMethod", "").startswith("Exact")
        qentry = atparams.get("FringeQuadEntrance", 0)
        qexit = atparams.get("FringeQuadExit", 0)
        self["edge_entry_model"] = self._edge_to_xs[exact][qentry > 0]
        self["edge_exit_model"] = self._edge_to_xs[exact][qexit > 0]

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
            passmethod = self._xsuite2at_model.get(self.get("model"), None)
        else:
            passmethod = Bend._xsuite2at_model.get(self.get("model"), None)
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
        hangle = abs(0.5 * self["angle"])
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

    @property
    def _l_anchor(self) -> float:
        return self.get("length_straight", 0.0)


class Cavity(XsElement):
    # Class variables
    _atClass = elt.RFCavity
    _at_integrator = "yoshida4"
    _xsuite2at_attr = XsElement._xsuite2at_attr | {
        "voltage": "Voltage",
        "frequency": "Frequency",
    }

    def _set_xs_lag(self, atparams: dict):
        freq = atparams["Frequency"]
        tl = atparams.get("TimeLag", 0.0)
        pl = 180.0 - 360.0 * tl * freq / cst.clight
        self["lag"] = pl

    def _params_to_at(self, **atparams) -> dict[str, Any]:
        atparams = super()._params_to_at(
            Voltage=0.0, Frequency=np.nan, Energy=0.0, **atparams
        )
        length = atparams.get("Length", 0.0)
        atparams["PassMethod"] = "IdentityPass" if length == 0.0 else "DriftPass"
        atparams["HarmNumber"] = self.get("harmonic", 0)
        if (phaselag := self.get("lag")) is not None:
            atparams["_phaselag"] = phaselag
        return atparams

    @classmethod
    def from_at(cls, **atparams):
        elem = super().from_at(**atparams)
        elem._set_xs_lag(atparams)
        return elem


class NotInAT:
    """Class for Xsuite elements without AT equivalent."""

    @classmethod
    def from_dict(
        cls,
        xsparams: dict[str, Any],
        name: str = "?",
        warn: Callable[[str], bool] = lambda clname: True,
    ) -> XsElement:
        origin = xsparams.get("__class__", "Unknown")
        xsclass = Marker if xsparams.get("length", 0.0) == 0.0 else Drift
        if warn(origin):
            repl = xsclass.__name__
            msg = f"Element {name!r}: unknown class {origin} replaced by {repl}"
            warnings.warn(AtWarning(msg), stacklevel=3)
        return xsclass.from_dict(xsparams, name=name, warn=warn)


class NotInXsuite:
    """Class for AT elements without Xsuite equivalent."""

    @classmethod
    def from_at(cls, **atparams) -> XsElement:
        name = atparams["FamName"]
        xsclass = Marker if atparams["Length"] == 0.0 else Drift
        msg = f"Element {name!r} is replaced by a {xsclass.__name__}"
        warnings.warn(AtWarning(msg), stacklevel=2)
        return xsclass.from_at(**atparams)


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
    elt.ThinMultipole: Multipole,
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

    def __init__(self, elements: dict[str, dict], **root):

        def warn(clname):
            # Warn only for the 1st element of each class
            nbu = unknown.setdefault(clname, 0)
            unknown[clname] += 1
            return nbu == 0

        unknown: dict[str, int] = {}  # counter of unknown classes

        self.elements = {
            k: XsElement.static_from_dict(el, name=k, warn=warn)
            for k, el in elements.items()
        }

        for k in self._line_keys:
            if (value := root.get(k)) is not None:
                setattr(self, k, value)

        if unknown:
            print("Unknown classes:")
            for cl, nb in unknown.items():
                print(f"    {cl}: {nb}")

    def to_dict(self) -> dict[str, Any]:
        """Create a Xsuite-compatible dictionary representation of the XsLine.

        Returns:
            dict[str, Any]: Dictionary representation of the XsLine.
        """
        return dict(
            (
                (k, v)
                for k in ["element_names", "particle_ref"]
                if (v := getattr(self, k)) is not None
            ),
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
                try:
                    for elem in self.elements[nm].to_at():
                        if isinstance(elem, elt.RFCavity):
                            cavities.append(elem)
                        cell_length += getattr(elem, "Length", 0.0)
                        yield elem
                except Exception as exc:  # noqa: PERF203
                    exc.args = (f"In element {nm!r}: {exc.args[0]}",)
                    raise

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
                tl = (180.0 - lag) * cst.clight / 360.0 / cav.Frequency
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

        def unique_names(names):
            """Generate unique names."""
            _, idx_inv, cnt = np.unique(names, return_inverse=True, return_counts=True)
            idx_list = np.split(np.argsort(idx_inv), np.cumsum(cnt[:-1]))
            for ia in idx_list:
                if len(ia) > 1:
                    for i, ii in enumerate(np.sort(ia)):
                        names[ii] = ".".join((names[ii], str(i + 1)))
            return names

        element_names = unique_names([el.FamName for el in ring])
        return cls(
            elements={
                nm: XsElement.static_from_at(el, match_model=match_model)
                for nm, el in zip(element_names, ring, strict=True)
            },
            element_names=element_names,
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

    See Also:
        :py:func:`.load_json`, which reads both AT and Xsuite .json files,

        :py:func:`.load_lattice` for a generic lattice-loading function,

        :py:meth:`.Lattice.load` for a generic lattice factory method.
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


def lattice_from_line(line: Line, use: str | None = None, **kwargs) -> Lattice:
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
