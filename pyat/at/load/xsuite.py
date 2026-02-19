from __future__ import annotations

__all__ = ["XsEnvironment", "load_xsuite", "save_xsuite"]

import json
import warnings
from collections import ChainMap
from collections.abc import Generator, Iterable
from math import pi, sqrt
from pathlib import Path
from typing import Any, ClassVar
import contextlib

import numpy as np

from ..lattice import AtWarning, ReferencePoint, constants as cst
from ..lattice import Lattice, Particle, elements as elt
from .madx import p_to_at

# Default particle: electrons at 1 GeV
_default_energy = 1.0e9
_default_mass = cst.e_mass
_default_charge = -1
_default_gamma = _default_energy / _default_mass
_default_beta = sqrt(1.0 - 1.0 / _default_gamma / _default_gamma)

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


class BasicElement(dict):
    """Base class for Xsuite elements."""

    # Class attributes
    _atClass: ClassVar[type[elt.Element]] = elt.Marker
    _xsuite2at_attr: ClassVar[dict[str, str]] = {"length": "Length"}
    _at2xsuite_attr: ClassVar[dict[str, str]] = {
        v: k for k, v in _xsuite2at_attr.items()
    }
    _at2xsuite_integrator: ClassVar[dict[str, str]] = {}
    _index_to_model: ClassVar[dict[int, str]] = {
        0: "adaptive",
        1: "full",
        2: "bend-kick-bend",
        3: "rot-kick-rot",
        4: "mat-kick-mat",
        5: "drift-kick-drift-exact",
        6: "drift-kick-drift-expanded",
    }

    def __init__(self, name: str, *args, origin: str = "Unknown", **kwargs):
        super().__init__(*args, **kwargs)
        self.setdefault("length", 0.0)
        self.name = name
        self.origin = origin
        self._xsuite2at_integrator = {
            v: k for k, v in self._at2xsuite_integrator.items()
        }

    def _params_to_at(self, **atparams) -> dict:
        """Translate Xsuite parameters to AT parameters."""
        atparams["FamName"] = self.name
        atparams["origin"] = self.origin
        for k, v in self.items():
            atk = self._xsuite2at_attr.get(k, None)
            if atk is not None:
                atparams[atk] = v
        atparams.update(self._integrator_to_at())
        return atparams

    def _class_to_at(self) -> type[elt.Element]:
        """Get the AT element class."""
        return self._atClass

    def _integrator_to_at(self) -> dict[str, Any]:
        """Get the AT passmethod."""
        try:
            # String model
            integrator = self.get("model")
        except KeyError:
            # Numeric model
            n_integr = self.get("_model", -1)
            integrator = self._index_to_model.get(n_integr, "unknown")
        at_integrator = self._xsuite2at_integrator.get(integrator, None)
        if at_integrator is not None:
            return {"PassMethod": at_integrator}
        else:
            return {}

    def to_at(self) -> elt.Element:
        """Generate the AT element."""
        atparams = self._params_to_at()
        cls = self._class_to_at()
        return cls.from_file(atparams)

    @classmethod
    def from_at(cls, FamName, **atparams):
        """Build a Xsuite element from an AT element."""
        xsparams = {}
        for k, v in atparams.items():
            xsk = cls._at2xsuite_attr.get(k, None)
            if xsk is not None:
                xsparams[xsk] = v
        passmethod = atparams.get("PassMethod")
        xs_integrator = cls._at2xsuite_integrator.get(passmethod, None)
        if xs_integrator is not None:
            xsparams["model"] = xs_integrator
        xsparams["__class__"] = cls.__name__

        transforms = atparams.get("transforms")
        if transforms is not None:
            xsparams.update(transforms)

        return cls(FamName, **xsparams)


class Marker(BasicElement):
    """Xsuite Marker element."""

    # Class variables
    _atClass: ClassVar[type[elt.Element]] = elt.Marker


class Drift(BasicElement):
    """Xsuite Drift element."""

    # Class variables
    _atClass: ClassVar[type[elt.Element]] = elt.Drift
    _at2xsuite_integrator = {
        "DriftPass": "expanded",
        "ExactDriftPass": "exact",
    }
    _index_to_model = {0: "adaptive", 1: "expanded", 2: "exact"}


class Multipole(BasicElement):
    """Base class for Xsuite magnetic elements."""

    # Class variables
    _atClass: ClassVar[type[elt.Element]] = elt.Multipole
    _at2xsuite_integrator = {
        "StrMPoleSymplectic4Pass": "mat-kick-mat",
        "ExactMultipolePass": "drift-kick-drift-exact",
        "DriftPass": "expanded",
        "ExactDriftPass": "exact",
    }
    _xsuite2at_attr = BasicElement._xsuite2at_attr
    _xsuite2at_attr.update(
        {
            "order": "MaxOrder",
            "num_multipole_kicks": "NumIntSteps",
        }
    )
    _at2xsuite_attr: ClassVar[dict[str, str]] = {
        v: k for k, v in _xsuite2at_attr.items()
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

        def xspoly(kmain: list[str], kerr: str) -> np.ndarray:
            # Get main components
            pmain = np.array([self.get(k, 0.0) for k in kmain])
            # Get errors
            perr = np.array(self.get(kerr, np.zeros(4)))
            # Pad to the same length
            lpoly = max(4, len(perr))
            perr = np.pad(perr, (0, lpoly - len(perr)))
            pmain = np.pad(pmain, (1, lpoly - 4))
            if length != 0.0:
                perr /= length
            # Convert to Xsuite
            return p_to_at(list(pmain + perr))

        length = self.get("length", 0.0)
        return {
            "PolynomB": xspoly(["k1", "k2", "k3"], "knl"),
            "PolynomA": xspoly(["k1s", "k2s", "k3s"], "ksl"),
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

    def _class_to_at(self) -> type[elt.Element]:
        if self.get("length", 0.0) == 0.0:
            return elt.ThinMultipole
        else:
            return super()._class_to_at()

    def _params_to_at(self, **atparams) -> dict[str, Any]:
        atparams = super()._params_to_at(**atparams)
        atparams.update(self._set_at_poly())
        atparams.update(self._set_at_fringe())
        atparams.update(self._set_at_transforms())
        atparams.update(self._integrator_to_at())
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
    # Class variables
    _atClass = elt.Bend
    _at2xsuite_integrator = {
        "BndMPoleSymplectic4Pass": "rot-kick-rot",
        "ExactSectorBendPass": "bend-kick-bend",
    }
    _xsuite2at_attr = Multipole._xsuite2at_attr
    _xsuite2at_attr.update(
        {
            "angle": "BendingAngle",
            "edge_entry_fint": "FringeInt1",
            "edge_exit_fint": "FringeInt2",
            "edge_entry_angle": "EntranceAngle",
            "edge_exit_angle": "ExitAngle",
        }
    )
    _at2xsuite_attr: ClassVar[dict[str, str]] = {
        v: k for k, v in _xsuite2at_attr.items()
    }

    def _set_at_fringe(self) -> dict[str, Any]:

        def fringe(kxs, kat):
            atparams = {}
            active = self.get("".join(("edge_", kxs, "_active")), 1)
            try:
                fringe_model = self.get("".join(("edge_", kxs, "_model")))
                fringe_model = _EDGE_MODEL_TO_INDEX[fringe_model]
            except KeyError:
                fringe_model = self.get("".join(("_edge_", kxs, "_model")), 0)
            if fringe_model == -1 or active == 0:
                atparams["".join(("FringeBend", kat))] = 0
            elif fringe_model == 1:
                atparams["".join(("FringeQuad", kat))] = 1
            return atparams

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
    def _params_to_at(self, **atparams) -> dict[str, Any]:
        atparams = super()._params_to_at(**atparams)
        hangle = abs(0.5 * self.get("angle"))
        atparams["EntranceAngle"] = self.get("edge_entry_angle", 0.0) + hangle
        atparams["ExitAngle"] = self.get("edge_exit_angle", 0.0) + hangle
        return atparams

    @classmethod
    def from_at(cls, **atparams):
        elem = super().from_at(**atparams)
        hangle = 0.5 * elem["angle"]
        elem["edge_entry_angle"] -= hangle
        elem["edge_exit_angle"] -= hangle
        elem["length_straight"] = elem.pop("length") * np.sinc(hangle / np.pi)
        return elem


class Cavity(BasicElement):
    # Class variables
    _atClass = elt.RFCavity
    _xsuite2at_attr = BasicElement._xsuite2at_attr | {
        "voltage": "Voltage",
        "frequency": "Frequency",
        "harmonic": "HarmNumber",
    }
    _at2xsuite_attr: ClassVar[dict[str, str]] = {
        v: k for k, v in _xsuite2at_attr.items()
    }

    def _set_xs_lag(self, atparams: dict):
        om = 2.0 * np.pi * atparams.get("Frequency")
        tl = atparams.get("TimeLag", 0.0)
        pl = tl * om / cst.clight * 360.0 + 180.0
        self["lag"] = pl

    def _params_to_at(self, **atparams) -> dict[str, Any]:
        atparams = super()._params_to_at(
            Voltage=0.0, Frequency=np.nan, HarmNumber=0, Energy=0.0, **atparams
        )
        if (phaselag := self.get("lag")) is not None:
            atparams["_phaselag"] = phaselag
        return atparams

    @classmethod
    def from_at(cls, HarmNumber=0, **atparams):
        # Discard HarmNumber
        elem = super().from_at(**atparams)
        elem._set_xs_lag(atparams)
        return elem


def Unknown(name: str, *args, origin: str = "Unknown", **kwargs):
    cls = Marker if kwargs.get("length", 0.0) == 0.0 else Drift
    msg = f"Element {name!r}: unknown class {origin} replaced by {cls.__name__}"
    warnings.warn(AtWarning(msg), stacklevel=2)
    return cls(name, *args, origin=origin, **kwargs)


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


class Dipole:
    @classmethod
    def from_at(cls, **atparams):
        if atparams["PassMethod"] in {
            "ExactRectangularBendPass",
            "ExactRectangularBendRadPass",
        }:
            return RBend.from_at(**atparams)
        else:
            return Bend.from_at(**atparams)


_at2xsclass: dict[type[elt.Element], type[BasicElement]] = {
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


class XsEnvironment(ChainMap):
    """Simulate an Xsuite Environment."""

    # Class attributes
    _env_keys: ClassVar[set[str]] = {"lines", "element_names", "particle_ref"}
    _line_keys: ClassVar[set[str]] = {"element_names", "particle_ref"}
    _default_particle: ClassVar[dict[str, float | list[float]]] = {
        "mass0": _default_mass,
        "gamma0": [_default_gamma],
        "beta0": [_default_beta],
        "q0": _default_charge
    }

    def __init__(self, elements=None, vars=None, lines=None, **kwargs):
        """Initialise an Xsuite Environment.

        Args:
            elements:   Dictionary of elements
            vars:       Dictionary of variables
            lines:      Dictionary of lines
            *:          All other arguments are set as attributes.
        """
        super().__init__(vars or {}, elements or {}, lines or {})
        self.vars = self.maps[0]
        self.elements = self.maps[1]
        self.lines = self.maps[2]
        for k, v in kwargs.items():
            setattr(self, k, v)

    @classmethod
    def from_dict(cls, root, **kwargs):
        """Create an Xsuite Environment from a dictionary.

        Args:
            root:       Dictionary representation of the environment.

        """

        def rr(name: str, elem_dict: dict[str, Any]):
            jcls = elem_dict.get("__class__", "Unknown")
            xtelement = _xsclass.get(jcls, Unknown)
            return xtelement(name, origin=jcls, **elem_dict)

        defclass = "Environment" if "lines" in root else "Line"
        if root.pop("__class__", defclass) == "Line":
            newline = {
                k: v for k in cls._line_keys if (v := root.pop(k, None)) is not None
            }
            root["lines"] = {"ring": newline}
        element_refs = root.pop("elements", {})
        elems = {k: rr(k, el) for k, el in element_refs.items()}
        items = {k: v for k in cls._env_keys if (v := root.pop(k, None)) is not None}
        vars = root.pop("_var_management_data", {}).pop("var_values", {})
        return cls(elements=elems, vars=vars, **items, **kwargs)

    @classmethod
    def from_json(cls, filename: str | Path):
        """Create an Xsuite Environment from a JSON file."""
        filename = Path(filename)
        with filename.open() as f:
            return cls.from_dict(json.load(f), in_file=str(filename.resolve()))

    @classmethod
    def from_at(cls, ring: Lattice, **kwargs):
        """Create an Xsuite environment from an AT lattice."""

        def name_gen(name: str, nmax: int = 10000) -> Generator[str, None, None]:
            """Generate a list of tentative names."""
            yield name
            for i in range(1, nmax):
                yield ".".join((name, str(i)))

        def check(name: str, elem: BasicElement) -> tuple[bool, bool]:
            """Check if the element is already in the elements dictionary."""
            if name in elements:
                return elem == elements[name], False
            return True, True

        def store_elem(elem: elt.Element):
            """Store the element in the dictionary with a unique name."""
            xsclass = _at2xsclass.get(type(elem), BasicElement)
            xselem = xsclass.from_at(**elem.to_file())
            for nm in name_gen(elem.FamName):
                valid, new = check(nm, xselem)
                if valid:
                    if new:
                        elements[nm] = xselem
                    return nm
            msg = f"Cannot store {elem.FamName}"
            raise NameError(msg)

        elements: dict[str, BasicElement] = {}
        line = [store_elem(el) for el in ring]
        lines = {ring.name or "ring": {"element_names": line}}
        return cls(elements=elements, vars={}, lines=lines, **kwargs)

    def lattice(self, **kwargs):
        """Create a lattice from the selected line of an Xsuite environment.

        Keyword Args:
            use (str):          Name of the sequence or line containing the desired
              lattice. Default: ``ring``
            name (str):         Name of the lattice. Default: Xsuite line name.
            particle(Particle): Circulating particle. Default: Xsuite reference
              particle
            energy (float):     Energy of the lattice [eV]. Default: from Xsuite
              reference particle
            periodicity(int):   Number of periods. Default: 1
            *:                  All other keywords will be set as Lattice attributes
        """

        def json_generator(params: dict[str, Any]):
            use = params.get("use")
            if use is None:
                lkeys = list(self.lines.keys())
                use = lkeys[0] if lkeys else "ring"
            params.setdefault("name", use)
            params.setdefault("energy", 1.0e9)
            params.setdefault("periodicity", 1)
            try:
                line = self.lines[use]
            except KeyError as exc:
                exc.args = (f"Unknown line: choose in {set(self.lines.keys())}",)
                raise

            element_names = line["element_names"]

            particle = line.get(
                "particle_ref", self.get("particle_ref", self._default_particle)
            )
            mass0 = float(particle["mass0"])
            gamma0 = float(particle["gamma0"][0])
            beta0 = float(particle["beta0"][0])
            energy = gamma0 * mass0

            beam_params = {
                "particle": Particle("", rest_energy=mass0, charge=particle["q0"]),
                "energy": energy,  # [eV]
                "periodicity": 1,
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

        in_file = getattr(self, "in_file", None)
        return Lattice(iterator=json_generator, in_file=in_file, **kwargs)

    def to_json(self, filename: str | Path | None = None, compact: bool = False):
        """Save an Xsuite environment into a JSON file.

        Parameters:
            filename:   File name of file path
            compact:    If :py:obj:`False` (default), the JSON file is pretty-printed
          with line feeds and indentation. Otherwise, the output is a single line.
        """
        indent = None if compact else 2
        root = {
            "elements": self.elements,
            "vars": self.vars,
            "lines": self.lines,
            "__class__": "Environment",
        }
        if filename is None:
            print(json.dumps(root, cls=_AtEncoder, indent=indent))
        else:
            filename = Path(filename)
            with filename.open("w") as jsonfile:
                json.dump(root, jsonfile, cls=_AtEncoder, indent=indent)

    to_at = lattice


def load_xsuite(filename: str | Path, **kwargs) -> Lattice:
    """Create a :py:class:`.Lattice`  from a Xsuite JSON file.

    Parameters:
        filename:           Name of a JSON file

    Keyword Args:
        use:                Name of the Xsuite line containing the desired lattice.
          Default: ``ring``
        *:                  All other keywords update the lattice properties

    Returns:
        lattice (Lattice):          New :py:class:`.Lattice` object
    """
    return XsEnvironment.from_json(filename).to_at(**kwargs)


def save_xsuite(
    lattice: Lattice, filename: str | Path | None, compact: bool = False
) -> None:
    """Save a :py:class:`.Lattice` as a Xsuite JSON file.

    Parameters:
        lattice:        Lattice description
        filename:       Name of the JSON file. Default: outputs on
          :py:obj:`sys.stdout`
        compact:        If :py:obj:`False` (default), the JSON file is pretty-printed
          with line feeds and indentation. Otherwise, the output is a single line.
    """
    XsEnvironment.from_at(lattice).to_json(filename, compact=compact)
