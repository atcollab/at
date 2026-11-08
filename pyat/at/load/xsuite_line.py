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

- ``match_model`` is False: ``model`` and ``integrator`` are set to their
  default value. This way, the default behaviour in AT is turned into the
  default behaviour in Xsuite,
- ``match_model`` is True (the default): a model matching at best the AT model is
  selected, the integrator is set to ``yoshida4``.

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
   |                 |*default*               |"adaptive"      |"drift-kick-drift-expanded"|
   +-----------------+------------------------+----------------+---------------------------+

Longitudinal motion
^^^^^^^^^^^^^^^^^^^^
Whatever the 6D status of the AT lattice, the Xsuite lattice will have synchrotron
radiation turned off and active RF cavities.

Ignored attributes
^^^^^^^^^^^^^^^^^^

The following AT attributes are ignored in writing Xsuite .json files:

- ``X0ref``: the conversion to the similar Xsuite attribute ``rbend_shift`` is
  approximate.
- ``RefDZ``: this defines the correction of path lengthening for straight-body bends.
  It has no equivalent in Xsuite.


.. _xsuite: https://xsuite.readthedocs.io
"""  # fmt: skip # noqa: E501

from __future__ import annotations

__all__ = [
    "XsLine",
    "load_xsuite",
    "save_xsuite",
]

import json
from math import sqrt
from pathlib import Path
from typing import Any, ClassVar
from collections.abc import Iterator
import contextlib
import warnings

import numpy as np

from ..lattice import AtWarning, constants as cst
from ..lattice import Lattice, Particle, elements as elt
from .xsuite_elements import XsElement


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

_transform_keys: set[str] = {"dx", "dy", "dz", "tilt"}


class _AtEncoder(json.JSONEncoder):
    """JSON encoder for specific AT types."""

    def default(self, o):
        if isinstance(o, np.ndarray):
            return o.tolist()
        elif isinstance(o, Particle):
            return o.to_dict()
        elif isinstance(o, np.integer):
            return int(o)
        elif isinstance(o, np.floating):
            return float(o)
        else:
            return super().default(o)


class _FilterWarnings(dict):
    """Summarises warnings for unknown classes."""

    def __init__(self, header: str):
        """
        Args:
            header:     Header to print before the summary
        """  # noqa: D415
        super().__init__()
        self.header = header

    def __call__(self, clname: str) -> bool:
        """Check if a warning should be issued for the given class name.

        Args:
            clname:     Class name of the element

        Returns:
            bool:       True if a warning should be issued, False otherwise
        """
        nbu = self.setdefault(clname, 0)
        self[clname] += 1
        return nbu == 0

    def __str__(self) -> str:
        """Return a string representation of the warnings."""
        flist = "\n".join(f"{clname:>20}: {nb}" for clname, nb in self.items())
        return f"{self.header}\n{flist}" if flist else ""

    def summarise(self) -> None:
        """Print a summary of the warnings."""
        if self:
            warnings.warn(AtWarning(str(self)), stacklevel=3)


def _separate_apertures(lattice: Lattice) -> Lattice:
    """Generate a Lattice with separate rectangular and elliptical apertures.

    Args:
        lattice:        AT lattice
    Returns:
        new_lattice:    New AT lattice with separate apertures
    """

    def create_aperture_marker(elem: elt.Element, attr: str) -> elt.Marker | None:
        if (ap := getattr(elem, attr, None)) is not None:
            at = {k: v for k in _transform_keys if (v := getattr(elem, k, 0.0)) != 0.0}
            at[attr] = ap
            return elt.Marker(elem.FamName + "_aperture", **at)
        return None

    def same_transforms(mk1: elt.Marker, mk2: elt.Marker) -> bool:
        at1 = {k: v for k in _transform_keys if (v := getattr(mk1, k, 0.0)) != 0.0}
        at2 = {k: v for k in _transform_keys if (v := getattr(mk2, k, 0.0)) != 0.0}
        return at1 == at2

    def same_apertures(mk1: elt.Marker, mk2: elt.Marker, aname: str) -> bool:
        return np.all(getattr(mk1, aname) == getattr(mk2, aname))

    def compare_markers(mk2: list[elt.Marker | None], aname: str) -> list[elt.Marker]:
        markers = [m for m in mk2 if m is not None]
        if (
            len(markers) == 2
            and same_apertures(markers[0], markers[1], aname)
            and same_transforms(markers[0], markers[1])
        ):
            return [markers[0]]
        else:
            return markers

    def element_generator(ring: Lattice) -> Iterator[elt.Element]:
        """Insert aperture markers in the lattice."""
        r_pending = None
        e_pending = None
        for el in ring:
            r_elem = None
            e_elem = None
            if not isinstance(el, elt.Marker):
                r_elem = create_aperture_marker(el, "RApertures")
                e_elem = create_aperture_marker(el, "EApertures")
            yield from compare_markers([r_pending, r_elem], "RApertures")
            yield from compare_markers([e_pending, e_elem], "EApertures")
            r_pending = None
            e_pending = None

            yield el

            if el.Length > 0.0:
                r_pending = r_elem
                e_pending = e_elem
        yield from compare_markers([r_pending], "RApertures")
        yield from compare_markers([e_pending], "EApertures")

    return Lattice(
        element_generator(lattice),
        name=lattice.name,
        particle=lattice.particle,
        energy=lattice.energy,
        periodicity=lattice.periodicity,
    )


class XsLine:
    # Class attributes
    _attr_names: ClassVar = ["element_names", "particle_ref", "name", "in_file"]

    # Instance attributes
    name: str = "line"
    element_names: list[str] = []
    particle_ref: dict[str, Any] = _default_particle
    in_file: str

    def __init__(self, xselems: dict[str, XsElement], /, **kwargs):
        """line object simulating a Xsuite Line.

        Args:
            xselems:    XsElements from which the line is built

        Keyword Args:
            element_names(Sequence[str]):  Name of elements constituting the line
            particle_ref:   Reference particle
            name(str):      Name of the line
            in_file(str):   Name of the source file
        """

        self.elements = xselems

        for k in self._attr_names:
            if (value := kwargs.get(k)) is not None:
                setattr(self, k, value)

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
                    atelem = self.elements[nm].to_at()
                    if isinstance(atelem, (elt.RFCavity, elt.CrabCavity)):
                        cavities.append(atelem)
                    cell_length += getattr(atelem, "Length", 0.0)
                    yield atelem
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
                cav.TimeLag = tl

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
            root:           Dictionary representation of the line or environment,
            use:            For environment, name of the sequence or line containing the
              desired lattice. Default: ``ring``.

        Keyword Args:
            name(str):      Name of the line (default: from JSON file)
            in_file(str):   Name of the source file

        Returns:
            xsline:     New :py:class:`XsLine`
        """
        notavailable = _FilterWarnings("Unknown classes:")
        defclass = "Environment" if "lines" in root else "Line"
        elements = {
            nm: XsElement.static_from_dict(el, name=nm, do_warn=notavailable)
            for nm, el in root["elements"].items()
        }
        if root.get("__class__", defclass) == "Line":
            line = root
        else:
            lines = root["lines"]
            lkeys = list(lines.keys())
            if use is None:
                use = lkeys[0] if len(lkeys) == 1 else "ring"
            if use not in lkeys:
                msg = (f"Selected '{use}' line does not exist, choose from: {lkeys}")
                raise KeyError(msg)
            line = lines[use]
            line.setdefault("name", use)
        line.update(kwargs)
        ln = cls(elements, **line)

        notavailable.summarise()
        return ln

    @classmethod
    def from_json(
        cls, filename: str | Path, use: str | None = None, **kwargs
    ) -> XsLine:
        """Create a XsLine from a JSON file.

        Args:
            filename:   Name or path of the JSON file.
            use:        For environment, name of the sequence or line containing the
              desired lattice. Default: ``ring``

        Keyword Args:
            name(str):  Name of the line (default: from JSON file)

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
    def from_at(cls, lattice: Lattice, match_model: bool = True) -> XsLine:
        """Create a XsLine from an AT :py:class:`.Lattice`.

        Args:
            lattice:        AT lattice
            match_model:    If :py:obj:`True`, set the Xsuite model and integrator
              matching at best the AT PassMethod. If :py:obj:`False`, the Xsuite
              default model and integrator will be used.

        Returns:
            xsline:     New :py:class:`XsLine`
        """

        def refpart(rng):
            prt = rng.particle
            if prt.name == "relativistic":
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

        ring = _separate_apertures(lattice)
        element_names = unique_names([el.FamName for el in ring])
        notavailable = _FilterWarnings("Unknown classes:")
        ln = cls(
            {
                nm: XsElement.static_from_at(
                    el, match_model=match_model, do_warn=notavailable
                )
                for nm, el in zip(element_names, ring, strict=True)
            },
            element_names=element_names,
            particle_ref=refpart(ring),
            name=ring.name,
        )
        notavailable.summarise()
        return ln


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
    match_model: bool = True,
    compact: bool = False,
) -> None:
    """Save a :py:class:`.Lattice` as a Xsuite JSON file.

    Parameters:
        lattice:        Lattice description
        filename:       Name of the JSON file. Default: outputs on
          :py:obj:`sys.stdout`
        match_model:    If :py:obj:`True`, set the Xsuite model matching at best
          the AT PassMethod. If :py:obj:`False`, the Xsuite default model will be used.
        compact:        If :py:obj:`False` (default), the JSON file is pretty-printed
          with line feeds and indentation. Otherwise, the output is a single line.
    """
    XsLine.from_at(lattice, match_model=match_model).to_json(filename, compact=compact)
