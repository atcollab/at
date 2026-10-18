"""
Conversion utilities for creating pyat elements.
"""

from __future__ import annotations

__all__ = [
    "RingParam",
    "element_from_dict",
    "find_class",
    "hide_under",
    "keep_attributes",
    "keep_elements",
    "protect",
    "restore",
    "show_under",
    "split_ignoring_parentheses",
]

import collections
import re
from typing import ClassVar
from warnings import warn, filterwarnings
from collections.abc import Generator

import numpy as np

from at.lattice import AtWarning
from at.lattice import elements as elt
from at.lattice import Lattice, Particle, Element, Drift, Marker

_plh = "placeholder"

filterwarnings("always", category=AtWarning, module=__name__)


def _particle(value) -> Particle:
    if isinstance(value, Particle):
        # Create from python: save_mat
        return value
    elif isinstance(value, dict):
        # Create from Matlab: load_mat
        return Particle(**value)
    else:
        return Particle(value)


def _warn(idx: int | None, message: str, elem: dict, stacklevel: int = 2) -> None:
    name = elem.get("FamName", "")
    location = f'"{name}"' if idx is None else f'{idx} ("{name}")'
    warning = f"\nIn element {location}: {message}\n{elem}\n"
    warn(AtWarning(warning), stacklevel=stacklevel)


class RingParam(elt.Element):
    """Private class for Matlab RingParam element.

    :meta private:
    """

    # noinspection PyProtectedMember
    _BUILD_ATTRIBUTES: ClassVar[list[str]] = [
        *elt.Element._BUILD_ATTRIBUTES,
        "Energy",
        "Periodicity",
    ]
    _conversions: ClassVar[dict] = dict(
        elt.Element._conversions,
        Energy=float,
        Periodicity=int,
        Particle=_particle,
        cell_harmnumber=float,
    )

    # noinspection PyPep8Naming
    def __init__(
        self,
        FamName: str,
        Energy: float,
        Periodicity: int,
        **kwargs,
    ):
        if not np.isnan(float(Energy)):
            kwargs["Energy"] = Energy
        kwargs.setdefault("PassMethod", "IdentityPass")
        super().__init__(FamName, Periodicity=Periodicity, **kwargs)


_alias_map: dict[str, type[Element]] = {
    "bend": elt.Dipole,
    "rbend": elt.Dipole,
    "sbend": elt.Dipole,
    "quad": elt.Quadrupole,
    "sext": elt.Sextupole,
    "rf": elt.RFCavity,
    "bpm": elt.Monitor,
    "ap": elt.Aperture,
    "ringparam": RingParam,
    "wig": elt.Wiggler,
    "matrix66": elt.M66,
    "M66": elt.M66,
}

# Map class names to Element classes
_CLASS_MAP: dict[str, type[Element]] = {
    cls.__name__.lower(): cls for cls in Element.subclasses()
}
_CLASS_MAP.update(_alias_map)

# Maps passmethods to Element classe names
_PASS_MAP: dict[str, str] = {
    "BendLinearPass": "Dipole",
    "BndMPoleSymplectic4RadPass": "Dipole",
    "BndMPoleSymplectic4Pass": "Dipole",
    "QuadLinearPass": "Quadrupole",
    "StrMPoleSymplectic4Pass": "Multipole",
    "StrMPoleSymplectic4RadPass": "Multipole",
    "CorrectorPass": "Corrector",
    "CavityPass": "RFCavity",
    "RFCavityPass": "RFCavity",
    "ThinMPolePass": "ThinMultipole",
    "Matrix66Pass": "M66",
    "AperturePass": "Aperture",
    "IdTablePass": "InsertionDeviceKickMap",
    "GWigSymplecticPass": "Wiggler",
    "VariableThinMPolePass": "VariableThinMultipole",
}

# Lattice attributes which must be dropped when writing a file
_drop_attrs: dict[str, str | None] = {
    "in_file": None,
    "use": None,
    "mat_key": None,
    "mat_file": None,  # Not used any more...
    "m_file": None,
    "repr_file": None,
}


def hide_under(key: str) -> str:
    """Hide leading underscores in attribute names.

    Leading underscores are not allowed in Matlab field names, so they are prepended
    with an x when saving to file.
    """
    return f"x{key}" if key.startswith("_") else key


def show_under(key: str) -> str:
    """Restore leading underscores in attribute names.
    """
    return key[1:] if key.startswith("x_") else key


def _hasattrs(kwargs: dict, *attributes) -> bool:
    """Checks the presence of keys in a :py:class:`dict`.

    Returns :py:obj:`True` if any of the ``attributes`` is in ``kwargs``

    Args:
        kwargs:     The dictionary of keyword arguments passed to the
          Element constructor.
        attributes: A list of strings, the attribute names to be checked.

    Returns:
        found (bool):   :py:obj:`True` if the element has any of the specified
          attributes.
    """
    return any(attribute in kwargs for attribute in attributes)


def keep_attributes(ring: Lattice):
    """Remove Lattice attributes which must not be saved on file."""
    return {k: v for k, v in ring.attrs.items() if _drop_attrs.get(k, k) is not None}


def keep_elements(ring: Lattice) -> Generator[Element, None, None]:
    """Remove the 'RingParam' Marker."""
    for elem in ring:
        if not (isinstance(elem, Marker) and getattr(elem, "tag", None) == "RingParam"):
            yield elem


def find_class(elem: dict, *, index: int | None = None, quiet: bool = False) -> str:
    """Deduce the class of an element from its attributes.

    `find_class` tries to deduce the class from "FamName", from "PassMethod",
    and finally form the element contents.

    Args:
        elem:           Dictionary of element attributes.
        quiet:          Suppress the warning for non-standard classes
        index:          Element index in the lattice

    Returns:
        element_class:  The guessed Class name
    """

    def low_order(key):
        polynom = np.array(elem[key], dtype=np.float64).reshape(-1)
        try:
            low = np.where(polynom != 0.0)[0][0]
        except IndexError:
            low = -1
        return low

    elname = elem.get("FamName", "")
    length = float(elem.get("Length", 0.0))
    pass_method = elem.get("PassMethod", "")

    if not quiet:  # try from element name
        _warn(index, "No class found, looking at the element name.", elem)
    cls = _CLASS_MAP.get(elname.lower())
    if cls is not None:
        return cls.__name__

    if not quiet:  # try from passmethod
        _warn(index, "No class found, looking at the element passmethod.", elem)
    classname = _PASS_MAP.get(pass_method)
    if classname is not None:
        return classname

    if not quiet:  # try from element contents
        _warn(index, "No class found, looking at the element contents.", elem)
    if _hasattrs(elem, "BendingAngle", "EntranceAngle", "ExitAngle"):
        return "Dipole"
    elif _hasattrs(elem, "Voltage", "Frequency", "HarmNumber"):
        return "RFCavity"
    elif _hasattrs(elem, "Periodicity"):
        return "RingParam"
    elif _hasattrs(elem, "Limits"):
        return "Aperture"
    elif _hasattrs(elem, "M66"):
        return "M66"
    elif _hasattrs(elem, "K"):
        return "Quadrupole"
    elif _hasattrs(elem, "PolynomB", "PolynomA"):
        loworder = low_order("PolynomB")
        if loworder == 1:
            return "Quadrupole"
        elif loworder == 2:
            return "Sextupole"
        elif loworder == 3:
            return "Octupole"
        elif pass_method.startswith("StrMPoleSymplectic4") or (length > 0):
            return "Multipole"
        else:
            return "ThinMultipole"
    elif _hasattrs(elem, "KickAngle"):
        return "Corrector"
    elif length > 0.0:
        return "Drift"
    elif _hasattrs(elem, "GCR"):
        return "Monitor"
    else:
        if not quiet:
            _warn(index, "No class found, using 'Marker'.", elem)
        return "Marker"


def element_from_dict(
    elem_dict: dict,
    index: int | None = None,
    quiet: bool = False,
) -> Element:
    """Builds an :py:class:`.Element` from a dictionary of attributes.

    Parameters:
        elem_dict:      Dictionary of element attributes
        index:          Element index in the lattice
        quiet:          Suppress the warning for undefined classes

    Returns:
        elem (Element): new :py:class:`.Element`
    """
    # Translate Matlab field names back to python
    dct = {show_under(k): v for k, v in elem_dict.items()}
    # Extract the class name, or try to deduce it from the element contents
    try:
        classname = dct.pop("Class")
    except KeyError:
        classname = find_class(elem_dict, index=index, quiet=quiet)
    # Get the class from the class name, or use Drift/Marker if not found
    try:
        cls = _CLASS_MAP[classname.lower()]
    except KeyError:
        if elem_dict.get("Length", 0.0) > 0:
            cls = Drift
            elem_dict["PassMethod"] = "DriftPass"
        else:
            cls = Marker
            elem_dict["PassMethod"] = "IdentityPass"
        msg = f"Class '{classname}' does not exist, using `{cls.__name__}` instead."
        _warn(index, msg, elem_dict, stacklevel=3)
    return cls.new_from_dict(**dct)


def split_ignoring_parentheses(
    string: str,
    delimiter: str = ",",
    fence: tuple[str, str] = ("\\(", "\\)"),
    maxsplit: int = -1,
) -> list[str]:
    """Split a string while keeping protected expressions intact.

    Example: "l=0,hom(4,0.0,0)" -> ["l=0", "hom(4,0.0,0)"]
    """
    substituted, matches = protect(string, fence=fence)
    parts = substituted.split(delimiter, maxsplit=maxsplit)
    return restore(matches, *parts)


def protect(
    string: str,
    fence: tuple[str, str] = ('"', '"'),
    *,
    placeholder: str = _plh,
):
    inf, outf = fence
    pattern = f"{inf}[^{inf}]*?{outf}"
    substituted = string[:]
    matches = collections.deque(re.finditer(pattern, string))
    for match in matches:
        substituted = substituted.replace(match.group(), placeholder, 1)
    return substituted, (placeholder, matches)


def restore(replmatch, *parts):
    def rep(part):
        while placeholder in part:
            next_match = matches.popleft()
            part = part.replace(placeholder, next_match.group(), 1)
        return part

    placeholder, matches = replmatch
    replaced_parts = [rep(part) for part in parts]
    assert not matches

    return replaced_parts
