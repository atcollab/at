"""
Conversion utilities for creating pyat elements
"""

from __future__ import annotations

__all__ = [
    "element_classes",
    "element_from_dict",
    "element_to_dict",
    "find_class",
    "keep_elements",
    "keep_attributes",
    "split_ignoring_parentheses",
    "RingParam",
    "protect",
    "restore",
]

import collections
import os
import re
import sysconfig
from typing import Any
from warnings import warn
from collections.abc import Callable, Generator

import numpy as np

from at import integrators
from at.lattice import AtWarning
from at.lattice import elements as elt
from at.lattice import Lattice, Particle, Element, Marker
from at.lattice import idtable_element

_ext_suffix = sysconfig.get_config_var("EXT_SUFFIX")
_plh = "placeholder"


def _no_encoder(v):
    """type encoding for .mat files"""
    return v


def _particle(value) -> Particle:
    if isinstance(value, Particle):
        # Create from python: save_mat
        return value
    elif isinstance(value, dict):
        # Create from Matlab: load_mat
        return Particle(**value)
    else:
        return Particle(value)


def _warn(index: int, message: str, elem_dict: dict) -> None:
    name = elem_dict.get("FamName", "")
    location = f'"{name}":\n' if index is None else f'{index} ("{name}"):\n'
    warning = "".join(("In element ", location, message, f"\n{elem_dict}\n"))
    warn(AtWarning(warning), stacklevel=2)


def element_classes() -> frozenset[type[Element]]:
    """Build a set of all Element subclasses"""

    # Misses class aliases (Bend, Matrix66)
    def subclasses_recursive(cl):
        direct = cl.__subclasses__()
        indirect = []
        for subclass in direct:
            indirect.extend(subclasses_recursive(subclass))
        return frozenset([cl] + direct + indirect)

    return subclasses_recursive(Element)


class RingParam(elt.Element):
    """Private class for Matlab RingParam element

    :meta private:
    """

    # noinspection PyProtectedMember
    _BUILD_ATTRIBUTES = elt.Element._BUILD_ATTRIBUTES + [
        "Energy",
        "Periodicity",
    ]
    _conversions = dict(
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


_alias_map = {
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
_CLASS_MAP = {cls.__name__.lower(): cls for cls in element_classes()}
_CLASS_MAP.update(_alias_map)

# Maps passmethods to Element classes
_PASS_MAP = {
    "BendLinearPass": elt.Dipole,
    "BndMPoleSymplectic4RadPass": elt.Dipole,
    "BndMPoleSymplectic4Pass": elt.Dipole,
    "QuadLinearPass": elt.Quadrupole,
    "StrMPoleSymplectic4Pass": elt.Multipole,
    "StrMPoleSymplectic4RadPass": elt.Multipole,
    "CorrectorPass": elt.Corrector,
    "CavityPass": elt.RFCavity,
    "RFCavityPass": elt.RFCavity,
    "ThinMPolePass": elt.ThinMultipole,
    "Matrix66Pass": elt.M66,
    "AperturePass": elt.Aperture,
    "IdTablePass": idtable_element.InsertionDeviceKickMap,
    "GWigSymplecticPass": elt.Wiggler,
}

# Maps python class name to Matlab class
# Default: element_class.__name__
_mat_class = {
    "Dipole": "Bend",
    "M66": "Matrix66",
}

# Lattice attributes which must be dropped when writing a file
_drop_attrs = {
    "in_file": None,
    "use": None,
    "mat_key": None,
    "mat_file": None,  # Not used anymore...
    "m_file": None,
    "repr_file": None,
}


def _hasattrs(kwargs: dict, *attributes) -> bool:
    """Checks the presence of keys in a :py:class:`dict`

    Returns :py:obj:`True` if any of the ``attributes`` is in ``kwargs``

    Args:
        kwargs:     The dictionary of keyword arguments passed to the
          Element constructor.
        attributes: A list of strings, the attribute names to be checked.

    Returns:
        found (bool):   :py:obj:`True` if the element has any of the specified
          attributes.
    """
    for attribute in attributes:
        if attribute in kwargs:
            return True
    return False


def keep_attributes(ring: Lattice):
    """Remove Lattice attributes which must not be saved on file"""
    return {k: v for k, v in ring.attrs.items() if _drop_attrs.get(k, k) is not None}


def keep_elements(ring: Lattice) -> Generator[Element, None, None]:
    """Remove the 'RingParam' Marker"""
    for elem in ring:
        if not (isinstance(elem, Marker) and getattr(elem, "tag", None) == "RingParam"):
            yield elem


def _from_contents(elem: dict) -> type[Element]:
    """Deduce the element class from its contents"""

    def low_order(key):
        polynom = np.array(elem[key], dtype=np.float64).reshape(-1)
        try:
            low = np.where(polynom != 0.0)[0][0]
        except IndexError:
            low = -1
        return low

    length = float(elem.get("Length", 0.0))
    pass_method = elem.get("PassMethod", "")
    if _hasattrs(
        elem, "FullGap", "FringeInt1", "FringeInt2", "gK", "EntranceAngle", "ExitAngle"
    ):
        return elt.Dipole
    elif _hasattrs(elem, "Voltage", "Frequency", "HarmNumber", "PhaseLag", "TimeLag"):
        return elt.RFCavity
    elif _hasattrs(elem, "Periodicity"):
        # noinspection PyProtectedMember
        return RingParam
    elif _hasattrs(elem, "Limits"):
        return elt.Aperture
    elif _hasattrs(elem, "M66"):
        return elt.M66
    elif _hasattrs(elem, "K"):
        return elt.Quadrupole
    elif _hasattrs(elem, "PolynomB", "PolynomA"):
        loworder = low_order("PolynomB")
        if loworder == 1:
            return elt.Quadrupole
        elif loworder == 2:
            return elt.Sextupole
        elif loworder == 3:
            return elt.Octupole
        elif pass_method.startswith("StrMPoleSymplectic4") or (length > 0):
            return elt.Multipole
        else:
            return elt.ThinMultipole
    elif _hasattrs(elem, "KickAngle"):
        return elt.Corrector
    elif length > 0.0:
        return elt.Drift
    elif _hasattrs(elem, "GCR"):
        return elt.Monitor
    elif pass_method == "IdentityPass":
        return elt.Marker
    else:
        return elt.Element


def find_class(
    elem_dict: dict, quiet: bool = False, index: int | None = None
) -> type(Element):
    """Deduce the class of an element from its attributes

    `find_class` looks first at the "Class" field, if existing. It then tries to deduce
    the class from "FamName", from "PassMethod", and finally form the element contents.

    Args:
        elem_dict:      The dictionary of keyword arguments passed to the
                        Element constructor.
        quiet:          Suppress the warning for non-standard classes
        index:          Element index in the lattice

    Returns:
        element_class:  The guessed Class name
    """

    def check_class(clname):
        if clname:
            _warn(index, f"Class '{clname}' does not exist.", elem_dict)

    def check_pass(passm):
        if not passm:
            _warn(index, "No PassMethod provided.", elem_dict)
        elif not passm.endswith("Pass"):
            message = (
                f"Invalid PassMethod '{passm}': "
                "provided pass methods should end in 'Pass'."
            )
            _warn(index, message, elem_dict)

    class_name = elem_dict.pop("Class", "")  # try from class name
    cls = _CLASS_MAP.get(class_name.lower(), None)
    if cls is not None:
        return cls
    elif not quiet:
        check_class(class_name)

    elname = elem_dict.get("FamName", "")  # try from element name
    cls = _CLASS_MAP.get(elname.lower(), None)
    if cls is not None:
        return cls

    pass_method = elem_dict.get("PassMethod", "")  # try from passmethod
    cls = _PASS_MAP.get(pass_method, None)
    if cls is not None:
        return cls
    elif not quiet:
        check_pass(pass_method)

    return _from_contents(elem_dict)  # look for contents


def element_from_dict(
    elem_dict: dict,
    index: int | None = None,
    check: bool = True,
    quiet: bool = False,
) -> Element:
    """Builds an :py:class:`.Element` from a dictionary of attributes

    Parameters:
        elem_dict:      Dictionary of element attributes
        index:          Element index
        check:          Check the compatibility of class and PassMethod
        quiet:          Suppress the warning for non-standard classes

    Returns:
        elem (Element): new :py:class:`.Element`
    """

    # noinspection PyShadowingNames
    def sanitise_class(index, cls, elem_dict):
        """Checks that the Class and PassMethod of the element are a valid
            combination. Some Classes and PassMethods are incompatible and
            would raise errors during calculation, so we send a
            warning here.

        Args:
            index:          element index
            cls:            Proposed class
            elem_dict:      The dictionary of keyword arguments passed to the
                            Element constructor.

        Raises:
            AttributeError: if the PassMethod and Class are incompatible.
        """
        class_name = cls.__name__
        pass_method = elem_dict.get("PassMethod")
        if pass_method is not None:
            pass_to_class = _PASS_MAP.get(pass_method)
            length = float(elem_dict.get("Length", 0.0))
            file_name = pass_method + _ext_suffix
            file_path = os.path.join(integrators.__path__[0], file_name)
            if not os.path.isfile(os.path.realpath(file_path)):
                message = f"PassMethod {pass_method} is missing {file_name}."
                _warn(index, message, elem_dict)
            elif (pass_method == "IdentityPass") and (length != 0.0):
                message = (
                    f"PassMethod {pass_method} is not compatible with length {length}."
                )
                _warn(index, message, elem_dict)
            elif pass_to_class is not None:
                if not issubclass(cls, pass_to_class):
                    message = (
                        f"PassMethod {pass_method} is not compatible "
                        f"with Class {class_name}."
                    )
                    _warn(index, message, elem_dict)

    cls = find_class(elem_dict, quiet=quiet, index=index)
    if check:
        sanitise_class(index, cls, elem_dict)
    # Remove mandatory attributes from the keyword arguments.
    # Create list rather than generator to ensure that elements are removed
    # from elem_dict.
    elem_args = [elem_dict.pop(attr, None) for attr in cls._BUILD_ATTRIBUTES]
    element = cls(*(arg for arg in elem_args if arg is not None), **elem_dict)
    return element


def element_to_dict(elem: Element, encoder: Callable[[Any], Any] = _no_encoder) -> dict:
    """Convert a :py:class:`.Element` to a :py:class:`dict`

    Parameters:
        elem:           :py:class:`.Element`
        encoder:        data converter

    Returns:
        dct (dict):     Dictionary of :py:class:`.Element` attributes
    """
    dct = {k: encoder(v) for k, v in elem.items()}
    class_name = elem.__class__.__name__
    dct["Class"] = _mat_class.get(class_name, class_name)
    return dct


def split_ignoring_parentheses(
    string: str,
    delimiter: str = ",",
    fence: tuple[str, str] = ("\\(", "\\)"),
    maxsplit: int = -1,
) -> list[str]:
    """Split a string while keeping protected expressions intact

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


Element.from_dict = staticmethod(element_from_dict)
Element.to_dict = element_to_dict
