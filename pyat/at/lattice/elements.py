"""
Module to define common elements used in AT.

Each element has a default PassMethod attribute for which it should have the
appropriate attributes.  If a different PassMethod is set, it is the caller's
responsibility to ensure that the appropriate attributes are present.
"""

from __future__ import annotations

import abc
import re
from abc import ABC
from collections.abc import Generator, Iterable
from copy import copy, deepcopy
from typing import Any, Optional

import numpy as np

# noinspection PyProtectedMember
from .variables import _nop

_zero6 = np.zeros(6)
_eye6 = np.eye(6, order="F")


def _array(value, shape=(-1,), dtype=np.float64):
    # Ensure proper ordering(F) and alignment(A) for "C" access in integrators
    return np.require(value, dtype=dtype, requirements=["F", "A"]).reshape(
        shape, order="F"
    )


def _array66(value):
    return _array(value, shape=(6, 6))


def _float(value) -> float:
    return float(value)


def _int(value, vmin: int | None = None, vmax: int | None = None) -> int:
    intv = int(value)
    if vmin is not None and intv < vmin:
        raise ValueError(f"Value must be greater of equal to {vmin}")
    if vmax is not None and intv > vmax:
        raise ValueError(f"Value must be smaller of equal to {vmax}")
    return intv


class LongtMotion(ABC):
    """Abstract Base class for all Element classes whose instances may modify
    the particle momentum

    Allows identifying elements potentially inducing longitudinal motion.

    Subclasses of :py:class:`LongtMotion` must provide two methods for
    enabling longitudinal motion:

    * ``_get_longt_motion(self)`` must return the activation state,
    * ``set_longt_motion(self, enable, new_pass=None, copy=False, **kwargs)``
      must enable or disable longitudinal motion.
    """

    @abc.abstractmethod
    def _get_longt_motion(self):
        return False

    # noinspection PyShadowingNames
    @abc.abstractmethod
    def set_longt_motion(self, enable, new_pass=None, copy=False, **kwargs):
        """Enable/Disable longitudinal motion

        Parameters:
            enable:     :py:obj:`True`: for enabling, :py:obj:`False` for
              disabling
            new_pass:   New PassMethod:

              * :py:obj:`None`: makes no change,
              * ``'auto'``: Uses the default conversion,
              * Anything else is used as the new PassMethod.
            copy:       If True, returns a modified copy of the element,
              otherwise modifies the element in-place
        """
        # noinspection PyUnresolvedReferences
        if new_pass is None or new_pass == self.PassMethod:
            return self if copy else None
        if copy:
            newelem = deepcopy(self)
            newelem.PassMethod = new_pass
            return newelem
        # noinspection PyAttributeOutsideInit
        self.PassMethod = new_pass


# noinspection PyUnresolvedReferences
class _DictLongtMotion(LongtMotion):
    # noinspection PyShadowingNames
    """Mixin class for elements implementing a 'default_pass' class attribute

    :py:class:`DictLongtMotion` provides:

    * a :py:meth:`set_longt_motion` method setting the PassMethod according
      to the ``default_pass`` dictionary.
    * a :py:obj:`.longt_motion` property set to :py:obj:`True` when the
      PassMethod is ``default_pass[True]``

    The class must have a ``default_pass`` class attribute, a dictionary
    such that:

    * ``default_pass[False]`` is the PassMethod when radiation is turned
      OFF,
    * ``default_pass[True]`` is the default PassMethod when radiation is
      turned ON.

    The :py:class:`DictLongtMotion` class must be set as the first base class.

    Example:

        >>> class QuantumDiffusion(_DictLongtMotion, Element):
        ...     default_pass = {False: "IdentityPass", True: "QuantDiffPass"}

        Defines a class such that :py:meth:`set_longt_motion` will select
        ``'IdentityPass'`` or ``'IdentityPass'``.
    """

    def _get_longt_motion(self):
        return self.PassMethod != self.default_pass[False]

    # noinspection PyShadowingNames
    def set_longt_motion(self, enable, new_pass=None, **kwargs):
        if new_pass == "auto":
            new_pass = self.default_pass[enable]
        return super().set_longt_motion(enable, new_pass=new_pass, **kwargs)


# noinspection PyUnresolvedReferences
class _Radiative(LongtMotion):
    # noinspection PyShadowingNames
    r"""Mixin class for radiating elements

    :py:class:`_Radiative` implements the mechanism for converting the pass
    methods of radiating elements. It provides:

    * a :py:meth:`set_longt_motion` method setting the PassMethod
      according to the following rule:

      * ``enable == True``: replace "\*Pass" by "\*RadPass"
      * ``enable == False``: replace "\*RadPass" by "\*Pass"
    * a :py:obj:`.longt_motion` property set to true when the PassMethod
      ends with "RadPass"

    The :py:class:`_Radiative` class must be set as the first base class.

    Example:
        >>> class Multipole(_Radiative, LongElement, ThinMultipole):

        Defines a class where :py:meth:`set_longt_motion` will convert the
        PassMethod according to the \*Pass or \*RadPass suffix.
    """

    def _get_longt_motion(self):
        return self.PassMethod.endswith(("RadPass", "QuantPass"))

    def _autopass(self, enable):
        if enable:
            root = self.PassMethod.replace("QuantPass", "Pass").replace(
                "RadPass", "Pass"
            )
            return "".join((root[:-4], "RadPass"))
        elif self.longt_motion:
            root = self.PassMethod.replace("QuantPass", "Pass").replace(
                "RadPass", "Pass"
            )
            return root
        else:
            return None

    # noinspection PyTypeChecker,PyShadowingNames
    def set_longt_motion(self, enable, new_pass=None, copy=False, **kwargs):
        if new_pass == "auto":
            new_pass = self._autopass(enable)
        if new_pass is None or new_pass == self.PassMethod:
            return self if copy else None
        if enable:

            def setpass(el):
                el.PassMethod = new_pass
                el.Energy = kwargs["energy"]

        else:

            def setpass(el):
                el.PassMethod = new_pass
                try:
                    del el.Energy
                except AttributeError:
                    pass

        if copy:
            newelem = deepcopy(self)
            setpass(newelem)
            return newelem
        setpass(self)


class Radiative(_Radiative):
    # noinspection PyShadowingNames
    r"""Mixin class for default radiating elements (:py:class:`.Dipole`,
    :py:class:`.Quadrupole`, :py:class:`.Wiggler`)

    :py:class:`Radiative` is a base class for the subset of radiative elements
    considered as the ones to be turned on by default: :py:class:`.Dipole`,
    :py:class:`.Quadrupole` and :py:class:`.Wiggler`, excluding the higher
    order multipoles.

    :py:class:`Radiative` inherits from :py:class:`_Radiative` and does not
    add any new functionality. Its purpose is to identify the default set of
    radiating elements.

    Example:
        >>> class Dipole(Radiative, Multipole):

        Defines a class belonging to the default radiating elements. It
        converts the PassMethod according to the "\*Pass" or "\*RadPass"
        suffix.
    """

    pass


class Collective(_DictLongtMotion):
    """Mixin class for elements representing collective effects

    Derived classes will automatically set the
    :py:attr:`~Element.is_collective` property when the element is active.

    The class must have a ``default_pass`` class attribute, a dictionary such
    that:

    * ``default_pass[False]`` is the PassMethod when collective effects
      are turned OFF,
    * ``default_pass[True]`` is the default PassMethod when collective effects
      are turned ON.

    The :py:class:`Collective` class must be set as the first base class.

    Example:
        >>> class WakeElement(Collective, Element):
        ...     default_pass = {False: "IdentityPass", True: "WakeFieldPass"}

        Defines a class where the :py:attr:`~Element.is_collective` property is
        handled
    """

    def _get_collective(self):
        # noinspection PyUnresolvedReferences
        return self.PassMethod != self.default_pass[False]

    @abc.abstractmethod
    def clear_history(self):
        pass


class Element:
    """Base class for AT elements"""

    _BUILD_ATTRIBUTES = ["FamName"]
    _conversions = {
        "FamName": str,
        "PassMethod": str,
        "Length": _float,
        "R1": _array66,
        "R2": _array66,
        "T1": lambda v: _array(v, (6,)),
        "T2": lambda v: _array(v, (6,)),
        "RApertures": lambda v: _array(v, (4,)),
        "EApertures": lambda v: _array(v, (2,)),
        "KickAngle": lambda v: _array(v, (2,)),
        "PolynomB": _array,
        "PolynomA": _array,
        "BendingAngle": _float,
        "MaxOrder": _int,
        "NumIntSteps": lambda v: _int(v, vmin=0),
        "Energy": _float,
    }

    _entrance_fields = ["T1", "R1"]
    _exit_fields = ["T2", "R2"]
    _no_swap = _entrance_fields + _exit_fields

    def __init__(self, family_name: str, **kwargs):
        """
        Parameters:
            family_name:    Name of the element

        All keywords will be set as attributes of the element
        """

        self.FamName = family_name
        self.Length = kwargs.pop("Length", 0.0)
        self.PassMethod = kwargs.pop("PassMethod", "IdentityPass")
        self.update(kwargs)

    def __setattr__(self, key, value):
        try:
            value = self._conversions.get(key, _nop)(value)
        except Exception as exc:
            exc.args = (f"In element {self.FamName}, parameter {key}: {exc}",)
            raise
        else:
            super().__setattr__(key, value)

    def __str__(self):
        return "\n".join(
            [self.__class__.__name__ + ":"]
            + [f"{k:>14}: {v!s}" for k, v in self.items()]
        )

    def __repr__(self):
        clsname, args, kwargs = self.definition
        keywords = [f"{arg!r}" for arg in args]
        keywords += [f"{k}={v!r}" for k, v in kwargs.items()]
        args = re.sub(r"\n\s*", " ", ", ".join(keywords))
        return f"{clsname}({args})"

    def equals(self, other) -> bool:
        """Whether an element is equivalent to another.

        This implementation was found to be too slow for the generic
        __eq__ method when comparing lattices.
        """
        return repr(self) == repr(other)

    def divide(self, frac) -> list[Element]:
        """split the element in len(frac) pieces whose length
        is frac[i]*self.Length

        Parameters:
            frac:           length of each slice expressed as a fraction of the
              initial length. ``sum(frac)`` may differ from 1.

        Returns:
            elem_list:  a list of elements equivalent to the original.

        Example:

            >>> Drift("dr", 0.5).divide([0.2, 0.6, 0.2])
            [Drift('dr', 0.1), Drift('dr', 0.3), Drift('dr', 0.1)]
        """
        # Bx default, the element is indivisible
        return [self]

    def swap_faces(self, copy=False):
        """Swap the faces of an element, alignment errors are ignored"""

        def swapattr(element, attro, attri):
            val = getattr(element, attri)
            delattr(element, attri)
            return attro, val

        if copy:
            el = self.copy()
        else:
            el = self
        # Remove and swap entrance and exit attributes
        fin = dict(
            swapattr(el, kout, kin)
            for kin, kout in zip(el._entrance_fields, el._exit_fields)
            if kin in vars(el) and kin not in el._no_swap
        )
        fout = dict(
            swapattr(el, kin, kout)
            for kin, kout in zip(el._entrance_fields, el._exit_fields)
            if kout in vars(el) and kout not in el._no_swap
        )
        # Apply swapped entrance and exit attributes
        for key, value in fin.items():
            setattr(el, key, value)
        for key, value in fout.items():
            setattr(el, key, value)
        return el if copy else None

    def update(self, *args, **kwargs):
        """
        update(**kwargs)
        update(mapping, **kwargs)
        update(iterable, **kwargs)
        Update the element attributes with the given arguments
        """
        attrs = dict(*args, **kwargs)
        for key, value in attrs.items():
            setattr(self, key, value)

    def copy(self) -> Element:
        """Return a shallow copy of the element"""
        return copy(self)

    def deepcopy(self) -> Element:
        """Return a deep copy of the element"""
        return deepcopy(self)

    @property
    def definition(self) -> tuple[str, tuple, dict]:
        """tuple (class_name, args, kwargs) defining the element"""
        attrs = dict(self.items())
        arguments = tuple(
            attrs.pop(k, getattr(self, k)) for k in self._BUILD_ATTRIBUTES
        )
        defelem = self.__class__(*arguments)
        keywords = {
            k: v
            for k, v in attrs.items()
            if not np.array_equal(v, getattr(defelem, k, None))
        }
        return self.__class__.__name__, arguments, keywords

    def items(self) -> Generator[tuple[str, Any], None, None]:
        """Iterates through the data members"""
        v = vars(self).copy()
        for k in ["FamName", "Length", "PassMethod"]:
            yield k, v.pop(k)
        for k, val in sorted(v.items()):
            yield k, val

    def is_compatible(self, other: Element) -> bool:
        """Checks if another :py:class:`Element` can be merged"""
        return False

    def merge(self, other) -> None:
        """Merge another element"""
        if not self.is_compatible(other):
            badname = getattr(other, "FamName", type(other))
            raise TypeError(f"Cannot merge {self.FamName} and {badname}")

    # noinspection PyMethodMayBeStatic
    def _get_longt_motion(self):
        return False

    # noinspection PyMethodMayBeStatic
    def _get_collective(self):
        return False

    @property
    def longt_motion(self) -> bool:
        """:py:obj:`True` if longitudinal motion is affected by the element"""
        return self._get_longt_motion()

    @property
    def is_collective(self) -> bool:
        """:py:obj:`True` if the element involves collective effects"""
        return self._get_collective()


class LongElement(Element):
    """Base class for long elements"""

    _BUILD_ATTRIBUTES = Element._BUILD_ATTRIBUTES + ["Length"]

    def __init__(self, family_name: str, length: float, *args, **kwargs):
        """
        Args:
            family_name:    Name of the element
            length:         Element length [m]

        Other arguments and keywords are given to the base class
        """
        kwargs.setdefault("Length", length)
        # Ancestor may be either Element or ThinMultipole
        # noinspection PyArgumentList
        super().__init__(family_name, *args, **kwargs)

    def _part(self, fr, sumfr):
        pp = self.copy()
        pp.Length = fr * self.Length
        if hasattr(self, "KickAngle"):
            pp.KickAngle = fr / sumfr * self.KickAngle
        return pp

    def divide(self, frac) -> list[Element]:
        def popattr(element, attr):
            val = getattr(element, attr)
            delattr(element, attr)
            return attr, val

        frac = np.asarray(frac, dtype=float)
        el = self.copy()
        # Remove entrance and exit attributes
        fin = dict(
            popattr(el, key) for key in vars(self) if key in self._entrance_fields
        )
        fout = dict(popattr(el, key) for key in vars(self) if key in self._exit_fields)
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
        for fname in ("RApertures", "EApertures"):
            if not compatible_field(fname):
                return False
        return True

    def merge(self, other) -> None:
        super().merge(other)
        self.Length += other.Length


class Marker(Element):
    """Marker element"""


class Monitor(Element):
    """Monitor element"""


class BeamMoments(Element):
    """Element to compute bunches mean and std"""

    def __init__(self, family_name: str, **kwargs):
        """
        Args:
            family_name:    Name of the element

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
        """Beam 6d standard deviation"""
        return self._stds

    @property
    def means(self):
        """Beam 6d center of mass"""
        return self._means


class SliceMoments(Element):
    """Element computing the mean and std of slices"""

    _BUILD_ATTRIBUTES = Element._BUILD_ATTRIBUTES + ["nslice"]
    _conversions = dict(Element._conversions, nslice=int)

    def __init__(self, family_name: str, nslice: int, **kwargs):
        """
        Args:
            family_name:    Name of the element
            nslice:         Number of slices

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
        """Slices x,y,dp standard deviation"""
        return self._stds.reshape((3, self._nbunch, self.nslice, self._dturns))

    @property
    def means(self):
        """Slices x,y,dp center of mass"""
        return self._means.reshape((3, self._nbunch, self.nslice, self._dturns))

    @property
    def spos(self):
        """Slices s position"""
        return self._spos.reshape((self._nbunch, self.nslice, self._dturns))

    @property
    def weights(self):
        """Slices weights in mA if beam current >0,
        otherwise fraction of total number of
        particles in the bunch
        """
        return self._weights.reshape((self._nbunch, self.nslice, self._dturns))

    @property
    def startturn(self):
        """Start turn of the acquisition"""
        return self._startturn

    @startturn.setter
    def startturn(self, value):
        if value < 0:
            raise ValueError("startturn must be greater or equal to 0")
        if value >= self._endturn:
            raise ValueError("startturn must be smaller than endturn")
        self._startturn = value

    @property
    def endturn(self):
        """End turn of the acquisition"""
        return self._endturn

    @endturn.setter
    def endturn(self, value):
        if value <= 0:
            raise ValueError("endturn must be greater than 0")
        if value <= self._startturn:
            raise ValueError("endturn must be greater than startturn")
        self._endturn = value


class Aperture(Element):
    """Transverse aperture element"""

    _BUILD_ATTRIBUTES = Element._BUILD_ATTRIBUTES + ["Limits"]
    _conversions = dict(Element._conversions, Limits=lambda v: _array(v, (4,)))

    def __init__(self, family_name, limits, **kwargs):
        """
        Args:
            family_name:    Name of the element
            limits:         (4,) array of physical aperture:
              [xmin, xmax, ymin, ymax]
        Default PassMethod: ``AperturePass``
        """
        kwargs.setdefault("PassMethod", "AperturePass")
        super().__init__(family_name, Limits=limits, **kwargs)


class LongtAperture(Element):
    """Longitudinal aperture element"""

    _BUILD_ATTRIBUTES = Element._BUILD_ATTRIBUTES + ["Limits"]
    _conversions = dict(Element._conversions, Limits=lambda v: _array(v, (4,)))

    def __init__(self, family_name, limits, **kwargs):
        """
        Args:
            family_name:    Name of the element
            limits:         (4,) array of physical aperture:
              [dpmin, dpmax, ctmin, ctmax]
        Default PassMethod: ``LongtAperturePass``
        """
        kwargs.setdefault("PassMethod", "LongtAperturePass")
        super().__init__(family_name, Limits=limits, **kwargs)


class Drift(LongElement):
    """Drift space element"""

    def __init__(self, family_name: str, length: float, **kwargs):
        """
        Args:
            family_name:    Name of the element
            length:         Element length [m]

        Default PassMethod: ``DriftPass``
        """
        kwargs.setdefault("PassMethod", "DriftPass")
        super().__init__(family_name, length, **kwargs)

    def insert(
        self, insert_list: Iterable[tuple[float, Element | None]]
    ) -> list[Element]:
        """insert elements inside a drift

        Arguments:
            insert_list: iterable, each item of insert_list is itself an
              iterable with 2 objects:

              1. the location where the center of the element
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
        frac, elements = zip(*insert_list)
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
    """Collimator element"""

    _BUILD_ATTRIBUTES = LongElement._BUILD_ATTRIBUTES + ["RApertures"]

    def __init__(self, family_name: str, length: float, limits, **kwargs):
        """
        Args:
            family_name:    Name of the element
            length:         Element length [m]
            limits:         (4,) array of physical aperture:
              [xmin, xmax, zmin, zmax] [m]

        Default PassMethod: ``DriftPass``
        """
        super().__init__(family_name, length, RApertures=limits, **kwargs)


class ThinMultipole(Element):
    """Thin multipole element"""

    _BUILD_ATTRIBUTES = Element._BUILD_ATTRIBUTES + ["PolynomA", "PolynomB"]

    def __init__(self, family_name: str, poly_a, poly_b, **kwargs):
        """
        Args:
            family_name:    Name of the element
            poly_a:         Array of skew multipole components
            poly_b:         Array of normal multipole components

        Keyword arguments:
            MaxOrder:       Number of desired multipoles. Default: highest
              index of non-zero polynomial coefficients
            FieldScaling:   Scaling factor applied to the magnetic field
              (*PolynomA* and *PolynomB*)

        Default PassMethod: ``ThinMPolePass``
        """

        def getpol(poly):
            nonzero = np.flatnonzero(poly != 0.0)
            return poly, len(poly), nonzero[-1] if len(nonzero) > 0 else -1

        def lengthen(poly, dl):
            if dl > 0:
                return np.concatenate((poly, np.zeros(dl)))
            else:
                return poly

        # PolynomA and PolynomB and convert to ParamArray
        prmpola = self._conversions["PolynomA"](kwargs.pop("PolynomA", poly_a))
        prmpolb = self._conversions["PolynomB"](kwargs.pop("PolynomB", poly_b))
        poly_a, len_a, ord_a = getpol(prmpola)
        poly_b, len_b, ord_b = getpol(prmpolb)
        deforder = max(getattr(self, "DefaultOrder", 0), ord_a, ord_b)
        # Remove MaxOrder
        maxorder = kwargs.pop("MaxOrder", deforder)
        kwargs.setdefault("PassMethod", "ThinMPolePass")
        super().__init__(family_name, **kwargs)
        # Set MaxOrder while PolynomA and PolynomB are not set yet
        super().__setattr__("MaxOrder", maxorder)
        # Adjust polynom lengths and set them
        len_ab = max(self.MaxOrder + 1, len_a, len_b)
        self.PolynomA = lengthen(prmpola, len_ab - len_a)
        self.PolynomB = lengthen(prmpolb, len_ab - len_b)

    def __setattr__(self, key, value):
        """Check the compatibility of MaxOrder, PolynomA and PolynomB"""
        polys = ("PolynomA", "PolynomB")
        if key in polys:
            lmin = self.MaxOrder
            if not len(value) > lmin:
                raise ValueError(f"Length of {key} must be larger than {lmin}")
        elif key == "MaxOrder":
            intval = int(value)
            lmax = min(len(getattr(self, k)) for k in polys)
            if not intval < lmax:
                raise ValueError(f"MaxOrder must be smaller than {lmax}")
        super().__setattr__(key, value)


class Multipole(_Radiative, LongElement, ThinMultipole):
    """Multipole element"""

    _BUILD_ATTRIBUTES = LongElement._BUILD_ATTRIBUTES + ["PolynomA", "PolynomB"]
    _conversions = dict(ThinMultipole._conversions, K=float, H=float)

    def __init__(self, family_name: str, length: float, poly_a, poly_b, **kwargs):
        """
        Args:
            family_name:    Name of the element
            length:         Element length [m]
            poly_a:         Array of skew multipole components
            poly_b:         Array of normal multipole components

        Keyword arguments:
            MaxOrder:       Number of desired multipoles. Default: highest
              index of non-zero polynomial coefficients
            NumIntSteps:    Number of integration steps (default: 10)
            KickAngle:      Correction deviation angles (H, V)
            FieldScaling:   Scaling factor applied to the magnetic field
              (*PolynomA* and *PolynomB*)

        Default PassMethod: ``StrMPoleSymplectic4Pass``
        """
        kwargs.setdefault("PassMethod", "StrMPoleSymplectic4Pass")
        kwargs.setdefault("NumIntSteps", 10)
        super().__init__(family_name, length, poly_a, poly_b, **kwargs)

    def is_compatible(self, other) -> bool:
        if super().is_compatible(other) and self.MaxOrder == other.MaxOrder:
            for i in range(self.MaxOrder + 1):
                if self.PolynomB[i] != other.PolynomB[i]:
                    return False
                if self.PolynomA[i] != other.PolynomA[i]:
                    return False
            return True
        else:
            return False

    # noinspection PyPep8Naming
    @property
    def K(self) -> float:
        """Focusing strength [mˆ-2]"""
        arr = self.PolynomB
        return 0.0 if len(arr) < 2 else float(arr[1])

    # noinspection PyPep8Naming
    @K.setter
    def K(self, strength: float):
        self.PolynomB[1] = strength

    # noinspection PyPep8Naming
    @property
    def H(self) -> float:
        """Sextupolar strength [mˆ-3]"""
        arr = self.PolynomB
        return 0.0 if len(arr) < 3 else float(arr[2])

    # noinspection PyPep8Naming
    @H.setter
    def H(self, strength):
        self.PolynomB[2] = strength


class Dipole(Radiative, Multipole):
    """Dipole element"""

    _BUILD_ATTRIBUTES = LongElement._BUILD_ATTRIBUTES + ["BendingAngle", "K"]
    _conversions = dict(
        Multipole._conversions,
        EntranceAngle=float,
        ExitAngle=float,
        FringeInt1=float,
        FringeInt2=float,
        FringeQuadEntrance=int,
        FringeQuadExit=int,
        FringeBendEntrance=int,
        FringeBendExit=int,
    )

    _entrance_fields = Multipole._entrance_fields + [
        "EntranceAngle",
        "FringeInt1",
        "FringeBendEntrance",
        "FringeQuadEntrance",
    ]
    _exit_fields = Multipole._exit_fields + [
        "ExitAngle",
        "FringeInt2",
        "FringeBendExit",
        "FringeQuadExit",
    ]

    DefaultOrder = 0

    def __init__(
        self,
        family_name: str,
        length: float,
        bending_angle: float | None = 0.0,
        k: float = 0.0,
        **kwargs,
    ):
        """
        Args:
            family_name:    Name of the element
            length:         Element length [m]
            bending_angle:  Bending angle [rd]
            k:              Focusing strength [m^-2]

        Keyword arguments:
            EntranceAngle=0.0:  entrance angle
            ExitAngle=0.0:      exit angle
            PolynomB:           straight multipoles
            PolynomA:           skew multipoles
            MaxOrder=0:         Number of desired multipoles
            NumIntSt=10:        Number of integration steps
            FullGap:            Magnet full gap
            FringeInt1:         Extension of the entrance fringe field
            FringeInt2:         Extension of the exit fringe field
            FringeBendEntrance: 1: legacy version Brown First Order (default)

              2: SOLEIL close to second order of Brown

              3: THOMX
            FringeBendExit:     See *FringeBendEntrance*
            FringeQuadEntrance: 0: no fringe field effect (default)

              1: Lee-Whiting's thin lens limit formula

              2: elegant-like
            FringeQuadExit:     See *FringeQuadEntrance*
            fringeIntM0:        Integrals for FringeQuad method 2
            fringeIntP0:
            KickAngle:          Correction deviation angles (H, V)
            FieldScaling:       Scaling factor applied to the magnetic field

        Available PassMethods: :ref:`BndMPoleSymplectic4Pass`,
        :ref:`BendLinearPass`, :ref:`ExactSectorBendPass`,
        :ref:`ExactRectangularBendPass`, :ref:`ExactRectBendPass`,
        BndStrMPoleSymplectic4Pass

        Default PassMethod: :ref:`BndMPoleSymplectic4Pass`
        """
        kwargs.setdefault("BendingAngle", bending_angle)
        kwargs.setdefault("EntranceAngle", 0.0)
        kwargs.setdefault("ExitAngle", 0.0)
        kwargs.setdefault("PassMethod", "BndMPoleSymplectic4Pass")
        super().__init__(family_name, length, [], [0.0, k], **kwargs)

    def items(self) -> Generator[tuple[str, Any], None, None]:
        yield from super().items()
        yield "K", self.K

    def _part(self, fr, sumfr):
        pp = super()._part(fr, sumfr)
        pp.BendingAngle = fr / sumfr * self.BendingAngle
        pp.EntranceAngle = 0.0
        pp.ExitAngle = 0.0
        return pp

    def is_compatible(self, other) -> bool:
        def invrho(dip: Dipole):
            return dip.BendingAngle / dip.Length

        return (
            super().is_compatible(other)
            and self.ExitAngle == -other.EntranceAngle
            and abs(invrho(self) - invrho(other)) <= 1.0e-6
        )

    def merge(self, other) -> None:
        super().merge(other)
        # noinspection PyAttributeOutsideInit
        self.ExitAngle = other.ExitAngle
        self.BendingAngle += other.BendingAngle


# Bend is a synonym of Dipole.
Bend = Dipole


class Quadrupole(Radiative, Multipole):
    """Quadrupole element"""

    _BUILD_ATTRIBUTES = LongElement._BUILD_ATTRIBUTES + ["K"]
    _conversions = dict(
        Multipole._conversions, FringeQuadEntrance=int, FringeQuadExit=int
    )

    _entrance_fields = Multipole._entrance_fields + ["FringeQuadEntrance"]
    _exit_fields = Multipole._exit_fields + ["FringeQuadExit"]

    DefaultOrder = 1

    def __init__(
        self, family_name: str, length: float, k: float | None = 0.0, **kwargs
    ):
        """Quadrupole(FamName, Length, Strength=0, **keywords)

        Args:
            family_name:    Name of the element
            length:         Element length [m]
            k:              Focusing strength [mˆ-2]

        Keyword Arguments:
            PolynomB:           straight multipoles
            PolynomA:           skew multipoles
            MaxOrder=1:         Number of desired multipoles
            NumIntSteps=10:     Number of integration steps
            FringeQuadEntrance: 0: no fringe field effect (default)

              1: Lee-Whiting's thin lens limit formula

              2: elegant-like
            FringeQuadExit:     See ``FringeQuadEntrance``
            fringeIntM0:        Integrals for FringeQuad method 2
            fringeIntP0:
            KickAngle:          Correction deviation angles (H, V)
            FieldScaling:       Scaling factor applied to the magnetic field
              (*PolynomA* and *PolynomB*)

        Default PassMethod: ``StrMPoleSymplectic4Pass``
        """
        kwargs.setdefault("PassMethod", "StrMPoleSymplectic4Pass")
        super().__init__(family_name, length, [], [0.0, k], **kwargs)

    def items(self) -> Generator[tuple[str, Any], None, None]:
        yield from super().items()
        yield "K", self.K


class Sextupole(Multipole):
    """Sextupole element"""

    _BUILD_ATTRIBUTES = LongElement._BUILD_ATTRIBUTES + ["H"]

    DefaultOrder = 2

    def __init__(
        self, family_name: str, length: float, h: float | None = 0.0, **kwargs
    ):
        """
        Args:
            family_name:    Name of the element
            length:         Element length [m]
            h:              strength [mˆ-3]

        Keyword Arguments:
            PolynomB:           straight multipoles
            PolynomA:           skew multipoles
            MaxOrder:           Number of desired multipoles
            NumIntSteps=10:     Number of integration steps
            KickAngle:          Correction deviation angles (H, V)
            FieldScaling:       Scaling factor applied to the magnetic field
              (*PolynomA* and *PolynomB*)

        Default PassMethod: ``StrMPoleSymplectic4Pass``
        """
        kwargs.setdefault("PassMethod", "StrMPoleSymplectic4Pass")
        super().__init__(family_name, length, [], [0.0, 0.0, h], **kwargs)

    def items(self) -> Generator[tuple[str, Any], None, None]:
        yield from super().items()
        yield "H", self.H


class Octupole(Multipole):
    """Octupole element, with no changes from multipole at present"""

    _BUILD_ATTRIBUTES = Multipole._BUILD_ATTRIBUTES

    DefaultOrder = 3


class RFCavity(LongtMotion, LongElement):
    """RF cavity element"""

    _BUILD_ATTRIBUTES = LongElement._BUILD_ATTRIBUTES + [
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
            energy:         ring energy [eV]

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
    """Linear (6, 6) transfer matrix"""

    _BUILD_ATTRIBUTES = Element._BUILD_ATTRIBUTES + ["M66"]
    _conversions = dict(Element._conversions, M66=_array66)

    def __init__(self, family_name: str, m66=None, **kwargs):
        """
        Args:
            family_name:    Name of the element
            m66:            Transfer matrix. Default: Identity matrix

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

        Optional Args:
            betax:         Horizontal beta function at element [m]
            betay:         Vertical beta function at element [m]
            emitx:         Horizontal equilibrium emittance [m.rad]
            emity:         Vertical equilibrium emittance [m.rad]
            espread:       Equilibrium energy spread
            taux:          Horizontal damping time [turns]
            tauy:          Vertical damping time [turns]
            tauz:          Longitudinal damping time [turns]

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
    """Simple radiation damping and energy loss"""

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

        Optional Args:
            taux:          Horizontal damping time [turns]
            tauy:          Vertical damping time [turns]
            tauz:          Longitudinal damping time [turns]
            U0:            Energy loss per turn [eV]

        Default PassMethod: ``SimpleRadiationRadPass``
        """
        assert taux >= 0.0, "taux must be greater than or equal to 0"
        if taux == 0.0:
            dampx = 1
        else:
            dampx = np.exp(-1 / taux)

        assert tauy >= 0.0, "tauy must be greater than or equal to 0"
        if tauy == 0.0:
            dampy = 1
        else:
            dampy = np.exp(-1 / tauy)

        assert tauz >= 0.0, "tauz must be greater than or equal to 0"
        if tauz == 0.0:
            dampz = 1
        else:
            dampz = np.exp(-1 / tauz)

        kwargs.setdefault("PassMethod", self.default_pass[True])
        kwargs.setdefault("U0", U0)
        kwargs.setdefault(
            "damp_mat_diag", np.array([dampx, dampx, dampy, dampy, dampz, dampz])
        )

        super().__init__(family_name, **kwargs)


class Corrector(LongElement):
    """Corrector element"""

    _BUILD_ATTRIBUTES = LongElement._BUILD_ATTRIBUTES + ["KickAngle"]

    def __init__(self, family_name: str, length: float, kick_angle, **kwargs):
        """
        Args:
            family_name:    Name of the element
            length:         Element length [m]
            KickAngle:      Correction deviation angles (H, V)

        Keyword Args:
            FieldScaling:   Scaling factor applied to the magnetic field
              (*KickAngle*)

        Default PassMethod: ``CorrectorPass``
        """
        kwargs.setdefault("PassMethod", "CorrectorPass")
        super().__init__(family_name, length, KickAngle=kick_angle, **kwargs)


class Wiggler(Radiative, LongElement):
    """Wiggler element

    See atwiggler.m
    """

    _BUILD_ATTRIBUTES = LongElement._BUILD_ATTRIBUTES + ["Lw", "Bmax"]
    _conversions = dict(
        Element._conversions,
        Lw=_float,
        Bmax=_float,
        Energy=_float,
        Bx=lambda v: _array(v, (6, -1)),
        By=lambda v: _array(v, (6, -1)),
        Nstep=int,
        Nmeth=int,
        NHharm=int,
        NVharm=int,
    )

    # noinspection PyPep8Naming
    def __init__(
        self,
        family_name: str,
        length: float,
        wiggle_period: float,
        b_max: float,
        energy: float = 0.0,
        *,
        Nstep: int | None = 5,
        Nmeth: int | None = 4,
        By=(1, 1, 0, 1, 1, 0),
        Bx=(),
        **kwargs,
    ):
        """
        Args:
            length:         total length of the wiggler
            wiggle_period:  length must be a multiple of this
            b_max:          peak wiggler field [Tesla]
            energy:         beam energy [eV]
            Nstep:          number of integration steps.
            Nmeth:          symplectic integration order: 2 or 4
            Bx:             harmonics for horizontal wiggler: (6, nHharm)
                              array-like object
            By:             harmonics for vertical wiggler (6, nHharm)
                              array-like object

        Default PassMethod: ``GWigSymplecticPass``
        """
        kwargs.setdefault("PassMethod", "GWigSymplecticPass")
        n_wiggles = length / wiggle_period
        if abs(round(n_wiggles) - n_wiggles) > 1e-6:
            raise ValueError(
                "Wiggler: length / wiggle_period is not an "
                f"integer. ({length}/{wiggle_period}={n_wiggles})"
            )
        super().__init__(
            family_name,
            length,
            Lw=wiggle_period,
            Bmax=b_max,
            Nstep=Nstep,
            Nmeth=Nmeth,
            By=By,
            Bx=Bx,
            Energy=energy,
            **kwargs,
        )

        for i, b in enumerate(self.By.T):
            dk = abs(b[3] ** 2 - b[4] ** 2 - b[2] ** 2) / abs(b[4])
            if dk > 1e-6:
                raise ValueError(f"Wiggler(H): kx^2 + kz^2 -ky^2 !=0, i = {i}")

        for i, b in enumerate(self.Bx.T):
            dk = abs(b[2] ** 2 - b[4] ** 2 - b[3] ** 2) / abs(b[4])
            if dk > 1e-6:
                raise ValueError(f"Wiggler(V): ky^2 + kz^2 -kx^2 !=0, i = {i}")

        self.NHharm = self.By.shape[1]
        self.NVharm = self.Bx.shape[1]

    def divide(self, frac) -> list[Element]:
        # A wiggler is indivisible
        return [self]


class QuantumDiffusion(_DictLongtMotion, Element):
    _BUILD_ATTRIBUTES = Element._BUILD_ATTRIBUTES + ["Lmatp"]
    default_pass = {False: "IdentityPass", True: "QuantDiffPass"}
    _conversions = dict(Element._conversions, Lmatp=_array66)

    def __init__(self, family_name: str, lmatp: np.ndarray, **kwargs):
        """Quantum diffusion element

        Args:
            family_name:    Name of the element
            lmatp      :    Diffusion matrix for generation (see
              :py:func:`.gen_quantdiff_elem`)

        Default PassMethod: ``QuantDiffPass``
        """
        kwargs.setdefault("PassMethod", self.default_pass[True])
        super().__init__(family_name, Lmatp=lmatp, **kwargs)


class EnergyLoss(_DictLongtMotion, Element):
    _BUILD_ATTRIBUTES = Element._BUILD_ATTRIBUTES + ["EnergyLoss"]
    _conversions = dict(Element._conversions, EnergyLoss=float)
    default_pass = {False: "IdentityPass", True: "EnergyLossRadPass"}

    def __init__(self, family_name: str, energy_loss: float, **kwargs):
        """Energy loss element

        the :py:class:`EnergyLoss` element is taken into account in
        :py:func:`.radiation_parameters`: it adds damping by contributing to the
        :math:`I_2` integral, thus reducing the equilibrium emittance. But it does not
        generate any diffusion. This makes sense only if the losses summarised in
        the element occur in non-dispersive locations.

        Args:
            family_name:    Name of the element
            energy_loss:    Energy loss [eV]

        """
        kwargs.setdefault("PassMethod", self.default_pass[False])
        super().__init__(family_name, EnergyLoss=energy_loss, **kwargs)


Radiative.register(EnergyLoss)


def build_class_map():  # Missing class aliases (Bend)
    global CLASS_MAP

    def subclasses_recursive(cl):
        direct = cl.__subclasses__()
        indirect = []
        for subclass in direct:
            indirect.extend(subclasses_recursive(subclass))
        return frozenset([cl] + direct + indirect)

    cls_list = subclasses_recursive(Element)
    CLASS_MAP = {cls.__name__: cls for cls in cls_list}


def get_class_map():
    return CLASS_MAP


# build_class_map()

CLASS_MAP = {
    k: v for k, v in locals().items() if isinstance(v, type) and issubclass(v, Element)
}
