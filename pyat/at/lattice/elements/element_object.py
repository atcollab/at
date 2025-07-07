"""Base Element object"""

from __future__ import annotations

__all__ = ["Element"]

import re
from collections.abc import Generator
from copy import copy, deepcopy
from typing import Any

import numpy as np

from .conversions import _array, _array66, _int, _float
# noinspection PyProtectedMember
from ..variables import _nop


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

    @classmethod
    def get_subclasses(cls) -> Generator[type[Element], None, None]:
        """Iterator over the subclasses of this element

        Because of multiple inheritance, some classes may appear several times
        """
        for subclass in cls.__subclasses__():
            yield from subclass.get_subclasses()
        yield cls

    def equals(self, other) -> bool:
        """Whether an element is equivalent to another.

        This implementation was found to be too slow for the generic
        __eq__ method when comparing lattices.
        """
        return repr(self) == repr(other)

    def divide(self, frac) -> list[Element]:
        # noinspection PyUnresolvedReferences
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
