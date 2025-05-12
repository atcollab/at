from __future__ import annotations

__all__ = ["Param", "ParamArray", "AttributeArray"]

from typing import Any
from collections.abc import Callable
import numpy as np
from .variables import Number, ParamBase, _Constant, _nop


class Param(ParamBase[Number]):
    """Standalone scalar parameter

    See :py:class:`.Variable` for a description of inherited methods
    """

    COUNTER_PREFIX = "param"

    _counter = 0

    def __init__(
        self,
        value: Number,
        *,
        name: str = "",
        conversion: Callable[[Any], Number] = _nop,
        bounds: tuple[Number, Number] = (-np.inf, np.inf),
        delta: Number = 1.0,
    ):
        """
        Args:
            value:      Initial value of the parameter
            name:       Name of the parameter
            conversion: data conversion function
            bounds:     Lower and upper bounds of the parameter value
            delta:      Initial variation step
        """
        super().__init__(
            _Constant(conversion(value)),
            name=name,
            conversion=conversion,
            bounds=bounds,
            delta=delta,
        )
        self._history.append(self._evaluator())

    def _getfun(self, ring=None):
        return self._evaluator()

    def _setfun(self, value, ring=None):
        self._evaluator = _Constant(self._conversion(value))

    def set_conversion(self, conversion: Callable[[Number], Number]):
        oldv = self._evaluator()
        super(Param, self).set_conversion(conversion)
        self._evaluate = _Constant(conversion(oldv))


class _SafeArray(np.ndarray):
    """Subclass of ndarray which forbids setting parameters as items"""

    def __setitem__(self, key, value):
        if isinstance(value, ParamBase):
            raise TypeError("Cannot set a parameter into an array")
        super().__setitem__(key, value)


def AttributeArray(value, shape=(-1,), dtype=float):
    v = np.asfortranarray(value).reshape(shape, order="F")
    if v.dtype == "O":
        return ParamArray(v, shape=shape, dtype=dtype)
    else:
        return v.astype(dtype, copy=False).view(_SafeArray)


class _PArray(np.ndarray):
    """Subclass of ndarray which reports to its parent ParamArray"""

    # This is the array obtained with an element get_attribute.
    # It is also the one used when setting an item of an array attribute.

    def __new__(cls, value, dtype=float):
        obj = np.array(value, dtype=dtype, order="F").view(cls)
        obj._parent = value
        return obj

    def __array_finalize__(self, obj):
        self._parent = getattr(obj, "_parent", None)

    def __setitem__(self, key, value):
        # report the value to the parent
        super().__setitem__(key, value)
        if self._parent is not None:
            self._parent[key] = value

    def __repr__(self):
        # Simulate a standard ndarray
        return repr(self.view(np.ndarray))


class ParamArray(np.ndarray):
    """Simulate a numpy array where items may be parametrised"""

    def __new__(cls, value, shape=(-1,), dtype=float):
        obj = np.asfortranarray(value, dtype="O").reshape(shape).view(cls)
        obj._value = _PArray(obj, dtype=dtype)
        return obj

    def __array_finalize__(self, obj):
        val = getattr(obj, "_value", None)
        if val is not None:
            self._value = _PArray(self, dtype=val.dtype)

    @property
    def value(self):
        self._value[:] = self
        return self._value

    def __repr__(self):
        return repr(self.value)

    def __str__(self):
        it = np.nditer(self, flags=["refs_ok"], order="C")
        contents = " ".join([str(el) for el in it])
        return f"[{contents}]"
