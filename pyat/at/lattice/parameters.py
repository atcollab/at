from __future__ import annotations

from typing import Any, Union
import abc
from collections.abc import Callable
from operator import add, sub, mul, truediv, pos, neg
import numpy as np
from .elements import Element
from .variables import VariableBase

Number = Union[int, float]


def _nop(value):
    return value


class _Evaluate(abc.ABC):
    @abc.abstractmethod
    def __call__(self): ...


class _Scalar(_Evaluate):
    __slots__ = "value"

    def __init__(self, value):
        if not isinstance(value, (int, float)):
            raise TypeError("The parameter value must be a scalar")
        self.value = value

    def __call__(self):
        return self.value


class _BinaryOp(_Evaluate):
    __slots__ = ["oper", "left", "right"]

    @staticmethod
    def _set_type(value):
        if isinstance(value, (int, float)):
            return _Scalar(value)
        elif isinstance(value, VariableBase):
            return value
        else:
            msg = "Param Operation not defined for type {0}".format(type(value))
            raise TypeError(msg)

    def __init__(self, oper, left, right):
        self.oper = oper
        self.right = self._set_type(right)
        self.left = self._set_type(left)

    def __call__(self):
        return self.oper(self.left.value, self.right.value)


class _UnaryOp(_Evaluate):
    __slots__ = ["oper", "param"]

    def __init__(self, oper, param):
        self.oper = oper
        self.param = param

    def __call__(self):
        return self.oper(self.param.value)


class ParamBase(VariableBase):
    """Read-only base class for parameters

    It is used for computed parameters and should not be instantiated
    otherwise. See :py:class:`.Variable` for a description of inherited
    methods
    """

    _counter = 0
    _prefix = "calc"

    def __init__(
        self,
        evaluate: _Evaluate,
        *,
        name: str = "",
        conversion: Callable[[Any], Number] = _nop,
        bounds: tuple[Number, Number] = (-np.inf, np.inf),
        delta: Number = 1.0,
    ):
        """

        Args:
            evaluate:   Evaluator function
            name:       Name of the parameter
            conversion: data conversion function
            bounds:     Lower and upper bounds of the parameter value
            delta:      Initial variation step
        """
        if not isinstance(evaluate, _Evaluate):
            raise TypeError("'Evaluate' must be an _Evaluate object")
        self._evaluate = evaluate
        self._conversion = conversion
        super(ParamBase, self).__init__(name=name, bounds=bounds, delta=delta)

    def _getfun(self, **kwargs):
        return self._conversion(self._evaluate())

    @property
    def _safe_value(self):
        return self._getfun()

    def set_conversion(self, conversion: Callable[[Number], Number]):
        """Set the data type. Called when a parameter is assigned to an
        :py:class:`.Element` attribute"""
        if conversion is not self._conversion:
            if self._conversion is _nop:
                self._conversion = conversion
            else:
                raise ValueError("Cannot change the data type of the parameter")

    def __add__(self, other):
        fun = _BinaryOp(add, self, other)
        return ParamBase(fun)

    def __radd__(self, other):
        return self.__add__(other)

    def __pos__(self):
        return ParamBase(_UnaryOp(pos, self))

    def __neg__(self):
        return ParamBase(_UnaryOp(neg, self))

    def __sub__(self, other):
        fun = _BinaryOp(sub, self, other)
        return ParamBase(fun)

    def __rsub__(self, other):
        fun = _BinaryOp(sub, other, self)
        return ParamBase(fun)

    def __mul__(self, other):
        fun = _BinaryOp(mul, self, other)
        return ParamBase(fun)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        fun = _BinaryOp(truediv, self, other)
        return ParamBase(fun)

    def __rtruediv__(self, other):
        fun = _BinaryOp(add, other, self)
        return ParamBase(fun)

    def __float__(self):
        return float(self._safe_value)

    def __int__(self):
        return int(self._safe_value)

    def __str__(self):
        return f"{self.__class__.__name__}({self._safe_value}, name={self.name!r})"

    def __repr__(self):
        return repr(self._safe_value)


class Param(ParamBase):
    """Standalone scalar parameter

    See :py:class:`.Variable` for a description of inherited methods
    """

    _counter = 0
    _prefix = "param"

    def __init__(
        self,
        value: Number,
        *,
        name: str = "",
        conversion: Callable[[Number], Number] = _nop,
        bounds: tuple[float, float] = (-np.inf, np.inf),
        delta: float = 1.0,
    ):
        """
        Args:
            value:      Initial value of the parameter
            name:       Name of the parameter
            conversion: data conversion function
            bounds:     Lower and upper bounds of the parameter value
            delta:      Initial variation step
        """
        super(Param, self).__init__(
            _Scalar(value), name=name, conversion=conversion, bounds=bounds, delta=delta
        )
        self._history.append(self._evaluate())

    def _getfun(self, ring=None):
        return self._evaluate()

    def _setfun(self, value, ring=None):
        self._evaluate = _Scalar(self._conversion(value))

    def set_conversion(self, conversion: Callable[[Number], Number]):
        oldv = self._evaluate()
        super(Param, self).set_conversion(conversion)
        self._evaluate = _Scalar(conversion(oldv))


class _PArray(np.ndarray):
    """Subclass of ndarray which reports to its parent ParamArray"""

    # This is the array obtained with an element get_attribute.
    # It is also the one used when setting an item of an array attribute.

    def __new__(cls, value, dtype=np.float64):
        obj = np.array(value, dtype=dtype, order="F").view(cls)
        obj._parent = value
        return obj

    def __array_finalize__(self, obj):
        self._parent = getattr(obj, "_parent", None)

    def __setitem__(self, key, value):
        super().__setitem__(key, value)
        if self._parent is not None:
            self._parent[key] = value

    def __repr__(self):
        # Simulate a standard ndarray
        return repr(self.view(np.ndarray))


class ParamArray(np.ndarray):
    """Simulate a numpy array where items may be parametrised"""

    def __new__(cls, value, shape=(-1,), dtype=np.float64):
        obj = np.asfortranarray(value, dtype=object).reshape(shape).view(cls)
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


def __setattr__(self, key, value):
    try:
        if isinstance(value, ParamBase):
            value.set_conversion(self._conversions.get(key, _nop))
        else:
            value = self._conversions.get(key, _nop)(value)
    except Exception as exc:
        exc.args = ('In element {0}, parameter {1}: {2}'.format(
            self.FamName, key, exc),)
        raise
    else:
        super(Element, self).__setattr__(key, value)


def __getattribute__(self, key):
    attr = super(Element, self).__getattribute__(key)
    if isinstance(attr, (ParamBase, ParamArray)):
        return attr.value
    else:
        return attr
