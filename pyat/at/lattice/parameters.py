from __future__ import annotations

__all__ = ["ParamBase", "Param", "ParamArray", "AttributeArray"]

from typing import Any, Generic
import abc
from collections.abc import Callable
from operator import add, sub, mul, truediv, pos, neg
import numpy as np
from .variables import Number, VariableBase, _nop


class _Evaluator(Generic[Number], abc.ABC):
    @abc.abstractmethod
    def __call__(self) -> Number: ...


class _Constant(_Evaluator[Number]):
    __slots__ = "value"

    def __init__(self, value):
        if not isinstance(value, (int, float)):
            raise TypeError("The parameter value must be a scalar")
        self.value = value

    def __call__(self) -> Number:
        return self.value


class _BinaryOperator(_Evaluator[Number]):
    __slots__ = ["operator", "left_operand", "right_operand"]

    @staticmethod
    def _convert_to_evaluator(value):
        if isinstance(value, (int, float)):
            return _Constant(value)
        elif isinstance(value, VariableBase):
            return value
        raise TypeError(f"Parameter operation not defined for type {type(value)}")

    def __init__(self, operator, left, right):
        self.operator = operator
        self.right_operand = self._convert_to_evaluator(right)
        self.left_operand = self._convert_to_evaluator(left)

    def __call__(self) -> Number:
        return self.operator(self.left_operand.value, self.right_operand.value)


class _UnaryOperator(_Evaluator[Number]):
    __slots__ = ["operator", "operand"]

    def __init__(self, operator, operand):
        self.operator = operator
        self.operand = operand

    def __call__(self) -> Number:
        return self.operator(self.operand.value)


class ParamBase(VariableBase[Number]):
    """Read-only base class for parameters

    It is used for computed parameters and should not be instantiated
    otherwise. See :py:class:`.Variable` for a description of inherited
    methods
    """

    COUNTER_PREFIX = "calc"

    _counter = 0
    _evaluator: _Evaluator[Number]
    _conversion: Callable[[Any], Number]

    def __init__(
        self,
        evaluator: _Evaluator[Number],
        *,
        name: str = "",
        conversion: Callable[[Any], Number] = _nop,
        bounds: tuple[Number, Number] = (-np.inf, np.inf),
        delta: Number = 1.0,
    ):
        """

        Args:
            evaluator:  Evaluator function
            name:       Name of the parameter
            conversion: data conversion function
            bounds:     Lower and upper bounds of the parameter value
            delta:      Initial variation step
        """
        if not isinstance(evaluator, _Evaluator):
            raise TypeError("'Evaluate' must be an _Evaluate object")
        self._evaluator = evaluator
        self._conversion = conversion
        super().__init__(name=name, bounds=bounds, delta=delta)

    def _getfun(self, **kwargs) -> Number:
        return self._conversion(self._evaluator())

    @property
    def _safe_value(self):
        return self._getfun()

    def set_conversion(self, conversion: Callable[[Any], Number]):
        """Set the data type. Called when a parameter is assigned to an
        :py:class:`.Element` attribute"""
        if conversion is not self._conversion:
            if self._conversion is _nop:
                self._conversion = conversion
            else:
                raise ValueError("Cannot change the data type of the parameter")

    def __add__(self, other):
        fun = _BinaryOperator(add, self, other)
        return ParamBase(fun)

    __radd__ = __add__

    def __pos__(self):
        return ParamBase(_UnaryOperator(pos, self))

    def __neg__(self):
        return ParamBase(_UnaryOperator(neg, self))

    def __sub__(self, other):
        fun = _BinaryOperator(sub, self, other)
        return ParamBase(fun)

    def __rsub__(self, other):
        fun = _BinaryOperator(sub, other, self)
        return ParamBase(fun)

    def __mul__(self, other):
        fun = _BinaryOperator(mul, self, other)
        return ParamBase(fun)

    __rmul__ = __mul__

    def __truediv__(self, other):
        fun = _BinaryOperator(truediv, self, other)
        return ParamBase(fun)

    def __rtruediv__(self, other):
        fun = _BinaryOperator(truediv, other, self)
        return ParamBase(fun)

    def __float__(self):
        return float(self._safe_value)

    def __int__(self):
        return int(self._safe_value)

    def __str__(self):
        return f"{self.__class__.__name__}({self._safe_value}, name={self.name!r})"

    def __repr__(self):
        return repr(self._safe_value)


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
            _Constant(value),
            name=name,
            conversion=conversion,
            bounds=bounds,
            delta=delta,
        )
        self._history.append(self._evaluator())

    def _getfun(self, **kwargs) -> Number:
        return self._evaluator()

    def _setfun(self, value, ring=None):
        self._evaluator = _Constant(self._conversion(value))

    def set_conversion(self, conversion: Callable[[Any], Number]):
        oldv = self._evaluator()
        super().set_conversion(conversion)
        self._evaluator = _Constant(conversion(oldv))


class _SafeArray(np.ndarray):
    """Subclass of ndarray which forbids setting parameters as items"""

    def __setitem__(self, key, value):
        if isinstance(value, ParamBase):
            raise TypeError("Cannot set a parameter into an array")
        super().__setitem__(key, value)


def AttributeArray(value, shape=(-1,), dtype=float):
    v = np.asfortranarray(value).reshape(shape, order="F")
    if v.dtype == np.dtype("O"):
        return ParamArray(v, shape=shape, dtype=dtype)
    else:
        return v.astype(dtype, copy=False).view(_SafeArray)


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
