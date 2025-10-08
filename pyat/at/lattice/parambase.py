from __future__ import annotations

__all__ = ["Combiner", "ParamDef", "ParamBase"]

import abc
from operator import add, sub, mul, truediv, neg
from collections.abc import Callable
from typing import Any, Generic, TypeVar

# Define a type variable for numeric types
Number = TypeVar("Number", int, float)


def _nop(value: Any) -> Any:
    """No-operation function that returns its input unchanged

    This function is used as a default conversion function in parameter classes.
    """
    return value


class _Evaluator(Generic[Number], abc.ABC):
    """Abstract base class for evaluators

    An evaluator is a callable object that returns a scalar value.
    """

    @abc.abstractmethod
    def __call__(self) -> Number:
        """Evaluate and return the value

        Returns:
            The evaluated value
        """
        ...


class _Constant(_Evaluator[Number]):
    """An evaluator that always returns a constant value"""

    __slots__ = "value"

    def __init__(self, value: Number):
        """Initialize a constant evaluator

        Args:
            value: The constant value to return

        Raises:
            TypeError: If the value is not a scalar (int or float)
        """
        if not isinstance(value, (int, float)):
            raise TypeError("The parameter value must be a scalar")
        self.value: Number = value

    def __call__(self) -> Number:
        return self.value


class _BinaryOperator(_Evaluator[Number]):
    __slots__ = ["operator", "left_operand", "right_operand"]

    @staticmethod
    def _convert_to_evaluator(value):
        """Convert a value to an evaluator"""
        if isinstance(value, (int, float)):
            return _Constant(value)
        elif isinstance(value, Combiner):
            return value
        raise TypeError(f"Parameter operation not defined for type {type(value)}")

    def __init__(self, operator, left, right) -> None:
        """Initialize a binary operator

        Args:
            operator: The operator function to apply
            left: The left operand of the operator
            right: The right operand of the operator
        """
        self.operator = operator
        self.right_operand = self._convert_to_evaluator(right)
        self.left_operand = self._convert_to_evaluator(left)

    def __call__(self) -> Number:
        return self.operator(self.left_operand.value, self.right_operand.value)


class _UnaryOperator(_Evaluator[Number]):
    __slots__ = ["operator", "operand"]

    def __init__(self, operator, operand) -> None:
        """Initialize a unary operator

        Args:
            operator: The operator function to apply
            operand: The operand to apply the operator to
        """
        self.operator = operator
        self.operand = operand

    def __call__(self) -> Number:
        return self.operator(self.operand.value)


class Combiner(abc.ABC):
    """Abstract base class for arithmetic combinations of parameters"""
    def __init__(self, name: str, **kwargs):
        self.name = name
        super().__init__(**kwargs)

    @staticmethod
    def _nm(obj, priority: int):
        """Return the parenthesised name of the object"""
        if isinstance(obj, ParamBase):
            return obj.name if obj._priority >= priority else f"({obj.name})"
        else:
            return str(obj)

    @property
    @abc.abstractmethod
    def _safe_value(self): ...

    def __add__(self, other):
        op = _BinaryOperator(add, self, other)
        name = "+".join((self._nm(self, 10), self._nm(other, 10)))
        return ParamBase(evaluator=op, name=name, priority=10)

    __radd__ = __add__

    def __pos__(self):
        return self

    def __neg__(self):
        name = "-" + self._nm(self, 20)
        op = _UnaryOperator(neg, self)
        return ParamBase(evaluator=op, name=name, priority=0)

    def __abs__(self):
        name = f"abs({self._nm(self, 0)})"
        op = _UnaryOperator(abs, self)
        return ParamBase(evaluator=op, name=name, priority=20)

    def __sub__(self, other):
        op = _BinaryOperator(sub, self, other)
        name = "-".join((self._nm(self, 10), self._nm(other, 10)))
        return ParamBase(evaluator=op, name=name, priority=10)

    def __rsub__(self, other):
        op = _BinaryOperator(sub, other, self)
        name = "-".join((self._nm(other, 10), self._nm(self, 10)))
        return ParamBase(evaluator=op, name=name, priority=10)

    def __mul__(self, other):
        op = _BinaryOperator(mul, self, other)
        name = "*".join((self._nm(self, 20), self._nm(other, 20)))
        return ParamBase(evaluator=op, name=name, priority=20)

    __rmul__ = __mul__

    def __truediv__(self, other):
        op = _BinaryOperator(truediv, self, other)
        name = "/".join((self._nm(self, 20), self._nm(other, 20)))
        return ParamBase(evaluator=op, name=name, priority=20)

    def __rtruediv__(self, other):
        op = _BinaryOperator(truediv, other, self)
        name = "/".join((self._nm(other, 20), self._nm(self, 20)))
        return ParamBase(evaluator=op, name=name, priority=20)

    def __float__(self):
        return float(self._safe_value)

    def __int__(self):
        return int(self._safe_value)


class ParamDef(abc.ABC):
    """Abstract base class for parameter definitions

    This class defines the interface for parameter objects that can be used
    as element attributes. It provides a *value* property and a method for converting
    values to the appropriate type.
    """

    def __init__(self, *, conversion: Callable[[Any], Any] | None = None):
        """
        Args:
            conversion: Function to convert values to the appropriate type
        """
        self._conversion = _nop if conversion is None else conversion

    def __copy__(self):
        # Parameters are not copied
        return self

    def __deepcopy__(self, memo):
        # Parameters are not deep-copied
        return self

    def set_conversion(self, conversion: Callable[[Any], Any]) -> None:
        """Set the data type conversion function

        This method is called when a parameter is assigned to an
        :py:class:`.Element` attribute. It can only be set once.

        Args:
            conversion: Function to convert values to the appropriate type

        Raises:
            ValueError: If attempting to change an already set conversion function
        """
        if conversion is not self._conversion:
            if self._conversion is _nop:
                self._conversion = conversion
            else:
                raise ValueError("Cannot change the data type of the parameter")

    @property
    @abc.abstractmethod
    def value(self) -> Any:
        """Current value of the parameter"""
        ...


class ParamBase(Combiner, ParamDef):
    """Read-only base class for parameters

    It is used for computed parameters and should not be instantiated
    otherwise. See :py:class:`.Variable` for a description of inherited
    methods
    """

    _evaluator: _Evaluator

    def __init__(
            self,
            evaluator: _Evaluator,
            *,
            priority: int = 20,
            **kwargs
    ) -> None:
        """

        Args:
            evaluator:  Evaluator function
            name:       Name of the parameter
            conversion: data conversion function
            priority:   priority of the operator
        """
        if not isinstance(evaluator, _Evaluator):
            raise TypeError("'Evaluate' must be an _Evaluate object")
        self._evaluator = evaluator
        self._priority = priority
        super().__init__(**kwargs)

    @property
    def value(self) -> Any:
        return self._conversion(self._evaluator())

    @property
    def _safe_value(self):
        return self._conversion(self._evaluator())

    def __str__(self):
        return self.name

    def __repr__(self):
        return repr(self._safe_value)
