r"""
Definition of :py:class:`Variable <.VariableBase>` objects used in matching and
response matrices.

See :ref:`example-notebooks` for examples of matching and response matrices.

Each :py:class:`Variable <.VariableBase>` has a scalar value.

.. rubric:: Class hierarchy

:py:class:`VariableBase`\ (name, bounds, delta)

- :py:class:`~.lattice_variables.ElementVariable`\ (elements, attrname, index, ...)
- :py:class:`~.lattice_variables.RefptsVariable`\ (refpts, attrname, index, ...)
- :py:class:`CustomVariable`\ (setfun, getfun, ...)

.. rubric:: VariableBase methods

:py:class:`VariableBase` provides the following methods:

- :py:meth:`~VariableBase.get`
- :py:meth:`~VariableBase.set`
- :py:meth:`~VariableBase.set_previous`
- :py:meth:`~VariableBase.reset`
- :py:meth:`~VariableBase.increment`
- :py:meth:`~VariableBase.step_up`
- :py:meth:`~VariableBase.step_down`

.. rubric:: VariableBase properties

:py:class:`.VariableBase` provides the following properties:

- :py:attr:`~VariableBase.initial_value`
- :py:attr:`~VariableBase.last_value`
- :py:attr:`~VariableBase.previous_value`
- :py:attr:`~VariableBase.history`

The :py:class:`VariableBase` abstract class may be used as a base class to define
custom variables (see examples). Typically, this consist in overloading the abstract
methods *_setfun* and *_getfun*

.. rubric:: Examples

Write a subclass of :py:class:`VariableBase` which varies two drift lengths so
that their sum is constant:

.. code-block:: python

    class ElementShifter(at.VariableBase):
        '''Varies the length of the elements identified by *ref1* and *ref2*
        keeping the sum of their lengths equal to *total_length*.

        If *total_length* is None, it is set to the initial total length
        '''

        def __init__(self, drift1, drift2, total_length=None, **kwargs):
            # store the 2 variable elements
            self.drift1 = drift1
            self.drift2 = drift2
            # store the initial total length
            if total_length is None:
                total_length = drift1.Length + drift2.Length
            self.length = total_length
            super().__init__(bounds=(0.0, total_length), **kwargs)

        def _setfun(self, value, **kwargs):
            self.drift1.Length = value
            self.drift2.Length = self.length - value

        def _getfun(self, **kwargs):
            return self.drift1.Length

And create a variable varying the length of drifts *DR_01* and *DR_01* and
keeping their sum constant:

.. code-block:: python

    drift1 = hmba_lattice["DR_01"]
    drift2 = hmba_lattice["DR_02"]
    var2 = ElementShifter(drift1, drift2, name="DR_01")

"""

from __future__ import annotations

__all__ = [
    "VariableBase",
    "ParamBase",
    "CustomVariable",
    "VariableList",
]

import abc
from operator import add, sub, mul, truediv, neg
from collections import deque
from collections.abc import Iterable, Sequence, Callable
from typing import Any, Generic, TypeVar

import numpy as np
import numpy.typing as npt

from .parser import ParamDef, _nop

# Define a type variable for numeric types
Number = TypeVar("Number", int, float)


def _name(obj, priority: int):
    """Return the parenthesised name of the object"""
    if isinstance(obj, ParamBase):
        return obj.name if obj._priority >= priority else f"({obj.name})"
    else:
        return str(obj)


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
        elif isinstance(value, VariableBase):
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


class VariableBase(Generic[Number], abc.ABC):
    """A Variable abstract base class

    Derived classes must implement the :py:meth:`~VariableBase._getfun` and
    :py:meth:`~VariableBase._setfun` methods
    """

    # Class constants
    DEFAULT_DELTA = 1.0
    COUNTER_PREFIX = "var"

    _counter = 0

    def __init__(
        self,
        *,
        name: str = "",
        bounds: tuple[Number, Number] | None = None,
        delta: Number = DEFAULT_DELTA,
        history_length: int | None = None,
        ring=None,
        **kwargs,
    ) -> None:
        """
        Parameters:
            name:       Name of the Variable
            bounds:     Lower and upper bounds of the variable value
            delta:      Initial variation step
            history_length: Maximum length of the history buffer. :py:obj:`None`
              means infinite
            ring:       provided to an attempt to get the initial value of the
              variable
        """
        self.name: str = self._generate_name(name)  #: Variable name
        if bounds is None:
            bounds = (None, None)
        self.bounds: tuple[Number | None, Number | None] = bounds  #: Variable bounds
        self.delta: Number = delta  #: Increment step
        #: Maximum length of the history buffer. :py:obj:`None` means infinite
        self.history_length = history_length
        self._initial = np.nan
        self._history: deque[Number] = deque([], self.history_length)
        try:
            self.get(ring=ring, initial=True)
        except ValueError:
            pass
        super().__init__(**kwargs)

    @classmethod
    def _generate_name(cls, name: str) -> str:
        """Generate unique name for variable"""
        cls._counter += 1
        return name if name else f"{cls.COUNTER_PREFIX}{cls._counter}"

    def _check_bounds(self, value: Number) -> None:
        """Verify value is within bounds"""
        min_val, max_val = self.bounds
        if min_val is not None and value < min_val:
            raise ValueError(f"Value {value} must be larger or equal to {min_val}")
        if max_val is not None and value > max_val:
            raise ValueError(f"Value {value} must be smaller or equal to {max_val}")

    # noinspection PyUnusedLocal
    def _setfun(self, value: Number, ring=None) -> None:
        classname = self.__class__.__name__
        raise TypeError(f"{classname!r} is read-only")

    @abc.abstractmethod
    def _getfun(self, **kwargs) -> Number: ...

    @property
    def history(self) -> list[Number]:
        """History of the values of the variable"""
        return list(self._history)

    @property
    def initial_value(self) -> Number:
        """Initial value of the variable

        Raises:
            IndexError: If no initial value has been set yet
        """
        if not np.isnan(self._initial):
            return self._initial
        else:
            raise IndexError(f"{self.name}: No initial value has been set yet")

    @property
    def last_value(self) -> Number:
        """Last value of the variable

        Raises:
            IndexError: If no value has been set yet
        """
        try:
            return self._history[-1]
        except IndexError as exc:
            exc.args = (f"{self.name}: No value has been set yet",)
            raise

    @property
    def previous_value(self) -> Number:
        """Value before the last one

        Raises:
            IndexError: If there are fewer than 2 values in the history
        """
        try:
            return self._history[-2]
        except IndexError as exc:
            exc.args = (f"{self.name}: history too short (need at least 2 values)",)
            raise

    def set(self, value: Number, ring=None) -> None:
        """Set the variable value

        Args:
            value:  New value to be applied on the variable
            ring:   Depending on the variable type, a :py:class:`.Lattice` argument
              may be necessary to set the variable.
        """
        self._check_bounds(value)
        self._setfun(value, ring=ring)
        if np.isnan(self._initial):
            self._initial = value
        self._history.append(value)

    def get(
        self, ring=None, *, initial: bool = False, check_bounds: bool = False
    ) -> Number:
        """Get the actual variable value

        Args:
            ring:   Depending on the variable type, a :py:class:`.Lattice` argument
              may be necessary to get the variable value.
            initial:    If :py:obj:`True`, clear the history and set the variable
              initial value
            check_bounds: If :py:obj:`True`, raise a ValueError if the value is out
              of bounds

        Returns:
            value:      Value of the variable
        """
        value = self._getfun(ring=ring)
        if initial or np.isnan(self._initial):
            self._initial = value
            self._history = deque([value], self.history_length)
        if check_bounds:
            self._check_bounds(value)
        return value

    value = property(get, set, doc="Actual value")

    @property
    def _safe_value(self):
        try:
            return self._history[-1]
        except IndexError:
            return np.nan

    def set_previous(self, ring=None) -> None:
        """Reset to the value before the last one

        Args:
            ring:   Depending on the variable type, a :py:class:`.Lattice` argument
              may be necessary to set the variable.

        Raises:
            IndexError: If there are fewer than 2 values in the history
        """
        if len(self._history) >= 2:
            self._history.pop()  # Remove the last value
            previous_value = self._history.pop()  # retrieve the previous value
            self.set(previous_value, ring=ring)
        else:
            raise IndexError(f"{self.name}: history too short (need at least 2 values)")

    def reset(self, ring=None) -> None:
        """Reset to the initial value and clear the history buffer

        Args:
            ring:   Depending on the variable type, a :py:class:`.Lattice` argument
              may be necessary to reset the variable.

        Raises:
            IndexError: If no initial value has been set yet
        """
        initial_value = self._initial
        if not np.isnan(initial_value):
            self._history = deque([], self.history_length)
            self.set(initial_value, ring=ring)
        else:
            raise IndexError(
                f"Cannot reset {self.name}: No initial value has been set yet"
            )

    def increment(self, incr: Number, ring=None) -> None:
        """Increment the variable value

        Args:
            incr:   Increment value
            ring:   Depending on the variable type, a :py:class:`.Lattice` argument
              may be necessary to increment the variable.
        """
        try:
            current_value = self.last_value
        except IndexError:
            # If no value has been set yet, get the initial value
            self.get(ring=ring, initial=True)
            current_value = self.last_value

        self.set(current_value + incr, ring=ring)

    def _step(self, step: Number, ring=None) -> None:
        """Apply a step relative to the initial value

        Args:
            step:   Step value to add to the initial value
            ring:   Depending on the variable type, a :py:class:`.Lattice` argument
              may be necessary to set the variable.
        """
        try:
            initial_value = self.initial_value
        except IndexError:
            # If no initial value has been set yet, get it
            self.get(ring=ring, initial=True)
            initial_value = self.initial_value

        self.set(initial_value + step, ring=ring)

    def step_up(self, ring=None) -> None:
        """Set to initial_value + delta

        Args:
            ring:   Depending on the variable type, a :py:class:`.Lattice` argument
              may be necessary to set the variable.
        """
        self._step(self.delta, ring=ring)

    def step_down(self, ring=None) -> None:
        """Set to initial_value - delta

        Args:
            ring:   Depending on the variable type, a :py:class:`.Lattice` argument
              may be necessary to set the variable.
        """
        self._step(-self.delta, ring=ring)

    @staticmethod
    def _header():
        return "\n{:>12s}{:>13s}{:>16s}{:>16s}\n".format(
            "Name", "Initial", "Final ", "Variation"
        )

    def _line(self):
        vnow = self._safe_value
        vini = self._initial
        return f"{self.name:>12s}{vini: 16e}{vnow: 16e}{vnow - vini: 16e}"

    def status(self):
        """Return a string describing the current status of the variable

        Returns:
            status: Variable description
        """
        return "\n".join((self._header(), self._line()))

    def __add__(self, other):
        fun = _BinaryOperator(add, self, other)
        name = "+".join((_name(self, 10), _name(other, 10)))
        return ParamBase(fun, name=name, priority=10)

    __radd__ = __add__

    def __pos__(self):
        return self

    def __neg__(self):
        name = "-" + _name(self, 20)
        return ParamBase(_UnaryOperator(neg, self), name=name, priority=0)

    def __abs__(self):
        name = f"abs({_name(self, 0)})"
        return ParamBase(_UnaryOperator(abs, self), name=name, priority=20)

    def __sub__(self, other):
        fun = _BinaryOperator(sub, self, other)
        name = "-".join((_name(self, 10), _name(other, 10)))
        return ParamBase(fun, name=name, priority=10)

    def __rsub__(self, other):
        fun = _BinaryOperator(sub, other, self)
        name = "-".join((_name(other, 10), _name(self, 10)))
        return ParamBase(fun, name=name, priority=10)

    def __mul__(self, other):
        fun = _BinaryOperator(mul, self, other)
        name = "*".join((_name(self, 20), _name(other, 20)))
        return ParamBase(fun, name=name, priority=20)

    __rmul__ = __mul__

    def __truediv__(self, other):
        fun = _BinaryOperator(truediv, self, other)
        name = "/".join((_name(self, 20), _name(other, 20)))
        return ParamBase(fun, name=name, priority=20)

    def __rtruediv__(self, other):
        fun = _BinaryOperator(truediv, other, self)
        name = "/".join((_name(other, 20), _name(self, 20)))
        return ParamBase(fun, name=name, priority=20)

    def __float__(self):
        return float(self._safe_value)

    def __int__(self):
        return int(self._safe_value)

    def __str__(self):
        return f"{self.__class__.__name__}({self._safe_value}, name={self.name!r})"

    def __repr__(self):
        return repr(self._safe_value)


class ParamBase(VariableBase[Number], ParamDef):
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
        bounds: tuple[Number, Number] | None = None,
        delta: Number = 1.0,
        priority: int = 20,
    ) -> None:
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
        self._priority = priority
        super().__init__(name=name, bounds=bounds, delta=delta)

    def _getfun(self, **kwargs) -> Number:
        return self._conversion(self._evaluator())

    @property
    def _safe_value(self):
        return self._getfun()


class CustomVariable(VariableBase[Number]):
    r"""A Variable with user-defined get and set functions

    This is a convenience function allowing user-defined *get* and *set*
    functions. But subclassing :py:class:`.Variable` should always be preferred
    for clarity and efficiency.

    """

    def __init__(
        self,
        setfun: Callable[..., None],
        getfun: Callable[..., Number],
        *args,
        name: str = "",
        bounds: tuple[Number, Number] | None = None,
        delta: Number = 1.0,
        history_length: int | None = None,
        ring=None,
        **kwargs,
    ) -> None:
        """
        Parameters:
            getfun:     Function for getting the variable value. Called as
              :pycode:`getfun(*args, ring=ring, **kwargs) -> Number`
            setfun:     Function for setting the variable value. Called as
              :pycode:`setfun(value: Number, *args, ring=ring, **kwargs): None`
            name:       Name of the Variable
            bounds:     Lower and upper bounds of the variable value
            delta:      Initial variation step
            *args:      Variable argument list transmitted to both the *getfun*
              and *setfun* functions. Such arguments can always be avoided by
              using :py:func:`~functools.partial` or callable class objects.

        Keyword Args:
            **kwargs:   Keyword arguments transmitted to both the *getfun*
              and *setfun* functions. Such arguments can always be avoided by
              using :py:func:`~functools.partial` or callable class objects.
        """
        self.getfun = getfun
        self.setfun = setfun
        self.args = args
        self.kwargs = kwargs
        super().__init__(
            name=name,
            bounds=bounds,
            delta=delta,
            history_length=history_length,
            ring=ring,
        )

    def _getfun(self, ring=None) -> Number:
        return self.getfun(*self.args, ring=ring, **self.kwargs)

    def _setfun(self, value: Number, ring=None):
        self.setfun(value, *self.args, ring=ring, **self.kwargs)


class VariableList(list):
    """Container for Variable objects

    :py:class:`VariableList` supports all :py:class:`list` methods, like
    appending, insertion or concatenation with the "+" operator.
    """

    def __getitem__(self, index):
        if isinstance(index, slice):
            return VariableList(super().__getitem__(index))
        else:
            return super().__getitem__(index)

    def get(self, ring=None, **kwargs) -> Sequence[float]:
        r"""Get the current values of Variables

        Args:
            ring:   Depending on the variable type, a :py:class:`.Lattice` argument
              may be necessary to set the variable.

        Keyword Args:
            initial:    If :py:obj:`True`, set the Variables'
              initial value
            check_bounds: If :py:obj:`True`, raise a ValueError if the value is out
              of bounds

        Returns:
            values:     1D array of values of all variables
        """
        return np.array([var.get(ring=ring, **kwargs) for var in self])

    def set(self, values: Iterable[float], ring=None) -> None:
        r"""Set the values of Variables

        Args:
            values:     Iterable of values
            ring:   Depending on the variable type, a :py:class:`.Lattice` argument
              may be necessary to set the variable.
        """
        for var, val in zip(self, values):
            var.set(val, ring=ring)

    def increment(self, increment: Iterable[float], ring=None) -> None:
        r"""Increment the values of Variables

        Args:
            increment:  Iterable of values
            ring:   Depending on the variable type, a :py:class:`.Lattice` argument
              may be necessary to increment the variable.
        """
        for var, incr in zip(self, increment):
            var.increment(incr, ring=ring)

    # noinspection PyProtectedMember
    def status(self, **kwargs) -> str:
        """String description of the variables"""
        values = "\n".join(var._line(**kwargs) for var in self)
        return "\n".join((VariableBase._header(), values))

    def __str__(self) -> str:
        return self.status()

    @property
    def deltas(self) -> npt.NDArray[Number]:
        """delta values of the variables"""
        return np.array([var.delta for var in self])

    @deltas.setter
    def deltas(self, value: Number | Sequence[Number]) -> None:
        deltas = np.broadcast_to(value, len(self))
        for var, delta in zip(self, deltas):
            var.delta = delta
