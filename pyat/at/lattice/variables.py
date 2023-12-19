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
import numpy as np
import abc
from numbers import Number
from operator import add, sub, mul, truediv, pos, neg
from collections import deque
from collections.abc import Iterable, Sequence, Callable
from typing import Any

__all__ = [
    "VariableBase",
    "CustomVariable",
    "ParamBase",
    "Param",
    "ParamArray",
    "VariableList",
]


def _nop(value):
    return value


def _default_array(value):
    return np.require(value, dtype=float, requirements=["F", "A"])


class _Evaluate(abc.ABC):
    @abc.abstractmethod
    def __call__(self): ...


class _Scalar(_Evaluate):
    __slots__ = "value"

    def __init__(self, value):
        if not isinstance(value, Number):
            raise TypeError("The parameter value must be a scalar")
        self.value = value

    def __call__(self):
        return self.value


class _BinaryOp(_Evaluate):
    __slots__ = ["oper", "left", "right"]

    @staticmethod
    def _set_type(value):
        if isinstance(value, Number):
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


class VariableBase(abc.ABC):
    """A Variable abstract base class

    Derived classes must implement the :py:meth:`~VariableBase._getfun` and
    :py:meth:`~VariableBase._getfun` methods
    """

    _counter = 0
    _prefix = "var"

    def __init__(
        self,
        *,
        name: str = "",
        bounds: tuple[Number, Number] = (-np.inf, np.inf),
        delta: Number = 1.0,
        history_length: int = None,
        ring=None,
    ):
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
        self.name = self._setname(name)  #: Variable name
        self.bounds = bounds  #: Variable bounds
        self.delta = delta  #: Increment step
        #: Maximum length of the history buffer. :py:obj:`None` means infinite
        self.history_length = history_length
        self._initial = np.nan
        self._history = deque([], self.history_length)
        try:
            self.get(ring=ring, initial=True)
        except ValueError:
            pass

    @classmethod
    def _setname(cls, name):
        cls._counter += 1
        if name:
            return name
        else:
            return f"{cls._prefix}{cls._counter}"

    # noinspection PyUnusedLocal
    def _setfun(self, value: Number, ring=None):
        classname = self.__class__.__name__
        raise TypeError(f"{classname!r} is read-only")

    @abc.abstractmethod
    def _getfun(self, ring=None) -> Number: ...

    @property
    def history(self) -> list[Number]:
        """History of the values of the variable"""
        return list(self._history)

    @property
    def initial_value(self) -> Number:
        """Initial value of the variable"""
        if not np.isnan(self._initial):
            return self._initial
        else:
            raise IndexError(f"{self.name}: No value has been set yet")

    @property
    def last_value(self) -> Number:
        """Last value of the variable"""
        try:
            return self._history[-1]
        except IndexError as exc:
            exc.args = (f"{self.name}: No value has been set yet",)
            raise

    @property
    def previous_value(self) -> Number:
        """Value before the last one"""
        try:
            return self._history[-2]
        except IndexError as exc:
            exc.args = (f"{self.name}: history too short",)
            raise

    def set(self, value: Number, ring=None) -> None:
        """Set the variable value

        Args:
            value:  New value to be applied on the variable
            ring:   Depending on the variable type, a :py:class:`.Lattice` argument
              may be necessary to set the variable.
        """
        if value < self.bounds[0] or value > self.bounds[1]:
            raise ValueError(f"set value must be in {self.bounds}")
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
            if value < self.bounds[0] or value > self.bounds[1]:
                raise ValueError(f"value out of {self.bounds}")
        return value

    value = property(get, set, doc="Actual value")

    @property
    def _safe_value(self):
        try:
            v = self._history[-1]
        except IndexError:
            v = np.nan
        return v

    def set_previous(self, ring=None) -> None:
        """Reset to the value before the last one

        Args:
            ring:   Depending on the variable type, a :py:class:`.Lattice` argument
              may be necessary to set the variable.
        """
        if len(self._history) >= 2:
            self._history.pop()  # Remove the last value
            value = self._history.pop()  # retrieve the previous value
            self.set(value, ring=ring)
        else:
            raise IndexError(f"{self.name}: history too short",)

    def reset(self, ring=None) -> None:
        """Reset to the initial value and clear the history buffer

        Args:
            ring:   Depending on the variable type, a :py:class:`.Lattice` argument
              may be necessary to reset the variable.
        """
        iniv = self._initial
        if not np.isnan(iniv):
            self._history = deque([], self.history_length)
            self.set(iniv, ring=ring)
        else:
            raise IndexError(f"reset {self.name}: No value has been set yet")

    def increment(self, incr: Number, ring=None) -> None:
        """Increment the variable value

        Args:
            incr:   Increment value
            ring:   Depending on the variable type, a :py:class:`.Lattice` argument
              may be necessary to increment the variable.
        """
        if self._initial is None:
            self.get(ring=ring, initial=True)
        self.set(self.last_value + incr, ring=ring)

    def _step(self, step: Number, ring=None) -> None:
        if self._initial is None:
            self.get(ring=ring, initial=True)
        self.set(self._initial + step, ring=ring)

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

        return "{:>12s}{: 16e}{: 16e}{: 16e}".format(
            self.name, vini, vnow, (vnow - vini)
        )

    def status(self):
        """Return a string describing the current status of the variable

        Returns:
            status: Variable description
        """
        return "\n".join((self._header(), self._line()))

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


class CustomVariable(VariableBase):
    r"""A Variable with user-defined get and set functions

    This is a convenience function allowing user-defined *get* and *set*
    functions. But subclassing :py:class:`.Variable` should always be preferred
    for clarity and efficiency.

    """

    def __init__(
        self,
        setfun: Callable,
        getfun: Callable,
        *args,
        name: str = "",
        bounds: tuple[Number, Number] = (-np.inf, np.inf),
        delta: Number = 1.0,
        history_length: int = None,
        ring=None,
        **kwargs,
    ):
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


class ParamBase(VariableBase):
    """Read-only base class for parameters

    It is used for computed parameters, and should not be instantiated
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
        bounds: tuple[float, float] = (-np.inf, np.inf),
        delta: float = 1.0,
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


class VariableList(list):
    """Container for Variable objects

    :py:class:`VariableList` supports all :py:class:`list` methods, like
    appending, insertion or concatenation with the "+" operator.
    """

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
    def deltas(self) -> Sequence[Number]:
        """delta values of the variables"""
        return np.array([var.delta for var in self])
