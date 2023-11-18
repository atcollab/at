r"""
Definition of :py:class:`.Variable` objects used in matching and
response matrices.

See :ref:`example-notebooks` for examples of matching and response matrices.

Each :py:class:`Variable` has a scalar value.

.. rubric:: Class hierarchy

:py:class:`Variable`\ (name, bounds, delta)

- :py:class:`.ParamBase`\ (...)

  - :py:class:`.Param`\ (value)
- :py:class:`~.element_variables.ElementVariable`\ (elements, attrname, index, ...)
- :py:class:`~.element_variables.RefptsVariable`\ (refpts, attrname, index, ...)
- :py:class:`CustomVariable`\ (setfun, getfun, ...)

.. rubric:: Variable methods

:py:class:`Variable` provides the following methods:

- :py:meth:`~Variable.get`
- :py:meth:`~Variable.set`
- :py:meth:`~Variable.set_previous`
- :py:meth:`~Variable.set_initial`
- :py:meth:`~Variable.increment`
- :py:meth:`~Variable.step_up`
- :py:meth:`~Variable.step_down`

.. rubric:: Variable properties

:py:class:`.Variable` provides the following properties:

- :py:attr:`~Variable.initial_value`
- :py:attr:`~Variable.last_value`
- :py:attr:`~Variable.previous_value`
- :py:attr:`~Variable.history`

The :py:class:`Variable` abstract class may be used as a base class to define
custom variables (see examples). Typically, this consist in overloading the abstract
methods *_setfun* and *_getfun*

.. rubric:: Examples

Write a subclass of :py:class:`Variable` which varies two drift lengths so
that their sum is constant:

.. code-block:: python

    class ElementShifter(at.Variable):
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
from collections.abc import Iterable, Sequence, Callable
from typing import Any

__all__ = [
    "Variable",
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
    def __call__(self):
        ...


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
        elif isinstance(value, Variable):
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


class Variable(abc.ABC):
    """A :py:class:`Variable` abstract base class

    Derived classes must implement the :py:meth:`_getfun` and
    :py:meth:`_getfun` methods
    """

    _counter = 0
    _prefix = "var"

    def __init__(
        self,
        *,
        name: str = "",
        bounds: tuple[Number, Number] = (-np.inf, np.inf),
        delta: Number = 1.0,
    ):
        """
        Parameters:
            name:       Name of the Variable
            bounds:     Lower and upper bounds of the variable value
            delta:      Initial variation step
        """
        self.name = self._setname(name)
        self.bounds = bounds
        self.delta = delta
        self._history = []

    @classmethod
    def _setname(cls, name):
        cls._counter += 1
        if name:
            return name
        else:
            return f"{cls._prefix}{cls._counter}"

    # noinspection PyUnusedLocal
    def _setfun(self, value: Number, **kwargs):
        classname = self.__class__.__name__
        raise TypeError(f"{classname!r} is read-only")

    @abc.abstractmethod
    def _getfun(self, **kwargs) -> Number:
        ...

    @property
    def history(self) -> list[Number]:
        """History of the values of the variable"""
        return self._history

    @property
    def initial_value(self) -> Number:
        """Initial value of the variable"""
        if len(self._history) > 0:
            return self._history[0]
        else:
            raise IndexError(f"{self.name}: No value has been set yet")

    @property
    def last_value(self) -> Number:
        """Last value of the variable"""
        if len(self._history) > 0:
            return self._history[-1]
        else:
            raise IndexError(f"{self.name}: No value has been set yet")

    @property
    def previous_value(self) -> Number:
        """Value before the last one"""
        if len(self._history) > 1:
            return self._history[-2]
        else:
            raise IndexError(f"{self.name}: history too short")

    def set(self, value: Number, **kwargs) -> None:
        """Set the variable value

        Args:
            value:  New value to be applied on the variable
        """
        if value < self.bounds[0] or value > self.bounds[1]:
            raise ValueError(f"set value must be in {self.bounds}")
        self._setfun(value, **kwargs)
        self._history.append(value)

    def get(self, initial=False, **kwargs) -> Number:
        """Get the actual variable value

        Args:
            initial:    If :py:obj:`True`, set the variable initial value

        Returns:
            value:      Value of the variable
        """
        value = self._getfun(**kwargs)
        if initial:
            self._history = [value]
        return value

    @property
    def value(self):
        return self.get()

    @value.setter
    def value(self, value: Number):
        self.set(value)

    def set_previous(self, **kwargs) -> None:
        """Reset to the value before the last one"""
        if len(self._history) > 1:
            self._history.pop()  # Remove the last value
            prev = self._history.pop()  # retrieve the previous value
            self.set(prev, **kwargs)
        else:
            raise IndexError(f"{self.name}: history too short")

    def set_initial(self, **kwargs) -> None:
        """Reset to the initial value"""
        if len(self._history) > 0:
            iniv = self._history[0]
            self._history = []
            self.set(iniv, **kwargs)
        else:
            raise IndexError(f"{self.name}: No value has been set yet")

    def increment(self, incr: Number, **kwargs) -> None:
        """Increment the variable value

        Args:
            incr:   Increment value
        """
        if len(self._history) == 0:
            self.get(initial=True, **kwargs)
        self.set(self.last_value + incr, **kwargs)

    def _step(self, step: Number, **kwargs) -> None:
        self.set(self.initial_value + step, **kwargs)

    def step_up(self, **kwargs) -> None:
        """Set to initial_value + delta"""
        self._step(self.delta, **kwargs)

    def step_down(self, **kwargs) -> None:
        """Set to initial_value - delta"""
        self._step(-self.delta, **kwargs)

    @staticmethod
    def _header():
        return "\n{:>12s}{:>13s}{:>16s}{:>16s}\n".format(
            "Name", "Initial", "Final ", "Variation"
        )

    def _line(self, ring=None):
        if ring is not None:
            vnow = self.get(ring)
            vini = self._history[0]
        elif len(self._history) > 0:
            vnow = self._history[-1]
            vini = self._history[0]
        else:
            vnow = vini = np.nan

        return "{:>12s}{: 16e}{: 16e}{: 16e}".format(
            self.name, vini, vnow, (vnow - vini)
        )

    def status(self, **kwargs):
        """Return a string describing the current status of the variable

        Returns:
            status: Variable description
        """
        return "\n".join((self._header(), self._line(**kwargs)))

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
        return float(self.value)

    def __int__(self):
        return int(self.value)

    def __str__(self):
        return f"{self.__class__.__name__}({self.value}, name={self.name!r})"

    def __repr__(self):
        return repr(self.value)


class CustomVariable(Variable):
    r"""A :py:class:`.Variable` with user-defined get and set functions

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
        super().__init__(name=name, bounds=bounds, delta=delta)
        self.getfun = getfun
        self.setfun = setfun
        self.args = args
        self.kwargs = kwargs

    def _getfun(self, ring=None) -> Number:
        return self.getfun(*self.args, ring=ring, **self.kwargs)

    def _setfun(self, value: Number, ring=None):
        self.setfun(value, *self.args, ring=ring, **self.kwargs)


class ParamBase(Variable):
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
        super(ParamBase, self).__init__(name=name, bounds=bounds, delta=delta)
        if not isinstance(evaluate, _Evaluate):
            raise TypeError("'Evaluate' must be an _Evaluate object")
        self._evaluate = evaluate
        self._conversion = conversion

    def _getfun(self, **kwargs):
        return self._conversion(self._evaluate())

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
    """Container for :py:class:`Variable` objects

    :py:class:`VariableList` supports all :py:class:`list` methods, like
    appending, insertion or concatenation with the "+" operator.
    """

    def get(self, initial=False, **kwargs) -> Sequence[float]:
        r"""Get the current :py:class:`Variable`\ s' values

        Args:
            initial:    If :py:obj:`True`, set the :py:class:`Variable`\ s'
              initial value

        Returns:
            values:     1D array of values of all variables
        """
        return np.array([var.get(initial=initial, **kwargs) for var in self])

    def set(self, values: Iterable[float], **kwargs) -> None:
        r"""Set the :py:class:`Variable`\ s' values

        Args:
            values:     Iterable of values
        """
        for var, val in zip(self, values):
            var.set(val, **kwargs)

    def increment(self, increment: Iterable[float], **kwargs) -> None:
        r"""Increment the :py:class:`Variable`\ s' values

        Args:
            increment:  Iterable of values
        """
        for var, incr in zip(self, increment):
            var.increment(incr, **kwargs)

    # noinspection PyProtectedMember
    def status(self, **kwargs) -> str:
        """String description of the variables"""
        values = "\n".join(var._line(**kwargs) for var in self)
        return "\n".join((Variable._header(), values))

    def __str__(self) -> str:
        return self.status()

    @property
    def deltas(self) -> Sequence[Number]:
        """delta values of the variables"""
        return np.array([var.delta for var in self])
