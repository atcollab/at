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
    "CustomVariable",
    "VariableList",
]

import abc
from collections import deque
from collections.abc import Iterable, Sequence, Callable
from typing import Union

import numpy as np
import numpy.typing as npt

Number = Union[int, float]


def _nop(value):
    return value


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
        self.name: str = self._setname(name)  #: Variable name
        self.bounds: tuple[Number, Number] = bounds  #: Variable bounds
        self.delta: Number = delta  #: Increment step
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
            raise IndexError(f"{self.name}: history too short")

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
