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
    :py:meth:`~VariableBase._setfun` methods
    """

    # Class constants
    DEFAULT_DELTA = 1.0
    _COUNTER_PREFIX = "var"

    _counter = 0

    def __init__(
        self,
        *args,
        name: str = "",
        bounds: tuple[Number | None, Number | None] | None = None,
        delta: Number = DEFAULT_DELTA,
        history_length: int | None = None,
        **kwargs,
    ) -> None:
        """
        Parameters:
            *args:      Positional arguments passed to the _setfun and _getfun methods
            name:       Name of the Variable
            bounds:     Lower and upper bounds of the variable value
            delta:      Initial variation step
            history_length: Maximum length of the history buffer. :py:obj:`None`
              means infinite

        Keyword Args:
            **kwargs:    Keyword arguments passed to the _setfun and _getfun methods
        """
        self.name: str = self._generate_name(name)  #: Variable name
        self.args = args
        self.kwargs = kwargs
        if bounds is None:
            bounds = (None, None)
        self._bounds: tuple[Number | None, Number | None] = bounds  #: Variable bounds
        self.delta: Number = delta  #: Increment step
        #: Maximum length of the history buffer. :py:obj:`None` means infinite
        self.history_length = history_length
        self._initial = np.nan
        self._history = deque([], self.history_length)
        try:
            self.get(initial=True)
        except ValueError:
            pass

    @classmethod
    def _generate_name(cls, name: str) -> str:
        """Generate unique name for variable"""
        cls._counter += 1
        return name if name else f"{cls._COUNTER_PREFIX}{cls._counter}"

    @property
    def bounds(self) -> tuple[float, float]:
        """Bounds of the variable"""
        vmin, vmax = self._bounds
        return -np.inf if vmin is None else vmin, np.inf if vmax is None else vmax

    def check_bounds(self, value: Number) -> None:
        """Check that a value is within the variable bounds

        Raises:
            ValueError: If the value is not within bounds
        """
        min_val, max_val = self._bounds
        if min_val is not None and value < min_val:
            raise ValueError(f"Value {value} must be larger or equal to {min_val}")
        if max_val is not None and value > max_val:
            raise ValueError(f"Value {value} must be smaller or equal to {max_val}")

    def _setfun(self, value: Number, *args, **kwargs) -> None:
        classname = self.__class__.__name__
        raise TypeError(f"{classname!r} is read-only")

    @abc.abstractmethod
    def _getfun(self, *args, **kwargs) -> Number: ...

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

    def set(self, value: Number, **setkw) -> None:
        """Set the variable value

        Args:
            value:      New value to be applied on the variable

        Keyword Args:
            ring:       Depending on the variable type, a :py:class:`.Lattice` argument
              may be necessary to set the variable.
            **setkw:    Other keyword arguments to be passed to the setfun function.
              They augment the keyword arguments given in the constructor.
        """
        self.check_bounds(value)
        kw = self.kwargs.copy()
        kw.update(setkw)
        self._setfun(value, *self.args, **kw)
        if np.isnan(self._initial):
            self._initial = value
        self._history.append(value)

    def get(
        self, *, initial: bool = False, check_bounds: bool = False, **getkw
    ) -> Number:
        """Get the actual variable value

        Args:
            initial:    If :py:obj:`True`, clear the history and set the variable
              initial value
            check_bounds: If :py:obj:`True`, raise a ValueError if the value is out
              of bounds

        Keyword Args:
            ring:       Depending on the variable type, a :py:class:`.Lattice` argument
              may be necessary to set the variable.
            **getkw:    Other keyword arguments to be passed to the setfun function.
              They augment the keyword arguments given in the constructor.

        Returns:
            value:      Value of the variable
        """
        kw = self.kwargs.copy()
        kw.update(getkw)
        value = self._getfun(*self.args, **kw)
        if initial or np.isnan(self._initial):
            self._initial = value
            self._history = deque([value], self.history_length)
        if check_bounds:
            self.check_bounds(value)
        return value

    value = property(get, set, doc="Actual value")

    @property
    def _safe_value(self):
        try:
            return self._history[-1]
        except IndexError:
            return np.nan

    def set_previous(self, **setkw) -> None:
        """Reset to the value before the last one

        Keyword Args:
            ring:       Depending on the variable type, a :py:class:`.Lattice` argument
              may be necessary to set the variable.
            **setkw:    Other keyword arguments to be passed to the setfun function.
              They augment the keyword arguments given in the constructor.

        Raises:
            IndexError: If there are fewer than 2 values in the history
        """
        if len(self._history) >= 2:
            self._history.pop()  # Remove the last value
            previous_value = self._history.pop()  # retrieve the previous value
            self.set(previous_value, **setkw)
        else:
            raise IndexError(f"{self.name}: history too short (need at least 2 values)")

    def reset(self, **setkw) -> None:
        """Reset to the initial value and clear the history buffer

        Keyword Args:
            ring:       Depending on the variable type, a :py:class:`.Lattice` argument
              may be necessary to set the variable.
            **setkw:    Other keyword arguments to be passed to the setfun function.
              They augment the keyword arguments given in the constructor.

        Raises:
            IndexError: If no initial value has been set yet
        """
        initial_value = self._initial
        if not np.isnan(initial_value):
            self._history = deque([], self.history_length)
            self.set(initial_value, **setkw)
        else:
            raise IndexError(
                f"Cannot reset {self.name}: No initial value has been set yet"
            )

    def increment(self, incr: Number, **setkw) -> None:
        """Increment the variable value

        Args:
            incr:   Increment value

        Keyword Args:
            ring:       Depending on the variable type, a :py:class:`.Lattice` argument
              may be necessary to set the variable.
            **setkw:    Other keyword arguments to be passed to the setfun function.
              They augment the keyword arguments given in the constructor.
        """
        try:
            current_value = self.last_value
        except IndexError:
            # If no value has been set yet, get the initial value
            self.get(initial=True, **setkw)
            current_value = self.last_value

        self.set(current_value + incr, **setkw)

    def _step(self, step: Number, **setkw) -> None:
        """Apply a step relative to the initial value"""
        try:
            initial_value = self.initial_value
        except IndexError:
            # If no initial value has been set yet, get it
            self.get(initial=True, **setkw)
            initial_value = self.initial_value

        self.set(initial_value + step, **setkw)

    def step_up(self, **setkw) -> None:
        """Set to initial_value + delta

        Keyword Args:
            ring:       Depending on the variable type, a :py:class:`.Lattice` argument
              may be necessary to set the variable.
            **setkw:    Other keyword arguments to be passed to the setfun function.
              They augment the keyword arguments given in the constructor.
        """
        self._step(self.delta, **setkw)

    def step_down(self, **setkw) -> None:
        """Set to initial_value - delta

        Keyword Args:
            ring:       Depending on the variable type, a :py:class:`.Lattice` argument
              may be necessary to set the variable.
            **setkw:    Other keyword arguments to be passed to the setfun function.
              They augment the keyword arguments given in the constructor.
        """
        self._step(-self.delta, **setkw)

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

    def __float__(self):
        return float(self._safe_value)

    def __int__(self):
        return int(self._safe_value)

    def __str__(self):
        return self.name
#       return f"{self.__class__.__name__}({self._safe_value}, name={self.name!r})"

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
        setfun: Callable[..., None],
        getfun: Callable[..., Number],
        *args,
        name: str = "",
        bounds: tuple[Number, Number] | None = None,
        delta: Number = 1.0,
        history_length: int | None = None,
        **kwargs,
    ) -> None:
        """
        Parameters:
            getfun:     Function for getting the variable value. Called as
              :pycode:`getfun(*args, **kwargs) -> Number`.
              *args* are the positional arguments given to the constructor, *kwargs*
              are the keyword arguments given to the constructor augmented with the
              keywords given to the :py:meth:`~.Variable.get` function.
            setfun:     Function for setting the variable value. Called as
              :pycode:`setfun(value: Number, *args, **kwargs): None`.
              *args* are the positional arguments given to the constructor, *kwargs*
              are the keyword arguments given to the constructor augmented with the
              keywords given to the :py:meth:`~.Variable.set` function.
            name:       Name of the Variable. If empty, a unique name is generated.
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
        super().__init__(
            *args,
            name=name,
            bounds=bounds,
            delta=delta,
            history_length=history_length,
            **kwargs,
        )

    def _getfun(self, *args, **kwargs) -> Number:
        return self.getfun(*args, **kwargs)

    def _setfun(self, value: Number, *args, **kwargs) -> None:
        self.setfun(value, *self.args, **self.kwargs)


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

    def reset(self, ring=None) -> None:
        """Reset to all variables their initial value and clear their history buffer

        Args:
            ring:   Depending on the variable type, a :py:class:`.Lattice` argument
              may be necessary to reset the variable.
        """
        for var in self:
            var.reset(ring=ring)

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
