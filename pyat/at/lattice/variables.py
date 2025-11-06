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
    "AttributeVariable",
    "CustomVariable",
    "ItemVariable",
    "VariableBase",
    "VariableList",
    "attr_",
    "membergetter",
]

import abc
from collections import deque
from collections.abc import (
    Iterable,
    Sequence,
    Callable,
    Generator,
    MutableMapping,
    MutableSequence,
)
from contextlib import suppress, contextmanager
from operator import itemgetter, attrgetter

import numpy as np
import numpy.typing as npt

Number = int | float


class attr_(str):
    """subclass of :py:class:`str` used to differentiate directory keys and
    attribute names.
    """

    __slots__ = []


class _AttributeAccessor:
    """Class object setting/getting an attribute of an object."""

    __slots__ = ["attrname", "obj"]

    def __init__(self, obj, attrname: str):
        self.obj = obj
        self.attrname = attrname

    def set(self, value: float):
        setattr(self.obj, self.attrname, value)

    def get(self) -> float:
        return getattr(self.obj, self.attrname)


class _ItemAccessor:
    """Class object setting/getting an item of an object."""

    __slots__ = ["key", "obj"]

    def __init__(self, obj: MutableMapping | MutableSequence, key):
        self.obj = obj
        self.key = key

    def set(self, value: float):
        self.obj[self.key] = value

    def get(self) -> float:
        return self.obj[self.key]


class membergetter:
    """Generalised attribute and item lookup.

    Callable object fetching attributes or items from its operand object. This
    generalises :py:func:`~operator.attrgetter` and :py:func:`~operator.itemgetter` and
    allows to extract elements deep in the object structure. For example:

    - With ``f1 = membergetter("key1", [2])``, then ``f1(obj)`` returns
      ``obj["key1"][2]``,
    - With ``f2 = membergetter("key2", attr_("attr1"), [key3])``, then ``f2(obj)``
      returns ``obj["key2"].attr1["key3"]``.
    """

    def __init__(self, *args):
        r"""
        Args:
            *args:      Sequence of dictionary keys, sequence indices or
              attribute names. A :py:class:`str` argument is interpreted as a
              dictionary key. Attribute names must be decorated with ``attr_(attrname)``
              to distinguish them from directory keys.

              - ``f1 = SetGet("key1")(obj)`` returns ``obj["key1"]``,
              - ``f2 = SetGet(attr("attr1"))(obj)`` returns ``obj.attr1``

        Example:
            >>> dct = {"a": 42.0, "b": [0.0, 1.0, 2.0, 3.0]}
            >>> f = membergetter("b", 1)
            >>> f(dct)
            1.0
        """

        def getter(key):
            return attrgetter(key) if isinstance(key, attr_) else itemgetter(key)

        self.key = args[-1]
        self.getters = [getter(key) for key in args]

    def __call__(self, obj):
        for getter in self.getters:
            obj = getter(obj)
        return obj

    def accessor(self, obj):
        """Return an accessor object.

        The returned object has *set* and *get* methods acting on the selected item of
        the object *obj*.

        Example:
            >>> dct = {"a": 42.0, "b": [0.0, 1.0, 2.0, 3.0]}
            >>> v2 = membergetter("b", 1)
            >>> accessor = v2.accessor(dct)
            >>> accessor.get()
            1.0
            >>> accessor.set(42.0)
            >>> dct
            {'a': 42.0, 'b': [0.0, 42.0, 2.0, 3.0]}
        """
        for getter in self.getters[:-1]:
            obj = getter(obj)
        if isinstance(self.key, attr_):
            return _AttributeAccessor(obj, self.key)
        else:
            return _ItemAccessor(obj, self.key)


class VariableBase(abc.ABC):
    """A Variable abstract base class.

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
            name:       Name of the Variable. If omitted or blank, a unique name is
              generated.
            bounds:     Lower and upper bounds of the variable value
            delta:      Initial variation step
            history_length: Maximum length of the history buffer. :py:obj:`None`
              means infinite.

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
        self._history: deque[Number] = deque([], self.history_length)
        with suppress(ValueError):
            self.get(initial=True)

    @classmethod
    def _generate_name(cls, name: str) -> str:
        """Generate unique name for variable."""
        cls._counter += 1
        return name if name else f"{cls._COUNTER_PREFIX}{cls._counter}"

    @property
    def bounds(self) -> tuple[float, float]:
        """Lower and upper bounds of the variable."""
        vmin, vmax = self._bounds
        return -np.inf if vmin is None else vmin, np.inf if vmax is None else vmax

    def check_bounds(self, value: Number) -> None:
        """Check that a value is within the variable bounds.

        Raises:
            ValueError: If the value is not within bounds
        """
        min_val, max_val = self._bounds
        if min_val is not None and value < min_val:
            msg = f"Value {value} must be larger or equal to {min_val}"
            raise ValueError(msg)
        if max_val is not None and value > max_val:
            msg = f"Value {value} must be smaller or equal to {max_val}"
            raise ValueError(msg)

    def _setfun(self, value: Number, *args, **kwargs) -> None:
        classname = self.__class__.__name__
        msg = f"{classname!r} is read-only"
        raise TypeError(msg)

    @abc.abstractmethod
    def _getfun(self, *args, **kwargs) -> Number: ...

    @property
    def history(self) -> list[Number]:
        """History of the values of the variable."""
        return list(self._history)

    @property
    def initial_value(self) -> Number:
        """Initial value of the variable.

        Raises:
            IndexError: If no initial value has been set yet
        """
        if not np.isnan(self._initial):
            return self._initial
        else:
            msg = f"{self.name}: No initial value has been set yet"
            raise IndexError(msg)

    @property
    def last_value(self) -> Number:
        """Last value of the variable.

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
        """Value before the last one.

        Raises:
            IndexError: If there are fewer than 2 values in the history
        """
        try:
            return self._history[-2]
        except IndexError as exc:
            exc.args = (f"{self.name}: history too short (need at least 2 values)",)
            raise

    def set(self, value: Number, **setkw) -> None:
        """Set the variable value.

        Args:
            value:      New value to be applied on the variable

        Keyword Args:
            ring:       Depending on the variable type, a :py:class:`.Lattice` argument
              may be necessary to set the variable.
            **setkw:    Other keyword arguments to be passed to the *setfun* function.
              They augment the keyword arguments given in the constructor.
        """
        self.check_bounds(value)
        self._setfun(value, *self.args, **(self.kwargs | setkw))
        if np.isnan(self._initial):
            self._initial = value
        self._history.append(value)

    def get(
        self, *, initial: bool = False, check_bounds: bool = False, **getkw
    ) -> Number:
        """Get the actual variable value.

        Args:
            initial:    If :py:obj:`True`, clear the history and set the variable
              initial value
            check_bounds: If :py:obj:`True`, raise a ValueError if the value is out
              of bounds

        Keyword Args:
            ring:       Depending on the variable type, a :py:class:`.Lattice` argument
              may be necessary to set the variable.
            **getkw:    Other keyword arguments to be passed to the *getfun* function.
              They augment the keyword arguments given in the constructor.

        Returns:
            value:      Value of the variable
        """
        value = self._getfun(*self.args, **(self.kwargs | getkw))
        if initial or np.isnan(self._initial):
            self._initial = value
            self._history = deque([value], self.history_length)
        if check_bounds:
            self.check_bounds(value)
        return value

    value = property(get, set, doc="Actual value")

    @property
    def _print_value(self):
        try:
            return self._history[-1]
        except IndexError:
            return np.nan

    def set_previous(self, **setkw) -> None:
        """Reset to the value before the last one.

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
            msg = f"{self.name}: history too short (need at least 2 values)"
            raise IndexError(msg)

    def reset(self, **setkw) -> None:
        """Reset to the initial value and clear the history buffer.

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
            msg = f"Cannot reset {self.name}: No initial value has been set yet"
            raise IndexError(msg)

    def increment(self, incr: Number, **setkw) -> None:
        """Increment the variable value.

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
        """Apply a step relative to the initial value."""
        try:
            initial_value = self.initial_value
        except IndexError:
            # If no initial value has been set yet, get it
            self.get(initial=True, **setkw)
            initial_value = self.initial_value

        self.set(initial_value + step, **setkw)

    def step_up(self, **setkw) -> None:
        """Set to initial_value + delta.

        Keyword Args:
            ring:       Depending on the variable type, a :py:class:`.Lattice` argument
              may be necessary to set the variable.
            **setkw:    Other keyword arguments to be passed to the setfun function.
              They augment the keyword arguments given in the constructor.
        """
        self._step(self.delta, **setkw)

    def step_down(self, **setkw) -> None:
        """Set to initial_value - delta.

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
        vnow = self._print_value
        vini = self._initial
        return f"{self.name:>12s}{vini: 16e}{vnow: 16e}{vnow - vini: 16e}"

    def status(self):
        """Return a string describing the current status of the variable.

        Returns:
            status: Variable description
        """
        return "\n".join((self._header(), self._line()))

    @contextmanager
    def restore(self, **setkw) -> Generator[None, None, None]:
        # noinspection PyUnresolvedReferences
        """Context manager that saves and restore a variable.

        The value of the :py:class:`Variable <VariableBase>` is initially saved, and
        then restored when exiting the context.

        Keyword Args:
            ring:       Depending on the variable type, a :py:class:`.Lattice` argument
              may be necessary to set the variable.
            **setkw:    Other keyword arguments to be passed to the setfun function.
              They augment the keyword arguments given in the constructor.

        Example:
            >>> var = AttributeVariable(ring, "energy")
            >>> with var.restore():
            ...     do_something(var)
        """
        # print(f"save {self.name}")
        v0 = self.get(**setkw)
        try:
            yield
        finally:
            # print(f"restore {self.name}")
            self.set(v0, **setkw)

    def __float__(self):
        return float(self._print_value)

    def __int__(self):
        return int(self._print_value)

    def __str__(self):
        return self.name

    def __repr__(self):
        return repr(self._print_value)


class ItemVariable(VariableBase):
    """A Variable controlling an item of a dictionary or a sequence."""

    def __init__(
        self, obj: MutableSequence | MutableMapping, key, *args, **kwargs
    ) -> None:
        # noinspection PyUnresolvedReferences
        """
        Args:
            obj:        Mapping or Sequence containing the variable value,
            key:        Index or attribute name of the variable.  A :py:class:`str`
              argument is interpreted as a dictionary key. Attribute names must be
              decorated with ``attr_(attrname)`` to distinguish them from directory
              keys.
            *args:      additional sequence of indices or attribute names allowing to
              extract elements deeper in the object structure.

        Keyword Args:
            name:       Name of the Variable. If empty, a unique name is generated.
            bounds:     Lower and upper bounds of the variable value
            delta:      Initial variation step
            history_length: Maximum length of the history buffer. :py:obj:`None`
              means infinite.

        Example:
            >>> dct = {"a": 42.0, "b": [0.0, 1.0, 2.0, 3.0]}
            >>> v1 = at.ItemVariable(dct, "a")
            >>> v1.value
            42.0

            *v1* points to ``dct["a"]``

            >>> v2 = at.ItemVariable(dct, "b", 1)
            >>> v2.value
            1.0

            *v2* points to ``dct["b"][1]``
        """
        super().__init__(membergetter(key, *args).accessor(obj), **kwargs)

    def _setfun(self, value, obj):
        obj.set(value)

    def _getfun(self, obj):
        return obj.get()


class AttributeVariable(VariableBase):
    """A Variable controlling an attribute of an object."""

    def __init__(self, obj, attrname: str, index: int | None = None, **kwargs):
        # noinspection PyUnresolvedReferences
        """
        Args:
            obj:        Object containing the variable value
            attrname:   attribute name of the variable
            index:      Index in the attribute array. Use :py:obj:`None` for
              scalar attributes.

        Keyword Args:
            name:       Name of the Variable. If empty, a unique name is generated.
            bounds:     Lower and upper bounds of the variable value
            delta:      Initial variation step
            history_length: Maximum length of the history buffer. :py:obj:`None`
              means infinite.

        Example:
            >>> ring = at.Lattice.load("hmba.mat")
            >>> v3 = at.AttributeVariable(ring, "energy")
            >>> v3.value
            6000000000.0

            *v3* points to the *"energy"* attribute of *ring*
        """
        args = (attr_(attrname),) if index is None else (attr_(attrname), index)
        super().__init__(membergetter(*args).accessor(obj), **kwargs)

    def _setfun(self, value, obj):
        obj.set(value)

    def _getfun(self, obj):
        return obj.get()


class CustomVariable(VariableBase):
    r"""A Variable with user-defined get and set functions.

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
    """Container for Variable objects.

    :py:class:`VariableList` supports all :py:class:`list` methods, like
    appending, insertion or concatenation with the "+" operator.
    """

    def __getitem__(self, index):
        if isinstance(index, slice):
            return VariableList(super().__getitem__(index))
        else:
            return super().__getitem__(index)

    def get(self, ring=None, **kwargs) -> Sequence[float]:
        r"""Get the current values of Variables.

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
        r"""Set the values of Variables.

        Args:
            values:     Iterable of values
            ring:   Depending on the variable type, a :py:class:`.Lattice` argument
              may be necessary to set the variable.
        """
        for var, val in zip(self, values, strict=False):
            var.set(val, ring=ring)

    def increment(self, increment: Iterable[float], ring=None) -> None:
        r"""Increment the values of Variables.

        Args:
            increment:  Iterable of values
            ring:   Depending on the variable type, a :py:class:`.Lattice` argument
              may be necessary to increment the variable.
        """
        for var, incr in zip(self, increment, strict=False):
            var.increment(incr, ring=ring)

    def reset(self, ring=None) -> None:
        """Reset to all variables their initial value and clear their history buffer.

        Args:
            ring:   Depending on the variable type, a :py:class:`.Lattice` argument
              may be necessary to reset the variable.
        """
        for var in self:
            var.reset(ring=ring)

    # noinspection PyProtectedMember
    def status(self, **kwargs) -> str:
        """String description of the variables."""
        values = "\n".join(var._line(**kwargs) for var in self)
        return "\n".join((VariableBase._header(), values))

    @contextmanager
    def restore(self, **kwargs):
        """Context manager that saves and restore the variable list.

        The value of the :py:class:`VariableList` is initially saved, and then restored
        when exiting the context.

        Keyword Args:
            ring:       Depending on the variable type, a :py:class:`.Lattice` argument
              may be necessary to set the variable.
        """
        # print("Saving variables")
        v0 = self.get(**kwargs)
        try:
            yield
        finally:
            # print("Restoring variables")
            self.set(v0, **kwargs)

    def __str__(self) -> str:
        return self.status()

    @property
    def deltas(self) -> npt.NDArray[Number]:
        """delta values of the variables."""
        return np.array([var.delta for var in self])

    @deltas.setter
    def deltas(self, value: Number | Sequence[Number]) -> None:
        deltas = np.broadcast_to(value, len(self))
        for var, delta in zip(self, deltas, strict=False):
            var.delta = delta
