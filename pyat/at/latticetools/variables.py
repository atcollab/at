r"""
Definition of :py:class:`.Variable` objects used in matching and
response matrices.

See :ref:`example-notebooks` for examples of matching and response matrices

Each :py:class:`Variable` has a scalar value

.. rubric:: Class hierarchy

:py:class:`Variable`\ (setfun, getfun, ...)
    :py:class:`ElementVariable`\ (refpts, attrname, index, ...)

:py:class:`Variable` provides the following methods:

- :py:meth:`~Variable.get`
- :py:meth:`~Variable.set`
- :py:meth:`~Variable.set_previous`
- :py:meth:`~Variable.set_initial`
- :py:meth:`~Variable.increment`
- :py:meth:`~Variable.step_up`
- :py:meth:`~Variable.step_down`

:py:class:`.Variable` provides the following properties:

- :py:attr:`~Variable.initial_value`
- :py:attr:`~Variable.last_value`
- :py:attr:`~Variable.previous_value`
- :py:attr:`~Variable.history`

An :py:class:`ElementVariable` is associated with an attribute of one or
several lattice elements,

The :py:class:`Variable` class may be used as a base class to define custom
variables (see examples).

:py:class:`VariableList`\ ( ...)
    A container based on :py:class:`list` which allows setting several
    :py:class:`Variable`\ s simultaneously. Used in
    :py:mod:`~at.matching`. :py:class:`Variable`\ s are set in order.

.. rubric:: Examples

Write a subclass of :py:class:`Variable` which varies two drift lengths so
that their sum is constant:

.. code-block:: python

    class ElementShifter(at.Variable):
        \"""Varies the length of the elements identified by *ref1* and *ref2*
        keeping the sum of their lengths equal to *total_length*.

        If *total_length* is None, it is set to the initial total length
        \"""
        def __init__(self, ring, ref1, ref2, total_length=None, **kwargs):
            # store indexes of the 2 variable elements
            self.ref1 = ring.get_uint32_index(ref1)[0]
            self.ref2 = ring.get_uint32_index(ref2)[0]
            # store the initial total length
            if total_length is None:
                total_length = ring[self.ref1].Length + ring[self.ref2].Length
            self.length = total_length
            super().__init__(self.setfun, self.getfun,
                         refpts=[self.ref1, self.ref2],
                         bounds=(0.0, total_length), **kwargs)

        def setfun(self, ring, value):
            ring[self.ref1].Length = value
            ring[self.ref2].Length = self.length - value

        def getfun(self, ring):
            return  ring[self.ref1].Length

And create a variable varying the length of drifts *DR_01* and *DR_01* and
keeping their sum constant:

.. code-block:: python

    var2 = ElementShifter(hmba_lattice, "DR_01", "DR_02", name="DR_01")

"""
from __future__ import annotations
from typing import Optional
import abc
from collections.abc import Iterable, Sequence
import numpy as np
from ..lattice import Lattice, Refpts, Variable

# Observables must be pickleable. For this, the set and get functions must be
# module-level functions. No inner, nested function are allowed. So nested
# functions are replaced be module-level callable class instances


class _Getf(object):
    def __init__(self, index):
        self.index = index

    def __call__(self, elem, attrname):
        return getattr(elem, attrname)[self.index]


class _Setf(object):
    def __init__(self, index: int):
        self.index = index

    def __call__(self, elem, attrname, value):
        getattr(elem, attrname)[self.index] = value


class LatticeVariable(abc.ABC):
    """A :py:class:`Variable` is a scalar value acting on a lattice
    """
    counter = 0

    def __init__(self,
                 name: str = '',
                 refpts: Refpts = None,
                 bounds: tuple[float, float] = (-np.inf, np.inf),
                 delta: float = 1.0,
                 fun_args: tuple = ()):
        """
        Parameters:
            name:       Name of the Variable
            refpts:     Identifies variable elements. It is used when a copy of
              the variable elements is needed
            bounds:     Lower and upper bounds of the variable value
            delta:      Initial variation step
        """
        self.name = name if name else self.newname()
        self.refpts = refpts
        self.bounds = bounds
        self.delta = delta
        self._history = []
        super(LatticeVariable, self).__init__(name=name,
                                              bounds=bounds,
                                              delta=delta,
                                              fun_args=fun_args,
                                              needs_ring=True)

    @abc.abstractmethod
    def setfun(self, ring: Lattice, value: float):
        ...

    @abc.abstractmethod
    def getfun(self, ring: Lattice) -> float:
        ...

    def set(self, ring: Lattice, value: float) -> None:
        """Set the variable value

        Args:
            ring:   Lattice description
            value:  New value to be applied on the variable
        """

        if value < self.bounds[0] or value > self.bounds[1]:
            raise ValueError(f"set value must be in {self.bounds}")
        self.setfun(ring, value)
        self._history.append(value)

    def get(self, ring: Lattice, initial=False) -> float:
        """Get the actual variable value

        Args:
            ring:   Lattice description
            initial:    If :py:obj:`True`, set the variable initial value

        Returns:
            value:      Value of the variable
        """
        value = self.getfun(ring)
        if initial:
            self._history = []
        h = self._history
        if len(h) == 0 or value != h[-1]:
            h.append(value)
        return value

    def set_previous(self, ring: Lattice) -> None:
        """Set to the value before the last one

        Args:
            ring:   Lattice description
        """
        self.set(ring, self.previous_value)

    def set_initial(self, ring: Lattice) -> None:
        """Set to initial value

        Args:
            ring:   Lattice description
        """
        self.set(ring, self.initial_value)

    def increment(self, ring: Lattice, incr: float) -> None:
        """Increment the variable value

        Args:
            ring:   Lattice description
            incr:   Increment value
        """
        if len(self._history) == 0:
            self.get(ring, initial=True)
        self.set(ring, self.last_value + incr)

    def _step(self, ring: Lattice, step: float) -> None:
        self.set(ring, self.initial_value + step)

    def step_up(self, ring: Lattice) -> None:
        """Set to initial_value + delta"""
        self._step(ring, self.delta)

    def step_down(self, ring: Lattice) -> None:
        """Set to initial_value - delta"""
        self._step(ring, -self.delta)

    def _line(self, ring: Lattice):
        vnow = self.get(ring)
        vini = self._history[0]

        return '{:>12s}{: 16e}{: 16e}{: 16e}'.format(
            self.name, vini, vnow, (vnow - vini))

    def status(self, ring: Lattice):
        """Return a string describing the current status of the variable

        Args:
            ring:   Lattice description. If :py:obj:`None`, the last archived
              value is returned. Otherwise, current values are extracted from
              *ring*

        Returns:
            status: Variable description
        """
        return "\n".join((self._header(), self._line(ring)))

    def __str__(self):
        return self.status()


class ElementVariable(LatticeVariable):
    r"""A :py:class:`Variable` associated with an attribute of
    :py:class:`.Lattice` elements.

    It can be:

    * a scalar attribute or
    * an element of an array attribute

    of one or several :py:class:`.Element`\ s of a lattice.
    """
    def __init__(self,
                 refpts: Refpts, attrname: str,
                 index: Optional[int] = None,
                 **kwargs):
        r"""
        Parameters:
            refpts:     Location of variable :py:class:`.Element`\ s
            attrname:   Attribute name
            index:      Index in the attribute array. Use :py:obj:`None` for
              scalar attributes

        Keyword Args:
            name (str):     Name of the Variable. Default: ``''``
            bounds (tuple[float, float]):   Lower and upper bounds of the
              variable value. Default: (-inf, inf)
            delta (float):  Step. Default: 1.0
        """
        if index is None:
            self._getf = getattr
            self._setf = setattr
        else:
            self._getf = _Getf(index)
            self._setf = _Setf(index)
        self.attrname = attrname
        super(ElementVariable, self).__init__(refpts=refpts, **kwargs)

    def setfun(self, ring: Lattice, value: float):
        for elem in ring.select(self.refpts):
            self._setf(elem, self.attrname, value)

    def getfun(self, ring: Lattice) -> float:
        values = np.array([self._getf(elem, self.attrname) for elem in
                           ring.select(self.refpts)])
        return np.average(values)


class VariableList(list):
    """Container for :py:class:`Variable` objects

    :py:class:`VariableList` supports all :py:class:`list` methods, like
    appending, insertion or concatenation with the "+" operator.
    """

    def get(self, ring: Lattice, initial=False) -> Sequence[float]:
        r"""Get the current :py:class:`Variable`\ s' values

        Args:
            ring:       Lattice description
            initial:    If :py:obj:`True`, set the :py:class:`Variable`\ s'
              initial value

        Returns:
            values:     1D array of values of all variables
        """
        return np.array([var.get(ring, initial=initial) if var.needs_ring
                         else var.get(initial=initial) for var in self])

    def set(self, ring: Lattice, values: Iterable[float]) -> None:
        r"""Set the :py:class:`Variable`\ s' values

        Args:
            ring:       Lattice description
            values:     Iterable of values
        """
        [var.set(ring, val) if var.needs_ring
         else var.set(val) for var, val in zip(self, values)]

    def increment(self, ring: Lattice, increment: Iterable[float]) -> None:
        r"""Increment the :py:class:`Variable`\ s' values

        Args:
            ring:       Lattice description
            increment:  Iterable of values
        """
        [var.increment(ring, incr) if self.needs_ring
        else var.increment(incr)
        for var, incr in zip(self, increment)]

    # noinspection PyProtectedMember
    def status(self, ring: Lattice = None) -> str:
        """String description of the variables"""
        values = "\n".join(var._line(ring) if var.needs_ring
                           else var._line() for var in self)
        return "\n".join((Variable._header(), values))

    def __str__(self) -> str:
        return self.status()

    @property
    def deltas(self) -> Sequence[float]:
        return np.array([var.delta for var in self])
