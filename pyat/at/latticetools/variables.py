"""
Definition of :py:class:`.Variable` objects used in matching and
response matrices
"""
from __future__ import annotations
from typing import Optional
from collections.abc import Callable, Iterable, Sequence
import numpy as np
from ..lattice import Lattice, Refpts

# Observables must be pickleable. For this, the set and get functions must be
# module-level functions. No inner, nested function are allowed. So nested
# functions are replaced be module-level callable class instances


class _getf(object):
    def __init__(self, index):
        self.index = index

    def __call__(self, elem, attrname):
        return getattr(elem, attrname)[self.index]


class _setf(object):
    def __init__(self, index: int):
        self.index = index

    def __call__(self, elem, attrname, value):
        getattr(elem, attrname)[self.index] = value


class _setfun(object):
    def __init__(self, refpts: Refpts, attrname: str, index: int):
        self.refpts = refpts
        self.attrname = attrname
        if index is None:
            self.fun = setattr
        else:
            self.fun = _setf(index)

    def __call__(self, ring: Lattice, value: float):
        for elem in ring.select(self.refpts):
            self.fun(elem, self.attrname, value)


class _getfun(object):
    def __init__(self, refpts: Refpts, attrname: str, index: int):
        self.refpts = refpts
        self.attrname = attrname
        if index is None:
            self.fun = getattr
        else:
            self.fun = _getf(index)

    def __call__(self, ring: Lattice):
        values = np.array([self.fun(elem, self.attrname) for elem in
                           ring.select(self.refpts)])
        return np.average(values)


class Variable(object):
    """A :py:class:`Variable` is a scalar value acting on a lattice
    """
    def __init__(self,
                 setfun: Callable[[Lattice, float, ...], None],
                 getfun: Callable[[Lattice, ...], float],
                 *args,
                 name: str = '',
                 bounds: tuple[float, float] = (-np.inf, np.inf),
                 delta: float = 1.0,
                 **kwargs):
        """
        Parameters:
            setfun:     User-defined function for setting the Variable. Called
              as:

              :pycode:`setfun(ring, value, *args, **kwargs)`

              where *value* is the scalar value to apply

              The positional and keyword parameters come from the
              :py:class:`Variable` initialisation
            getfun:     User-defined function for retrieving the actual value
              of the variable: Called as:

              :pycode:`value = getfun(ring, *args, **kwargs)`

              The positional and keyword parameters come from the
              :py:class:`Variable` initialisation
            name:       Name of the Variable. Default: ``''``
            bounds:     Lower and upper bounds of the variable value
            delta:      Step
            *args:      Positional arguments transmitted to ``setfun`` and
              ``getfun`` functions
            **kwargs:   Keyword arguments transmitted to ``setfun``and
              ``getfun`` functions
        """
        self.setfun = setfun
        self.getfun = getfun
        self.name = name
        self.bounds = bounds
        self.delta = delta
        self.args = args
        self.kwargs = kwargs
        self._history = []

    @property
    def history(self):
        """History of the values of the variable"""
        return self._history

    @property
    def initial_value(self):
        if len(self._history) > 0:
            return self._history[0]
        else:
            mess = "{}: Initial value has not been defined"
            raise ValueError(mess.format(self.name))

    @property
    def last_value(self):
        """Last value of the variable"""
        if len(self._history) > 0:
            return self._history[-1]
        else:
            mess = "{}: No value has been set yet"
            raise ValueError(mess.format(self.name))

    @property
    def previous_value(self):
        """Value before the last one"""
        if len(self._history) > 1:
            return self._history[-2]
        else:
            mess = "{}: history too short"
            raise ValueError(mess.format(self.name))

    def set(self, ring: Lattice, value: float):
        """Set the variable value

        Args:
            ring:   Lattice description
            value:  New value to be applied on the variable
        """
        self._history.append(value)
        self.setfun(ring, value, *self.args, **self.kwargs)

    def get(self, ring: Lattice, initial=False) -> float:
        """Get the actual variable value

        Args:
            ring:   Lattice description
            initial:    If :py:obj:`True`, set the variable initial value

        Returns:
            value:      Value of the variable
        """
        value = self.getfun(ring, *self.args, **self.kwargs)
        if initial:
            self._history = []
        h = self._history
        if len(h) == 0 or value != h[-1]:
            h.append(value)
        return value

    def set_previous(self, ring: Lattice):
        """Set to the value before the last one"""
        try:
            _ = self._history.pop()
            previous = self._history.pop()
        except IndexError:
            mess = "{}: cannot find a previous value: history too short"
            raise IndexError(mess.format(self.name))
        else:
            self.set(ring, previous)

    def set_initial(self, ring: Lattice):
        """Set to initial"""
        self.set(ring, self.initial_value)

    def increment(self, ring: Lattice, incr: float):
        if len(self._history) == 0:
            self.get(ring, initial=True)
        self.set(ring, self.last_value + incr)

    def _step(self, ring: Lattice, step: float):
        self.set(ring, self.initial_value + step)

    def step_up(self, ring: Lattice):
        """Set to initial + delta"""
        self._step(ring, self.delta)

    def step_down(self, ring: Lattice):
        """Set to initial - delta"""
        self._step(ring, -self.delta)

    @staticmethod
    def _header():
        return '\n{:>12s}{:>13s}{:>16s}{:>16s}\n'.format(
            'Name', 'Initial', 'Final ', 'Variation')

    def _line(self, ring: Lattice = None):
        if ring is not None:
            vnow = self.get(ring)
            vini = self._history[0]
        elif len(self._history) > 0:
            vnow = self._history[-1]
            vini = self._history[0]
        else:
            vnow = vini = np.nan

        return '{:>12s}{: 16e}{: 16e}{: 16e}'.format(
            self.name, vini, vnow, (vnow - vini))

    def status(self, ring: Lattice = None):
        """Return a string describing the current status of the observable

        Args:
            ring:   Lattice description. If :py:obj:`None`, the last archived
              value is returned. Otherwise, new values are extracted from
              *ring*

        Returns:
            status: Variable description
        """
        return "\n".join((self._header(), self._line(ring)))

    def __str__(self):
        return self.status()


class ElementVariable(Variable):
    r"""An :py:class:`ElementVariable` associated with an attribute of
    :py:class:`.Lattice` elements.

    It can be:

    * a scalar attribute or
    * an element of an array attribute

    of one or several :py:class:`.Element`\ (s) of a lattice.
    """

    def __init__(self,
                 refpts: Refpts, attrname: str,
                 index: Optional[int] = None,
                 **kwargs):
        """
        Parameters:
            refpts:     Location of variable :py:class:`.Element` (s)
            attrname:   Attribute name
            index:      Index in the attribute array. Use :py:obj:`None` for
              scalar attributes

        Keyword Args:
            name:       Name of the Variable; Default: ``''``
            bounds:     Lower and upper bounds of the variable value
            delta:      Step
        """
        setfun = _setfun(refpts, attrname, index)
        getfun = _getfun(refpts, attrname, index)

        self.refpts = refpts
        super(ElementVariable, self).__init__(setfun, getfun, **kwargs)


class VariableList(list):

    def get(self, ring: Lattice, initial=False) -> Sequence[float]:
        return np.array([var.get(ring, initial=initial) for var in self])

    def set(self, ring: Lattice, values: Iterable[float]) -> None:
        for var, val in zip(self, values):
            var.set(ring, val)

    def increment(self, ring: Lattice, values: Iterable[float]):
        for var, val in zip(self, values):
            var.increment(ring, val)

    # noinspection PyProtectedMember
    def status(self, ring: Lattice = None) -> str:
        values = "\n".join(var._line(ring) for var in self)
        return "\n".join((Variable._header(), values))

    def __str__(self) -> str:
        return self.status()

    @property
    def deltas(self) -> Sequence[float]:
        return np.array([var.delta for var in self])
