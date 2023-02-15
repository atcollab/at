"""
Classes for response matrix observables and variables
"""
import sys
from typing import Optional
if sys.version_info.minor < 9:
    from typing import Callable, Tuple
else:
    from collections.abc import Callable
    Tuple = tuple
import numpy as np
from ..lattice import Lattice, Refpts


class Variable(object):
    """A :py:class:`Variable` is a scalar value acting on a lattice
    """
    def __init__(self,
                 setfun: Callable[..., None],
                 getfun: Callable[..., float],
                 *args,
                 name: str = '',
                 bounds: Tuple[float, float] = (-np.inf, np.inf),
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
            self._history = [value]
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
    def header():
        return '\n{:>12s}{:>13s}{:>16s}{:>16s}\n'.format(
            'Name', 'Initial', 'Final ', 'Variation')

    def status(self, ring: Lattice):
        vnow = self.get(ring)
        vini = self.initial_value
        return '{:>12s}{: 16e}{: 16e}{: 16e}'.format(
            self.name, vini, vnow, (vnow - vini))


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
        setf, getf = self._access(index)

        def setfun(rng: Lattice, value):
            for elem in rng.select(refpts):
                setf(elem, attrname, value)

        def getfun(rng: Lattice):
            values = np.array([getf(elem, attrname) for elem in
                               rng.select(refpts)])
            return np.average(values)

        self.refpts = refpts
        super(ElementVariable, self).__init__(setfun, getfun, **kwargs)

    @staticmethod
    def _access(index):
        """Access to element attributes"""
        if index is None:
            setf = setattr
            getf = getattr
        else:
            def setf(elem, attrname, value):
                getattr(elem, attrname)[index] = value

            def getf(elem, attrname):
                return getattr(elem, attrname)[index]
        return setf, getf


class VariableList(list):

    def increment(self, ring: Lattice, values):
        for var, val in zip(self, values):
            var.increment(ring, val)

    @property
    def delta(self):
        return np.array([v.delta for v in self])
