from __future__ import annotations
import numpy
import abc
from numbers import Number
from operator import add, sub, mul, truediv


class Variable(abc.ABC):
    """A :py:class:`Variable base class`
    """
    counter = 0

    def __init__(self,
                 name: str = '',
                 bounds: tuple[float, float] = (-numpy.inf, numpy.inf),
                 delta: float = 1.0,
                 fun_args: tuple = (),
                 needs_ring: bool = False):
        """
        Parameters:
            name:       Name of the Variable
            bounds:     Lower and upper bounds of the variable value
            delta:      Initial variation step
        """
        self.name = name if name else self.newname()
        self.bounds = bounds
        self.delta = delta
        self._history = []
        self.args = fun_args
        self.needs_ring = needs_ring

    @staticmethod
    def newname():
        Variable.counter = Variable.counter+1
        return f"var{Variable.counter}"

    @abc.abstractmethod
    def setfun(self, value: float, *args):
        ...

    @abc.abstractmethod
    def getfun(self, *args) -> float:
        ...

    @property
    def history(self) -> list[float]:
        """History of the values of the variable"""
        return self._history

    @property
    def initial_value(self) -> float:
        """Initial value of the variable"""
        if len(self._history) > 0:
            return self._history[0]
        else:
            raise IndexError(f"{self.name}: No value has been set yet")

    @property
    def last_value(self) -> float:
        """Last value of the variable"""
        if len(self._history) > 0:
            return self._history[-1]
        else:
            raise IndexError(f"{self.name}: No value has been set yet")

    @property
    def previous_value(self) -> float:
        """Value before the last one"""
        if len(self._history) > 1:
            return self._history[-2]
        else:
            raise IndexError(f"{self.name}: history too short")

    def set(self, value: float) -> None:
        """Set the variable value

        Args:
            value:  New value to be applied on the variable
        """
        if value < self.bounds[0] or value > self.bounds[1]:
            raise ValueError(f"set value must be in {self.bounds}")
        self.setfun(value, *self.args)
        self._history.append(value)

    def get(self, initial=False) -> float:
        """Get the actual variable value

        Args:
            initial:    If :py:obj:`True`, set the variable initial value

        Returns:
            value:      Value of the variable
        """
        value = self.getfun(*self.args)
        if initial:
            self._history = []
        h = self._history
        if len(h) == 0 or value != h[-1]:
            h.append(value)
        return value

    def set_previous(self) -> None:
        """Set to the value before the last one
        """
        self.set(self.previous_value)

    def set_initial(self) -> None:
        """Set to initial value
        """
        self.set(self.initial_value)

    def increment(self, incr: float) -> None:
        """Increment the variable value
            incr:   Increment value
        """
        if len(self._history) == 0:
            self.get(initial=True)
        self.set(self.last_value + incr)

    def _step(self, step: float) -> None:
        self.set(self.initial_value + step)

    def step_up(self) -> None:
        """Set to initial_value + delta"""
        self._step(self.delta)

    def step_down(self) -> None:
        """Set to initial_value - delta"""
        self._step(-self.delta)

    @staticmethod
    def _header():
        return '\n{:>12s}{:>13s}{:>16s}{:>16s}\n'.format(
            'Name', 'Initial', 'Final ', 'Variation')

    def _line(self):
        if len(self._history) > 0:
            vnow = self._history[-1]
            vini = self._history[0]
        else:
            vnow = vini = numpy.nan

        return '{:>12s}{: 16e}{: 16e}{: 16e}'.format(
            self.name, vini, vnow, (vnow - vini))

    def status(self):
        """Return a string describing the current status of the variable

        Returns:
            status: Variable description
        """
        return "\n".join((self._header(), self._line()))

    def __str__(self):
        return self.status()


class _Numeric(object):
    __slots__ = 'value'

    def __init__(self, value):
        self.value = value

    def __call__(self):
        return self.value


class _OpParam(object):
    __slots__ = ['oper', 'left', 'right']

    @staticmethod
    def _set_type(value):
        if isinstance(value, Number):
            return _Numeric(value)
        elif isinstance(value, Param):
            return value
        else:
            msg = ("Param Operation not defined for "
                   "type {0}".format(type(value)))
            raise TypeError(msg)

    def __init__(self, oper, left, right):
        self.oper = oper
        self.right = self._set_type(right)
        self.left = self._set_type(left)

    def __call__(self):
        return self.oper(self.left.value, self.right.value)


class Param(Variable):

    def __init__(self, value=numpy.NaN,
                 name: str = '',
                 bounds: tuple[float, float] = (-numpy.inf, numpy.inf),
                 delta: float = 1.0):
        self.setvalue(value)
        super(Param, self).__init__(name=name, bounds=bounds, delta=delta)

    def setfun(self, value):
        self.value = value

    def getfun(self):
        return float(self.value)

    def status(self, _, **kwargs):
        return "\n".join((self._header(), self._line()))

    def set(self, _, value: float) -> None:
        if value < self.bounds[0] or value > self.bounds[1]:
            raise ValueError(f"set value must be in {self.bounds}")
        self.setfun(value, *self.args)
        self._history.append(value)

    def get(self, _, initial=False) -> float:
        value = self.getfun(*self.args)
        if initial:
            self._history = []
        h = self._history
        if len(h) == 0 or value != h[-1]:
            h.append(value)
        return value

    def setvalue(self, param):
        if isinstance(param, Number):
            self._fun = _Numeric(param)
        else:
            self._fun = param

    @property
    def value(self):
        return self._fun()

    @value.setter
    def value(self, value):
        self.setvalue(value)

    def __add__(self, other):
        fun = _OpParam(add, self, other)
        return Param(fun)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        fun = _OpParam(sub, self, other)
        return Param(fun)

    def __rsub__(self, other):
        fun = _OpParam(sub, other, self)
        return Param(fun)

    def __mul__(self, other):
        fun = _OpParam(mul, self, other)
        return Param(fun)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        fun = _OpParam(truediv, self, other)
        return Param(fun)

    def __rtruediv__(self, other):
        fun = _OpParam(add, other, self)
        return Param(fun)

    def __repr__(self):
        return repr(self.value)

    def __str__(self):
        return str(self.value)

    def __float__(self):
        return float(self.value)

    def __int__(self):
        return int(self.value)


class ParamArray(list):

    def __init__(self, value, shape=-1, dtype=numpy.float64):
        if not isinstance(value, list):
            value = list(value)
        self.shape = shape
        self.dtype = dtype
        self._list = value
        self._value = self._build_array()
        super(ParamArray, self).__init__(self._list)

    def _build_array(self):
        arr = numpy.require(self._list, dtype=self.dtype,
                            requirements=['F', 'A'])
        return arr.reshape(self.shape, order='F')

    @property
    def value(self):
        self._value[:] = self._list[:]
        return self._value

    @value.setter
    def value(self, value):
        if not isinstance(value, list):
            value = list(value)
        self._list = value
        self._value = self._build_array()

    def __repr__(self):
        return repr(self.value)

    def __str__(self):
        return str(self.value)
