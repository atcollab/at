from __future__ import annotations

__all__ = ["AttributeArray", "Param", "ParamArray"]

from collections.abc import Callable
from typing import Any
import weakref

import numpy as np
import numpy.typing as npt

from .parambase import Operand, ParamBase, _Constant, ParamDef, _nop
from .variables import VariableBase, Number

_ACCEPTED = ParamDef


class Param(ParamBase, VariableBase):
    """Standalone scalar parameter.

    See :py:class:`.VariableBase` for a description of inherited methods
    """

    _COUNTER_PREFIX = "param"

    _counter = 0

    def __init__(
        self,
        value: Number,
        *,
        name: str = "",
        conversion: Callable[[Any], Any] = _nop,
        bounds: tuple[Number, Number] | None = None,
        delta: Number = 1.0,
    ) -> None:
        """
        Args:
            value:      Initial value of the parameter
            name:       Name of the parameter. If omitted or blank, a unique name is
              generated.
            conversion: data conversion function
            bounds:     Lower and upper bounds of the parameter value
            delta:      Initial variation step.
        """
        value = conversion(value)
        super().__init__(
            _Constant(value),
            name=name,
            conversion=conversion,
            bounds=bounds,
            delta=delta,
        )

    def _getfun(self, **_) -> Number:
        return self._evaluator()

    def _setfun(self, value: Number, **_) -> None:
        self._evaluator = _Constant(self._conversion(value))

    def fast_value(self) -> Number:
        return self._evaluator()

    @property
    def value(self) -> Number:
        return self.fast_value()

    @value.setter
    def value(self, value: Number):
        self.set(value)

    def set_conversion(self, conversion: Callable[[Any], Any]) -> None:
        oldv = self._evaluator()
        super().set_conversion(conversion)
        self._evaluator = _Constant(conversion(oldv))


class _SafeArray(np.ndarray):
    """Subclass of ndarray which forbids setting parameters as items.

    This array type is used for element attributes that should not contain
    parameters. It raises a TypeError if a parameter is assigned to any element.
    """

    def __setitem__(self, key: Any, value: Any) -> None:
        """Set an item in the array, preventing parameter assignment."""
        if isinstance(value, Operand):
            msg = "Cannot set a parameter into an array"
            raise TypeError(msg)
        super().__setitem__(key, value)

    def __repr__(self) -> str:
        # Simulate a standard ndarray
        return repr(self.view(np.ndarray))


def AttributeArray(
    value: Any, shape: tuple[int, ...] = (-1,), dtype: npt.DTypeLike = float
) -> np.ndarray:
    """Create an array of attributes, which may contain parameters.

    This function creates either a :py:class:`ParamArray` (if the input contains
    parameters) or a :py:class:`_SafeArray` (if the input contains only regular values).

    Args:
        value: Input array or sequence
        shape: Shape of the output array
        dtype: Data type of the output array

    Returns:
        Either a ParamArray or a SafeArray depending on the input
    """

    v = np.asfortranarray(value).reshape(shape, order="F")
    if any(isinstance(el, _ACCEPTED) for el in v.flat):
        return ParamArray(v, shape=shape, dtype=dtype)
    else:
        return v.astype(dtype, copy=False).view(_SafeArray)


class ParamArray(np.ndarray):
    """Simulate a numpy array where items may be parameterised.

    This class allows creating arrays that can contain Parameter objects. It provides
    a `value` property that returns a numeric array with the current values of all
    parameters. Changes to the numeric array are propagated back to the parameters.

    This is primarily used for element attributes that can be parameterised.
    """

    class ValueArray(np.ndarray):
        """Subclass of ndarray which reports changes back to its parent ParamArray.

        This array is used as the value property of :py:class:`ParamArray`. When items
        in this array are modified, the changes are propagated back to the
        :py:class:`ParamArray` parent.

        This is the array obtained with an element get_attribute.
        It is also the one used when setting an item of an array attribute.
        """

        def __new__(cls, parent: ParamArray, dtype: npt.DTypeLike = float):
            """
            Args:
                parent: The parent ParamArray
                dtype: Data type of the array.

            Returns:
                A new ValueArray instance
            """
            obj = np.array(parent, dtype=dtype).view(cls)
            obj._parent = weakref.proxy(parent)
            return obj

        def __array_finalize__(self, obj: Any) -> None:
            pass

        def __setitem__(self, key: Any, value: Any) -> None:
            super().__setitem__(key, value)
            # report the value to the parent
            if self._parent is not None:
                self._parent[key] = value

        def __repr__(self) -> str:
            # Simulate a standard ndarray
            return repr(self.view(np.ndarray))

    def __new__(
        cls, value: Any, shape: tuple[int, ...] = (-1,), dtype: npt.DTypeLike = float
    ):
        obj = np.asfortranarray(value, dtype="O").reshape(shape).view(cls)
        obj._dtype = dtype
        return obj

    def __array_finalize__(self, obj: Any) -> None:
        dtype = getattr(obj, "_dtype", float)
        self._dtype = dtype
        self._value = ParamArray.ValueArray(self, dtype=dtype)

    def __setstate__(self, state) -> None:
        """Rebuild the value array after unpickling."""
        # Numpy uses the __reduce__ protocol to pickle objects, so there is no control
        # over the creation of the object. We need to rebuild the value array manually.
        super().__setstate__(state)
        self._value = ParamArray.ValueArray(self, dtype=self._dtype)

    def fast_value(self):
        """Numeric array with the current values of all the parameters.

        This property returns a :py:class:`ValueArray` that contains the numeric values
        of all parameters in the array. Changes to this array are propagated back to
        the :py:class:`ParamArray` parent.

        Returns:
            A numeric array with the current parameter values
        """
        # Update the numeric array with current parameter values
        # We use np.nditer to iterate over both arrays simultaneously
        # This is necessary because self may contain parameter objects
        # that need to be evaluated to get their current values
        with np.nditer(
            (self._value, self),
            ["external_loop", "refs_ok", "zerosize_ok"],
            [["writeonly"], ["readonly"]],
        ) as it:
            for x, y in it:
                x[...] = y[...]
        return self._value

    @property
    def value(self) -> np.ndarray:
        return self.fast_value()

    def __repr__(self) -> str:
        return repr(self.value)

    def __str__(self) -> str:
        return np.array2string(self, formatter={"all": str})
