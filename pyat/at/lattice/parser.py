from __future__ import annotations

__all__ = ["StrParser", "ParamDef", "StrParameter"]

import abc
from collections.abc import Callable
from typing import Any


def _nop(value: Any) -> Any:
    """No-operation function that returns its input unchanged

    This function is used as a default conversion function in parameter classes.

    Args:
        value: Any value

    Returns:
        The input value unchanged
    """
    return value


class StrParser(abc.ABC):
    """Abstract base class for string expression parsers

    This class defines the interface for parsers that can evaluate string expressions.
    Concrete implementations should provide the evaluate method.
    """
    @abc.abstractmethod
    def evaluate(self, expr: str) -> Any:
        """Evaluate a string expression in the context of this parser

        Args:
            expr: The string expression to evaluate

        Returns:
            The result of evaluating the expression
        """
        ...


class ParamDef(abc.ABC):
    """Abstract base class for parameter definitions

    This class defines the interface for parameter objects that can be used
    as element attributes. It provides methods for getting and setting values,
    as well as for converting values to the appropriate type.
    """
    def __init__(self, *, conversion: Callable[[Any], Any] = _nop):
        """Initialise a parameter definition

        Args:
            conversion: Function to convert values to the appropriate type
        """
        self._conversion = conversion

    @abc.abstractmethod
    def get(self, **kwargs) -> Any:
        """Get the current value of the parameter

        Args:
            **kwargs: Additional arguments for specific implementations

        Returns:
            The current value of the parameter
        """
        ...

    # noinspection PyUnusedLocal
    def set(self, value: Any) -> None:
        """Set the value of the parameter

        By default, parameters are read-only. Subclasses should override this
        method to provide write access.

        Args:
            value: The new value for the parameter

        Raises:
            TypeError: Always raised by this default implementation
        """
        classname = self.__class__.__name__
        raise TypeError(f"{classname!r} is read-only")

    value = property(get, set, doc="Actual value of the parameter")

    def set_conversion(self, conversion: Callable[[Any], Any]) -> None:
        """Set the data type conversion function

        This method is called when a parameter is assigned to an
        :py:class:`.Element` attribute. It can only be set once.

        Args:
            conversion: Function to convert values to the appropriate type

        Raises:
            ValueError: If attempting to change an already set conversion function
        """
        if conversion is not self._conversion:
            if self._conversion is _nop:
                self._conversion = conversion
            else:
                raise ValueError("Cannot change the data type of the parameter")


class StrParameter(ParamDef):
    """String expression parameter

    A StrParameter represents an expression which can be evaluated in the context
    of a StrParser. This is typically used for MAD-X style parameters where
    expressions can reference other parameters.
    """

    def __init__(
        self, parser: StrParser, expr: str, conversion: Callable[[Any], Any] = _nop
    ):
        """Initialise a string parameter

        Args:
            parser: StrParser instance defining the context for evaluation
            expr: Expression to be evaluated
            conversion: Function to convert the evaluated result to the appropriate type

        The expression may contain parameter names, arithmetic operators and
        mathematical functions, depending on the capabilities of the parser.
        """
        self.expr = expr
        self.parser = parser
        super().__init__(conversion=conversion)

    def __str__(self):
        return self.expr

    def __repr__(self):
        return f"{self.value}"

    def __float__(self):
        return float(self.value)

    def __int__(self):
        return int(self.value)

    def __add__(self, other):
        return self.__class__(self.parser, f"({self.expr})+({other})")

    __radd__ = __add__

    def __sub__(self, other):
        return self.__class__(self.parser, f"({self.expr})-({other})")

    def __rsub__(self, other):
        return self.__class__(self.parser, f"({other})-({self.expr})")

    def __mul__(self, other):
        return self.__class__(self.parser, f"({self.expr})*({other})")

    __rmul__ = __mul__

    def __truediv__(self, other):
        return self.__class__(self.parser, f"({self.expr})/({other})")

    def __rtruediv__(self, other):
        return self.__class__(self.parser, f"({other})/({self.expr})")

    def __pow__(self, other):
        return self.__class__(self.parser, f"({self.expr})**({other})")

    def __rpow__(self, other):
        return self.__class__(self.parser, f"({other})**({self.expr})")

    def __neg__(self):
        return self.__class__(self.parser, f"-({self.expr})")

    def __pos__(self):
        return self.__class__(self.parser, f"({self.expr})")

    def get(self, **kwargs) -> Any:
        """Get the current value of the parameter

        Evaluates the expression using the parser and applies the conversion function.

        Args:
            **kwargs: Additional arguments (not used in this implementation)

        Returns:
            The evaluated and converted value of the expression
        """
        return self._conversion(self.parser.evaluate(self.expr))

    @property
    def value(self) -> Any:
        """Current value of the parameter"""
        return self.get()
