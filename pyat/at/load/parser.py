from __future__ import annotations

__all__ = ["StrParser", "StrParameter"]

import abc
from collections.abc import Callable
from typing import Any
from ..lattice import ParamDef, _nop


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

    @abc.abstractmethod
    def _check_constant(self, expr: str) -> Any:
        """Check if an expression is constant

        This method attempts to evaluate the expression in a context where no variables
        are defined. If the evaluation succeeds, the expression is considered constant
        and the evaluated value is returned. If the evaluation fails with a NameError,
        the expression is considered non-constant (i.e., it depends on variables).

        Args:
            expr: The string expression to evaluate

        Returns:
            The result of evaluating the expression if it's constant

        Raises:
            NameError: If the expression contains variables
        """
        ...


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
        if str(other) == "0.0":
            return self
        if self.expr == "0.0":
            return other
        return self.__class__(self.parser, f"({self.expr})+({other})")

    __radd__ = __add__

    def __sub__(self, other):
        if str(other) == "0.0":
            return self
        if self.expr == "0.0":
            return -other
        return self.__class__(self.parser, f"({self.expr})-({other})")

    def __rsub__(self, other):
        return self.__class__(self.parser, f"({other})-({self.expr})")

    def __mul__(self, other):
        if str(other) == "0.0":
            return 0.0
        if self.expr == "1.0":
            return other
        if str(other) == "1.0":
            return self
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

    def __abs__(self):
        return self.__class__(self.parser, f"abs({self.expr})")

    def __gt__(self, other):
        return self.value > other

    def __ge__(self, other):
        return self.value >= other

    def __lt__(self, other):
        return self.value < other

    def __le__(self, other):
        return self.value <= other

    @classmethod
    def parameter(cls, parser, expr: str):
        try:
            val = parser._check_constant(expr)
        except NameError:
            return cls(parser, expr)
        else:
            return val

    @property
    def value(self) -> Any:
        """Current value of the parameter"""
        return self._conversion(self.parser.evaluate(self.expr))
