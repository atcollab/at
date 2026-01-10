from __future__ import annotations

__all__ = ["StrParameter", "StrParser"]

import abc
from typing import Any


class StrParser(abc.ABC):
    """Abstract base class for string expression parsers.

    This class defines the interface for parsers that can evaluate string expressions.
    Concrete implementations should provide the evaluate method.
    """

    @abc.abstractmethod
    def evaluate(self, expr: str) -> Any:
        """Evaluate a string expression in the context of this parser.

        Args:
            expr: The string expression to evaluate

        Returns:
            The result of evaluating the expression
        """
        ...

    @abc.abstractmethod
    def _check_constant(self, expr: str) -> Any:
        """Check if an expression is constant.

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


class StrParameter:
    """MAD parameter.

    A MAD parameter is an expression which can be evaluated n the context
    of a MAD parser
    """

    def __init__(self, parser: StrParser, expr: str):
        """Initialise a string parameter.

        Args:
            parser: StrParser instance defining the context for evaluation
            expr: Expression to be evaluated

        The expression may contain parameter names, arithmetic operators and
        mathematical functions, depending on the capabilities of the parser.
        """
        self.expr = expr
        self.parser = parser

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
        """Current value of the parameter."""
        return self.parser.evaluate(self.expr)

    def __str__(self):
        return self.expr

    def __repr__(self):
        return repr(self.value)

    def __float__(self):
        return float(self.value)

    def __int__(self):
        return int(self.value)

    def __add__(self, other):
        return float(self) + float(other)

    def __radd__(self, other):
        return float(other) + float(self)

    def __mul__(self, other):
        return float(self) * float(other)

    def __rmul__(self, other):
        return float(other) * float(self)

    def __sub__(self, other):
        return float(self) - float(other)

    def __rsub__(self, other):
        return float(other) - float(self)

    def __truediv__(self, other):
        return float(self) / float(other)

    def __rtruediv__(self, other):
        return float(other) / float(self)

    def __pow__(self, other):
        return pow(float(self), other)

    def __rpow__(self, other):
        return pow(float(other), float(self))

    def __neg__(self):
        return -float(self)

    def __pos__(self):
        return +float(self)
