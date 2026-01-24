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


class StrParameter:
    """String expression parameter.

    A StrParameter represents an expression which can be evaluated in the context
    of a StrParser. This is typically used for MAD-X style parameters where
    expressions can reference other parameters.
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
