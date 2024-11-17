"""RPN interpreter. From Onel Harrison:
https://onelharrison.medium.com/watch-building-a-reverse-polish-notation-rpn-evaluator-in-python-75b1c910fab6
"""
from __future__ import annotations

__all__ = ["Rpn"]

import operator as op
from math import sqrt, sin, cos, pi, pow
from typing import Any, Union

Number = Union[int, float, str]


# noinspection PyUnusedLocal
def _pop(x):
    return []


def _swap(x, y):
    return [y, x]


supported_operators = {
    "+": (op.add, 2),
    "-": (op.sub, 2),
    "*": (op.mul, 2),
    "/": (op.truediv, 2),
    "sqrt": (sqrt, 1),
    "pi": (pi, 0),
    "sin": (sin, 1),
    "cos": (cos, 1),
    "pow": (pow, 2),
    "pop": (_pop, 1),
    "swap": (_swap, 2),
}


class Rpn:

    def __init__(self):
        self.stack = []
        self.store = False
        self.operators = supported_operators.copy()

    def _mpop(self, n: int = 1) -> list[Any]:
        """Pops and returns `n` items from the stack."""
        try:
            return [self.stack.pop() for _ in range(n)]
        except IndexError:
            raise SyntaxError("RPN: Malformed expression: empty stack") from None

    @staticmethod
    def _to_num(x: Any) -> Number:
        """Converts a value to its appropriate numeric type"""
        try:
            n = float(x)
        except ValueError:
            return x
        else:
            return int(n) if n.is_integer() else n

    def _consume_token(self, token: str) -> None:
        """Consumes a token given the current stack and returns the updated stack"""
        if token == "sto":
            self.store = True
            return
        elif self.store:
            self.operators[token] = (self.stack[-1], 0)
            self.store = False
            return
        try:
            oper, narg = self.operators[token]
        except KeyError:
            result = self._to_num(token)
        else:
            if narg == 0:
                result = oper
            else:
                result = oper(*reversed(self._mpop(narg)))

        if isinstance(result, list):
            self.stack.extend(result)
        else:
            self.stack.append(result)

    def input(self, expr: str) -> None:
        for token in expr.split():
            self._consume_token(token)

    def evaluate(self, expr: str) -> Number:
        """Evaluate an RPN expression and return the result"""
        for token in expr.split():
            self._consume_token(token)

        result, = self._mpop(1)
        if self.stack:
            raise SyntaxError("RPN: Found extra tokens")
        return result
