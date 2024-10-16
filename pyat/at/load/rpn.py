"""RPN interpreter. From Onel Harrison:
https://onelharrison.medium.com/watch-building-a-reverse-polish-notation-rpn-evaluator-in-python-75b1c910fab6
"""
from __future__ import annotations

__all__ = ["evaluate"]

import operator as op
from typing import Any, Union

supported_operators = {"+": op.add, "-": op.sub, "*": op.mul, "/": op.truediv}

Number = Union[int, float]


def tokenize(expr: str) -> list[str]:
    """Breaks expression `expr` into a list of tokens"""
    return expr.split()


def mpop(stack: list[Any], n: int = 1) -> list[Any]:
    """Pops and returns `n` items from a stack. Mutates `stack`"""
    return [stack.pop() for _ in range(n)]


def to_num(x: Any) -> Number:
    """Converts a value to its appropriate numeric type"""
    n = float(x)
    return int(n) if n.is_integer() else n


def consume_token(token: str, stack: list[Number]) -> list[Number]:
    """Consumes a token given the current stack and returns the updated stack"""
    if token in supported_operators:
        try:
            num1, num2 = mpop(stack, 2)
        except IndexError:
            raise SyntaxError("RPN: Malformed expression") from None

        result = supported_operators[token](num2, num1)
        return [*stack, result]
    else:
        try:
            return [*stack, to_num(token)]
        except ValueError:
            raise SyntaxError("RPN: Unsupported token '%s'", token) from None


def get_result_from_stack(stack: list[Number]) -> Number:
    """Gets the result from `stack`"""
    result, *rest = mpop(stack, 1)
    if rest:
        raise SyntaxError("RPN: Found extra tokens")
    return result


def evaluate(expr: str) -> Number:
    """Evaluate an RPN expression and return the result"""
    stack: list = []

    for token in tokenize(expr):
        stack = consume_token(token, stack)

    return get_result_from_stack(stack)
#
#
# def evaluate_v2(tokens: list[str]) -> Number:
#     """Evaluates a tokenized expression and returns the result"""
#
#     def _evaluate(tokens: list[str], stack: list) -> Number:
#         if not tokens:
#             return get_result_from_stack(stack)
#
#         stack = consume_token(tokens[0], stack)
#
#         return _evaluate(tokens[1:], stack)
#
#     return _evaluate(tokens, [])
