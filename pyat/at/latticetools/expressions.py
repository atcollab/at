from __future__ import annotations
import re
from collections.abc import Callable, Generator
from typing import Optional
from numbers import Number
from operator import add, sub, mul, truediv, pow
from .variables import Variable
from ..lattice import Lattice, Refpts


class _Oper(object):
    def __init__(self, fun: Callable[[float, float], float],
                 precedence: int, right_associative: bool):
        self.fun = fun
        self.pre = precedence
        self.right = right_associative


_opers = {
    '^': _Oper(pow, 4, True),
    '*': _Oper(mul, 3, False),
    '/': _Oper(truediv, 3, False),
    '+': _Oper(add, 2, False),
    '-': _Oper(sub, 2, False),
}


class _Setf(object):
    def __init__(self, index: int):
        self.index = index

    def __call__(self, elem, attrname, value):
        getattr(elem, attrname)[self.index] = value


def _parse(fmt: str, *args: Variable) -> Generator:
    """Split a formula in tokens
    constants, variables, operators, parenthesis
    """
    inpt = fmt.strip()
    while len(inpt) > 0:
        if inpt[0] == '(':
            yield '('
            inpt = inpt[1:].lstrip()
            continue
        if inpt[0] == ')':
            yield ')'
            inpt = inpt[1:].lstrip()
            continue
        if inpt[:2] == '**':
            yield _opers['^']
            inpt = inpt[2:]
            continue
        mtch = re.match(r'(?P<v>(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)', inpt)
        if mtch:
            # number
            yield float(mtch.group('v'))
            inpt = inpt[mtch.span()[1]:].lstrip()
            continue
        mtch = re.match(r'(?P<v>[*+-/^]|\*\*)', inpt)
        if mtch:
            # operator
            yield _opers[mtch.group('v')]
            inpt = inpt[mtch.span()[1]:].lstrip()
            continue
        mtch = re.match(r'\{(?P<v>\d+)}', inpt)
        if mtch:
            # variable
            idx = int(mtch.group('v'))
            if len(args) < idx+1:
                raise ValueError(f"Missing argument {idx}")
            var = args[idx]
            if not isinstance(var, Variable):
                raise ValueError(f"{var!r} is not a Variable")
            yield var
            inpt = inpt[mtch.span()[1]:].lstrip()
            continue
        raise(ValueError(f"Cannot parse {inpt[0]}"))


def shunting_yard(fmt: str, *args: Variable) -> Generator:
    """Convert a formula to RPN"""
    opstack = []
    for token in _parse(fmt, *args):
        if isinstance(token, Number) or isinstance(token, Variable):
            yield token
            continue
        if isinstance(token, _Oper):
            pre1 = token.pre
            while len(opstack) > 0:
                top = opstack[-1]
                if not isinstance(top, _Oper):
                    break
                pre2 = top.pre
                if pre1 > pre2:
                    break
                if (pre1 == pre2) and not token.right:
                    break
                yield opstack.pop()
            opstack.append(token)
        if token == '(':
            opstack.append(token)
        if token == ')':
            try:
                while isinstance(opstack[-1], _Oper):
                    yield opstack.pop()
            except IndexError:
                raise ValueError("Unmatched right parenthesis")
            if len(opstack) == 0 or opstack[-1] != '(':
                raise ValueError("Unmatched right parenthesis")
            opstack.pop()  # discard left parenthesis
            if len(opstack) > 0 and isinstance(opstack[-1], Callable):
                yield opstack.pop()
    while len(opstack) > 0:
        el = opstack.pop()
        if el == '(':
            raise ValueError("Unmatched left parenthesis")
        yield el


class Expression(object):
    def __init__(self, fmt: str, *args: Variable):
        self.rpnstack = list(shunting_yard(fmt, *args))
        self.value = 0.0
        Variable.expressions.add(self)

    def __del__(self):
        Variable.expressions.remove(self)

    def evaluate(self):
        """Evaluate an expression in RPN form"""
        stack = []
        for elem in self.rpnstack:
            if isinstance(elem, Number):
                stack.append(elem)
            elif isinstance(elem, Variable):
                stack.append(elem.last_value)
            else:
                right = stack.pop()
                left = stack.pop()
                stack.append(elem.fun(left, right))
        self.value = stack.pop()


class ElementExpression(object):
    def __init__(self,
                 expr: Expression,
                 refpts: Refpts, attrname: str,
                 index: Optional[int] = None):
        if index is None:
            self._setf = setattr
        else:
            self._setf = _Setf(index)
        self.expr = expr
        self.refpts = refpts
        self.attrname = attrname
        Variable.dependents.add(self)

    def __del__(self):
        Variable.dependents.remove(self)

    def evaluate(self, ring: Lattice):
        value = self.expr.value
        for elem in ring.select(self.refpts):
            self._setf(elem, self.attrname, value)
