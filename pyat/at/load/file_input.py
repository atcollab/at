"""Generic text parsers for conversion of lattices in different formats to AT"""

from __future__ import annotations

__all__ = ["AnyDescr", "ElementDescr", "SequenceDescr", "BaseParser", "UnorderedParser"]

from os import getcwd
from os.path import join, normpath, dirname
import re
import abc
from collections.abc import Sequence, Callable, Iterable, Generator, Mapping
from typing import Union, Optional

from .utils import split_ignoring_parentheses, protect, restore
from ..lattice import Element


def _clean_expr(expr):
    """Replace "." by "_" """

    def repl(match):
        return match.group().replace(".", "_")

    expr, matches = protect(expr, fence=('"', '"'))
    expr = re.sub(r"[a-z][\w.]*", repl, expr)  # Replace "." by "_"
    expr = expr.replace("->", ".")  # for attribute access
    (expr,) = restore(matches, expr)
    return expr


def _default_arg_parser(parser: BaseParser, argstr: str):
    """
    Evaluate a keyword argument of a command and return the pair (key, value)
    """
    try:
        key, value = split_ignoring_parentheses(argstr, delimiter="=")
    except ValueError:  # Positional argument -> boolean flag
        if argstr.startswith("-"):
            v = False
            key = argstr[1:]
        else:
            key = argstr
            v = True
    else:  # Keyword argument
        try:
            v = parser.evaluate(value)
        except (NameError, SyntaxError):
            # noinspection PyProtectedMember
            if key in parser._soft_eval:
                v = value
            else:
                raise
    return key, v


class AnyDescr(abc.ABC):
    """Base class for source object descriptors"""

    def __init__(self, *args, **kwargs):
        self.name = kwargs.pop("name", self.__class__.__name__)
        self.inverse = kwargs.pop("inverse", False)
        super().__init__(*args, **kwargs)

    def __neg__(self):
        return self.inversed(copy=True)

    def __call__(self, *args, copy: bool = True, **kwargs) -> AnyDescr:
        """Create a copy of the element with updated fields"""
        if copy:
            b = dict((key, kwargs.pop(key, value)) for key, value in vars(self).items())
            return type(self)(self, *args, **b, **kwargs)
        else:
            self.update(*args, **kwargs)
            return self

    def update(self, *args, **kwargs):
        for key, value in vars(self).items():
            setattr(self, key, kwargs.pop(key, value))
        # if isinstance(self, ElementDescr):
        if isinstance(self, Mapping):
            super().update(*args, **kwargs)
        else:
            for key, value in kwargs:
                setattr(self, key, value)

    def inversed(self, copy=False):
        """Return a reversed element or line"""
        instance = self() if copy else self
        instance.inverse = not self.inverse
        return instance

    @abc.abstractmethod
    def expand(self, parser: BaseParser) -> Generator[Element, None, None]:
        """Iterator on the possibly inverted sequence of elements"""
        pass


class ElementDescr(AnyDescr, dict):
    """Simple representation of an element as a dict"""

    def __init__(self, *args, **kwargs):
        kwargs.pop("copy", False)
        source = kwargs.pop("source", self.__class__.__name__)
        super().__init__(*args, source=source, **kwargs)

    def __getattr__(self, item):
        # Allows accessing items using the attribute access syntax
        return self[item]

    def __repr__(self):
        descr = super().copy()
        cls = descr.pop("source")
        keywords = [f"name={self.name}"]
        keywords += [f"{k}={v!r}" for k, v in descr.items()]
        return f"{cls}({', '.join(keywords)})"

    @staticmethod
    def convert(name: str, *args, **params) -> list[Element]:
        """Generate the AT element, Most be overloaded for each specific element"""
        return []

    # noinspection PyUnusedLocal
    def expand(self, parser: BaseParser) -> Generator[Element, None, None]:
        """Iterator on the possibly inverted sequence of elements"""
        try:
            elems = self.convert(self.name, **self)
        except Exception as exc:
            exc.args += (f"{self}",)
            raise

        if self.inverse:
            for elem in reversed(elems):
                yield elem.swap_faces(copy=True)
        else:
            yield from elems

    def _length(self) -> float:
        """Element length"""
        return self.get("l", 0.0)

    length = property(_length)


class SequenceDescr(AnyDescr, list, abc.ABC):
    """Simple representation of a sequence of elements as a list"""

    @property
    def length(self) -> float:
        return getattr(self, "l", 0.0)


class BaseParser(dict):
    """Generic file parser

    Analyses files with the following MAD-like format:

    ``variable = value``

    ``label : command [,attribute=value] [,attribute=value]...``

    The parser builds a database of all the defined objects
    """

    _soft_eval = set()
    _argument_parser = {}

    def __init__(
        self,
        env: dict,
        *args,
        delimiter: Optional[str] = None,
        continuation: str = "\\",
        linecomment: Union[str, Sequence[str], None] = "#",
        blockcomment: Optional[tuple[str, str]] = None,
        endfile: Optional[str] = None,
        **kwargs,
    ):
        """
        Args:
            env: global namespace used for evaluating commands
            delimiter: command delimiter
            continuation: command continuation character
            linecomment: Line comment character
            blockcomment: Block comment delimiter
            endfile: "End of input" marker
            *args: dict initializer
            **kwargs: dict initializer
        """
        if not (linecomment is None or isinstance(linecomment, Sequence)):
            linecomment = (linecomment,)
        self.delimiter = delimiter
        self.continuation = continuation
        self.linecomment = linecomment
        if blockcomment is None:
            self.begcomment = self.endcomment = None
        else:
            self.begcomment, self.endcomment = blockcomment
        self.endfile = endfile
        self.env = env
        self.bases = [getcwd()]
        self.kwargs = kwargs
        super().__init__(*args, **kwargs)

    def clear(self):
        """Clean the database"""
        super().clear()
        self.update(self.kwargs)

    # noinspection PyUnusedLocal
    def evaluate(self, item, no_global: bool = False):
        """Evaluate an expression using *self* as local namespace"""
        return eval(_clean_expr(item), self.env, self)

    def _eval_cmd(self, cmdname: str, no_global: bool = False):
        """Evaluate a command"""
        cname = _clean_expr(cmdname)
        cmd = self.get(cname, None)
        if cmd is not None:
            return cmd
        else:
            cmd = self.env.get(cname, None)
        if cmd is None:
            raise KeyError(cname)
        elif no_global:
            raise TypeError(f"{cmdname!r} is not allowed in this context")
        else:
            return cmd

    @staticmethod
    def _reason(exc: Exception) -> str:
        """Extract the element name from the exception"""
        if isinstance(exc, KeyError):
            return exc.args[0]
        elif isinstance(exc, NameError):
            return re.search(r"'(\w*)'", exc.args[0])[1]
        elif isinstance(exc, TypeError):
            idx = re.search(r"name=(\w*)", exc.args[-1])
            if idx is None:
                idx = re.search(r"'(\w*)'", exc.args[0])
            if idx is None:
                return "TypeError"
            else:
                return idx[1]
        else:
            return type(exc).__name__

    def _assign(self, label: str, key: str, value: str):
        """Variable assignment"""
        return self.evaluate(value)

    def _raw_command(
        self,
        label: Optional[str],
        cmdname: str,
        *argnames: str,
        no_global: bool = False,
        **kwargs,
    ):
        """Command execution"""
        argparser = self._argument_parser.get(cmdname, _default_arg_parser)
        kwargs.update(argparser(self, arg) for arg in argnames)
        if label is None:
            kwargs.setdefault("copy", False)
        else:
            kwargs.setdefault("name", label)
        func = self._eval_cmd(cmdname, no_global=no_global)
        return func(**kwargs)

    def _command(self, *args, **kwargs):
        return self._raw_command(*args, **kwargs)

    def _format_statement(self, line: str) -> str:
        """Reformat the input line

        Overload this method for specific languages"""
        line, matches = protect(line, fence=('"', '"'))  # protect the quoted parts
        line = "".join(line.split()).lower()  # Remove all spaces, lower
        (line,) = restore(matches, line)
        return line

    def _statement(self, line: str) -> bool:
        """Split a statement in 'label: command'"""
        fmtline = self._format_statement(line)
        if self.endfile is not None and fmtline.startswith(self.endfile):
            return False
        label, *cmd = fmtline.split(":", maxsplit=1)
        if cmd:  # label
            fmtline = cmd[0]
        else:  # no label
            label = None

        arguments = split_ignoring_parentheses(fmtline)
        self._decode(label, *arguments)
        return True

    def _decode(self, label: str, cmdname: str, *argnames: str) -> None:
        """Execute the split statement"""
        key, *vals = cmdname.split("=")
        if vals:
            result = self._assign(label, key, vals[0])
        else:
            key = label
            result = self._command(label, cmdname, *argnames)
        if not (key is None or result is None):
            self[_clean_expr(key)] = result

    def _finalise(self) -> None:
        """Called at the end of processing"""
        pass

    def _analyse(self, key: str) -> None:
        """Print info on failed statements"""
        print()
        print(f"{key!r} is not defined\n")

    def expand(self, key: str) -> Generator[Element, None, None]:
        """iterator over AT objects generated by a source object"""
        try:
            v = self[key]
            if isinstance(v, AnyDescr):
                yield from v.expand(self)
            else:
                yield v
        except Exception as exc:
            print(f"{type(exc).__name__}: {exc.args[0]}")
            for arg in exc.args[1:]:
                print(arg)
            self._analyse(self._reason(exc))
            raise

    def parse_lines(
        self,
        lines: Iterable[str],
    ) -> None:
        """Process input lines and fill the database

        Args:
            lines: Iterable of input lines
        """
        buffer = []
        in_comment = False
        ok = True
        for line_number, contents in enumerate(lines):

            # Handle block comments
            if in_comment:
                parts = contents.split(sep=self.endcomment, maxsplit=1)
                if len(parts) == 1:
                    in_comment = True
                    continue
                else:
                    in_comment = False
                    contents = parts[1]

            # Handle comments
            if self.linecomment is not None:
                for linecom in self.linecomment:
                    contents, *cmt = contents.split(sep=linecom, maxsplit=1)
            if self.begcomment is not None:
                contents, *cmt = contents.split(sep=self.begcomment, maxsplit=1)
                if len(cmt) > 0:
                    parts = cmt[0].split(sep=self.endcomment, maxsplit=1)
                    if len(parts) == 1:
                        in_comment = True
                        buffer.append(contents)
                        continue
                    else:
                        in_comment = False
                        contents += parts[1]

            # Handle continuations
            c1 = contents.rstrip()
            if c1.endswith(self.continuation):
                buffer.append(c1[: -len(self.continuation)])
                continue
            else:
                buffer.append(contents)
            contents = "".join(buffer)
            buffer = []

            # Split line
            if self.delimiter is not None:
                cmd = contents.split(sep=self.delimiter)
                commands = cmd[:-1]
                buffer.append(cmd[-1])
                if not commands:
                    continue
            else:
                commands = [contents]

            # Process commands
            for cmd in commands:
                cmd = cmd.strip()
                if not cmd:
                    continue
                try:
                    ok = self._statement(cmd)
                except Exception as exc:
                    message = f"Line {line_number} {cmd!r}, {exc}"
                    raise type(exc)(message) from exc
                if not ok:
                    break

            if not ok:
                break

        self._finalise()

    def parse_files(
        self,
        *filenames: str,
        prolog: Union[None, int, Callable[..., None]] = None,
        epilog: Optional[Callable[..., None]] = None,
    ) -> None:
        """Process files and fill the database

        Args:
            *filenames: Files to process
            prolog:
            epilog:
        """
        for fn in filenames:
            fn = normpath(join(self.bases[-1], fn))
            self.bases.append(dirname(fn))
            print("Processing", fn)
            try:
                with open(fn, "rt") as f:
                    if callable(prolog):
                        prolog(f)
                    elif isinstance(prolog, int):
                        for _i in range(prolog):
                            next(f)

                    self.parse_lines(f)

                    if callable(epilog):
                        epilog(f)
            finally:
                print("End", fn)
                self.bases.pop()


class UnorderedParser(BaseParser):
    """parser allowing definitions in any order

    This is done by storing the failed statements in a queue and iteratively trying
    to execute them after all imput statements have been processed, until the number
    of failures is constant (hopefully zero)
    """

    def __init__(self, *args, **kwargs):
        """
        Args:
            env: global namespace
            delimiter: command delimiter
            continuation: command continuation character
            linecomment: Line comment character
            blockcomment: Block comment delimiter
            endfile: End of input marker
            *args: dict initializer
            **kwargs: dict initializer
        """
        super().__init__(*args, **kwargs)
        self.delayed = []
        self.missing = set()

    def clear(self):
        super().clear()
        self.delayed = []
        self.missing = set()

    def _assign(self, label: str, key: str, value: str):
        """Store the failing assignments in self.delayed for later evaluation"""
        try:
            result = self.evaluate(value)
        except NameError as exc:  # store the failing assignment
            self.delayed.append((label, "=".join((key, value))))
            self.missing.add(self._reason(exc))
            result = None
        return result

    def _command(self, label: Optional[str], cmdname: str, *argnames: str, **kwargs):
        """Store failing commands in self.delayed for later evaluation"""
        try:
            result = super()._command(label, cmdname, *argnames, **kwargs)
        except (KeyError, NameError) as exc:  # store the failing assignment
            # NameError for evaluate, KeyError for eval_cmd
            self.delayed.append((label, cmdname, *argnames))
            self.missing.add(self._reason(exc))
            result = None
        return result

    def _finalise(self) -> None:
        """Loop on evaluation of the pending statements"""
        nend = len(self.delayed)
        while nend > 0:
            statements = self.delayed
            self.delayed = []
            self.missing = set()
            nstart = nend
            for args in statements:
                self._decode(*args)
            nend = len(self.delayed)
            if nend == nstart:
                break

    def _analyse(self, key: str) -> None:
        """Print the chain of failing commands"""

        def lookup(item):
            """Search for an object in the pending statements"""
            for label, command, *args in self.delayed:
                if label == item:
                    return (label, command, *args)
                elif label is None and command == item:
                    return (label, command, *args)
            return None

        def decode(label, cmdname, *argnames):
            """Try evaluating the command to identify the error"""
            _, *vals = cmdname.split("=")
            if vals:
                self.evaluate(vals[0])
            else:
                self._raw_command(label, cmdname, *argnames)

        cmd = lookup(key)
        print()
        while cmd is not None:
            try:
                decode(*cmd)
            except Exception as exc:
                newkey = self._reason(exc)
            else:  # Should not happen since the command is in the pending list
                newkey = "NoException"
            print(f"{key!r} depends on {newkey!r}: \"{cmd[0]} : {', '.join(cmd[1:])}\"")
            key = newkey
            cmd = lookup(key)
        print(f"{key!r} is not defined\n")
