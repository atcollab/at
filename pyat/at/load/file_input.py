"""Generic text parsers for conversion of lattices in different formats to AT"""

from __future__ import annotations

__all__ = [
    "AnyDescr",
    "ElementDescr",
    "SequenceDescr",
    "BaseParser",
    "UnorderedParser",
    "DictNoDot",
]

from os import getcwd
from os.path import join, normpath, dirname
import re
import abc
from itertools import repeat
from collections.abc import Callable, Iterable, Generator, Mapping
from typing import Union, Optional

from .utils import split_ignoring_parentheses, protect, restore
from ..lattice import Element, Lattice, params_filter

_dot = re.compile(r"\.?[a-z][\w.]*")  # An identifier: starts with a letter
_singlequoted = re.compile(r"'([\w.]*)'")
_named = re.compile(r"name=([\w.]*)")
_doublequoted = re.compile(r'^"(.*)"$')


def _default_arg_parser(parser: BaseParser, argstr: str):
    """
    Evaluate a keyword argument of a command and return the pair (key, value)
    """
    try:
        key, value = split_ignoring_parentheses(argstr, delimiter="=")
    except ValueError:  # Positional argument -> boolean flag
        argstr = argstr.lower()
        if argstr.startswith("-"):
            v = False
            key = argstr[1:]
        else:
            key = argstr
            v = True
    else:  # Keyword argument
        key = key.lower().replace("from", "frm")
        # noinspection PyProtectedMember
        if key in parser._str_arguments:
            v = _doublequoted.sub(r"\1", value)
        else:
            v = parser.evaluate(value)
    return key, v


class DictNoDot(dict):

    @staticmethod
    def _no_dot(expr):
        def repl(match):
            return match.group().replace(".", "_")

        return _dot.sub(repl, expr.lower())

    def __setitem__(self, key, value):
        super().__setitem__(self._no_dot(key), value)

    def __getitem__(self, key):
        return super().__getitem__(self._no_dot(key))

    def get(self, key, *args):
        return super().get(self._no_dot(key), *args)


class AnyDescr(abc.ABC):
    """Base class for source object descriptors"""

    def __init__(self, *args, **kwargs):
        self.name = kwargs.pop("name", self.__class__.__name__)
        self.inverse = kwargs.pop("inverse", False)
        # kwargs.setdefault("madclass", self.__class__.__name__)
        super().__init__(*args, **kwargs)

    def __neg__(self):
        return self.inverted(copy=True)

    def __call__(self, *args, copy: bool = True, **kwargs) -> Optional[AnyDescr]:
        """Create a copy of the element with updated fields"""
        if copy:
            b = dict((key, kwargs.pop(key, value)) for key, value in vars(self).items())
            b.update(kwargs)
            # b.update(madclass=self.name)
            return type(self)(self, *args, **b)
        else:
            self.update(*args, **kwargs)
            return None

    def update(self, *args, **kwargs):
        # Update attributes
        for key, value in vars(self).items():
            setattr(self, key, kwargs.pop(key, value))
        if isinstance(self, Mapping):
            # Update mapping
            super().update(*args, **kwargs)
        else:
            # Add new attributes
            for key, value in kwargs:
                setattr(self, key, value)

    def inverted(self, copy=False):
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
        kwargs.setdefault("madtype", self.__class__.__name__)
        super().__init__(*args, **kwargs)

    def __getattr__(self, item):
        # Allows accessing items using the attribute access syntax
        return self[item]

    def __rmul__(self, other):
        """Element repetition"""
        return list(repeat(self, other))

    def __repr__(self):
        descr = super().copy()
        cls = descr.pop("madtype")
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

    def __repr__(self):
        str = super().__repr__()
        return f"{self.__class__.__name__}({str})"

    @property
    def length(self) -> float:
        return getattr(self, "l", 0.0)


class BaseParser(DictNoDot):
    """Generic file parser

    Analyses files with the following MAD-like format:

    ``variable = value``

    ``label : command [,attribute=value] [,attribute=value]...``

    The parser builds a database of all the defined objects
    """

    _str_arguments = set()
    _argument_parser = {}

    def __init__(
        self,
        env: dict,
        *args,
        delimiter: Optional[str] = None,
        continuation: str = "\\",
        linecomment: Union[str, tuple[str], None] = "#",
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
        if isinstance(linecomment, tuple):

            def line_comment(line):
                for linecom in linecomment:
                    line, *_ = line.split(sep=linecom, maxsplit=1)
                return line

        else:
            if linecomment is None:

                def line_comment(line):
                    return line

            else:

                def line_comment(line):
                    line, *_ = line.split(sep=linecomment, maxsplit=1)
                    return line

        if blockcomment is None:
            # noinspection PyUnusedLocal
            def handle_comments(buffer, line, in_comment):
                buffer.append(line_comment(line))
                return False, ""

        else:

            def handle_comments(buffer, line, in_comment):
                if in_comment:
                    *rest, line = line.split(sep=endcomment, maxsplit=1)
                    in_comment = len(rest) <= 0
                    return in_comment, "" if in_comment > 0 else line
                else:
                    line = line_comment(line)
                    contents, *rest = line.split(sep=begcomment, maxsplit=1)
                    buffer.append(contents)
                    in_comment = len(rest) > 0
                    return in_comment, rest[0] if in_comment else ""

            begcomment, endcomment = blockcomment

        self.skip_comments = handle_comments
        self.delimiter = delimiter
        self.continuation = continuation
        self.endfile = endfile
        self.env = env
        self.bases = [getcwd()]
        self.kwargs = kwargs

        super().__init__(*args, **kwargs)

    def clear(self):
        """Clean the database"""
        super().clear()
        self.update(self.kwargs)

    def evaluate(self, expr):
        """Evaluate an expression using *self* as local namespace"""
        return eval(expr, self.env, self)

    def _eval_cmd(self, cmdname: str, no_global: bool = False) -> Callable:
        """Evaluate a command"""
        cmd: Optional[Callable] = self.get(cmdname, None)
        if cmd is not None:
            return cmd
        else:
            cmdname = cmdname.lower()
            cmd = self.env.get(cmdname, None)
        if cmd is None:
            raise KeyError(cmdname)
        elif no_global:
            raise TypeError(f"{cmdname!r} is not allowed in this context")
        else:
            return cmd

    @staticmethod
    def _reason(exc: Exception) -> str:
        """Extract the element name from the exception"""
        if isinstance(exc, KeyError):  # Undefined element, attribute
            return exc.args[0]
        elif isinstance(exc, NameError):  # refpos missing
            return _singlequoted.search(exc.args[0])[1]
        elif isinstance(exc, TypeError):
            idx = _named.search(exc.args[-1])  # Missing pos. arg.
            if idx is None:
                idx = _singlequoted.search(exc.args[0])  # Not allowed in seq.
            if idx is None:
                return "TypeError"
            else:
                return idx[1]
        elif isinstance(exc, ValueError):  # overlap
            print(exc.args[0])
            return _singlequoted.search(exc.args[0])[1]
        else:
            return type(exc).__name__

    def _argparser(self, command: str, *args: str):
        argparser = self._argument_parser.get(command.lower(), _default_arg_parser)
        return (argparser(self, arg) for arg in args)

    def _assign(self, label: str, key: str, value: str):
        """Variable assignment"""
        return key, self.evaluate(value)

    def _raw_command(
        self,
        label: Optional[str],
        cmdname: str,
        *argnames: str,
        no_global: bool = False,
        **kwargs,
    ):
        """Command execution"""
        func = self._eval_cmd(cmdname, no_global=no_global)
        kwargs.update(self._argparser(cmdname, *argnames))
        if label is None:
            kwargs.setdefault("copy", False)
        else:
            kwargs.setdefault("name", label)
        return func(**kwargs)

    def _command(self, *args, **kwargs):
        return self._raw_command(*args, **kwargs)

    def _format_statement(self, line: str) -> str:
        """Reformat the input line

        Overload this method for specific languages"""
        line, matches = protect(line, fence=('"', '"'))  # protect the quoted parts
        line = "".join(line.split())  # Remove all spaces
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
        left, *right = cmdname.split("=")
        if right:
            label, result = self._assign(label, left, right[0])
        else:
            result = self._command(label, cmdname, *argnames)
        if not (label is None or result is None):
            self[label] = result

    def _finalise(self, final: bool = True) -> None:
        """Called at the end of processing"""
        pass

    @staticmethod
    def _command_str(_, label, cmdname, *argnames):
        string = ", ".join((cmdname, *argnames))
        if label is not None:
            string = " : ".join((label, string))
        return string

    def _analyse(self, key: str) -> None:
        """Print info on failed statements"""
        print(f"\n{key!r} is not defined\n")

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

    def _generator(self, params):
        """Generate AT elements for the Lattice constructor"""
        use = params.setdefault("use", "ring")
        params.setdefault("name", use)

        # Iterate from the elements
        yield from self.expand(use)

    def lattice(self, use="ring", **kwargs):
        """Create a lattice from the selected sequence

        Parameters:
            use:                Name of the MADX sequence or line containing the desired
              lattice. Default: ``ring``

        Keyword Args:
            name (str):         Name of the lattice. Default: MADX sequence name.
            particle(Particle): Circulating particle. Default: from MADX
            energy (float):     Energy of the lattice [eV], Default: from MADX
            periodicity(int):   Number of periods. Default: 1
            *:                  All other keywords will be set as Lattice attributes
        """
        return Lattice(self._generator, iterator=params_filter, use=use, **kwargs)

    def parse_lines(
        self,
        lines: Iterable[str],
        final: bool = True,
    ) -> None:
        """Process input lines and fill the database

        Args:
            lines: Iterable of input lines
            final: If :py:obj:`True`, signals that the undefined variables may be set
              to the defalt value
        """
        buffer = []
        in_comment: bool = False
        ok: bool = True
        for line_number, contents in enumerate(lines):

            # Handle comments
            while contents:
                in_comment, contents = self.skip_comments(buffer, contents, in_comment)

            if buffer:
                contents = "".join(buffer).strip()
                buffer = []
                if not contents:
                    continue
            else:
                continue

            # print(repr(contents))
            # Handle delimiters
            if self.delimiter is None:
                statements = []
                last = contents
            else:
                *statements, last = contents.split(sep=self.delimiter)

            # Handle continuation
            if self.continuation is None:
                buffer.append(last)
            else:
                idc = last.find(self.continuation)
                if idc >= 0:
                    buffer.append(last[:idc])
                else:
                    statements.append(last)

            # Process statements
            for stmnt in statements:
                stmnt = stmnt.strip()
                if not stmnt:
                    continue
                try:
                    ok = self._statement(stmnt)
                except Exception as exc:
                    message = f"Line {line_number} {stmnt!r}, {exc}"
                    raise type(exc)(message) from exc
                if not ok:
                    break

            if not ok:
                break

        self._finalise(final=final)

    def parse_files(
        self,
        *filenames: str,
        final: bool = True,
        prolog: Union[None, int, Callable[..., None]] = None,
        epilog: Optional[Callable[..., None]] = None,
    ) -> None:
        """Process files and fill the database

        Args:
            *filenames: Files to process
            final: If :py:obj:`True`, signals that the undefined variables may be set
              to the default value
            prolog:
            epilog:
        """
        last = len(filenames) - 1
        for nf, fn in enumerate(filenames):
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

                    self.parse_lines(f, final=final and (nf == last))

                    if callable(epilog):
                        epilog(f)
            finally:
                print("End", fn)
                self.bases.pop()


class UnorderedParser(BaseParser):
    """parser allowing definitions in any order

    This is done by storing the failed statements in a queue and iteratively trying
    to execute them after all input statements have been processed, until the number
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

    def clear(self):
        super().clear()
        self.delayed = []

    def _postpone(self, reason, label, cmdname, *argnames):
        """Store failing commands in self.delayed for later evaluation"""
        self.delayed.append((reason, label, cmdname, *argnames))

    def _decode(self, label: str, cmdname: str, *argnames: str) -> None:
        """Postpone failing commands"""
        try:
            super()._decode(label, cmdname, *argnames)
        except (KeyError, NameError) as exc:  # store the failing assignment
            self._postpone(self._reason(exc), label, cmdname, *argnames)

    def _finalise(self, final: bool = True) -> None:
        """Loop on evaluation of the pending statements"""
        nend = len(self.delayed)
        while nend > 0:
            statements = self.delayed
            self.delayed = []
            nstart = nend
            for reason, *args in statements:
                self._decode(*args)
            nend = len(self.delayed)
            if nend == nstart:
                break

    def _lookup(self, item: str):
        """Search for an object in the pending statements"""
        for reason, label, *args in self.delayed:
            if label is not None and label.lower() == item:
                return (reason, label, *args)
        return None

    def _missing(self, verbose: bool = False):
        miss = set()
        for cmd in self.delayed:
            reason = cmd[0]
            if reason == cmd[2].lower():
                if verbose:
                    print(
                        f"Unknown command {cmd[2]!r} ignored: "
                        f"{self._command_str(*cmd)!r}"
                    )
                continue
            while cmd is not None:
                reason = cmd[0]
                cmd = self._lookup(reason)
            miss.add(reason)
        return miss

    def _analyse(self, key: str) -> None:
        """Print the chain of failing commands"""
        cmd = self._lookup(key.lower())
        reason = key
        while cmd is not None:
            reason = cmd[0]
            print(f"\n{key!r} depends on {reason!r}: {self._command_str(*cmd)!r}")
            key = reason
            cmd = self._lookup(reason)
        print(f"\n{reason!r} is not defined\n")

    missing = property(_missing, doc="Set of missing definitions")
