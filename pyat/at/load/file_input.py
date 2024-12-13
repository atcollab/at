"""Generic text parsers for conversion of lattices in different formats to AT"""

from __future__ import annotations

__all__ = [
    "set_argparser",
    "skip_class",
    "ignore_class",
    "skip_names",
    "ignore_names",
    "AnyDescr",
    "ElementDescr",
    "SequenceDescr",
    "BaseParser",
    "UnorderedParser",
    "LowerCaseParser",
    "UpperCaseParser",
    "DictNoDot",
]

from os import getcwd
from os.path import join, abspath, normpath, dirname
import re
from itertools import repeat, count
from collections.abc import Callable, Iterable, Generator, Mapping, Sequence

import numpy as np

from .utils import split_ignoring_parentheses, protect, restore
from ..lattice import Lattice, elements as elt, params_filter

_dot = re.compile(r'("?)(\.?[a-zA-Z_][\w.:]*)\1')  # look for MAD identifiers
_colon = re.compile(r":(?!=)")  # split on :=

# For decoding error messages:
_singlequoted = re.compile(r"'([\w.]*)'")  # look for single-quoted items
_named = re.compile(r"name=([\w.]*)")  # look for 'name=MADid' items


def set_argparser(argparser):
    """Decorator which adds an "argparser" attribute to a function"""

    def decorator(func):
        func.argparser = argparser
        return func

    return decorator


def skip_class(classname: str, baseclass: type[ElementDescr], **kwargs):
    """Generate a class for skipped elements.

    Args:
        classname: Name of the generated class
        baseclass: Base class, must be a subclass of :py:class:`ElementDescr`
        **kwargs: dictionary of additional attributes and methods. See :py:func:`type`.

    Returns:
        cls: Element class, skipped when generating AT elements
    """

    def init(self, *args, **kwargs):
        baseclass.__init__(self, *args, **kwargs)
        # if type(self) not in self.mentioned:
        type1 = self.__class__.__name__
        print(f"Element {self.name} ({type1}) is ignored.")
        self._mentioned.add(type(self))

    kwargs.update(__init__=init)
    return type(classname, (baseclass,), kwargs)


def ignore_class(classname: str, baseclass: type[ElementDescr], **kwargs):
    """Generate a class for ignored elements.

    Args:
        classname: Name of the generated class
        baseclass: Base class, must be a subclass of :py:class:`ElementDescr`
        **kwargs: dictionary of additional attributes and methods. See :py:func:`type`.

    Returns:
        cls: Element class, converted to :py:class:`.Marker` or :py:class:`.Drift`
          when generating AT elements
    """

    def init(self, *args, **kwargs):
        baseclass.__init__(self, *args, **kwargs)
        # if type(self) not in self.mentioned:
        if not args:  # not for replication
            type1 = self.__class__.__name__
            type2 = "Marker" if self.get("l", 0.0) == 0.0 else "Drift"
            print(f"Element {self.name} ({type1}) is replaced by a {type2}.")
            self._mentioned.add(type(self))

    def to_at(self, l=0.0, origin="", **params):  # noqa: E741
        if l == 0.0:
            return [elt.Marker(self.name, origin=origin, **self.meval(params))]
        else:
            return [elt.Drift(self.name, l, origin=origin, **self.meval(params))]

    kwargs.update(__init__=init, to_at=to_at)
    return type(classname, (baseclass,), kwargs)


def ignore_names(
    namespace: dict, baseclass: type[ElementDescr], classnames: Iterable[str]
) -> None:
    """Add classes for ignored elements in the given namespace.

    Args:
        namespace:
        baseclass:  Base class, must be a subclass of :py:class:`ElementDescr`
        classnames: Class names of the generated classes
    """
    namespace.update((nm, ignore_class(nm, baseclass)) for nm in classnames)


def skip_names(
    namespace: dict, baseclass: type[ElementDescr], classnames: Iterable[str]
) -> None:
    """Add classes for ignored elements in the given namespace.

    Args:
        namespace:
        baseclass:  Base class, must be a subclass of :py:class:`ElementDescr`
        classnames: Class names of the generated classes
    """
    namespace.update((nm, skip_class(nm, baseclass)) for nm in classnames)


class DictNoDot(dict):
    @classmethod
    def _defkey(cls, expr: str, quoted: bool) -> str:
        """substitutions to get a valid python identifier"""
        # Using classmethod to allow using super() in subclasses
        return expr.replace(".", "_").replace(":", "_")

    @classmethod
    def _gen_key(cls, expr: str) -> str:
        """Generate a dict key"""
        if expr and expr[0] == '"':
            return cls._defkey(expr[1:-1], True)
        else:
            return cls._defkey(expr, False)

    @classmethod
    def _gen_expr(cls, expr) -> str:
        """Generate a valid python expression"""

        def repl(match):
            return cls._defkey(match.group(2), match.group(1))

        return _dot.sub(repl, expr)

    def __init__(self, *args, verbose: bool = False, **kwargs):
        super().__init__(*args, **kwargs)
        self.verbose = verbose

    def __setitem__(self, key, value):
        super().__setitem__(self._gen_key(key), value)

    def __getitem__(self, key):
        return super().__getitem__(self._gen_key(key))

    def get(self, key, *args):
        return super().get(self._gen_key(key), *args)

    def _print(self, *args, **kwargs):
        if self.verbose:
            print(*args, **kwargs)


class AnyDescr:
    """Base class for source object descriptors"""

    str_attr = ()
    "list of names of str attributes"
    bool_attr = ()
    "list of names of bool attributes"
    pos_args = ()
    "list of names of positional arguments"

    @classmethod
    def argparser(cls, parser, argcount, argstr):
        """Specialised argument parser"""
        return parser._argparser(
            argcount, argstr, bool_attr=cls.bool_attr, str_attr=cls.str_attr
        )

    def __init__(self, *args, **kwargs):
        self.name = kwargs.pop("name", self.__class__.__name__)
        self.inverse = kwargs.pop("inverse", False)
        super().__init__(*args, **kwargs)

    def __neg__(self):
        return self.inverted(copy=True)

    def __call__(self, *args, copy: bool = True, **kwargs) -> AnyDescr | None:
        """Create a copy of the element with updated fields"""
        if copy:
            b = {key: kwargs.pop(key, value) for key, value in vars(self).items()}
            b.update(kwargs)
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
            # noinspection PyUnresolvedReferences
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

    def expand(self, parser: BaseParser) -> Generator[elt.Element, None, None]:
        """Iterator on the generated AT elements"""
        yield from ()


class ElementDescr(AnyDescr, dict):
    """Simple representation of an element as a :py:class:`dict`"""

    _mentioned = set()
    at2mad = {"Length": "L"}
    array_fmt = str.maketrans("[]()", "{}{}")
    bool_fmt = None

    def __init__(self, *args, origin=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.origin = origin or self.__class__.__name__

    def __getattr__(self, item):
        # Allows accessing items using the attribute access syntax
        # tried after lookup for a real attribute failed
        try:
            return self[item]
        except KeyError as exc:
            # necessary for getattr to return the default value for a missing attribute
            name = self.__class__.__name__
            raise AttributeError(f"{name!r} object has no {item!r} attribute") from exc

    def __rmul__(self, other):
        """Element repetition"""
        return list(repeat(self, other))

    def __repr__(self):
        keywords = [f"name={self.name}"]
        keywords += [f"{k}={v!r}" for k, v in self.items()]
        return f"{self.__class__.__name__}({', '.join(keywords)})"

    def __str__(self):
        attrs = [f"{key}={self.attr_format(value)}" for key, value in self.items()]
        return ", ".join([self.__class__.__name__.upper().ljust(10)] + attrs)

    @staticmethod
    def attr_format(value):
        if isinstance(value, bool):
            return ElementDescr.bool_fmt[value]
        elif isinstance(value, np.ndarray):
            return np.array2string(value, separator=", ").translate(
                ElementDescr.array_fmt
            )
        elif isinstance(value, Sequence):
            return str(value).translate(ElementDescr.array_fmt)
        else:
            return str(value)

    @classmethod
    def from_at(cls, kwargs):
        def translate(attributes):
            for at, mad in cls.at2mad.items():
                v = attributes.pop(at, None)
                if v is not None:
                    yield mad, v

        params = {"name": kwargs.pop("FamName", "?")}
        params.update(translate(kwargs))
        return cls(**params)

    def to_at(self, *args, **params) -> list[elt.Element]:
        """Generate the AT element. Must be overloaded for each specific element"""
        return []

    # noinspection PyUnusedLocal
    def expand(self, parser: BaseParser) -> Generator[elt.Element, None, None]:
        """Iterator on the generated AT elements"""
        try:
            elems = self.to_at(**self)
        except Exception as exc:
            exc.args += (f"{self}",)
            raise

        if self.inverse:
            for elem in reversed(elems):
                yield elem.swap_faces(copy=True)
        else:
            yield from elems

    @property
    def length(self) -> float:
        """Element length"""
        return self.get("l", 0.0)

    @staticmethod
    def meval(params: dict):
        """Evaluation of superfluous parameters"""
        # return params
        return {}


class SequenceDescr(AnyDescr, list):
    """Simple representation of a sequence of elements as a :py:class:`list`"""

    def __repr__(self):
        string = super().__repr__()
        return f"{self.__class__.__name__}({string})"

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

    _delimiter: str | None = None
    _continuation: str = "\\"
    _linecomment: str | tuple[str] | None = "#"
    _blockcomment: tuple[str, str] | None = None
    _endfile: str | None = None
    _undef_key: str = "missing"

    def __init__(
        self,
        env: dict,
        *args,
        strict: bool = True,
        always_force: bool = True,
        **kwargs,
    ):
        """
        Args:
            env: global namespace used for evaluating commands
            verbose: If True, print detail on the processing
            strict: If :py:obj:`False`, assign 0 to undefined variables
            *args: dict initializer
            **kwargs: dict initializer
        """
        linecomment = self._linecomment
        blockcomment = self._blockcomment
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
                line = line_comment(line)
                if line:
                    buffer.append(line_comment(line))
                    return False, ""
                else:
                    # Special case to avoid that empty lines break the continuation
                    return False, None

        else:

            def handle_comments(buffer, line, in_comment):
                if in_comment:
                    *rest, line = line.split(sep=endcomment, maxsplit=1)
                    in_comment = len(rest) <= 0
                    return in_comment, "" if in_comment > 0 else line
                else:
                    line = line_comment(line)
                    if line:
                        contents, *rest = line.split(sep=begcomment, maxsplit=1)
                        buffer.append(contents)
                        in_comment = len(rest) > 0
                        return in_comment, rest[0] if in_comment else ""
                    else:
                        # Special case to avoid that empty lines break the continuation
                        return False, None

            begcomment, endcomment = blockcomment

        self.skip_comments = handle_comments
        self.env = env
        self.bases = [getcwd()]
        self.kwargs = kwargs
        self.strict = strict
        self.always_force = always_force
        self.force = always_force

        super().__init__(*args, **kwargs)

        if not strict:
            self[self._undef_key] = 0
        self.postponed = []
        self.in_file = []

    def clear(self):
        """Clear the database: remove all parameters and objects"""
        super().clear()
        self.update(self.kwargs)
        if not self.strict:
            self[self._undef_key] = 0
        self.postponed = []
        self.in_file = []

    def _format_command(self, expr: str) -> str:
        """Format a command for evaluation

        Overload this method for specific languages"""
        return expr

    def evaluate(self, expr: str):
        """Evaluate the right side of an expression

        Args:
            expr: expression to evaluate

        Returns:
            value: evaluated expression
        """
        expr = self._format_command(self._gen_expr(expr))
        default_value = self.get(self._undef_key)
        if self.force and default_value is not None:
            for _loop in range(5):
                try:
                    return eval(expr, self.env, self)
                except NameError as exc:
                    var = self._reason(exc)
                    self._print(f"Set {var!r} to {default_value} ({_loop})")
                    self[var] = default_value
        return eval(expr, self.env, self)

    def _eval_cmd(self, cmdname: str, no_global: bool = False) -> Callable:
        """Evaluate a command"""
        cmd: Callable | None = self.get(cmdname, None)
        if cmd is not None:
            return cmd
        elif no_global:
            raise TypeError(f"{cmdname!r} is not allowed in this context")
        else:
            return self.env[self._gen_key(cmdname)]

    @staticmethod
    def _reason(exc: Exception) -> str | None:
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
            names = _singlequoted.search(exc.args[0])
            return names[-1]
        else:
            return None

    def _argparser(
        self,
        argcount: int,
        argstr: str,
        *,
        bool_attr: tuple[str] = (),
        str_attr: tuple[str] = (),
        pos_args: tuple[str] = (),
    ):
        """Evaluate the value of a command argument and return the pair (key, value)"""

        def arg_value(k, v):
            if k in str_attr:
                return v[1:-1] if v[0] == '"' else v
            else:
                return self.evaluate(v)

        key, *value = split_ignoring_parentheses(
            argstr, delimiter="=", fence=('"', '"'), maxsplit=1
        )
        if value:  # Keyword argument
            return key, arg_value(key, value[0])
        else:
            ok = argstr[0] != "-"
            key = argstr if ok else argstr[1:]
            if key in bool_attr:  # boolean flag
                return key, ok
            else:  # positional parameter
                try:
                    key = pos_args[argcount]
                except IndexError:
                    print(f"Unexpected positional argument '{argstr}' ignored")
                    return None
                return key, arg_value(key, argstr)

    def _assign(self, label: str, key: str, value: str):
        """Variable assignment"""
        return key, self.evaluate(value)

    def _raw_command(
        self,
        label: str | None,
        cmdname: str,
        *args: str,
        no_global: bool = False,
        **kwargs,
    ):
        """Command execution"""
        # Get the command
        cmd = self._eval_cmd(cmdname, no_global=no_global)
        # Evaluate the arguments
        argp = getattr(cmd, "argparser", None)
        if argp:
            ags = (argp(self, n, arg) for n, arg in enumerate(args) if arg)
        else:
            ags = (self._argparser(n, arg) for n, arg in enumerate(args) if arg)
        kwargs.update(arg for arg in ags if arg is not None)
        if label is None:
            kwargs.setdefault("copy", False)
        else:
            kwargs.setdefault("name", label.replace('"', ""))
        # Execute the command
        return cmd(**kwargs)

    def _command(self, *args: str, **kwargs):
        return self._raw_command(*args, **kwargs)

    def _decode(self, label: str | None, cmdname: str, *args: str) -> None:
        """Execute the split statement"""
        left, *right = cmdname.split("=")
        try:
            if right:
                label, result = self._assign(label, left, right[0])
            else:
                result = self._command(label, cmdname, *args)
            if not (label is None or result is None):
                self[label] = result
        except (KeyError, NameError) as exc:  # store the failing assignment
            self._fallback(exc, label, cmdname, *args)

    def _fallback(self, exc: Exception, lbl: str | None, cmd: str, *args: str) -> None:
        """Store failing commands in self.postponed for later evaluation"""
        self.postponed.append((self._reason(exc), lbl, cmd, *args))

    def _format_statement(self, line: str) -> str:
        """Reformat the input line"""
        return line.replace(" ", "")  # Remove all spaces

    def _statement(self, line: str) -> bool:
        # protect quoted items. Make sure placeholder cannot be modified
        line, match1 = protect(line, fence=('"', '"'), placeholder="_0_")

        if self._endfile is not None and line.startswith(self._endfile):
            return False

        *left, right = _colon.split(line, maxsplit=1)
        right = self._format_statement(right)

        # protect nested parentheses
        cpt = count()
        rpt = []
        while "(" in right:
            plh = str(next(cpt)).join(("_par", "_"))
            right, match = protect(right, fence=("\\(", "\\)"), placeholder=plh)
            rpt.append(match)

        cmdargs = right.split(",")

        # restore parentheses
        for match in reversed(rpt):
            cmdargs = restore(match, *cmdargs)

        # restore quoted items
        b = restore(match1, *left, *cmdargs)

        # handle label
        if left:
            id = 1
            label = b[0].replace(" ", "")
        else:
            id = 0
            label = None

        # ignore MAD qualifiers
        while b[id] in ["const", "int", "real"]:
            id += 1

        # process statement
        self._decode(label, *b[id:])
        return True

    def _finalise(self, final: bool = True) -> None:
        """Called at the end of processing"""
        if final:
            undefined = self._missing(verbose=self.verbose)
            self._print(f"{len(undefined)} missing definitions.")

    @property
    def sequences(self):
        """List of available sequences or lines"""
        return [k for k, v in self.items() if isinstance(v, SequenceDescr)]

    @staticmethod
    def _command_str(label: str | None, cmdname: str, *argnames: str):
        string = ", ".join((cmdname, *argnames))
        if label is not None:
            string = " : ".join((label, string))
        return string

    def _lookup(self, item: str):
        """Search for an object in the pending statements"""
        for reason, label, *args in self.postponed:
            if label is not None and self._gen_key(label) == item:
                return (reason, label, *args)
        return None

    def _missing(self, verbose=False):
        """Return the set of missing definitions"""
        miss = set()
        for cmd in self.postponed:
            reason = cmd[0]
            if reason == self._gen_key(cmd[2]):
                if verbose:
                    cmdstr = self._command_str(*cmd[1:])
                    self._print(f"Command {cmdstr!r} ignored")
                continue
            while cmd is not None:
                reason = cmd[0]
                cmd = self._lookup(reason)
            miss.add(reason)
        return miss

    missing = property(_missing, doc="Set of missing definitions")

    def _analyse(self, key: str) -> None:
        """Print the chain of failing commands"""
        if isinstance(key, str):
            reason = self._gen_key(key)
            cmd = self._lookup(reason)
            while cmd is not None:
                reason, *cmdargs = cmd
                cmdstr = self._command_str(*cmdargs)
                print(f"\n{key!r} depends on {reason!r}: {cmdstr!r}")
                key = reason
                cmd = self._lookup(reason)
            print(f"\n{reason!r} is not defined\n")

    @property
    def ignored(self):
        """Set of ignored commands"""
        ignd = set()
        for cmd in self.postponed:
            if cmd[0] == self._gen_key(cmd[2]):
                ignd.add(cmd[2])
        return ignd

    def expand(self, key: str) -> Generator[elt.Element, None, None]:
        """iterator over AT objects generated by a source object"""
        self.force = True
        try:
            v = self[key]
            if isinstance(v, AnyDescr):
                yield from v.expand(self)
            else:
                yield v
        except Exception as exc:
            strargs = (arg for arg in exc.args if isinstance(arg, str))
            print(f"{type(exc).__name__}: {': '.join(strargs)}")
            self._analyse(self._reason(exc))
            raise
        finally:
            self.force = self.always_force

    def _generator(self, params):
        """Generate AT elements for the Lattice constructor"""
        use = params.setdefault("use", "ring")
        params.setdefault("name", use)
        params.setdefault("energy", 1.0e9)
        params.setdefault("periodicity", 1)

        # Iterate from the elements
        yield from self.expand(use)

    def lattice(self, use="ring", **kwargs):
        """Create a lattice from the selected sequence

        Parameters:
            use:                Name of the sequence or line containing the desired
              lattice. Default: ``ring``

        Keyword Args:
            name (str):         Name of the lattice. Default: sequence name.
            particle(Particle): Circulating particle. Default: Particle("relativistic")
            energy (float):     Energy of the lattice [eV]. Default: 1.0 GeV
            periodicity(int):   Number of periods. Default: 1
            *:                  All other keywords will be set as Lattice attributes
        """
        return Lattice(
            self._generator,
            iterator=params_filter,
            in_file=self.in_file or None,
            use=use,
            **kwargs,
        )

    def parse_lines(
        self,
        lines: Iterable[str],
        *,
        final: bool = True,
        **kwargs,
    ) -> None:
        """Process input lines and fill the database

        Args:
            lines: Iterable of input lines
            final: If :py:obj:`True`, signals that the undefined variables may be set
              to the default value
            **kwargs:   Initial variable definitions
        """
        self.update(**kwargs)
        buffer = []
        in_comment: bool = False
        ok: bool = True
        for line_number, contents in enumerate(lines):
            # Handle comments
            while contents:
                in_comment, contents = self.skip_comments(buffer, contents, in_comment)
            if contents is None:
                # Special case to avoid that empty lines break the continuation
                continue

            if buffer:
                contents = "".join(buffer).strip()
                buffer = []
                if not contents:
                    continue
            else:
                continue

            # Handle delimiters
            if self._delimiter is None:
                statements = []
                last = contents
            else:
                *statements, last = contents.split(sep=self._delimiter)

            # Handle continuation
            if self._continuation is None:
                buffer.append(last)
            else:
                idc = last.find(self._continuation)
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
                    print(message)
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
        prolog: None | int | Callable[..., None] = None,
        epilog: Callable[..., None] | None = None,
        **kwargs,
    ) -> None:
        """Process files and fill the database

        Args:
            *filenames: Files to process
            final: If :py:obj:`True`, signals that the undefined variables may be set
              to the default value
            prolog:
            epilog:
            **kwargs:   Initial variable definitions
        """
        self.update(**kwargs)
        filenames = tuple(abspath(file) for file in filenames)
        self.in_file.extend(filenames)
        last = len(filenames) - 1
        ElementDescr._mentioned.clear()
        for nf, fn in enumerate(filenames):
            fn = normpath(join(self.bases[-1], fn))
            self.bases.append(dirname(fn))
            print("Processing", fn)
            try:
                with open(fn) as f:
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
    """Parser allowing definitions in any order

    This is done by storing the failed statements in a queue and iteratively trying
    to execute them after all input statements have been processed, until the number
    of failures is constant (hopefully zero)
    """

    def __init__(self, env: dict, *args, **kwargs):
        """
        Args:
            env: global namespace
            verbose: If True, print detail on the processing
            *args: dict initializer
            **kwargs: dict initializer
        """
        super().__init__(env, *args, always_force=False, **kwargs)

    def _finalise(self, final: bool = True) -> None:
        """Loop on evaluation of the pending statements"""

        def replay():
            nend = len(self.postponed)
            while nend > 0:
                self._print(f"Delayed evaluation of {nend} statements")
                nstart = nend
                statements = self.postponed
                self.postponed = []
                for _reason, *args in statements:
                    self._decode(*args)
                nend = len(self.postponed)
                if nend == nstart:
                    break

        # at the end of each file: try again previously failing commands
        replay()

        # After the last file: initialize the remaining undefined variables
        default_value = self.get(self._undef_key)
        if final:
            undefined = self._missing(verbose=self.verbose)
            self._print(f"{len(undefined)} missing definitions.")
            if undefined and default_value is not None:
                for var in undefined:
                    self._print(f"Set {var} to {default_value}")
                    self[var] = default_value
                # last trial
                replay()


class LowerCaseParser(BaseParser):
    """Case independent parser"""

    @classmethod
    def _defkey(cls, expr: str, quoted: bool) -> str:
        """substitutions to get a valid python identifier"""
        expr = super()._defkey(expr, quoted)
        return expr if quoted else expr.lower()

    def _format_statement(self, line: str) -> str:
        """Reformat the input line"""
        return super()._format_statement(line.lower())


class UpperCaseParser(BaseParser):
    """Case independent parser"""

    @classmethod
    def _defkey(cls, expr: str, quoted: bool) -> str:
        """substitutions to get a valid python identifier"""
        expr = super()._defkey(expr, quoted)
        return expr if quoted else expr.upper()

    def _format_statement(self, line: str) -> str:
        """Reformat the input line"""
        return super()._format_statement(line.upper())
