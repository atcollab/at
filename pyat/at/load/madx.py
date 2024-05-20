"""Load a lattice from a MADX file (.seq)."""

from __future__ import annotations

__all__ = ["MadxParser", "load_madx"]

import functools
import warnings

# functions known by MADX
from math import pi, e, sqrt, exp, log, log10, sin, cos, tan  # noqa: F401
from math import asin, acos, atan, sinh, cosh, tanh, erf, erfc  # noqa: F401
from os.path import abspath
from typing import Optional
from itertools import chain
from collections.abc import Sequence, Generator

import numpy as np

# constants known by MADX
from scipy.constants import c as clight, hbar as _hb, e as qelect
from scipy.constants import physical_constants as _cst

from .file_input import UnorderedParser, AnyDescr, ElementDescr, SequenceDescr
from ..lattice import Lattice, Particle, elements as elt, tilt_elem
from .utils import split_ignoring_parentheses, protect, restore


# Constants known by MADX
true = True
false = False
twopi = 2 * pi
degrad = 180.0 / pi
raddeg = pi / 180.0
emass = 1.0e-03 * _cst["electron mass energy equivalent in MeV"][0]  # [GeV]
pmass = 1.0e-03 * _cst["proton mass energy equivalent in MeV"][0]  # [GeV]
nmass = 1.0e-03 * _cst["neutron mass energy equivalent in MeV"][0]  # [GeV]
umass = 1.0e-03 * _cst["atomic mass constant energy equivalent in MeV"][0]  # [GeV]
mumass = 1.0e-03 * _cst["muon mass energy equivalent in MeV"][0]  # [GeV]
hbar = _hb / qelect * 1.0e-09  # [GeV.s]
erad = _cst["classical electron radius"][0]  # [m]
prad = erad * emass / pmass  # [m]


def sinc(x: float) -> float:
    return sin(x) / x


# -------------------
#  Utility functions
# -------------------


def _value_arg_parser(parser: UnorderedParser, argstr: str):
    """
    Evaluate a keyword argument of a command and return the pair (key, value)
    """
    return argstr, parser.evaluate(argstr)


def set_tilt(func):
    """Decorator which tilts the decorated AT element"""

    @functools.wraps(func)
    def wrapper(name, *args, tilt=None, **kwargs):
        elems = func(name, *args, **kwargs)
        if tilt is not None:
            tilt_elem(elems[0], tilt)
        return elems

    return wrapper


def polyn(a: Sequence[float]) -> np.ndarray:
    """Convert polynomials from MADX to AT"""

    def ref(n: int, t: float):
        nonlocal f
        v = t / f
        f *= n + 1
        return v

    a = np.ravel(a)
    f = 1.0
    return np.array([ref(n, t) for n, t in enumerate(a)], dtype=float)


# ------------------------------
#  Base class for MAD-X elements
# ------------------------------


class _MadElement(ElementDescr):
    """Description of MADX elements"""

    def __init__(self, *args, at=0.0, frm=None, **kwargs):
        self.at = at
        self.frm = frm
        super().__init__(*args, **kwargs)

    def limits(self, refer):
        """return [entrance, exit] coordinate"""
        half_length = 0.5 * self.length
        offset = self.at + refer * half_length
        if self.frm is not None:
            offset += self.frm.at
        return np.array([-half_length, half_length]) + offset


# ------------------------------
#  MAD-X elements
# ------------------------------


# noinspection PyPep8Naming
class drift(_MadElement):
    @staticmethod
    def convert(name: str, l, **params):  # noqa: E741
        return [elt.Drift(name, l, **params)]


# noinspection PyPep8Naming
class marker(_MadElement):
    @staticmethod
    def convert(name, **params):
        return [elt.Marker(name, **params)]


# noinspection PyPep8Naming
class quadrupole(_MadElement):
    @staticmethod
    @set_tilt
    def convert(name, l, k1=0.0, k1s=None, **params):  # noqa: E741
        if k1s is not None:
            params["PolynomA"] = [0.0, k1s]
        return [elt.Quadrupole(name, l, k1, **params)]


# noinspection PyPep8Naming
class sextupole(_MadElement):
    @staticmethod
    @set_tilt
    def convert(name, l, k2=0.0, k2s=None, **params):  # noqa: E741
        if k2s is not None:
            params["PolynomA"] = [0.0, 0.0, k2s / 2.0]
        return [elt.Sextupole(name, l, k2 / 2.0, **params)]


# noinspection PyPep8Naming
class octupole(_MadElement):
    @staticmethod
    @set_tilt
    def convert(name, l, k3=0.0, k3s=0.0, **params):  # noqa: E741
        poly_b = [0.0, 0.0, 0.0, k3 / 6.0]
        poly_a = [0.0, 0.0, 0.0, k3s / 6.0]
        return [elt.Multipole(name, l, poly_a, poly_b, **params)]


# noinspection PyPep8Naming
class multipole(_MadElement):
    @staticmethod
    @set_tilt
    def convert(name, knl, ksl=(), **params):  # noqa: E741
        params.pop("l", None)
        return [elt.ThinMultipole(name, polyn(ksl), polyn(knl), **params)]


# noinspection PyPep8Naming
class sbend(_MadElement):
    @staticmethod
    @set_tilt
    def convert(
        name,
        l,
        angle,
        e1=0.0,
        e2=0.0,
        k1=0.0,
        k2=None,
        hgap=None,
        fint=0.0,
        **params,  # noqa: E741
    ):  # noqa: E741
        if hgap is not None:
            fintx = params.pop("fintx", fint)
            params.update(FullGap=2.0 * hgap, FringeInt1=fint, FringeInt2=fintx)
        if k2 is not None:
            params["PolynomB"] = [0.0, k1, k2 / 2.0]
        return [
            elt.Dipole(
                name,
                l,
                angle,
                k1,
                EntranceAngle=e1,
                ExitAngle=e2,
                **params,
            )
        ]


# noinspection PyPep8Naming
class rbend(_MadElement):
    @staticmethod
    @set_tilt
    def convert(name, l, angle, e1=0.0, e2=0.0, **params):  # noqa: E741
        hangle = 0.5 * angle
        arclength = l * hangle / sin(hangle)
        return sbend.convert(
            name, arclength, angle, e1=hangle - e1, e2=hangle - e2, **params
        )

    def _length(self):
        hangle = 0.5 * self.angle
        return self["l"] * hangle / sin(hangle)


# noinspection PyPep8Naming
class kicker(_MadElement):
    @staticmethod
    @set_tilt
    def convert(name, l, hkick=0.0, vkick=0.0, **params):  # noqa: E741
        kicks = np.array([hkick, vkick])
        return [elt.Corrector(name, l, kicks, **params)]


# noinspection PyPep8Naming
class hkicker(_MadElement):
    @staticmethod
    @set_tilt
    def convert(name, l, kick=0.0, **params):  # noqa: E741
        return kicker.convert(name, l, hkick=kick, **params)


# noinspection PyPep8Naming
class vkicker(_MadElement):
    @staticmethod
    @set_tilt
    def convert(name, l, kick=0.0, **params):  # noqa: E741
        return kicker.convert(name, l, vkick=kick, **params)


# noinspection PyPep8Naming
class rfcavity(_MadElement):
    @staticmethod
    def convert(
        name, l=0.0, volt=0.0, freq=0.0, lag=0.0, harmon=0, **params  # noqa: E741
    ):
        cavity = elt.RFCavity(
            name,
            l,
            1.0e6 * volt,
            1.0e6 * freq,
            harmon,
            0.0,
            PassMethod="IdentityPass" if l == 0.0 else "DriftPass",
            **params,
        )
        # del cavity.energy
        return [cavity]


# noinspection PyPep8Naming
class monitor(_MadElement):
    @staticmethod
    def convert(name, l=0.0, **params):  # noqa: E741
        if l == 0.0:
            return [elt.Monitor(name, **params)]
        else:
            hl = 0.5 * l
            return [
                elt.Drift(name, hl, source="monitor"),
                elt.Monitor(name, **params),
                elt.Drift(name, hl, source="monitor"),
            ]


# noinspection PyPep8Naming
class hmonitor(monitor):
    pass


# noinspection PyPep8Naming
class vmonitor(monitor):
    pass


class _Ignored(_MadElement):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        type1 = self["source"].title()
        type2 = "Marker" if kwargs.get("l", 0.0) == 0.0 else "Drift"
        print(f"MADX type {type1} is ignored and replaced by a {type2}.")

    @staticmethod
    def convert(name, l=0.0, **params):  # noqa: E741
        if l == 0.0:
            return [elt.Marker(name, **params)]
        else:
            return [elt.Drift(name, l, **params)]


# noinspection PyPep8Naming
class solenoid(_Ignored):
    pass


# noinspection PyPep8Naming
class rfmultipole(_Ignored):
    pass


# noinspection PyPep8Naming
class crabcavity(_Ignored):
    pass


# noinspection PyPep8Naming
class elseparator(_Ignored):
    pass


def value(**kwargs):
    kwargs.pop("copy", False)
    for key, v in kwargs.items():
        print(f"{key}: {v}")


# ------------------------------
#  MAD-X LIST element
# ------------------------------


class _List(SequenceDescr):
    """Descriptor for the MADX LIST"""

    def __add__(self, other):
        return type(self)(chain(self, other))

    def __mul__(self, other):
        def repeat(n):
            for _i in range(n):
                yield from self

        return type(self)(repeat(other))

    def __rmul__(self, other):
        return self.__mul__(other)

    def expand(self, parser: MadxParser) -> Generator[elt.Element, None, None]:
        if self.inverse:
            for elemstr in reversed(self):
                elem = parser.evaluate(elemstr)
                if isinstance(elem, AnyDescr):  # Element or List
                    yield from (-elem).expand(parser)
                elif isinstance(elem, Sequence):  # Other sequence (tuple)
                    for el in reversed(elem):
                        yield from (-el).expand(parser)
        else:
            for elemstr in self:
                elem = parser.evaluate(elemstr)
                if isinstance(elem, AnyDescr):  # Element or List
                    yield from elem.expand(parser)
                elif isinstance(elem, Sequence):  # Other sequence (tuple)
                    for el in elem:
                        yield from el.expand(parser)


# ------------------------------
#  MAD-X SEQUENCE
# ------------------------------


# noinspection PyPep8Naming
class _Sequence(SequenceDescr):
    """Descriptor for the MADX SEQUENCE"""

    _offset = {"entry": 1.0, "centre": 0.0, "exit": -1.0}

    def __init__(
        self,
        *args,
        l=0,
        refer="centre",
        at=0.0,
        frm=None,
        valid=True,
        **kwargs,  # noqa: E741
    ):
        self.l = l  # noqa: E741
        try:
            self.refer = self._offset[refer]
        except KeyError as exc:
            raise ValueError(f"REFER must be in {set(self._offset.keys())}") from exc
        self.at = at
        self.frm = frm
        self.valid = valid
        super().__init__(*args, **kwargs)

    def __call__(self, *args, **kwargs):
        # Never make a copy
        kwargs.pop("copy", None)
        return super().__call__(*args, copy=False, **kwargs)

    def limits(self, refer):
        half_length = 0.5 * self.length
        offset = self.at + refer * half_length
        if self.frm is not None:
            offset += self.frm.at
        return np.array([-half_length, half_length]) + offset

    def expand(self, parser: MadxParser) -> Generator[elt.Element, None, None]:
        def insert_drift(dl, el):
            if dl > 1.0e-12:
                yield from drift(name="filler", l=dl).expand(parser)
            elif dl < -1.0e-12:
                elemtype = type(el).__name__.upper()
                raise ValueError(f"{elemtype}({el.name!r}) is overlapping")

        if not self.valid:
            raise NameError(f"name {self.name!r} is not defined")

        end = 0.0
        elem = None
        for fname, *anames in self:
            try:
                # noinspection PyProtectedMember
                elem = parser._raw_command(
                    None, fname, *anames, no_global=True, copy=True
                )
                entry, ext = elem.limits(self.refer)
                yield from insert_drift(entry - end, elem)
                end = ext
                yield from elem.expand(parser)
            except Exception as exc:
                mess = (f"In sequence {self.name!r}: \"{fname}, {', '.join(anames)}\"",)
                exc.args += mess
                raise
        try:
            yield from insert_drift(self.length - end, elem)
        except Exception as exc:  # Cannot happen if no element in sequence
            mess = (f"Unexpected error in sequence {self.name!r}",)
            exc.args += mess
            raise


class _BeamDescr(ElementDescr):

    def __init__(
        self,
        callfun,
        *args,
        particle="positron",
        mass=emass,  # GeV
        charge=1.0,
        energy=1.0,  # GeV
        bcurrent=0.0,
        **kwargs,
    ):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
            atparticle = Particle(particle, rest_energy=1.0e09 * mass, charge=charge)
        mass = 1.0e-09 * atparticle.rest_energy
        charge = atparticle.charge
        gamma = energy / mass
        beta = sqrt(1.0 - 1.0 / gamma / gamma)
        pc = beta * energy  # GeV
        kwargs.setdefault("gamma", gamma)
        kwargs.setdefault("beta", beta)
        kwargs.setdefault("pc", pc)
        kwargs.setdefault("brho", 1.0e09 * pc / abs(charge) / clight)
        super().__init__(
            *args,
            source="beam",
            atparticle=atparticle,
            mass=mass,
            charge=charge,
            energy=energy,
            bcurrent=bcurrent,
            **kwargs,
        )
        self.callfun = callfun

    @property
    def particle(self):
        return self["atparticle"].name

    def __call__(self, *args, **kwargs):
        return self.callfun(*args, **kwargs)

    @staticmethod
    def convert(name, *args, **params):
        return []

    def expand(self, parser: MadxParser) -> dict:
        res = dict(
            particle=self["atparticle"],
            energy=1.0e09 * self.energy,
            beam_current=self.bcurrent,
            periodicity=1,
            radiate=self.get("radiate", False),
        )
        return res


class MadxParser(UnorderedParser):
    # noinspection PyUnresolvedReferences
    """MAD-X specific parser

    The parser is a subclass of :py:class:`dict` and is database containing all the
    MAD-X variables.

    Example:
        Parse a 1st file:

        >>> parser = at.MadxParser()
        >>> parser.parse_file("file1")

        Parse another file:

        >>> parser.parse_file("file2")

        Get the variable "vkick"

        >>> parser["vkick"]
        0.003

        Define a new variable:

        >>> parser["hkick"] = -0.0024

        Get the "qf1" element

        >>> parser["qf1"]
        quadrupole(name=qf1, l=1.0, k1=0.5, tilt=0.001)

        Generate an AT :py:class:`.Lattice` from the "ring" sequence

        >>> ring = parser.lattice(use="ring")  # generate an AT Lattice
    """

    _soft_eval = {"file", "refer"}
    _argument_parser = {"value": _value_arg_parser, "show": _value_arg_parser}

    def __init__(self):
        """"""
        super().__init__(
            globals(),
            delimiter=";",
            linecomment=("!", "//"),
            blockcomment=("/*", "*/"),
            endfile="return",
            call=self._call_cmd,
            beam=self._beam_cmd,
            sequence=_Sequence,
            centre="centre",
            entry="entry",
            exit="exit",
        )
        self.current_sequence = None
        self._beam_cmd()

    def clear(self):
        super().clear()
        self.current_sequence = None
        self._beam_cmd()

    # noinspection PyUnusedLocal
    def _call_cmd(self, file=None, name=None, copy=False):
        """Implement the CALL MAD-X command"""
        self.parse_files(file)

    def _beam_cmd(
        self,
        *args,
        sequence=None,
        **kwargs,
    ):
        """Implement the BEAM MAD-X statement"""
        if isinstance(sequence, _Sequence):
            name = f"beam%{sequence.name}"
        elif isinstance(sequence, str):
            name = f"beam%{sequence}"
        else:
            name = "beam"
        beam = _BeamDescr(
            self._beam_cmd,
            *args,
            **kwargs,
        )
        self[name] = beam

    def _format_statement(self, line: str) -> str:
        line, matches = protect(line, fence=('"', '"'))
        line = "".join(line.split()).lower()  # Remove all spaces, lower
        line = line.replace("from", "frm")  # Replace from by frm
        line = line.replace("{", "(").replace("}", ")")
        line = line.replace(":=", "=")  # since we evaluate only once
        (line,) = restore(matches, line)
        return line

    def _command(self, label: Optional[str], cmdname: str, *argnames: str, **kwargs):
        res = None
        if self.current_sequence is None:
            try:
                res = self._raw_command(label, cmdname, *argnames, **kwargs)
            except (KeyError, NameError) as exc:
                if cmdname == "sequence":
                    # if sequence creation failed, create a dummy sequence anyway
                    cmdname = label
                    res = self._raw_command(label, cmdname, "valid=False", **kwargs)
                    argnames += ("valid=True",)
                    self.current_sequence = res
                self.delayed.append((label, cmdname, *argnames))
                self.missing.add(self._reason(exc))
            finally:
                if isinstance(res, _Sequence):
                    self.current_sequence = res
        else:
            if cmdname == "endsequence":
                self.current_sequence = None
            else:
                if label is not None:
                    try:
                        res = self._raw_command(label, cmdname, *argnames, **kwargs)
                    except (NameError, KeyError) as exc:
                        self.delayed.append((label, cmdname, *argnames))
                        self.missing.add(self._reason(exc))
                    finally:
                        self.current_sequence.append((label, *argnames))
                else:
                    self.current_sequence.append((cmdname, *argnames))
        return res

    def _assign(self, label: str, key: str, value: str):
        if key == "list":
            members = split_ignoring_parentheses(value[1:-1])
            self[label] = _List(members, name=label)
            return None
        else:
            return super()._assign(label, key, value)

    def _get_beam(self, key: str):
        try:
            beam = self[f"beam%{key}"]
        except KeyError:
            beam = self["beam"]
        return beam

    def lattice(self, use: str = "cell", **kwargs):
        """Create a lattice from the selected sequence

        Parameters:
            use:                Name of the MADX sequence or line containing the desired
              lattice. Default: ``cell``

        Keyword Args:
            name (str):         Name of the lattice. Default: MADX sequence name.
            energy (float):     Energy of the lattice [eV]
            periodicity(int):   Number of periods. Default: 1
            *:                  All other keywords will be set as Lattice attributes
        """

        def gener(params, *args):
            params.setdefault("use", use)
            params.setdefault("name", use)

            # generate the Lattice attributes from the BEAM object
            beam = self._get_beam(use).expand(self)
            rad = beam.pop("radiate")
            if args:
                args[0].update(radiate=rad)
            for key, value in beam.items():
                params.setdefault(key, value)

            # Iterate from the elements
            yield from self.expand(use)

        options = {}
        lat = Lattice(iterator=gener, **kwargs)
        if options.get("radiate", False):
            lat.enable_6d(copy=False)
        return lat


def load_madx(*files: str, use: str = "cell", verbose=False, **kwargs) -> Lattice:
    """Create a :py:class:`.Lattice`  from MAD-X files

    Parameters:
        files:              Names of one or several MAD-X files
        use:                Name of the MADX sequence or line containing the desired
          lattice. Default: ``cell``
        verbose:            Print details on processing

    Keyword Args:
        name (str):         Name of the lattice. Default: MADX sequence name.
        energy (float):     Energy of the lattice [eV]
        periodicity(int):   Number of periods. Default: 1
        *:                  All other keywords will be set as Lattice attributes

    Returns:
        lattice (Lattice):  New :py:class:`.Lattice` object

    See Also:
        :py:func:`.load_lattice` for a generic lattice-loading function.
    """
    parser = MadxParser()
    absfiles = tuple(abspath(file) for file in files)
    kwargs.setdefault("in_file", absfiles)
    parser.parse_files(*absfiles)
    nmiss = len(parser.missing)
    if nmiss > 0:
        if verbose:
            print(f"\nMissing definitions: {parser.missing}\n")
        else:
            print(f'\n{nmiss} missing definitions\nUse "verbose=True" to see them.\n')

    return parser.lattice(use=use)
