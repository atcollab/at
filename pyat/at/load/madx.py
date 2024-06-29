"""Load a lattice from a MADX file (.seq)."""

from __future__ import annotations

__all__ = ["MadxParser", "load_madx"]

import functools
import warnings

# functions known by MAD-X
from math import pi, e, sqrt, exp, log, log10, sin, cos, tan  # noqa: F401
from math import asin, acos, atan, sinh, cosh, tanh, erf, erfc  # noqa: F401
from os.path import abspath
from typing import Optional
from itertools import chain
from collections.abc import Sequence, Generator, Iterable

import numpy as np

# constants known by MAD-X
from scipy.constants import c as clight, hbar as _hb, e as qelect
from scipy.constants import physical_constants as _cst

from . import register_format
from .file_input import UnorderedParser, AnyDescr, ElementDescr, SequenceDescr
from .utils import protect, restore
from ..lattice import Lattice, Particle, Filter, elements as elt, tilt_elem

_default_beam = dict(
    particle="positron",
    energy=1.0,  # GeV
    bcurrent=0.0,
    kbunch=1,
    radiate=False,
)

# Constants known by MAD-X
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
    """Cardinal sine function known by MAD-X"""
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

    def limits(self, parser: MadxParser, offset: float, refer: float):
        half_length = 0.5 * self.length
        offset = offset + refer * half_length + self.at
        if self.frm is not None:
            offset += parser[self.frm].at
        return np.array([-half_length, half_length]) + offset


# ------------------------------
#  MAD-X classes
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
    def convert(name, knl=(), ksl=(), **params):
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
        **params,
    ):
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
            name, arclength, angle, e1=hangle + e1, e2=hangle + e2, **params
        )

    def _length(self):
        hangle = 0.5 * self.angle
        return self["l"] * hangle / sin(hangle)


# noinspection PyPep8Naming
class kicker(_MadElement):
    @staticmethod
    @set_tilt
    def convert(name, l=0.0, hkick=0.0, vkick=0.0, **params):  # noqa: E741
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
        name, l=0.0, volt=0.0, freq=np.nan, lag=0.0, harmon=0, **params  # noqa: E741
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
                elt.Drift(name, hl, madtype="monitor"),
                elt.Monitor(name, **params),
                elt.Drift(name, hl, madtype="monitor"),
            ]


# noinspection PyPep8Naming
class hmonitor(monitor):
    pass


# noinspection PyPep8Naming
class vmonitor(monitor):
    pass


# noinspection PyPep8Naming
class instrument(monitor):
    pass


class _Ignored(_MadElement):

    report = True

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        type1 = self["madtype"].title()
        type2 = "Marker" if self.length == 0.0 else "Drift"
        if self.report:
            self.__class__.report = False
            print(f"MADX type {type1} is ignored and replaced by a {type2}.")

    @staticmethod
    def convert(name, l=0.0, **params):  # noqa: E741
        if l == 0.0:
            return [elt.Marker(name, **params)]
        else:
            return [elt.Drift(name, l, **params)]


# noinspection PyPep8Naming
class solenoid(_Ignored):
    report = True
    pass


# noinspection PyPep8Naming
class rfmultipole(_Ignored):
    report = True
    pass


# noinspection PyPep8Naming
class crabcavity(_Ignored):
    report = True
    pass


# noinspection PyPep8Naming
class elseparator(_Ignored):
    report = True
    pass


# noinspection PyPep8Naming
class collimator(_Ignored):
    report = True
    pass


# noinspection PyPep8Naming
class tkicker(_Ignored):
    report = True
    pass


def value(**kwargs):
    kwargs.pop("copy", False)
    for key, v in kwargs.items():
        print(f"{key}: {v}")


class _Line(SequenceDescr):
    """Descriptor for the MADX LINE"""

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
            for elem in reversed(self):
                if isinstance(elem, AnyDescr):  # Element or List
                    yield from (-elem).expand(parser)
                elif isinstance(elem, Sequence):  # Other sequence (tuple)
                    for el in reversed(elem):
                        yield from (-el).expand(parser)
        else:
            for elem in self:
                if isinstance(elem, AnyDescr):  # Element or List
                    yield from elem.expand(parser)
                elif isinstance(elem, Sequence):  # Other sequence (tuple, list)
                    for el in elem:
                        yield from el.expand(parser)


# noinspection PyPep8Naming
class _Sequence(SequenceDescr):
    """Descriptor for the MADX SEQUENCE"""

    _offset = {"entry": 1.0, "centre": 0.0, "exit": -1.0}

    def __init__(
        self,
        *args,
        l: float = 0,
        refer: str = "centre",
        refpos: Optional[str] = None,
        at: float = 0.0,
        frm: Optional[str] = None,
        valid: int = 1,
        **kwargs,  # noqa: E741
    ):
        self.l = l  # noqa: E741
        try:
            self.refer = self._offset[refer.lower()]
        except KeyError as exc:
            raise ValueError(f"REFER must be in {set(self._offset.keys())}") from exc
        self.refpos = refpos
        self.at = at
        self.frm = frm
        self.valid = bool(valid)
        super().__init__(*args, **kwargs)

    def __call__(self, *args, **kwargs):
        # Never make a copy
        copy = kwargs.pop("copy", None)
        super().__call__(*args, copy=False, **kwargs)
        return self if copy else None

    def origin(self, parser, refer, refpos):
        if refpos is None:
            orig = 0.5 * (refer - 1.0) * self.length
        else:
            orig = None
            for fname, *anames in self:
                if fname == refpos:
                    # noinspection PyProtectedMember
                    elem = parser._raw_command(
                        None, fname, *anames, no_global=True, copy=True
                    )
                    orig = -elem.at
                    break
            if orig is None:
                raise NameError(
                    f"REFPOS {refpos!r} is not in the sequence {self.name!r}"
                )
        if self.frm is not None:
            orig += parser[self.frm].at
        return orig

    def flatten(self, parser, offset: float = 0.0, refer: float = 1.0, refpos=None):
        offset = offset + self.origin(parser, refer, refpos) + self.at
        for fname, *anames in self:
            try:
                # noinspection PyProtectedMember
                elem = parser._raw_command(
                    None, fname, *anames, no_global=True, copy=True
                )
            except Exception as exc:
                mess = (f"In sequence {self.name!r}: \"{fname}, {', '.join(anames)}\"",)
                exc.args += mess
                raise
            if isinstance(elem, _Sequence):
                yield from elem.flatten(parser, offset, self.refer, elem.refpos)
            elif isinstance(elem, _MadElement):
                yield elem.limits(parser, offset, self.refer), elem

    def expand(self, parser: MadxParser) -> Generator[elt.Element, None, None]:
        def insert_drift(dl, el):
            if dl > 1.0e-10:
                yield from drift(name="filler", l=dl).expand(parser)
            elif dl < -1.0e-10:
                elemtype = type(el).__name__.upper()
                raise ValueError(f"{elemtype}({el.name!r}) is overlapping by {-dl} m")

        end = 0.0
        elem = None
        self.at = 0.0
        for (entry, ext), elem in self.flatten(parser):
            yield from insert_drift(entry - end, elem)
            end = ext
            yield from elem.expand(parser)

        try:
            yield from insert_drift(self.length - end, elem)
        except Exception as exc:  # Cannot happen if no element in sequence
            mess = (f"Unexpected error in sequence {self.name!r}",)
            exc.args += mess
            raise


class _BeamDescr(ElementDescr):
    """Descriptor for the MAD-X BEAM object"""

    @staticmethod
    def convert(name, *args, **params):
        return []

    def expand(self, parser: MadxParser) -> dict:

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
            atparticle = Particle(
                self["particle"],
                rest_energy=self.get("mass", emass),
                charge=self.get("charge", 1),
            )
        mass = 1.0e-09 * atparticle.rest_energy  # [GeV]
        charge = atparticle.charge

        # Respect the precedence of MAD-X
        if "energy" in self:
            energy = self["energy"]
            gamma = energy / mass
            # betagamma = sqrt(gamma * gamma - 1.0)
        elif "pc" in self:
            pc = self["pc"]
            betagamma = pc / mass
            gamma = sqrt(betagamma * betagamma + 1.0)
        elif "gamma" in self:
            gamma = self["gamma"]
            # betagamma = sqrt(gamma * gamma - 1.0)
        elif "beta" in self:
            beta = self["beta"]
            gamma = 1.0 / sqrt(1.0 - beta * beta)
            # betagamma = beta * gamma
        elif "brho" in self:
            brho = self["brho"]
            betagamma = 1.0e-9 * brho * clight * abs(charge) / mass
            gamma = sqrt(betagamma * betagamma + 1.0)
        else:
            gamma = 1.0 / mass

        energy = gamma * mass  # [GeV]
        # pc = betagamma * mass  # [GeV]
        # beta = betagamma / gamma
        # brho = 1.0e09 * pc / abs(charge) / clight

        return dict(
            particle=atparticle,
            energy=1.0e09 * energy,  # [eV]
            beam_current=self["bcurrent"],
            nbunch=self["kbunch"],
            periodicity=1,
            radiate=self["radiate"],
        )


class _MadParser(UnorderedParser):
    """Common class for both MAD8 anf MAD-X parsers"""

    _str_arguments = {"file", "refer", "refpos", "sequence", "frm"}
    _argument_parser = {"value": _value_arg_parser, "show": _value_arg_parser}

    def __init__(self, *args, strict: bool = True, **kwargs):
        """"""
        if not strict:
            kwargs.update(none=0.0)
        super().__init__(
            *args,
            delimiter=";",
            linecomment=("!", "//"),
            endfile="return",
            call=self._call_cmd,
            beam=self._beam_cmd,
            sequence=_Sequence,
            centre="centre",
            entry="entry",
            exit="exit",
            **kwargs,
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
        self.parse_files(file, final=False)

    def _beam_cmd(self, sequence=None, **kwargs):
        """Implement the BEAM MAD-X statement"""
        name = "beam%" if sequence is None else f"beam%{sequence}"
        beamobj = self.get(name, None)
        if beamobj is None:
            beamobj = _BeamDescr(_default_beam)
            self[name] = beamobj

        beamobj.update(**kwargs)

    def _command(self, label: Optional[str], cmdname: str, *argnames: str, **kwargs):
        res = None
        if self.current_sequence is None:
            try:
                res = self._raw_command(label, cmdname, *argnames, **kwargs)
            except (KeyError, NameError) as exc:
                if cmdname.lower() == "sequence":
                    # if sequence creation failed, create a dummy sequence anyway
                    res = self._raw_command(label, cmdname, "valid=0", **kwargs)
                    # But store the command for later update
                    reason = self._reason(exc)
                    self._postpone(reason, None, label, "valid=1", *argnames)
                else:
                    raise
            finally:
                if isinstance(res, _Sequence):
                    self.current_sequence = res
        else:
            if cmdname.lower() == "endsequence":
                self.current_sequence = None
            else:
                if label is not None:
                    try:
                        res = self._raw_command(label, cmdname, *argnames, **kwargs)
                    finally:
                        self.current_sequence.append((label, *argnames))
                else:
                    self.current_sequence.append((cmdname, *argnames))
        return res

    def _format_statement(self, line: str) -> str:
        line, matches = protect(line, fence=('"', '"'))
        line = "".join(line.split())  # Remove all spaces
        line = line.replace("{", "(").replace("}", ")")
        line = line.replace(":=", "=")  # since we evaluate only once
        (line,) = restore(matches, line)
        return line

    def _assign(self, label: str, key: str, val: str):
        if key.lower() == "line":
            return label, _Line(self.evaluate(val), name=label)
        else:
            return super()._assign(label, key, val)

    def _finalise(self, final: bool = True) -> None:
        super()._finalise(final=final)
        if final:
            try:
                default_value = self["none"]
                for var in self._missing(verbose=True):
                    self[var] = default_value
                super()._finalise()
            except KeyError:
                pass

    def _get_beam(self, key: str):
        """Get the beam object for a given sequence"""
        try:
            beam = self[f"beam%{key}"]
        except KeyError:
            beam = self["beam%"]
        return beam

    def lattice(self, use: str = "ring", **kwargs):
        """Create a lattice from the selected sequence

        Parameters:
            use:                Name of the MAD sequence or line containing the desired
              lattice. Default: ``ring``

        Keyword Args:
            name (str):         Name of the lattice. Default: MAD sequence name.
            particle(Particle): Circulating particle. Default: from MAD
            energy (float):     Energy of the lattice [eV], Default: from MAD
            periodicity(int):   Number of periods. Default: 1
            *:                  All other keywords will be set as Lattice attributes
        """

        def mad_filter(params: dict, elems: Filter, *args) -> Iterable[elt.Element]:
            def beta() -> float:
                rest_energy = params["particle"].rest_energy
                if rest_energy == 0.0:
                    return 1.0
                else:
                    gamma = float(params["energy"] / rest_energy)
                    return sqrt(1.0 - 1.0 / gamma / gamma)

            use = params.setdefault("use", "ring")
            # generate the Lattice attributes from the BEAM object
            beam = self._get_beam(use).expand(self)
            for key, val in beam.items():
                params.setdefault(key, val)

            cavities = []
            cell_length = 0

            for elem in elems(params, *args):
                if isinstance(elem, elt.RFCavity):
                    cavities.append(elem)
                cell_length += getattr(elem, "Length", 0.0)
                yield elem

            params["_length"] = cell_length
            rev = beta() * clight / cell_length
            for cav in cavities:
                if np.isnan(cav.Frequency):
                    cav.Frequency = rev * cav.HarmNumber
            if cavities:
                cavities.sort(key=lambda el: el.Frequency)
                c0 = cavities[0]
                params["_harmnumber"] = getattr(c0, "HarmNumber", np.nan)

        part = kwargs.get("particle", None)
        if isinstance(part, str):
            kwargs["particle"] = Particle(part)
        lat = Lattice(self._generator, iterator=mad_filter, use=use, **kwargs)
        try:
            radiate = lat.radiate
        except AttributeError:
            radiate = False
        else:
            del lat.radiate
        if radiate:
            lat.enable_6d(copy=False)
        return lat

    @property
    def sequences(self):
        """List of available sequences or lines"""
        return [k for k, v in self.items() if isinstance(v, (_Line, _Sequence))]


class MadxParser(_MadParser):
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

    def __init__(self, **kwargs):
        """"""
        super().__init__(
            globals(), continuation=None, blockcomment=("/*", "*/"), **kwargs
        )

    def evaluate(self, expr):
        """Evaluate an expression using *self* as local namespace"""
        expr = self._no_dot(expr)  # Replace "." by "_", lower case
        expr = expr.replace("->", ".")  # Attribute access
        expr = expr.replace("^", "**")  # Exponentiation
        return super().evaluate(expr)


def load_madx(*files: str, use: str = "ring", strict: bool = True, **kwargs) -> Lattice:
    """Create a :py:class:`.Lattice`  from MAD-X files

    - The *energy* and *particle* of the generated lattice are taken from the MAD-X
      ``BEAM`` object, using the MAD-X default parameters: positrons at 1 Gev.
      These parameters are overloaded by the value given in the *energy* and
      *particle* keyword arguments.
    - The radiation state is given by the ``RADIATE`` flag of the ``BEAM`` object,
      using the AT defaults: RF cavities active, synchrotron radiation in dipoles and
      quadrupoles.
    - Long elements are split according to the default AT value for *NumIntSteps* (10).

    Parameters:
        files:              Names of one or several MAD8 files
        use:                Name of the MAD8 sequence or line containing the desired
          lattice. Default: ``ring``

    Parameters:
        files:              Names of one or several MAD-X files
        strict:             If :py:obj:`False`, assign 0 to undefined variables
        use:                Name of the MADX sequence or line containing the desired
          lattice. Default: ``ring``

    Keyword Args:
        name (str):         Name of the lattice. Default: MADX sequence name.
        particle(Particle): Circulating particle. Default: from MADX
        energy (float):     Energy of the lattice [eV]. Default: from MADX
        periodicity(int):   Number of periods. Default: 1
        *:                  All other keywords will be set as Lattice attributes

    Returns:
        lattice (Lattice):  New :py:class:`.Lattice` object

    See Also:
        :py:func:`.load_lattice` for a generic lattice-loading function.
    """
    parser = MadxParser(strict=strict)
    absfiles = tuple(abspath(file) for file in files)
    kwargs.setdefault("in_file", absfiles)
    parser.parse_files(*absfiles)
    return parser.lattice(use=use, **kwargs)


register_format(
    ".seq", load_madx, descr="MAD-X lattice description. See :py:func:`.load_madx`."
)
