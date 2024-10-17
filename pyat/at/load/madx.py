"""Load a lattice from a MADX file (.seq)."""

from __future__ import annotations

__all__ = ["MadxParser", "load_madx"]

import functools
import warnings

# functions known by MAD-X
from math import pi, e, sqrt, exp, log, log10, sin, cos, tan  # noqa: F401
from math import asin, acos, atan, sinh, cosh, tanh, erf, erfc  # noqa: F401
from os.path import abspath
from itertools import chain
from collections.abc import Sequence, Generator, Iterable
import re

import numpy as np

# constants known by MAD-X
from scipy.constants import c as clight, hbar as _hb, e as qelect
from scipy.constants import physical_constants as _cst

from . import register_format
from .file_input import AnyDescr, ElementDescr, SequenceDescr, BaseParser
from .file_input import CaseIndependentParser, UnorderedParser, DeferredParser
from .file_input import set_argparser, ignore_names
from ..lattice import Lattice, Particle, Filter, elements as elt, tilt_elem

_kconst = re.compile("^ *const +")
_kreal = re.compile("^ *real +")
_kint = re.compile("^ *int +")
_separator = re.compile(r"(?<=[\w.)])\s+(?=[\w.(])")

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
    if abs(x) < 1e-10:
        return x
    else:
        return sin(x) / x


# -------------------
#  Utility functions
# -------------------


def set_tilt(func):
    """Decorator which tilts the decorated AT element"""

    @functools.wraps(func)
    def wrapper(name, *args, tilt=None, **kwargs):
        elems = func(name, *args, **kwargs)
        if tilt is not None:
            tilt_elem(elems[0], tilt)
        return elems

    return wrapper


def polyn(a: float | Sequence[float]) -> np.ndarray:
    """Convert polynomials from MADX to AT"""

    def ref(n: int, t: float):
        nonlocal f
        v = t / f
        f *= n + 1
        return v

    f = 1.0
    if not isinstance(a, Sequence):
        # In case of a single element, we have a scalar instead of a tuple
        a = (a,)
    return np.array([ref(n, t) for n, t in enumerate(a)], dtype=float)


# noinspection PyUnusedLocal
def _keyparser(parser, argcount, argstr):
    """Return the pair key, value for the given 'key' argument"""
    return argstr, parser._evaluate(argstr)


# ------------------------------
#  Base class for MAD-X elements
# ------------------------------


class _MadElement(ElementDescr):
    """Description of MADX elements"""

    str_attr = {"refer", "refpos", "sequence", "from"}

    def __init__(self, *args, at=0.0, **kwargs):
        self.at = at
        # Cannot use "from" as argument or attribute
        setattr(self, "from", kwargs.pop("from", None))
        # kwargs.pop("copy", False)
        super().__init__(*args, **kwargs)

    def limits(self, parser: MadxParser, offset: float, refer: float):
        half_length = 0.5 * self.length
        offset = offset + refer * half_length + self.at
        frm = getattr(self, "from")
        if frm is not None:
            offset += parser[frm].at
        return np.array([-half_length, half_length]) + offset


# ------------------------------
#  MAD-X classes
# ------------------------------


# noinspection PyPep8Naming
class drift(_MadElement):
    @staticmethod
    def convert(name: str, l=0.0, **params):  # noqa: E741
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
        l,  # noqa: E741
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
        arclength = l * sinc(hangle)
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
        kicks = np.array([hkick, vkick], dtype=float)
        return [elt.Corrector(name, l, kicks, **params)]


# noinspection PyPep8Naming
class hkicker(_MadElement):
    @staticmethod
    @set_tilt
    def convert(name, l=0.0, kick=0.0, **params):  # noqa: E741
        return kicker.convert(name, l=l, hkick=kick, **params)


# noinspection PyPep8Naming
class vkicker(_MadElement):
    @staticmethod
    @set_tilt
    def convert(name, l=0.0, kick=0.0, **params):  # noqa: E741
        return kicker.convert(name, l=l, vkick=kick, **params)


# noinspection PyPep8Naming
class rfcavity(_MadElement):
    @staticmethod
    def convert(
        name,
        l=0.0,  # noqa: E741
        volt=0.0,
        freq=np.nan,
        lag=0.0,
        harmon=0,
        **params,
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
                elt.Drift(name, hl, origin="monitor"),
                elt.Monitor(name, **params),
                elt.Drift(name, hl, origin="monitor"),
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


ignore_names(
    globals(),
    _MadElement,
    ["solenoid", "rfmultipole", "crabcavity", "elseparator", "collimator", "tkicker"],
)


@set_argparser(_keyparser)
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

    def expand(self, parser: BaseParser) -> Generator[elt.Element, None, None]:
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
        l: float = 0,  # noqa: E741
        refer: str = "centre",
        refpos: str | None = None,
        at: float = 0.0,
        valid: int = 1,
        **kwargs,
    ):
        self.l = l  # noqa: E741
        try:
            self.refer = self._offset[refer]
        except KeyError as exc:
            raise ValueError(f"REFER must be in {set(self._offset.keys())}") from exc
        self.refpos = refpos
        self.at = at
        # Cannot use "from" as argument or attribute name:
        setattr(self, "from", kwargs.pop("from", None))
        self.valid = bool(valid)
        super().__init__(*args, **kwargs)

    def __call__(self, *args, copy: bool = True, **kwargs):
        # Never make a copy
        super().__call__(*args, copy=False, **kwargs)
        return self if copy else None

    def reference(self, parser, refer, refpos):
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
        frm = getattr(self, "from")
        if frm is not None:
            orig += parser[frm].at
        return orig

    def flatten(self, parser, offset: float = 0.0, refer: float = 1.0, refpos=None):
        offset = offset + self.reference(parser, refer, refpos) + self.at
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
        elem = self
        self.at = 0.0
        for (entry, ext), elem in self.flatten(parser):
            yield from insert_drift(entry - end, elem)
            end = ext
            yield from elem.expand(parser)

        yield from insert_drift(self.length - end, elem)  # Final drift


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
            energy = float(self["energy"])  # force evaluation
            gamma = energy / mass
            # betagamma = sqrt(gamma * gamma - 1.0)
        elif "pc" in self:
            pc = float(self["pc"])  # force evaluation
            betagamma = pc / mass
            gamma = sqrt(betagamma * betagamma + 1.0)
        elif "gamma" in self:
            gamma = float(self["gamma"])  # force evaluation
            # betagamma = sqrt(gamma * gamma - 1.0)
        elif "beta" in self:
            beta = float(self["beta"])  # force evaluation
            gamma = 1.0 / sqrt(1.0 - beta * beta)
            # betagamma = beta * gamma
        elif "brho" in self:
            brho = float(self["brho"])  # force evaluation
            betagamma = 1.0e-9 * brho * clight * abs(charge) / mass
            gamma = sqrt(betagamma * betagamma + 1.0)
        else:
            gamma = 1.0 / mass

        # convert from array scalar to float
        energy = float(gamma * mass)  # [GeV]
        # pc = betagamma * mass  # [GeV]
        # beta = betagamma / gamma
        # brho = 1.0e09 * pc / abs(charge) / clight

        return {
            "particle": atparticle,
            "energy": 1.0e09 * energy,  # [eV]
            "beam_current": float(self["bcurrent"]),  # force deferred evaluation
            "nbunch": self["kbunch"],
            "periodicity": 1,
            "radiate": self["radiate"],
        }


class _Call:
    """Implement the CALL Mad command"""

    @staticmethod
    def argparser(parser, argcount, argstr):
        return parser._argparser(
            argcount, argstr, pos_args=("file",), str_attr=("file",)
        )

    def __init__(self, parser: _MadParser):
        self.parser = parser

    def __call__(self, file=None, name=None, copy=False):
        self.parser.parse_files(file, final=False)


class _Beam:
    """Implement the BEAM Mad command"""

    default_beam = {
        "particle": "positron",
        "energy": 1.0,  # GeV
        "bcurrent": 0.0,
        "kbunch": 1,
        "radiate": False,
    }

    @staticmethod
    def argparser(parser, argcount, argstr):
        return parser._argparser(
            argcount,
            argstr,
            str_attr=("particle", "sequence"),
            bool_attr=("radiate", "bunched"),
        )

    def __init__(self, parser: _MadParser):
        self.parser = parser

    def __call__(self, sequence=None, **kwargs):
        """create a :py:class:`_BeamDescr` object and store it as 'beam%sequence'"""
        name = "beam%" if sequence is None else f"beam%{sequence}"
        beamobj = self.parser.get(name, None)
        if beamobj is None:
            beamobj = _BeamDescr(self.default_beam)
            self.parser[name] = beamobj

        for k, v in kwargs.items():
            beamobj[k] = v


class _MadParser(CaseIndependentParser, DeferredParser, UnorderedParser):
    """Common class for both MAD8 anf MAD-X parsers"""

    _str_arguments = {"file", "refer", "refpos", "sequence", "from"}

    def __init__(self, env: dict, **kwargs):
        """Common behaviour for MAD-X and MAD8

        Args:
            env: global namespace used for evaluating commands
            verbose:    If :py:obj:`True`, print details on the processing
            strict: If :py:obj:`False`, assign 0 to undefined variables
            **kwargs: Initial variable definitions
        """
        super().__init__(
            env,
            delimiter=";",
            linecomment=("!", "//"),
            endfile="return",
            call=_Call(self),
            beam=_Beam(self),
            sequence=_Sequence,
            centre="centre",
            entry="entry",
            exit="exit",
            undef_key="none",
            **kwargs,
        )
        self.current_sequence = None
        self["beam"]()

    def clear(self):
        super().clear()
        self.current_sequence = None
        self["beam"]()

    def _assign(self, label: str, key: str, val: str):
        # Special treatment of "line=(...)" assignments
        if key == "line":
            val = val.replace(")", ",)")  # For tuples with a single item
            return label, _Line(self._evaluate(val), name=label)
        else:
            return super()._assign(label, key, val)

    def _command(self, label: str | None, cmdname: str, *argnames: str, **kwargs):
        # Special treatment of SEQUENCE definitions
        res = None
        if cmdname == "endsequence":
            self.current_sequence = None
            return None
        if self.current_sequence is None:
            try:
                res = super()._command(label, cmdname, *argnames, **kwargs)
            except (KeyError, NameError) as exc:
                if cmdname == "sequence":
                    # if sequence creation failed, create a dummy sequence anyway
                    res = super()._command(label, cmdname, "valid=0", **kwargs)
                    # But store the command for later update
                    self._fallback(exc, None, label, "valid=1", *argnames)
                else:
                    raise
            finally:
                if isinstance(res, _Sequence):
                    self.current_sequence = res
        else:
            if label is None:
                self.current_sequence.append((cmdname, *argnames))
            else:
                try:
                    res = super()._command(label, cmdname, *argnames, **kwargs)
                finally:
                    self.current_sequence.append((label, *argnames))
        return res

    def _format_statement(self, line: str) -> str:
        line = super()._format_statement(line)
        # Remove the MAD const, real, int qualifiers
        line = _kint.sub("", _kreal.sub("", _kconst.sub("", line)))
        # Accept space as separator (after removing qualifiers)
        line = _separator.sub(",", line)
        # turn curly braces into parentheses (MAD arrays)
        line = line.replace("{", "(").replace("}", ")")
        return line

    def _get_beam(self, key: str):
        """Get the beam object for a given sequence"""
        try:
            beam = self[f"beam%{self._gen_key(key)}"]
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

            # generate the Lattice attributes from the BEAM object
            beam = self._get_beam(params["use"]).expand(self)
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

            # Set the frequency of cavities in which it is not specified
            for cav in cavities:
                if np.isnan(cav.Frequency):
                    cav.Frequency = rev * cav.HarmNumber

            # Set the lattice harmonic number
            if cavities:
                cavities.sort(key=lambda el: el.Frequency)
                c0 = cavities[0]
                params["_cell_harmnumber"] = getattr(c0, "HarmNumber", np.nan)

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
        # noinspection PyUnboundLocalVariable
        if radiate:
            lat.enable_6d(copy=False)
        return lat


class MadxParser(_MadParser):
    # noinspection PyUnresolvedReferences
    """MAD-X specific parser

    The parser is a subclass of :py:class:`dict` and is database containing all the
    MAD-X parameters and objects.

    Example:
        Parse a 1st file:

        >>> parser = at.MadxParser()
        >>> parser.parse_file("file1")

        Parse another file. This adds the new definitions to the database:

        >>> parser.parse_file("file2")

        Look at the *rf_on* variable

        >>> parser["rf_on"]
        0

        Get the list of available sequences/lines:

        >>> parser.sequences
        ['arca',
         'arca_inj',
         'arcb_inj',
         'low_emit_ring',
         'arc_inj',
         'low_emit_ring_inj']

        Generate an AT :py:class:`.Lattice` from the *low_emit_ring* sequence

        >>> ring1 = parser.lattice(use="low_emit_ring")

        Change the value of a variable:

        >>> parser["rf_on"] = 1

        Generate a new AT :py:class:`.Lattice` taking into account the new variables:

        >>> ring2 = parser.lattice(use="low_emit_ring")

        Generate an AT :py:class:`.Lattice` from another sequence:

        >>> arca = parser.lattice(use="ring")
    """

    def __init__(self, **kwargs):
        """
        Args:
            strict:     If :py:obj:`False`, assign 0 to undefined variables
            verbose:    If :py:obj:`True`, print details on the processing
            **kwargs:   Initial variable definitions
        """
        super().__init__(
            globals(),
            continuation=None,
            blockcomment=("/*", "*/"),
            **kwargs,
        )

    def _format_command(self, expr: str) -> str:
        """Format a command for evaluation"""
        expr = expr.replace("->", ".")  # Attribute access: VAR->ATTR
        expr = expr.replace("^", "**")  # Exponentiation
        return super()._format_command(expr)


def load_madx(
    *files: str, use: str = "ring", strict: bool = True, verbose: bool = False, **kwargs
) -> Lattice:
    """Create a :py:class:`.Lattice` from MAD-X files

    - The *energy* and *particle* of the generated lattice are taken from the MAD-X
      ``BEAM`` object, using the MAD-X default parameters: positrons at 1 Gev.
      These parameters are overloaded by the value given in the *energy* and
      *particle* keyword arguments.
    - The radiation state is given by the ``RADIATE`` flag of the ``BEAM`` object,
      using the AT defaults: RF cavities active, synchrotron radiation in dipoles and
      quadrupoles.
    - Long elements are split according to the default AT value of *NumIntSteps* (10).

    Parameters:
        files:              Names of one or several MAD-X files
        strict:             If :py:obj:`False`, assign 0 to undefined variables
        use:                Name of the MADX sequence or line containing the desired
          lattice. Default: ``ring``
        verbose:            If :py:obj:`True`, print details on the processing

    Keyword Args:
        name (str):         Name of the lattice. Default: MADX sequence name.
        particle(Particle): Circulating particle. Default: from MADX
        energy (float):     Energy of the lattice [eV]. Default: from MADX
        periodicity(int):   Number of periods. Default: 1
        *:                  Other keywords will be used as initial variable definitions

    Returns:
        lattice (Lattice):  New :py:class:`.Lattice` object

    See Also:
        :py:func:`.load_lattice` for a generic lattice-loading function.
    """
    parser = MadxParser(strict=strict, verbose=verbose)
    absfiles = tuple(abspath(file) for file in files)
    params = {
        key: kwargs.pop(key)
        for key in ("name", "particle", "energy", "periodicity")
        if key in kwargs
    }
    parser.parse_files(*absfiles, **kwargs)
    return parser.lattice(use=use, in_file=absfiles, **params)


register_format(
    ".seq", load_madx, descr="MAD-X lattice description. See :py:func:`.load_madx`."
)
