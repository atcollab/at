# noinspection PyUnresolvedReferences
r"""Using `MAD-X`_ files with PyAT
==================================

PyAT can read lattice descriptions in Mad-X format (.seq files), and can export
lattices in MAD-X format.

However, because of intrinsic differences between PyAT and MAD-X, some
incompatibilities must be taken into account.

1. Translation losses
---------------------

6-D motion
^^^^^^^^^^^^^^
While AT allows setting 6-D motion and synchrotron radiation on individual elements,
MAD has a global setting for the whole lattice. When reading a MAD-X lattice without
radiation, 6-D motion in the resulting AT lattice is turned off, including in RF
cavities. If the MAD-X lattice has radiation on, 6-D motion is activated on the AT
lattice according to default settings (RF cavities active, radiation in dipoles,
quadrupoles and wigglers).

Combined function magnets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

MAD has little support for combined function magnets. When exporting a lattice in MAD
format, the main field component for each magnet class in kept but other components
disappear, except in the few cases handled by MAD (quadrupole and sextupole
components in ``SBEND`` or ``RBEND``, skew quad component in ``QUADRUPOLE``, see the
MAD user's reference manual for more).

Multipoles
^^^^^^^^^^^^^^

MAD has no thick multipoles. Multipoles in MAD format (``MULTIPOLE`` element) are
interpreted as :py:class:`.ThinMultipole` elements. In the other direction, an AT
:py:class:`.Multipole` is converted to a thin ``MULTIPOLE`` surrounded by two
drift spaces.

MAD elements absent from AT
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Many MAD elements have no equivalent in AT. They are replaced by
:py:class:`.Marker` or :py:class:`.Drift` elements, depending on their length.

Incompatible attributes
^^^^^^^^^^^^^^^^^^^^^^^^^^^
Some AT element attributes have no MAD equivalent, and vice versa.

When exporting to a MAD-X file:

- `NumIntSteps`, `MaxOrder` are discarded,
- `FringeBendEntrance`, `FringeBendExit`, `FringeQuadEntrance`, `FringeQuadExit` are
  discarded,
- `R1`, `R2`, `T1`, `T2` are discarded,
- `RApertures`, `EApertures` are discarded.

When reading a MAD-X file:

- `TILT` is interpreted and converted to `R1` and `R2` attributes,
- `KTAP` is interpreted and converted to `FieldScaling`.

.. _using-mad-x-files:

2. Reading MAD-X files
----------------------

Short way
^^^^^^^^^
The :py:func:`load_madx` function created a :py:class:`.Lattice` from a ``LINE`` or
``SEQUENCE`` in a MAD-X file.

>>> ring = at.load_madx("lattice.seq", use="PS")

>>> print(ring)
Lattice(<1295 elements>, name='PS', energy=1000000000.0, particle=Particle('positron'),
periodicity=1, harmonic_number=0, beam_current=0.0, nbunch=1, use='PS')

Detailed way
^^^^^^^^^^^^
:py:class:`MadxParser` creates an empty database which can be populated with the
elements of a MAD-X file.

>>> parser = at.MadxParser()

The :py:meth:`MadxParser.parse_files` method reads one or several MAD-X files and
populates the parser

>>> parser.parse_files("lattice.seq", use="PS")

The parser can be examined and modified using the standard python syntax:

>>> parser["pr_bht91_f"]
sbend(name=PR.BHT91.F, l=2.1975925, angle='angle.prbhf', k1='k1prbhf+k1prpfwf-k1prf8l')

>>> parser["angle.prbhf"]
0.03135884818

>>> parser["angle.prbhf"] = 0.032

All MAD parameters can be interactively modified and their last value will be taken
into account when generating a PyAT lattice.

The :py:meth:`MadxParser.lattice` method creates a :py:class:`.Lattice` from a ``LINE``
or ``SEQUENCE`` of the parser:

>>> ring = parser.lattice(use="ps")

>>> print(ring)
Lattice(<1295 elements>, name='PS', energy=1000000000.0, particle=Particle('positron'),
periodicity=1, harmonic_number=0, beam_current=0.0, nbunch=1, use='PS')

3. Exporting to MAD-X files
---------------------------
Exporting a PyAT lattice to a MAD-X files produces a single MAD ``SEQUENCE`` of
``LINE``.

See :py:func:`save_madx` for usage.

4. Functions and classes
------------------------

.. _mad-x: https://mad.web.cern.ch/mad/webguide/manual.html
"""

from __future__ import annotations

__all__ = ["MadParameter", "MadxParser", "load_madx", "save_madx"]

import functools
import warnings

# functions known by MAD-X
from math import pi, e, sqrt, exp, log, log10, sin, cos, tan  # noqa: F401
from math import asin, acos, atan, sinh, cosh, tanh, erf, erfc  # noqa: F401
from itertools import chain
from collections.abc import Sequence, Generator, Iterable
import re

import numpy as np

# constants known by MAD-X
from scipy.constants import c as clight, hbar as _hb, e as qelect
from scipy.constants import physical_constants as _cst

from .allfiles import register_format
from .utils import split_ignoring_parentheses, protect, restore
from .file_input import AnyDescr, ElementDescr, SequenceDescr, BaseParser
from .file_input import LowerCaseParser, UnorderedParser
from .file_input import set_argparser, ignore_names
from .file_output import Exporter
from ..lattice import Lattice, Particle, elements as elt, tilt_elem

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


class MadParameter(str):
    """MAD parameter

    A MAD parameter is an expression which can be evaluated n the context
    of a MAD parser
    """

    def __new__(cls, parser, expr):
        return super().__new__(cls, expr)

    # noinspection PyUnusedLocal
    def __init__(self, parser: _MadParser, expr: str):
        """Args:
            parser: MadParser instance defining the context for evaluation
            expr:   expression to be evaluated

        The expression may contain MAD parameter names, arithmetic operators and
        mathematical functions known by MAD
        """
        self.parser = parser

    def __float__(self):
        return float(self.parser.evaluate(self))

    def __int__(self):
        return int(self.parser.evaluate(self))

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

    def evaluate(self):
        return self.parser.evaluate(self)


def sinc(x: float) -> float:
    """Cardinal sine function known by MAD-X"""
    if abs(x) < 1e-10:
        return x
    else:
        return sin(x) / x


# -------------------
#  Utility functions
# -------------------


def mad_element(func):
    """Decorator for AT elements"""

    @functools.wraps(func)
    def wrapper(self, *args, tilt=0.0, ktap=0.0, **kwargs):
        elems = func(self, *args, **kwargs)
        for el in elems:
            tilt = float(tilt)  # MadParameter conversion
            if tilt != 0.0:
                tilt_elem(el, tilt)
            ktap = float(ktap)  # MadParameter conversion
            if ktap != 0.0:
                el.Scaling = 1.0 + ktap
            el.origin = self.origin
        return elems

    return wrapper


def poly_to_mad(x: Iterable[float], factor: float = 1.0) -> Generator[float]:
    """Convert polynomials from AT to MAD"""
    f = 1.0
    for n, vx in enumerate(x):
        yield factor * float(vx * f)
        f *= n + 1


def poly_from_mad(x: Iterable[float], factor: float = 1.0) -> Generator[float]:
    """Convert polynomials from MAD to AT"""
    f = 1.0
    for n, vx in enumerate(x):
        yield factor * float(vx / f)
        f *= n + 1


def p_to_at(a: float | Sequence[float]) -> np.ndarray:
    """Convert polynomials from MADX to AT"""
    if not isinstance(a, Sequence):
        # In case of a single element, we have a scalar instead of a tuple
        a = (a,)
    return np.fromiter(poly_from_mad(a), dtype=float)


def p_dict(keys: Iterable[str], a: Iterable[float]) -> dict[str, float]:
    """return K1, K2... from an AT Polynom"""
    return {k: v for k, v in zip(keys, poly_to_mad(a)) if k and (v != 0.0)}


def p_list(a: Iterable[float], factor: float = 1.0):
    """Return a Polynom list"""
    return list(poly_to_mad(a, factor=factor))


# noinspection PyUnusedLocal
def _keyparser(parser, argcount, argstr):
    """Return the pair key, value for the given 'key' argument"""
    return argstr, parser.evaluate(argstr)


# ------------------------------
#  Base class for MAD-X elements
# ------------------------------


class _MadElement(ElementDescr):
    """Description of MADX elements"""

    str_attr = {"apertype", "refer", "refpos", "sequence", "from"}

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

    def meval(self, params: dict):
        """Evaluation of superfluous parameters"""

        def mpeval(v):
            if isinstance(v, MadParameter):
                return v.evaluate()
            elif isinstance(v, str):
                return v
            elif isinstance(v, Sequence):
                return tuple(mpeval(a) for a in v)
            else:
                return v

        # return {k: mpeval(v) for k, v in params.items()}
        return {}


# ------------------------------
#  MAD-X element classes
# ------------------------------


# noinspection PyPep8Naming
class drift(_MadElement):
    @mad_element
    def to_at(self, l=0.0, **params):  # noqa: E741
        return [elt.Drift(self.name, l, **self.meval(params))]


# noinspection PyPep8Naming
class marker(_MadElement):
    at2mad = {}

    @mad_element
    def to_at(self, **params):
        return [elt.Marker(self.name, **self.meval(params))]


# noinspection PyPep8Naming
class quadrupole(_MadElement):
    @mad_element
    def to_at(self, l, k1=0.0, k1s=0.0, **params):  # noqa: E741
        atparams = {}
        k1s = float(k1s)  # MadParameter conversion  # MadParameter conversion
        if k1s != 0.0:
            atparams["PolynomA"] = [0.0, k1s]
        return [elt.Quadrupole(self.name, l, k1, **atparams, **self.meval(params))]

    @classmethod
    def from_at(cls, kwargs):
        el = super().from_at(kwargs)
        el.update(p_dict(["K0", "K1"], kwargs.pop("PolynomB", ())))
        el.update(p_dict(["K0S", "K1S"], kwargs.pop("PolynomA", ())))
        return el


# noinspection PyPep8Naming
class sextupole(_MadElement):
    @mad_element
    def to_at(self, l, k2=0.0, k2s=0.0, **params):  # noqa: E741
        atparams = {}
        k2s = float(k2s)  # MadParameter conversion
        if k2s != 0.0:
            atparams["PolynomA"] = [0.0, 0.0, k2s / 2.0]
        return [elt.Sextupole(self.name, l, k2 / 2.0, **atparams, **self.meval(params))]

    @classmethod
    def from_at(cls, kwargs):
        el = super().from_at(kwargs)
        el.update(p_dict(["K0", "K1", "K2"], kwargs.pop("PolynomB", ())))
        el.update(p_dict(["K0S", "K1S", "K2S"], kwargs.pop("PolynomA", ())))
        return el


# noinspection PyPep8Naming
class octupole(_MadElement):
    @mad_element
    def to_at(self, l, k3=0.0, k3s=0.0, **params):  # noqa: E741
        poly_b = [0.0, 0.0, 0.0, k3 / 6.0]
        poly_a = [0.0, 0.0, 0.0, k3s / 6.0]
        return [elt.Octupole(self.name, l, poly_a, poly_b, **self.meval(params))]

    @classmethod
    def from_at(cls, kwargs):
        el = super().from_at(kwargs)
        el.update(p_dict(["K0", "K1", "K2", "K3"], kwargs.pop("PolynomB", ())))
        el.update(p_dict(["K0S", "K1S", "K2S", "K3S"], kwargs.pop("PolynomA", ())))
        return el


# noinspection PyPep8Naming
class multipole(_MadElement):
    at2mad = {}

    @mad_element
    def to_at(self, knl=(), ksl=(), **params):
        params.pop("l", None)
        return [
            elt.ThinMultipole(
                self.name, p_to_at(ksl), p_to_at(knl), **self.meval(params)
            )
        ]

    @classmethod
    def from_at(cls, kwargs, factor=1.0):
        el = super().from_at(kwargs)
        el["KNL"] = p_list(kwargs.pop("PolynomB", ()), factor=factor)
        el["KSL"] = p_list(kwargs.pop("PolynomA", ()), factor=factor)
        return el


# noinspection PyPep8Naming
class sbend(_MadElement):
    at2mad = {
        "Length": "L",
        "BendingAngle": "ANGLE",
        "EntranceAngle": "E1",
        "ExitAngle": "E2",
    }

    @mad_element
    def to_at(
        self,
        l,  # noqa: E741
        angle,
        e1=0.0,
        e2=0.0,
        k1=0.0,
        k2=0.0,
        k1s=0.0,
        hgap=None,
        fint=0.0,
        **params,
    ):
        atparams = {}
        if hgap is not None:
            fintx = params.pop("fintx", fint)
            atparams.update(FullGap=2.0 * hgap, FringeInt1=fint, FringeInt2=fintx)
        k2 = float(k2)  # MadParameter conversion
        if k2 != 0.0:
            atparams["PolynomB"] = [0.0, k1, k2 / 2.0]
        k1s = float(k1s)  # MadParameter conversion
        if k1s != 0.0:
            atparams["PolynomA"] = [0.0, k1s]

        return [
            elt.Dipole(
                self.name,
                l,
                angle,
                k1,
                EntranceAngle=e1,
                ExitAngle=e2,
                **atparams,
                **self.meval(params),
            )
        ]

    @classmethod
    def from_at(cls, kwargs):
        el = super().from_at(kwargs)
        el.update(p_dict(["K0", "K1", "K2"], kwargs.pop("PolynomB", ())))
        el.update(p_dict(["K0S", "K1S", "K2S"], kwargs.pop("PolynomA", ())))
        return el


# noinspection PyPep8Naming
class rbend(sbend):
    @mad_element
    def to_at(self, l, angle, e1=0.0, e2=0.0, **params):  # noqa: E741
        hangle = abs(0.5 * angle)
        arclength = l / sinc(hangle)
        return super().to_at(arclength, angle, e1=hangle + e1, e2=hangle + e2, **params)

    @property
    def length(self):
        """Element length"""
        hangle = 0.5 * self["angle"]
        return self["l"] / sinc(hangle)


# noinspection PyPep8Naming
class kicker(_MadElement):
    @mad_element
    def to_at(self, l=0.0, hkick=0.0, vkick=0.0, **params):  # noqa: E741
        kicks = np.array([hkick, vkick], dtype=float)
        return [elt.Corrector(self.name, l, kicks, **self.meval(params))]

    @classmethod
    def from_at(cls, kwargs):
        el = super().from_at(kwargs)
        kicks = kwargs.pop("KickAngle", (0.0, 0.0))
        el["HKICK"] = kicks[0]
        el["VKICK"] = kicks[1]
        return el


# noinspection PyPep8Naming
class hkicker(kicker):
    @mad_element
    def to_at(self, l=0.0, kick=0.0, **params):  # noqa: E741
        return super().to_at(l=l, hkick=kick, **params)


# noinspection PyPep8Naming
class vkicker(kicker):
    @mad_element
    def to_at(self, l=0.0, kick=0.0, **params):  # noqa: E741
        return super().to_at(l=l, vkick=kick, **params)


# noinspection PyPep8Naming
class rfcavity(_MadElement):
    @mad_element
    def to_at(
        self,
        l=0.0,  # noqa: E741
        volt=0.0,
        freq=np.nan,
        lag=0.0,
        harmon=0,
        **params,
    ):
        cavity = elt.RFCavity(
            self.name,
            l,
            1.0e6 * volt,
            1.0e6 * freq,
            harmon,
            0.0,
            PassMethod="IdentityPass" if float(l) == 0.0 else "DriftPass",
            **self.meval(params),
        )
        return [cavity]

    @classmethod
    def from_at(cls, kwargs):
        el = super().from_at(kwargs)
        el["VOLT"] = 1.0e-6 * kwargs.pop("Voltage")
        el["FREQ"] = 1.0e-6 * kwargs.pop("Frequency")
        return el


# noinspection PyPep8Naming
class monitor(_MadElement):
    @mad_element
    def to_at(self, l=0.0, **params):  # noqa: E741
        hl = 0.5 * l  # MadParameter conversion
        if hl == 0.0:
            return [elt.Monitor(self.name, **self.meval(params))]
        else:
            drname = self.name + ".1"
            return [
                elt.Drift(drname, hl),
                elt.Monitor(self.name, **self.meval(params)),
                elt.Drift(drname, hl),
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
    """VALUE command"""
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
            nonlocal drcounter
            if dl > 1.0e-5:
                yield from drift(name=f"drift_{drcounter}", l=dl).expand(parser)
                drcounter += 1
            elif dl < -1.0e-5:
                elemtype = type(el).__name__.upper()
                raise ValueError(f"{elemtype}({el.name!r}) is overlapping by {-dl} m")

        drcounter = 0
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
    def to_at(name, *args, **params):
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


class _MadParser(LowerCaseParser, UnorderedParser):
    """Common class for both MAD8 anf MAD-X parsers"""

    _delimiter = ";"
    _linecomment = ("!", "//")
    _endfile = "return"
    _undef_key = "none"

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
            call=_Call(self),
            beam=_Beam(self),
            sequence=_Sequence,
            centre="centre",
            entry="entry",
            exit="exit",
            **kwargs,
        )
        self.current_sequence = None
        self["beam"]()

    def clear(self):
        super().clear()
        self.current_sequence = None
        self["beam"]()

    def _assign_deferred(self, value: str):
        """Deferred assignment"""
        if value[0] == "(" and value[-1] == ")":
            # Array variable: convert to tuple
            value, matches = protect(value[1:-1], fence=(r"\(", r"\)"))
            return tuple(
                MadParameter(self, v) for v in restore(matches, *value.split(","))
            )
        else:
            # Scalar variable
            return MadParameter(self, value)

    def _argparser(self, argcount, argstr: str, **kwargs):
        key, *value = split_ignoring_parentheses(
            argstr, delimiter=":=", fence=('"', '"'), maxsplit=1
        )
        if value:
            return key, self._assign_deferred(value[0])
        else:
            return super()._argparser(argcount, argstr, **kwargs)

    def _assign(self, label: str, key: str, val: str):
        # Special treatment of "line=(...)" assignments
        if key == "line":
            val = val.replace(")", ",)")  # For tuples with a single item
            return label, _Line(self.evaluate(val), name=label)
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

    def _decode(self, label: str, cmdname: str, *argnames: str) -> None:
        left, *right = cmdname.split(":=")
        if right:
            self[left] = self._assign_deferred(right[0])
        else:
            super()._decode(label, cmdname, *argnames)

    def _format_statement(self, line: str) -> str:
        line = _separator.sub(",", line)  # Replace space separators with commas
        # turn curly braces into parentheses. Must be done before splitting arguments
        line = line.replace("{", "(").replace("}", ")")
        return super()._format_statement(line)

    def _get_beam(self, key: str):
        """Get the beam object for a given sequence"""
        try:
            beam = self[f"beam%{self._gen_key(key)}"]
        except KeyError:
            beam = self["beam%"]
        return beam

    def _generator(self, params: dict) -> Iterable[elt.Element]:
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

        for elem in super()._generator(params):
            if isinstance(elem, elt.RFCavity):
                cavities.append(elem)
            cell_length += getattr(elem, "Length", 0.0)
            yield elem

        rev = beta() * clight / cell_length

        # Set the frequency of cavities in which it is not specified
        for cav in cavities:
            if np.isnan(cav.Frequency):
                cav.Frequency = rev * cav.HarmNumber
            elif cav.HarmNumber == 0:
                cav.HarmNumber = round(cav.Frequency / rev)

    def lattice(self, use: str = "ring", **kwargs):
        """Create a lattice from the selected sequence

        Parameters:
            use:                Name of the MAD sequence or line containing the desired
              lattice. Default: ``ring``

        Keyword Args:
            name (str):         Name of the lattice. Default: MAD sequence name.
            particle(Particle): Circulating particle. Default: from MAD
            energy (float):     Energy of the lattice [eV]. Default: from MAD
            periodicity(int):   Number of periods. Default: 1
            *:                  All other keywords will be set as Lattice attributes
        """
        part = kwargs.get("particle", None)
        if isinstance(part, str):
            kwargs["particle"] = Particle(part)

        lat = super().lattice(use=use, **kwargs)

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

    _continuation = None
    _blockcomment = ("/*", "*/")

    def __init__(self, **kwargs):
        """
        Args:
            strict:     If :py:obj:`False`, assign 0 to undefined variables
            verbose:    If :py:obj:`True`, print details on the processing
            **kwargs:   Initial variable definitions
        """
        super().__init__(
            globals(),
            **kwargs,
        )

    def _format_command(self, expr: str) -> str:
        """Format a command for evaluation"""
        expr = expr.replace("->", ".")  # Attribute access: VAR->ATTR
        expr = expr.replace("^", "**")  # Exponentiation
        return super()._format_command(expr)


def load_madx(
    *files: str,
    use: str = "ring",
    strict: bool = True,
    verbose: bool = False,
    **kwargs,
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
    params = {
        key: kwargs.pop(key)
        for key in ("name", "particle", "energy", "periodicity")
        if key in kwargs
    }
    parser.parse_files(*files, **kwargs)
    return parser.lattice(use=use, **params)


def longmultipole(kwargs):
    length = kwargs.get("Length", 0.0)
    if length == 0.0:
        return multipole.from_at(kwargs)
    else:
        drname = kwargs["FamName"] + ".1"
        dr1 = drift.from_at({"FamName": drname, "Length": 0.5 * length})
        dr2 = drift.from_at({"FamName": drname, "Length": 0.5 * length})
        return [dr1, multipole.from_at(kwargs, factor=length), dr2]


def ignore(kwargs):
    length = kwargs.get("Length", 0.0)
    if length == 0.0:
        print(f"{kwargs['name']} is replaced by a marker")
        return marker.from_at(kwargs)
    else:
        print(f"{kwargs['name']} is replaced by a drift")
        return drift.from_at(kwargs)


_AT2MAD = {
    elt.Quadrupole: quadrupole.from_at,
    elt.Sextupole: sextupole.from_at,
    elt.Octupole: octupole.from_at,
    elt.ThinMultipole: multipole.from_at,
    elt.Multipole: longmultipole,
    elt.RFCavity: rfcavity.from_at,
    elt.Drift: drift.from_at,
    elt.Bend: sbend.from_at,
    elt.Marker: marker.from_at,
    elt.Monitor: monitor.from_at,
    elt.Corrector: kicker.from_at,
}


class _MadExporter(Exporter):
    use_line = False

    def generate_madelems(
        self, eltype: type[elt.Element], elemdict: dict
    ) -> ElementDescr | list[ElementDescr]:
        return _AT2MAD.get(eltype, ignore)(elemdict)

    def print_beam(self, file):
        part = str(self.particle)
        if part == "relativistic":
            part = "electron"
        data = {
            "ENERGY": 1.0e-9 * self.energy,
            "PARTICLE": part.upper(),
            "RADIATE": self.is_6d,
        }
        attrs = [f"{k}={ElementDescr.attr_format(v)}" for k, v in data.items()]
        line = ", ".join(["BEAM".ljust(10)] + attrs)
        print(f"{line}{self.delimiter}", file=file)


class _MadxExporter(_MadExporter):
    delimiter = ";"
    continuation = ""
    bool_fmt = {False: "FALSE", True: "TRUE"}


def save_madx(ring: Lattice, filename: str | None = None, **kwargs):
    """Save a :py:class:`.Lattice` as a MAD-X file

    Args:
        ring:   lattice
        filename: file to be created. If None, write to sys.stdout

    Keyword Args:
        use (str | None): name of the created SEQUENCE of LINE.
          Default: name of the PyAT lattice
        use_line (bool):  If True, use a MAD "LINE" format. Otherwise, use
          a MAD "SEQUENCE"
    """
    exporter = _MadxExporter(ring, **kwargs)
    exporter.export(filename)


register_format(
    ".seq",
    load_madx,
    save_madx,
    descr="MAD-X lattice description. See :py:func:`.load_madx`.",
)
