"""Load a lattice from an Elegant file"""

from __future__ import annotations

__all__ = ["ElegantParser", "load_elegant", "save_elegant"]

import sys
import functools
from math import sqrt, factorial
from os.path import abspath
from collections.abc import Iterable
import warnings

import numpy as np
from scipy.constants import c as clight

from .allfiles import register_format
from .file_input import ElementDescr, BaseParser
from .file_input import skip_names, ignore_names, ignore_class
from .file_output import translate
from ..lattice import Particle, Lattice, Filter, elements as elt, tilt_elem, shift_elem

# noinspection PyProtectedMember
from .madx import sinc, _Line, p_dict, p_list
from . import rpn


# -------------------
#  Utility functions
# -------------------


def elegant_element(func):
    """Decorator which tilts and shifts the decorated AT element"""

    @functools.wraps(func)
    def wrapper(
        *args,
        origin="",
        tilt=0.0,
        dx=0.0,
        dy=0.0,
        n_slices=None,
        n_kicks=None,
        **kwargs,
    ):
        if n_kicks is not None:  # Deprecated parameter
            kwargs["NumIntSteps"] = int(n_kicks / 4)
        if n_slices is not None:
            kwargs["NumIntSteps"] = n_slices
        elems = func(*args, **kwargs)
        for el in elems:
            if tilt != 0.0:
                tilt_elem(el, tilt)
            if not (dx == 0.0 and dy == 0.0):
                shift_elem(el, deltax=dx, deltaz=dy)
            el.origin = origin
        return elems

    return wrapper


class ElegantVar(str):
    def __new__(cls, expr):
        return super().__new__(cls, expr)

    def __float__(self):
        return float(rpn.evaluate(self))

    def __int__(self):
        return int(rpn.evaluate(self))


# ------------------------------
#  Base class for Elegant elements
# ------------------------------


class _ElegantElement(ElementDescr):
    """Description of MADX elements"""

    str_attr = {"filename", "group", "mode"}


# ------------------------------
#  Elegant classes
# ------------------------------


# noinspection PyPep8Naming
class DRIF(_ElegantElement):
    @elegant_element
    def to_at(self, l=0.0, **params):  # noqa: E741
        return [elt.Drift(self.name, l, **params)]


# noinspection PyPep8Naming
class MARK(_ElegantElement):
    at2mad = {}

    @elegant_element
    def to_at(self, **params):
        return [elt.Marker(self.name, **params)]


# noinspection PyPep8Naming
class KQUAD(_ElegantElement):
    @elegant_element
    def to_at(self, l, k1=0.0, **params):  # noqa: E741
        return [elt.Quadrupole(self.name, l, k1, **params)]

    @classmethod
    def from_at(cls, kwargs):
        el = super().from_at(kwargs)
        el.update(p_dict(["", "K1"], kwargs.pop("PolynomB", ())))
        return el


# noinspection PyPep8Naming
class KSEXT(_ElegantElement):
    @elegant_element
    def to_at(self, l, k2=0.0, **params):  # noqa: E741
        return [elt.Sextupole(self.name, l, k2 / 2.0, **params)]

    @classmethod
    def from_at(cls, kwargs):
        el = super().from_at(kwargs)
        el.update(p_dict(["", "K1", "K2"], kwargs.pop("PolynomB", ())))
        return el


# noinspection PyPep8Naming
class KQUSE(_ElegantElement):
    @elegant_element
    def to_at(self, l, k1=0.0, k2=0.0, **params):  # noqa: E741
        poly_b = [0.0, k1, k2 / 2.0]
        return [elt.Multipole(self.name, l, [], poly_b, **params)]


# noinspection PyPep8Naming
class KOCT(_ElegantElement):
    @elegant_element
    def to_at(self, l, k3=0.0, **params):  # noqa: E741
        poly_b = [0.0, 0.0, 0.0, k3 / 6.0]
        poly_a = [0.0, 0.0, 0.0, 0.0]
        return [elt.Octupole(self.name, l, poly_a, poly_b, **params)]

    @classmethod
    def from_at(cls, kwargs):
        el = super().from_at(kwargs)
        el.update(p_dict(["", "", "", "K3"], kwargs.pop("PolynomB", ())))
        return el


class MULT(_ElegantElement):
    at2mad = {"Length": "L", "Order": "ORDER", "Value": "KNL"}

    @elegant_element
    def to_at(self, l=0, knl=0.0, order=1, **params):  # noqa: E741
        poly_a = np.zeros(order + 1)
        poly_b = np.zeros(order + 1)
        if l == 0.0:
            poly_b[order] = knl / factorial(order)
            return [elt.ThinMultipole(self.name, poly_a, poly_b, **params)]
        else:
            poly_b[order] = knl / factorial(order) / l
            return [elt.Multipole(self.name, l, poly_a, poly_b, **params)]


# noinspection PyPep8Naming
class CSBEND(_ElegantElement):
    at2mad = {
        "Length": "L",
        "BendingAngle": "ANGLE",
        "EntranceAngle": "E1",
        "ExitAngle": "E2",
    }

    @elegant_element
    def to_at(
        self,
        l,  # noqa: E741
        angle,
        e1=0.0,
        e2=0.0,
        k1=0.0,
        k2=0.0,
        k3=0.0,
        k4=0.0,
        hgap=None,
        fint=0.0,
        **params,
    ):
        if hgap is not None:
            params.update(FullGap=2.0 * hgap, FringeInt1=fint, FringeInt2=fint)
        if k2 != 0.0 or k3 != 0.0 or k4 != 0.0:
            params["PolynomB"] = [0.0, k1, k2 / 2.0, k3 / 6.0, k4 / 24.0]
        return [
            elt.Dipole(
                self.name,
                l,
                angle,
                k1,
                EntranceAngle=e1,
                ExitAngle=e2,
                **params,
            )
        ]

    @classmethod
    def from_at(cls, kwargs):
        el = super().from_at(kwargs)
        el.update(p_dict(["K0", "K1", "K2", "K3"], kwargs.pop("PolynomB", ())))
        return el


# noinspection PyPep8Naming
class RBEN(CSBEND):
    @elegant_element
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
class KICKER(_ElegantElement):
    @elegant_element
    def to_at(self, l=0.0, hkick=0.0, vkick=0.0, **params):  # noqa: E741
        kicks = np.array([hkick, vkick], dtype=float)
        return [elt.Corrector(self.name, l, kicks, **params)]

    @classmethod
    def from_at(cls, kwargs):
        el = super().from_at(kwargs)
        kicks = kwargs.pop("KickAngle", (0.0, 0.0))
        el["HKICK"] = kicks[0]
        el["VKICK"] = kicks[1]
        return el


# noinspection PyPep8Naming
class HKICK(KICKER):
    @elegant_element
    def to_at(self, l=0.0, kick=0.0, **params):  # noqa: E741
        return super().to_at(l=l, hkick=kick, **params)


# noinspection PyPep8Naming
class VKICK(KICKER):
    @elegant_element
    def to_at(self, l=0.0, kick=0.0, **params):  # noqa: E741
        return super().to_at(l=l, vkick=kick, **params)


# noinspection PyPep8Naming
class RFCA(_ElegantElement):
    at2mad = {"Length": "L", "Voltage": "VOLT", "Frequency": "FREQ"}

    @elegant_element
    def to_at(
        self,
        l=0.0,  # noqa: E741
        volt=0.0,
        freq=np.nan,
        **params,
    ):
        cavity = elt.RFCavity(
            self.name,
            l,
            volt,
            freq,
            0,
            0.0,
            PassMethod="IdentityPass" if l == 0.0 else "DriftPass",
            **params,
        )
        return [cavity]


# noinspection PyPep8Naming
class MONI(_ElegantElement):
    @elegant_element
    def to_at(self, l=0.0, **params):  # noqa: E741
        if l == 0.0:
            return [elt.Monitor(self.name, **params)]
        else:
            hl = 0.5 * l
            return [
                elt.Drift(self.name, hl, origin="MONI"),
                elt.Monitor(self.name, **params),
                elt.Drift(self.name, hl, origin="MONI"),
            ]


# noinspection PyPep8Naming
class HMON(MONI):
    pass


# noinspection PyPep8Naming
class VMON(MONI):
    pass


def multipole(kwargs):
    """AT ThinMultipole or Multipole converted to Elegant MULT"""

    def singlemul(o, v):
        mn = ".".join((name, str(o)))
        return MULT.from_at({"FamName": mn, "Length": length, "Order": o, "Value": v})

    name = kwargs.pop("FamName")
    length = kwargs.pop("Length", 0.0)
    poly_b = p_list(kwargs.pop("PolynomB", ()))
    if length == 0.0:
        return [singlemul(order, v) for order, v in enumerate(poly_b) if v != 0.0]
    else:
        for order, v in enumerate(poly_b):
            if v != 0.0:
                return singlemul(order, length * v)


def ignore(kwargs):
    """AT element ignored in Elegant: convert to marker or drift"""
    length = kwargs.get("Length", 0.0)
    if length == 0.0:
        print(f"{kwargs['name']} is replaced by a marker")
        return MARK.from_at(kwargs)
    else:
        print(f"{kwargs['name']} is replaced by a drift")
        return DRIF.from_at(kwargs)


SOLENOID = ignore_class("SOLENOID", _ElegantElement)

skip_names(
    globals(),
    _ElegantElement,
    [
        "MAXAMP",
        "CHARGE",
        "RECIRC",
        "MALIGN",
        "SREFFECTS",
        "PFILTER",
        "ENERGY",
        "SCATTER",
        "WATCH",
        "WAKE",
    ],
)


ignore_names(globals(), _ElegantElement, ["SCRAPER", "ECOL", "RCOL", "CSRDRIF"])

EDRIFT = DRIFT = DRIF
QUAD = QUADRUPOLE = KQUAD
CSRCSBEN = CSBEN = SBEN = SBEND = CSBEND
CRBEN = CRBEND = RBEND = RBEN
SEXT = SEXTUPOLE = KSEXT
OCTU = OCTUPOLE = KOCT
HKICKER = HKICK
VKICKER = VKICK
MONITOR = MONI
HMONITOR = HMON
VMONITOR = VMON
RFCW = RFCA
SOLE = SOLENOID
MULTIPOLE = MULT


class ElegantParser(BaseParser):
    # noinspection PyUnresolvedReferences
    """Elegant parser

    The parser is a subclass of :py:class:`dict` and is a database containing all the
    Elegant objects.

    Example:
        Parse a file:

        >>> parser = at.ElegantParser()
        >>> parser.parse_file("file1")

        Look at the "qf1" element

        >>> parser["QF1"]
        QUAD(name=QF1, l=1.0, k1=0.5, tilt=0.001)

        Generate an AT :py:class:`.Lattice` from the "RING" sequence

        >>> ring = parser.lattice(use="RING")  # generate an AT Lattice
    """

    def __init__(self, **kwargs):
        """
        Args:
            verbose:    If :py:obj:`True`, print details on the processing
            **kwargs:   Initial variable definitions
        """
        super().__init__(
            globals(),
            continuation="&",
            linecomment="!",
            endfile="RETURN",
            **kwargs,
        )

    @classmethod
    def _defkey(cls, expr: str, quoted: bool) -> str:
        """substitutions to get a valid python identifier"""
        expr = super()._defkey(expr, quoted)
        return expr if quoted else expr.upper()

    def _format_statement(self, line: str) -> str:
        """Reformat the input line"""
        return line.upper()

    def _assign(self, label: str, key: str, val: str):
        # Special treatment of "line=(...)" commands
        if key == "LINE":
            val = val.replace(")", ",)")  # For tuples with a single item
            return label, _Line(self._evaluate(val), name=label)
        else:
            return super()._assign(label, key, val)

    def _command(self, label, cmdname, *args: str, **kwargs):
        # Special treatment of label #INCLUDE
        if label == "#INCLUDE":
            file = cmdname[1:-1] if cmdname[0] == '"' else cmdname
            self.parse_files(file, final=False)
        else:
            return super()._command(label, cmdname, *args, **kwargs)

    def _argparser(
        self,
        argcount: int,
        argstr: str,
        *,
        bool_attr: tuple[str] = (),
        str_attr: tuple[str] = (),
        pos_args: tuple[str] = (),
    ):
        """Evaluate a command argument and return a pair (key, value)"""

        def arg_value(k, v):
            if k in str_attr:
                return v[1:-1] if v[0] == '"' else v
            elif v[0] == '"':
                return ElegantVar(v[1:-1])  # Deferred RPN evalution
                # return rpn.evaluate(v[1:-1])  # Immediate RPN evaluation
            else:
                return self._evaluate(v)

        key, *value = argstr.split(sep="=", maxsplit=1)
        if value:  # Keyword argument
            key = key.lower()
            return key, arg_value(key, value[0])
        else:
            try:
                key = pos_args[argcount]
            except IndexError:
                print(f"Unexpected positional argument '{argstr}' ignored")
                return None
            return key, arg_value(key, argstr)

    def lattice(self, use: str = "RING", **kwargs):
        """Create a lattice from the selected sequence

        - Elegant lattice files do not specify the beam energy.
          :py:class:`ElegantParser` sets it by default to 1.0 GeV. Use the *energy*
          keyword to set it to the desired value.
        - Long elements are split according to the default AT value of *NumIntSteps*
          (10) unless *N_SLICES* is specified in the Elegant element definition.

        Parameters:
            use:                Name of the Elegant LINE describing the desired
              lattice. Default: ``RING``

        Keyword Args:
            name (str):         Name of the lattice. Default: line name.
            particle(Particle): Circulating particle. Default: Particle("relativistic")
            energy (float):     Energy of the lattice [eV], Default: 1.0E9
            periodicity(int):   Number of periods. Default: 1
            *:                  All other keywords will be set as Lattice attributes
        """

        def elegant_filter(params: dict, elems: Filter, *args) -> Iterable[elt.Element]:
            def beta() -> float:
                rest_energy = params["particle"].rest_energy
                if rest_energy == 0.0:
                    return 1.0
                else:
                    gamma = float(params["energy"] / rest_energy)
                    return sqrt(1.0 - 1.0 / gamma / gamma)

            cavities = []
            cell_length = 0

            for elem in elems(params, *args):
                if isinstance(elem, elt.RFCavity):
                    cavities.append(elem)
                cell_length += getattr(elem, "Length", 0.0)
                yield elem

            params["_length"] = cell_length
            rev = beta() * clight / cell_length

            # Set the cavities' Emergy and harmonic number
            if cavities:
                cavities.sort(key=lambda el: el.Frequency)
                for cav in cavities:
                    cav.Energy = params["energy"]
                    cav.HarmNumber = cav.Frequency / rev
                params["_cell_harmnumber"] = getattr(cavities[0], "HarmNumber", np.nan)

        kwargs.setdefault("energy", 1.0e9)
        part = kwargs.setdefault("particle", "relativistic")
        if isinstance(part, str):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=UserWarning)
                kwargs["particle"] = Particle(part)
        return Lattice(self._generator, iterator=elegant_filter, use=use, **kwargs)


def load_elegant(
    *files: str, use: str = "RING", verbose: bool = False, **kwargs
) -> Lattice:
    """Create a :py:class:`.Lattice` from Elegant lattice files

    - Elegant lattice files do not specify the beam energy. :py:class:`ElegantParser`
      sets it by default to 1.0 GeV. Use the *energy* keyword to set it to the
      desired value.
    - Long elements are split according to the default AT value of *NumIntSteps* (10)
      unless *N_SLICES* is specified in the Elegant element definition.

    Parameters:
        files:              Names of one or several Elegant lattice description files
        use:                Name of the Elegant LINE describing the desired
          lattice. Default: ``RING``
        verbose:            If :py:obj:`True`, print details on the processing

    Keyword Args:
        name (str):         Name of the lattice. Default: Elegant sequence name
        particle(Particle): Circulating particle. Default: Particle("relativistic")
        energy (float):     Energy of the lattice [eV]. Default: 1.0E9
        periodicity(int):   Number of periods. Default: 1
        *:                  Other keywords will be used as Lattice attributes

    Returns:
        lattice (Lattice):  New :py:class:`.Lattice` object

    See Also:
        :py:func:`.load_lattice` for a generic lattice-loading function.
    """
    parser = ElegantParser(verbose=verbose)
    absfiles = tuple(abspath(file) for file in files)
    parser.parse_files(*absfiles)
    return parser.lattice(use=use, in_file=absfiles, **kwargs)


_AT2EL = {
    elt.Quadrupole: KQUAD.from_at,
    elt.Sextupole: KSEXT.from_at,
    elt.Octupole: KOCT.from_at,
    elt.ThinMultipole: multipole,
    elt.Multipole: multipole,
    elt.RFCavity: RFCA.from_at,
    elt.Drift: DRIF.from_at,
    elt.Bend: SBEN.from_at,
    elt.Marker: MARK.from_at,
    elt.Monitor: MONI.from_at,
    elt.Corrector: KICKER.from_at,
}


def at2elegant(attype):
    return _AT2EL.get(attype, ignore)


def save_elegant(ring: Lattice, filename: str | None = None):
    kwargs = {
        "delimiter": "",
        "continuation": "&",
        "bool_fmt": {False: ".FALSE.", True: ".TRUE."},
        "use_line": True,
        "beam_descr": None,
    }
    if filename is None:
        translate(at2elegant, ring, file=sys.stdout, **kwargs)
    else:
        with open(filename, "w") as mfile:
            translate(at2elegant, ring, file=mfile, **kwargs)


register_format(
    ".lte",
    load_elegant,
    save_elegant,
    descr="Elegant lattice description. See :py:func:`.load_elegant`.",
)
