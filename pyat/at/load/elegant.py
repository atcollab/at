"""Load a lattice from an Elegant file"""

from __future__ import annotations

__all__ = ["ElegantParser", "load_elegant"]

import functools
from math import sqrt, sin, factorial
from os.path import abspath
from collections.abc import Iterable
import warnings

import numpy as np
from scipy.constants import c as clight

from ..lattice import Particle, Lattice, Filter, elements as elt, tilt_elem, shift_elem
from .file_input import ElementDescr, BaseParser
from .file_input import skip_names, ignore_names, ignore_class
from .madx import sinc, _Line
from .rpn import evaluate


# -------------------
#  Utility functions
# -------------------


def misalign(func):
    """Decorator which tilts and shifts the decorated AT element"""

    @functools.wraps(func)
    def wrapper(
        name, *args, tilt=None, dx=0.0, dy=0.0, n_slices=None, n_kicks=None, **kwargs
    ):
        if n_kicks is not None:  # Deprecated parameter
            kwargs["NumIntSteps"] = int(n_kicks / 4)
        if n_slices is not None:
            kwargs["NumIntSteps"] = n_slices
        elems = func(name, *args, **kwargs)
        if tilt is not None:
            tilt_elem(elems[0], tilt)
        if not (dx == 0.0 and dy == 0.0):
            shift_elem(elems[0], deltax=dx, deltaz=dy)
        return elems

    return wrapper


class ElegantVar(str):
    def __new__(cls, parser, expr):
        return super().__new__(cls, expr)

    # noinspection PyUnusedLocal
    def __init__(self, parser, expr):
        self.parser = parser

    def __float__(self):
        return float(evaluate(self))

    def __int__(self):
        return int(evaluate(self))


# ------------------------------
#  Elegant classes
# ------------------------------


# noinspection PyPep8Naming
class DRIF(ElementDescr):
    @staticmethod
    def convert(name: str, l=0.0, **params):  # noqa: E741
        return [elt.Drift(name, l, **params)]


# noinspection PyPep8Naming
class MARK(ElementDescr):
    @staticmethod
    def convert(name, **params):
        return [elt.Marker(name, **params)]


# noinspection PyPep8Naming
class QUAD(ElementDescr):
    @staticmethod
    @misalign
    def convert(name, l, k1=0.0, **params):  # noqa: E741
        return [elt.Quadrupole(name, l, k1, **params)]


# noinspection PyPep8Naming
class SEXT(ElementDescr):
    @staticmethod
    @misalign
    def convert(name, l, k2=0.0, order=2, **params):  # noqa: E741
        return [elt.Sextupole(name, l, k2 / 2.0, **params)]


# noinspection PyPep8Naming
class OCTU(ElementDescr):
    @staticmethod
    @misalign
    def convert(name, l, k3=0.0, **params):  # noqa: E741
        poly_b = [0.0, 0.0, 0.0, k3 / 6.0]
        poly_a = [0.0, 0.0, 0.0, 0.0]
        return [elt.Multipole(name, l, poly_a, poly_b, **params)]


class MULT(ElementDescr):
    @staticmethod
    @misalign
    def convert(name, l=0, knl=0.0, order=1, **params):  # noqa: E741
        poly_a = np.zeros(order + 1)
        poly_b = np.zeros(order + 1)
        poly_b[order] = knl / factorial(order)
        return [elt.Multipole(name, l, poly_a, poly_b, **params)]


# noinspection PyPep8Naming
class SBEN(ElementDescr):
    @staticmethod
    @misalign
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
            params.update(FullGap=2.0 * hgap, FringeInt1=fint, FringeInt2=fint)
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
class RBEN(ElementDescr):
    @staticmethod
    @misalign
    def convert(name, l, angle, e1=0.0, e2=0.0, **params):  # noqa: E741
        hangle = 0.5 * angle
        arclength = l * sinc(hangle)
        return SBEN.convert(
            name, arclength, angle, e1=hangle + e1, e2=hangle + e2, **params
        )

    def _length(self):
        hangle = 0.5 * self.angle
        return self["l"] * hangle / sin(hangle)


# noinspection PyPep8Naming
class KICKER(ElementDescr):
    @staticmethod
    @misalign
    def convert(name, l=0.0, hkick=0.0, vkick=0.0, **params):  # noqa: E741
        kicks = np.array([hkick, vkick], dtype=float)
        return [elt.Corrector(name, l, kicks, **params)]


# noinspection PyPep8Naming
class HKICK(ElementDescr):
    @staticmethod
    @misalign
    def convert(name, l=0.0, kick=0.0, **params):  # noqa: E741
        return KICKER.convert(name, l=l, hkick=kick, **params)


# noinspection PyPep8Naming
class VKICK(ElementDescr):
    @staticmethod
    @misalign
    def convert(name, l=0.0, kick=0.0, **params):  # noqa: E741
        return KICKER.convert(name, l=l, vkick=kick, **params)


# noinspection PyPep8Naming
class RFCA(ElementDescr):
    @staticmethod
    def convert(
        name,
        l=0.0,  # noqa: E741
        volt=0.0,
        freq=np.nan,
        **params,
    ):
        cavity = elt.RFCavity(
            name,
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
class MONI(ElementDescr):
    @staticmethod
    def convert(name, l=0.0, **params):  # noqa: E741
        if l == 0.0:
            return [elt.Monitor(name, **params)]
        else:
            hl = 0.5 * l
            return [
                elt.Drift(name, hl, origin="MONI"),
                elt.Monitor(name, **params),
                elt.Drift(name, hl, origin="MONI"),
            ]


# noinspection PyPep8Naming
class HMON(MONI):
    pass


# noinspection PyPep8Naming
class VMON(MONI):
    pass


SOLENOID = ignore_class("SOLENOID", ElementDescr)

skip_names(
    globals(),
    ElementDescr,
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


ignore_names(globals(), ElementDescr, ["SCRAPER", "ECOL", "RCOL", "CSRDRIF"])

EDRIFT = DRIFT = DRIF
KQUAD = QUADRUPOLE = QUAD
CSRCSBEN = CSBEN = CSBEND = SBEND = SBEN
CRBEN = CRBEND = RBEND = RBEN
KSEXT = SEXTUPOLE = SEXT
KOCT = OCTUPOLE = OCTU
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

        def arg_value(v):
            if v[0] == '"':
                return ElegantVar(self, v[1:-1])
            else:
                return self._evaluate(v)

        key, *value = argstr.split(sep="=", maxsplit=1)
        if value:  # Keyword argument
            return key.lower(), arg_value(value[0])
        else:
            try:
                key = pos_args[argcount]
            except IndexError as exc:
                exc.args = ("too many positional arguments",)
                raise
            return key, arg_value(argstr)

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
    """Create a :py:class:`.Lattice` from Elegant files

    - Long elements are split according to the default AT value of *NumIntSteps* (10)
      unless *N_SLICES* is specified in the Elegant element definition.

    Parameters:
        files:              Names of one or several Elegant lattice description files
        use:                Name of the MADX sequence or line containing the desired
          lattice. Default: ``ring``
        verbose:            If :py:obj:`True`, print details on the processing

    Keyword Args:
        name (str):         Name of the lattice. Default: Elegant sequence name
        particle(Particle): Circulating particle. Default: 'relativistic'
        energy (float):     Energy of the lattice [eV]. Default: 0.0
        periodicity(int):   Number of periods. Default: 1
        *:                  Other keywords will be used as initial variable definitions

    Returns:
        lattice (Lattice):  New :py:class:`.Lattice` object

    See Also:
        :py:func:`.load_lattice` for a generic lattice-loading function.
    """
    parser = ElegantParser(verbose=verbose)
    absfiles = tuple(abspath(file) for file in files)
    params = {
        key: kwargs.pop(key)
        for key in ("name", "particle", "energy", "periodicity")
        if key in kwargs
    }
    parser.parse_files(*absfiles, **kwargs)
    return parser.lattice(use=use, in_file=absfiles, **params)
