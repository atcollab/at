"""Load a lattice from a MAD8 file."""

from __future__ import annotations

__all__ = ["Mad8Parser", "load_mad8"]

from os.path import abspath

# functions known by MAD-8
from math import pi, e, sqrt, exp, log, sin, cos, tan  # noqa: F401
from math import asin  # noqa: F401

# constants known by MAD-8
from scipy.constants import c as clight, e as qelect  # noqa: F401
from scipy.constants import physical_constants as _cst

from ..lattice import Lattice
# noinspection PyProtectedMember
from .madx import _MadParser

# Commands known by MAD8
# noinspection PyProtectedMember
from .madx import (  # noqa: F401
    drift,
    marker,
    quadrupole,
    sextupole,
    octupole,
    multipole,
    sbend,
    rbend,
    kicker,
    hkicker,
    vkicker,
    rfcavity,
    monitor,
    hmonitor,
    vmonitor,
    solenoid,
)

# Constants known by MAD-8
_true_ = True
_yes_ = True
_t_ = True
_on_ = True
_false_ = False
_no_ = False
_f_ = False
_off_ = False
twopi = 2 * pi
degrad = 180.0 / pi
raddeg = pi / 180.0
emass = 1.0e-03 * _cst["electron mass energy equivalent in MeV"][0]  # [GeV]
pmass = 1.0e-03 * _cst["proton mass energy equivalent in MeV"][0]  # [GeV]


class Mad8Parser(_MadParser):
    # noinspection PyUnresolvedReferences
    """MAD-X specific parser

    The parser is a subclass of :py:class:`dict` and is database containing all the
    MAD-X variables.

    Example:
        Parse a 1st file:

        >>> parser = at.Mad8Parser()
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

    def __init__(self):
        """"""
        super().__init__(globals(), continuation="&")


def load_mad8(*files: str, use: str = "ring", **kwargs) -> Lattice:
    """Create a :py:class:`.Lattice`  from MAD8 files

    Parameters:
        files:              Names of one or several MAD8 files
        use:                Name of the MAD8 sequence or line containing the desired
          lattice. Default: ``ring``

    Keyword Args:
        name (str):         Name of the lattice. Default: MAD8 sequence name.
        particle(Particle): Circulating particle. Default: from MAD8
        energy (float):     Energy of the lattice [eV]. Default: from MAD8
        periodicity(int):   Number of periods. Default: 1
        *:                  All other keywords will be set as Lattice attributes

    Returns:
        lattice (Lattice):  New :py:class:`.Lattice` object

    See Also:
        :py:func:`.load_lattice` for a generic lattice-loading function.
    """
    parser = Mad8Parser()
    absfiles = tuple(abspath(file) for file in files)
    kwargs.setdefault("in_file", absfiles)
    parser.parse_files(*absfiles)
    return parser.lattice(use=use)
