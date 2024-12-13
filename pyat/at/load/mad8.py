r"""Using `MAD8`_ files with PyAT
=================================

Using MAD8 files is similar to using MAD-X files: see
:py:mod:`the MAD-X documentation <at.load.madx>`

.. _mad8: https://mad8.web.cern.ch/user/mad.html
"""

from __future__ import annotations

__all__ = ["Mad8Parser", "load_mad8", "save_mad8"]

# functions known by MAD-8
from math import pi, e, sqrt, exp, log, sin, cos, tan  # noqa: F401
from math import asin  # noqa: F401
import re

# constants known by MAD-8
from scipy.constants import c as clight, e as qelect  # noqa: F401
from scipy.constants import physical_constants as _cst

from ..lattice import Lattice

from .file_input import ignore_names

# noinspection PyProtectedMember
from .madx import _MadElement, _MadParser, _MadExporter

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
    value,
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

_attr = re.compile(r"\[([a-zA-Z_][\w.:]*)]")  # Identifier enclosed in square brackets

ignore_names(
    globals(),
    _MadElement,
    ["solenoid", "rfmultipole", "crabcavity", "elseparator", "collimator", "tkicker"],
)


class Mad8Parser(_MadParser):
    # noinspection PyUnresolvedReferences
    r"""MAD-X specific parser

    The parser is a subclass of :py:class:`dict` and is database containing all the
    MAD-X variables.

    Example:
        Parse a 1\ :sup:`st` file:

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

    _continuation = "&"
    _blockcomment = ("comment", "endcomment")

    def __init__(self, **kwargs):
        """
        Args:
            strict:     If :py:obj:`False`, assign 0 to undefined variables
            verbose:    If :py:obj:`True`, print details on the processing
            **kwargs:   Initial variable definitions
        """
        super().__init__(globals(), **kwargs)

    def _format_command(self, expr: str) -> str:
        """Evaluate an expression using *self* as local namespace"""
        expr = _attr.sub(r".\1", expr)  # Attribute access: VAR[ATTR]
        expr = expr.replace("^", "**")  # Exponentiation
        return super()._format_command(expr)


def load_mad8(
    *files: str,
    use: str = "ring",
    strict: bool = True,
    verbose: bool = False,
    **kwargs,
) -> Lattice:
    """Create a :py:class:`.Lattice` from MAD8 files

    - The *energy* and *particle* of the generated lattice are taken from the MAD8
      ``BEAM`` object, using the MAD8 default parameters: positrons at 1 Gev.
      These parameters are overloaded by the value given in the *energy* and
      *particle* keyword arguments.
    - The radiation state is given by the ``RADIATE`` flag of the ``BEAM`` object,
      using the AT defaults: RF cavities active, synchrotron radiation in dipoles and
      quadrupoles.
    - Long elements are split according to the default AT value for *NumIntSteps* (10).

    Parameters:
        files:              Names of one or several MAD8 files
        strict:             If :py:obj:`False`, assign 0 to undefined variables
        use:                Name of the MAD8 sequence or line containing the desired
          lattice. Default: ``ring``
        verbose:            If :py:obj:`True`, print details on the processing

    Keyword Args:
        name (str):         Name of the lattice. Default: MAD8 sequence name.
        particle(Particle): Circulating particle. Default: from MAD8
        energy (float):     Energy of the lattice [eV]. Default: from MAD8
        periodicity(int):   Number of periods. Default: 1
        *:                  Other keywords will be used as initial variable definitions

    Returns:
        lattice (Lattice):  New :py:class:`.Lattice` object
    """
    parser = Mad8Parser(strict=strict, verbose=verbose)
    params = {
        key: kwargs.pop(key)
        for key in ("name", "particle", "energy", "periodicity")
        if key in kwargs
    }
    parser.parse_files(*files, **kwargs)
    return parser.lattice(use=use, **params)


class _Mad8Exporter(_MadExporter):
    delimiter = ""
    continuation = "&"
    bool_fmt = {False: ".FALSE.", True: ".TRUE."}


def save_mad8(ring: Lattice, filename: str | None = None, **kwargs):
    """Save a :py:class:`.Lattice` as a MAD8 file

    Args:
        ring:   lattice
        filename: file to be created. If None, write to sys.stdout

    Keyword Args:
        use (str | None): name of the created SEQUENCE of LINE.
          Default: name of the PyAT lattice
        use_line (bool):  If True, use a MAD "LINE" format. Otherwise, use
          a MAD "SEQUENCE"
    """
    exporter = _Mad8Exporter(ring, **kwargs)
    exporter.export(filename)
