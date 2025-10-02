r"""Using `MAD8`_ files with PyAT
=================================

Using MAD8 files is similar to using MAD-X files: see
:py:mod:`the MAD-X documentation <at.load.madx>`

.. _mad8: https://mad8.web.cern.ch/user/mad.html
"""

from __future__ import annotations

__all__ = ["Mad8Parser", "load_mad8", "save_mad8"]

# functions known by MAD-8
from math import pi, e, sqrt, exp, log, sin, cos, tan
from math import asin, atan2
import re

# constants known by MAD-8
from scipy.constants import c as clight, e as qelect
from scipy.constants import physical_constants as _cst

from ..lattice import Lattice, elements as elt

from .file_input import ignore_class

# noinspection PyProtectedMember
from .madx import _MadElement, _MadParser, _MadExporter
from .madx import mad_element, p_to_at, poly_to_mad

# Commands known by MAD8
# noinspection PyProtectedMember
from .madx import (
    drift,
    marker,
    quadrupole,
    sextupole,
    octupole,
    longmultipole,
    sbend,
    rbend,
    kicker,
    hkicker,
    vkicker,
    rfcavity,
    monitor,
    hmonitor,
    vmonitor,
    _Sequence,
    _value,
)

twopi = 2 * pi
degrad = 180.0 / pi
raddeg = pi / 180.0
emass = 1.0e-03 * _cst["electron mass energy equivalent in MeV"][0]  # [GeV]
pmass = 1.0e-03 * _cst["proton mass energy equivalent in MeV"][0]  # [GeV]

_attr = re.compile(r"\[([a-zA-Z_][\w.:]*)]")  # Identifier enclosed in square brackets


_mad8_env = {
    # Constants
    "_true_": True,
    "_yes_": True,
    "_t_": True,
    "_on_": True,
    "_false_": False,
    "_no_": False,
    "_f_": False,
    "_off_": False,
    "pi": pi,
    "e": e,
    "abs": abs,
    "sqrt": sqrt,
    "exp": exp,
    "log": log,
    "sin": sin,
    "cos": cos,
    "tan": tan,
    "asin": asin,
    "twopi": 2 * pi,
    "degrad": 180.0 / pi,
    "raddeg": pi / 180.0,
    "emass": emass,  # [GeV]
    "pmass": pmass,  # [GeV]
    "clight": clight,
    "qelect": qelect,
    "centre": "centre",
    "entry": "entry",
    "exit": "exit",
    # Elements
    "drift": drift,
    "marker": marker,
    "quadrupole": quadrupole,
    "sextupole": sextupole,
    "octupole": octupole,
    "longmultipole": longmultipole,
    "sbend": sbend,
    "rbend": rbend,
    "kicker": kicker,
    "hkicker": hkicker,
    "vkicker": vkicker,
    "rfcavity": rfcavity,
    "monitor": monitor,
    "hmonitor": hmonitor,
    "vmonitor": vmonitor,
    "sequence": _Sequence,
    # Commands
    "value": _value,
    "__builtins__": {},
}


_ignore_names = [
    "solenoid",
    "rfmultipole",
    "crabcavity",
    "elseparator",
    "collimator",
    "tkicker",
]

_mad8_env.update(
    (name, ignore_class(name, _MadElement, module=__name__)) for name in _ignore_names
)


# noinspection PyPep8Naming
class multipole(_MadElement):
    at2mad = {}
    klist = ["k0l", "k1l", "k2l", "k3l", "k4l", "k5l", "k6l", "k7l", "k8l", "k9l"]
    tlist = ["t0", "t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8", "t9"]

    @mad_element
    def to_at(self, **params):
        polya = [0.0]*10
        polyb = [0.0]*10
        params.pop("l", None)
        for k in list(params.keys()):
            try:
                order = self.klist.index(k)
            except ValueError:
                pass
            else:
                strength = params.pop(k)
                angle = (order + 1) * params.get(self.tlist[order], 0.0)
                polyb[order] = strength * cos(-angle)
                polya[order] = strength * sin(-angle)
        return [
            elt.ThinMultipole(
                self.name, p_to_at(polya), p_to_at(polyb), **self.meval(params)
            )
        ]

    @classmethod
    def from_at(cls, kwargs, factor=1.0):
        el = super().from_at(kwargs)
        maxorder = kwargs.pop("MaxOrder", -1) + 1
        nlist = poly_to_mad(kwargs.pop("PolynomB", ())[:maxorder+1], factor=factor)
        slist = poly_to_mad(kwargs.pop("PolynomA", ())[:maxorder+1], factor=factor)
        for order, (va, vb) in enumerate(zip(slist, nlist)):
            if va != 0.0 or vb != 0.0:
                el[multipole.klist[order]] = sqrt(va * va + vb * vb)
                tilt = -atan2(va, vb) / (order + 1)
                if tilt != 0.0:
                    el[multipole.tlist[order]] = tilt
        return el


_mad8_env.update(multipole=multipole)


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

    def __init__(self, *filenames: str, **kwargs):
        """
        Args:
            *filenames: files to be read at initialisation
            strict:     If :py:obj:`False`, assign 0 to undefined variables
            verbose:    If :py:obj:`True`, print details on the processing
            **kwargs:   Initial variable definitions
        """
        super().__init__(_mad8_env, *filenames, **kwargs)

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
    parameterised: bool = False,
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
    return parser.lattice(use=use, parameterised=parameterised, **params)


class _Mad8Exporter(_MadExporter):
    delimiter = ""
    continuation = "&"
    bool_fmt = {False: ".FALSE.", True: ".TRUE."}

    AT2MAD = {
        elt.Quadrupole: quadrupole,
        elt.Sextupole: sextupole,
        elt.Octupole: octupole,
        elt.ThinMultipole: multipole,
        elt.Multipole: longmultipole,
        elt.RFCavity: rfcavity,
        elt.Drift: drift,
        elt.Bend: sbend,
        elt.Marker: marker,
        elt.Monitor: monitor,
        elt.Corrector: kicker,
    }


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
