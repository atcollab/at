"""Conversion from AT to xtrack.

These functions are available if xsuite is installed.
"""
__all__ = [
    "element_from_xtrack",
    "element_to_xtrack",
    "lattice_from_line",
    "line_from_lattice",
]

import xtrack as xt
from . import xsuite as xs
from ..lattice import Lattice, Element


_xtclass: dict[type[xs.XsElement], type[xt.BeamElement]] = {
    xs.Marker: xt.Marker,
    xs.Drift: xt.Drift,
    xs.Cavity: xt.Cavity,
    xs.Quadrupole: xt.Quadrupole,
    xs.Sextupole: xt.Sextupole,
    xs.Octupole: xt.Octupole,
    xs.Multipole: xt.Multipole,
    xs.Bend: xt.Bend,
    xs.RBend: xt.RBend,
}


def element_from_xtrack(xtelem: xt.BeamElement) -> Element:
    """Create an :py:class:`.xtrack.BeamElement` from an AT :py:class:`.Element`.

        Returns:
            atelem:   new :py:class:`.Element` object
    """
    return xs.XsElement.static_from_dict(xtelem.to_dict()).to_at()


def element_to_xtrack(atelem: Element, match_model: bool = True) -> xt.BeamElement:
    """Create an AT :py:class:`.Element` from an :py:class:`.xtrack.BeamElement`.

        Args:
            atelem:         AT :py:class:`.Element`
            match_model:    If :py:obj:`True`, set the Xsuite model and integrator
              matching at best the AT PassMethod. Otherwise, use Xsuite defaults,

        Returns:
            xtelem:  new :py:class:`.xtrack.BeamElement` object
    """
    xselem = xs.XsElement.static_from_at(atelem, match_model=match_model)
    xtclass = _xtclass[type(xselem)]
    return xtclass.from_dict(xselem.to_dict())


def lattice_from_line(line: xt.Line, use: str | None = None, **kwargs) -> Lattice:
    """Create a :py:class:`.Lattice` object from a Xsuite :py:class:`.Line` or
    :py:class:`.Environment`.

    Args:
        line: Xsuite line or environment
        use (str):  Name of the line containing the desired lattice.
          Default: ``ring``.

    Keyword Args:
        name (str):         Name of the lattice. Default: Xsuite line name.
        particle(Particle): Circulating particle. Default: Xsuite reference
          particle
        energy (float):     Energy of the lattice [eV]. Default: from Xsuite
          reference particle
        periodicity(int):   Number of periods. Default: 1
        *:                  All other keywords will be set as Lattice attributes

    Returns:
        lattice:    New :py:class:`.Lattice` object.
    """
    return xs.XsLine.from_xsuite(line, use=use).to_at(**kwargs)


def line_from_lattice(ring: Lattice, match_model: bool = True) -> xt.Line:
    """Create a Xsuite :py:class:`.Line` from an AT lattice.

    Args:
        ring:           AT lattice
        match_model:    if :py:obj:`True`, set the Xsuite model and integrator
          matching at best the AT PassMethod. if :py:obj:`False`, use Xsuite defaults.

    Returns:
        line:           new :py:class:`.Line` object.
    """
    return xt.Line.from_dict(xs.XsLine.from_at(ring, match_model=match_model).to_dict())


Element.to_xtrack = element_to_xtrack
Element.from_xtrack = staticmethod(element_from_xtrack)
Lattice.to_xtrack = line_from_lattice
Lattice.from_xtrack = staticmethod(lattice_from_line)
