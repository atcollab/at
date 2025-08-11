"""Text representation of a python AT lattice with each element represented by
its :py:func:`repr` string
"""

from __future__ import annotations

__all__ = ["load_repr", "save_repr"]

import sys
from os.path import abspath
from typing import Optional

import numpy as np

# imports necessary in 'globals()' for 'eval'
from numpy import array, uint8, nan as NaN  # noqa: F401

from at.lattice import Lattice, Element

# imports necessary in' globals()' for 'eval'
from at.lattice import Particle  # noqa: F401
from at.load import register_format
from at.load.utils import element_classes, keep_attributes, keep_elements

# Map class names to Element classes
_CLASS_MAP = dict((cls.__name__, cls) for cls in element_classes())


def _element_from_string(elem_string: str) -> Element:
    """Builds an :py:class:`.Element` from its python :py:func:`repr` string

    Parameters:
        elem_string:    String representation of an :py:class:`.Element`

    Returns:
        elem (Element): new :py:class:`.Element`
    """
    return eval(elem_string, globals(), _CLASS_MAP)


def load_repr(filename: str, **kwargs) -> Lattice:
    """Create a :py:class:`.Lattice`  from a text repr-file

    Parameters:
        filename:           Name of a '.m' file

    Keyword Args:
        name (str):         Name of the lattice. Default: taken from
          the lattice
        energy (float):     Energy of the lattice [eV]. Default: taken
          from the lattice elements
        periodicity(int):   Number of periods. Default: taken from the
          elements, or 1
        *:                  All other keywords will be set as Lattice
          attributes

    Returns:
        lattice (Lattice):          New :py:class:`.Lattice` object

    See Also:
        :py:func:`.load_lattice` for a generic lattice-loading function.
    """

    def elem_iterator(params, repr_file):
        with open(params.setdefault("in_file", repr_file), "rt") as file:
            # the 1st line is the dictionary of saved lattice parameters
            for k, v in eval(next(file)).items():
                params.setdefault(k, v)
            for line in file:
                yield _element_from_string(line.strip())

    return Lattice(abspath(filename), iterator=elem_iterator, **kwargs)


def save_repr(ring: Lattice, filename: Optional[str] = None) -> None:
    """Save a :py:class:`.Lattice` as a repr-file

    Parameters:
        ring:           Lattice description
        filename:       Name of the '.repr' file. Default: outputs on
          :py:obj:`sys.stdout`

    See Also:
        :py:func:`.save_lattice` for a generic lattice-saving function.
    """

    def save(file):
        print(repr(keep_attributes(ring)), file=file)
        for elem in keep_elements(ring):
            print(repr(elem), file=file)

    # Set options to print the full representation of float variables
    with np.printoptions(formatter={"float_kind": repr}):
        if filename is None:
            save(sys.stdout)
        else:
            with open(filename, "wt") as reprfile:
                save(reprfile)


register_format(
    ".repr",
    load_repr,
    save_repr,
    descr=("Text representation of a python AT Lattice. "
           "See :py:func:`.load_repr`."),
)
