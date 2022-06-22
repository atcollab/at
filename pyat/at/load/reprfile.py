"""Text representation of a python AT lattice with each element represented by
its :py:func:`repr` string
"""
from __future__ import print_function
import sys
from os.path import abspath
from typing import Optional
import numpy
from at.lattice import Lattice
from at.load import register_format
from at.load.utils import element_from_string
# imports necessary in' globals()' for 'eval'
# noinspection PyUnresolvedReferences
from at.lattice import Particle
from numpy import array
__all__ = ['load_repr', 'save_repr']


def load_repr(filename: str, **kwargs) -> Lattice:
    """Create a :py:class:`.Lattice`  from a text repr-file

    Parameters:
        filename:       name of a '.m' file

    Keyword Args:
        name (Optional[str]):       Name of the lattice. Default: taken from
          the lattice
        energy (Optional[float]):   Energy of the lattice [eV]. Default: taken
          from the lattice elements
        periodicity(Optional[int]): Number of periods. Default: taken from the
          elements, or 1
        *:                          All other keywords will be set as Lattice
          attributes

    Returns:
        lattice (Lattice):          New :py:class:`.Lattice` object
    """
    def elem_iterator(params, repr_file):
        with open(params.setdefault('repr_file', repr_file), 'rt') as file:
            # the 1st line is the dictionary of saved lattice parameters
            for k, v in eval(next(file)).items():
                params.setdefault(k, v)
            for line in file:
                yield element_from_string(line.strip())

    return Lattice(abspath(filename), iterator=elem_iterator, **kwargs)


def save_repr(ring: Lattice, filename: Optional[str] = None) -> None:
    """Save a :py:class:`.Lattice` as a repr-file

    Parameters:
        ring:           Lattice description
        filename:       Name of the '.repr' file. Default: outputs on
          :py:obj:`sys.stdout`
    """
    def save(file):
        # print(repr(dict((k, v) for k, v in vars(ring).items()
        #                 if not k.startswith('_'))), file=file)
        print(repr(ring.attrs), file=file)
        for elem in ring:
            print(repr(elem), file=file)

    # Save the current options
    opts = numpy.get_printoptions()
    # Set options to print the full representation of float variables
    numpy.set_printoptions(formatter={'float_kind': repr})
    if filename is None:
        save(sys.stdout)
    else:
        with open(filename, 'wt') as reprfile:
            save(reprfile)
    # Restore the current options
    numpy.set_printoptions(**opts)


register_format('.repr', load_repr, save_repr,
                descr='Text representation of a python AT Lattice')
