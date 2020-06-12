"""Text representation of a python AT lattice with each element represented by
its 'repr' string
"""
from __future__ import print_function
import sys
from os.path import abspath
import numpy
from at.lattice import Lattice
from at.load import register_format
from at.load.utils import element_from_string

__all__ = ['load_repr', 'save_repr']


def load_repr(filename, **kwargs):
    """Load a Lattice object stored as a repr-file

    PARAMETERS
        filename        name of the '.repr' file

    KEYWORDS
        name            Name of the lattice (default: taken from the file)
        energy          Energy of the lattice (default: taken from the file)
        periodicity     Number of periods (default: taken from the file)
        *               all other keywords will be set as Lattice attributes

    OUTPUT
        Lattice object
    """
    def elem_iterator(params, repr_file):
        with open(params.setdefault('repr_file', repr_file), 'rt') as file:
            # the 1st line is the dictionary of saved lattice parameters
            for k, v in eval(next(file)).items():
                params.setdefault(k, v)
            for line in file:
                yield element_from_string(line.strip())

    return Lattice(abspath(filename), iterator=elem_iterator, **kwargs)


def save_repr(ring, filename=None):
    """Save a Lattice object as a repr-file

    PARAMETERS
        ring            Lattice object
        filename=None   name of the '.repr' file. Default: output on sys.stdout
    """
    def save(file):
        print(repr(dict((k, v) for k, v in vars(ring).items()
                        if not k.startswith('_'))), file=file)
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
