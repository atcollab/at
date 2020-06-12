"""Generic function to save and load python AT lattices. The format is
determined by the file extension
"""
import os.path
from at.lattice import Lattice

__all__ = ['load_lattice', 'save_lattice', 'register_format']

_load_extension = {}
_save_extension = {}


def load_lattice(filepath, **kwargs):
    """Load a Lattice object from a file

    The file format is indicated by the filepath extension.

    PARAMETERS
        filepath        name of the file

    KEYWORDS
        name            Name of the lattice
                        (default: taken from the file, or '')
        energy          Energy of the lattice
                        (default: taken from the file)
        periodicity     Number of periods
                        (default: taken from the file, or 1)
        *               all other keywords will be set as Lattice attributes

    MAT-FILE SPECIFIC KEYWORDS
        mat_key         name of the Matlab variable containing the lattice.
                        Default: Matlab variable name if there is only one,
                        otherwise 'RING'
        check=True      if False, skip the coherence tests
        quiet=False     If True, suppress the warning for non-standard classes
        keep_all=False  if True, keep RingParam elements as Markers

    Known extensions are:
    """
    _, ext = os.path.splitext(filepath)
    try:
        return _load_extension[ext.lower()](filepath, **kwargs)
    except KeyError:
        print("Could not load lattice file with extension {}.".format(ext))


def save_lattice(ring, filepath, **kwargs):
    """Save a Lattice object

    The file format is indicated by the filepath extension.

    PARAMETERS
        ring            Lattice object
        filepath        name of the file

    MAT-FILE SPECIFIC KEYWORDS
        mat_key='RING'  Name of the Matlab variable

    Known extensions are:
    """
    _, ext = os.path.splitext(filepath)
    try:
        return _save_extension[ext.lower()](ring, filepath, **kwargs)
    except KeyError:
        print("Could not save lattice file with extension {}.".format(ext))


def register_format(extension, load_func=None, save_func=None, descr=''):
    """Register format-specific processing functions

    PARAMETERS
        extension       File extension string
        load_func       load function (default: None)
        save_func       save_lattice function (default: None)
        descr           File type description
    """
    if load_func is not None:
        _load_extension[extension] = load_func
        load_lattice.__doc__ += '\n    {0:<10}\t{1}'.format(extension, descr)
    if save_func is not None:
        _save_extension[extension] = save_func
        save_lattice.__doc__ += '\n    {0:<10}\t{1}'.format(extension, descr)


Lattice.load = staticmethod(load_lattice)
Lattice.save = save_lattice
