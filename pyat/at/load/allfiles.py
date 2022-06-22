"""Generic function to save and load python AT lattices. The format is
determined by the file extension
"""
import os.path
from typing import Optional
from at.lattice import Lattice

__all__ = ['load_lattice', 'save_lattice', 'register_format']

_load_extension = {}
_save_extension = {}


def load_lattice(filepath: str, **kwargs) -> Lattice:
    """Load a Lattice object from a file

The file format is indicated by the filepath extension.

Parameters:
    filepath:        name of the file

Keyword Args:
    name (Optional[str]):           Name of the lattice
      . Default: taken from the file, or ``''``
    energy (Optional[float]):       Energy of the lattice
                    (default: taken from the file)
    periodicity (Optional[int]):    Number of periods
                    (default: taken from the file, or 1)
    *:               all other keywords will be set as Lattice attributes

Specific keywords for .mat files

Keyword Args:
    mat_key (Optional[str]):        Name of the Matlab variable containing
      the lattice. Default: Matlab variable name if there is only one,
      otherwise ``'RING'``
    check (Optional[bool]):     If Default: :py:obj:`False`, skip the
      coherence tests. Default: :py:obj:`True`
    quiet (Optional[bool]):     If :py:obj:`True`, suppress the warning for
      non-standard classes. Default: :py:obj:`False`
    keep_all (Optional[bool]):  Keep Matlab RingParam elements as Markers.
      Default: :py:obj:`False`

Returns:
    lattice (Lattice):          New :py:class:`.Lattice` object

.. Admonition:: Known extensions are:
    """
    _, ext = os.path.splitext(filepath)
    try:
        load_func = _load_extension[ext.lower()]
    except KeyError:
        print("File load failed: unknow extension {}.".format(ext))
    else:
        return load_func(filepath, **kwargs)


def save_lattice(ring: Lattice, filepath: str, **kwargs):
    """Save a Lattice object

The file format is indicated by the filepath extension.

Parameters:
    ring:           Lattice description
    filepath:       Name of the file

Specific keywords for .mat files

Keyword Args:
    mat_key (Optional[str]):        Name of the Matlab variable containing
      the lattice. Default: ``'RING'``

.. Admonition:: Known extensions are:
    """
    _, ext = os.path.splitext(filepath)
    try:
        save_func = _save_extension[ext.lower()]
    except KeyError:
        print("File save failed: unknow extension {}.".format(ext))
    else:
        return save_func(ring, filepath, **kwargs)


def register_format(extension: str, load_func=None, save_func=None,
                    descr: Optional[str] = ''):
    """Register format-specific processing functions

    Parameters:
        extension:      File extension string
        load_func:      load function (default: None)
        save_func:      save_lattice function (default: None)
        descr:          File type description
    """
    if load_func is not None:
        _load_extension[extension] = load_func
        load_lattice.__doc__ += '\n    {0:<10}'\
                                '\n        {1}\n'.format(extension, descr)
    if save_func is not None:
        _save_extension[extension] = save_func
        save_lattice.__doc__ += '\n    {0:<10}'\
                                '\n        {1}\n'.format(extension, descr)


Lattice.load = staticmethod(load_lattice)
Lattice.save = save_lattice
