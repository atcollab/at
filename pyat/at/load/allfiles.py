"""Generic function to save and load python AT lattices. The format is
determined by the file extension
"""
import os.path
from at.lattice import Lattice

__all__ = ['load_lattice', 'save_lattice', 'register_format']

_load_extension = {}
_save_extension = {}


def load_lattice(filepath: str, **kwargs) -> Lattice:
    """Load a Lattice object from a file

The file format is indicated by the filepath extension. The file name is stored in
the *in_file* Lattice attribute. The selected variable, if relevant, is stored
in the *use* Lattice attribute.

Parameters:
    filepath:           Name of the file

Keyword Args:
    use (str):          Name of the variable containing the desired lattice.
      Default: if there is a single variable, use it, otherwise select ``"RING"``
    name (str):         Name of the lattice.
      Default: taken from the file, or ``""``
    energy (float):     Energy of the lattice
      (default: taken from the file)
    periodicity (int):  Number of periods
      (default: taken from the file, or 1)
    *:                  All other keywords will be set as :py:class:`.Lattice`
      attributes

Check the format-specific function for specific keyword arguments.

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
    ring:               Lattice description
    filepath:           Name of the file

Specific keywords for .mat files

Keyword Args:
    mat_key (str):      Name of the Matlab variable containing the lattice.
      Default: ``'RING'``

See Also:
    :py:func:`.save_mat`, :py:func:`.save_m`, :py:func:`.save_repr`

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
                    descr: str = ''):
    """Register format-specific processing functions

    Parameters:
        extension:      File extension string.
        load_func:      load function. Default: :py:obj:`None`
        save_func:      save_lattice function Default: :py:obj:`None`
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
