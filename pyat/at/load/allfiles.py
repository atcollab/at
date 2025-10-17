"""Generic function to save and load python AT lattices. The format is
determined by the file extension.
"""

from __future__ import annotations

__all__ = ["load_lattice", "register_format", "save_lattice"]

from pathlib import Path
from collections.abc import Callable

from at.lattice import Lattice

_load_extension: dict[str, Callable[..., Lattice]] = {}
_save_extension: dict[str, Callable[..., None]] = {}


def load_lattice(filepath: str | Path, **kwargs) -> Lattice:
    """Load a Lattice object from a file.

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

    Returns:
        lattice (Lattice):          New :py:class:`.Lattice` object

    Check the format-specific function for specific keyword arguments:

    .. Admonition:: Known extensions are:
    """
    filepath = Path(filepath)
    ext = filepath.suffix
    try:
        load_func = _load_extension[ext.lower()]
    except KeyError:
        print(f"File load failed: unknow extension {ext}.")
    else:
        return load_func(filepath, **kwargs)


def save_lattice(ring: Lattice, filepath: str | Path, **kwargs) -> None:
    """Save a Lattice object.

    The file format is indicated by the filepath extension.

    Parameters:
        ring:               Lattice description
        filepath:           Name of the file

    Check the format-specific function for specific keyword arguments:

    .. Admonition:: Known extensions are:
    """
    filepath = Path(filepath)
    ext = filepath.suffix
    try:
        save_func = _save_extension[ext.lower()]
    except KeyError:
        print(f"File save failed: unknow extension {ext}.")
    else:
        return save_func(ring, filepath, **kwargs)


def register_format(
    extension: str,
    load_func: Callable[..., Lattice] | None = None,
    save_func: Callable[..., None] | None = None,
    descr: str = "",
):
    """Register format-specific processing functions.

    Parameters:
        extension:      File extension string.
        load_func:      load function.
        save_func:      save function.
        descr:          File type description.
    """
    if load_func is not None:
        _load_extension[extension] = load_func
        load_lattice.__doc__ += f"\n        {extension:<10}\n            {descr}\n"
    if save_func is not None:
        _save_extension[extension] = save_func
        save_lattice.__doc__ += f"\n        {extension:<10}\n            {descr}\n"


Lattice.load = staticmethod(load_lattice)
Lattice.save = save_lattice
