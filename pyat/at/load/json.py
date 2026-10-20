"""
Handling of JSON files.
"""

from __future__ import annotations

__all__ = ["load_json", "save_json"]

from pathlib import Path
import json
from typing import Any

import numpy as np

from .allfiles import register_format
from .utils import keep_elements, keep_attributes, element_from_dict
from ..lattice import Element, Lattice, Particle
from .xsuite import XsLine
from .._version import __version_tuple__


class _AtEncoder(json.JSONEncoder):
    """JSON encoder for specific AT types."""

    def default(self, obj):
        if isinstance(obj, Element):
            return obj.to_file()
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, Particle):
            return obj.to_dict()
        elif isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        else:
            return super().default(obj)


def save_json(
    ring: Lattice, filename: str | Path | None = None, compact: bool = False
) -> None:
    """Save a :py:class:`.Lattice` as a JSON file.

    Parameters:
        ring:           Lattice description
        filename:       Name of the JSON file. Default: outputs on
          :py:obj:`sys.stdout`
        compact:        If :py:obj:`False` (default), the JSON file is pretty-printed
          with line feeds and indentation. Otherwise, the output is a single line.

    See Also:
        :py:func:`.save_lattice` for a generic lattice-saving function.
        :py:meth:`.Lattice.save` for a generic lattice-saving method.
    """
    indent = None if compact else 2
    data = {
        "atjson": 1,
        "at_version": ".".join(str(n) for n in __version_tuple__[:3]),
        "elements": list(keep_elements(ring)),
        "properties": keep_attributes(ring),
    }
    if filename is None:
        print(json.dumps(data, cls=_AtEncoder, indent=indent))
    else:
        filename = Path(filename)
        with filename.open("w") as jsonfile:
            json.dump(data, jsonfile, cls=_AtEncoder, indent=indent)


def _load_at(root: dict[str, Any], **kwargs) -> Lattice:
    """Create a :py:class:`.Lattice`  from a JSON file.

    Parameters:
        filename:           Name of a JSON file

    Keyword Args:
        *:                  All keywords update the lattice properties

    Returns:
        lattice (Lattice):          New :py:class:`.Lattice` object

    See Also:
        :py:meth:`.Lattice.load` for a generic lattice-loading method.
    """

    def json_generator(params: dict[str, Any], data: dict[str, Any]):
        # Check the file signature - For later use
        try:
            # noinspection PyUnusedLocal
            atjson = data["atjson"]
        except KeyError:
            atjson = 1  # noqa: F841
        # Get elements
        elements = data["elements"]
        # Get lattice properties
        try:
            properties = data["properties"]
        except KeyError:
            properties = {}
        particle_dict = properties.pop("particle", {})
        params.setdefault("particle", Particle(**particle_dict))

        for k, v in properties.items():
            params.setdefault(k, v)
        for idx, elem_dict in enumerate(elements):
            yield element_from_dict(elem_dict, index=idx, check=False)

    return Lattice(root, iterator=json_generator, **kwargs)


def load_json(
    filename: str | Path,
    use: str | None = None,
    from_at: bool = False,
    from_xsuite: bool = False,
    **kwargs,
) -> Lattice:
    """Create a :py:class:`.Lattice`  from a AT or Xsuite JSON file.

    The kind of file is derived from its contents. In case of ambiguity, the file
    kind can be explicitly selected with the *from_at* or *from_xsuite* keywords.

    Parameters:
        filename:           Name of a JSON file
        from_at:            Force the selection of AT json file
        from_xsuite:        Force the selection of Xsuite json file
        use:                Line name, for Xsuite files containing several lines.
          Default: ``ring``

    Keyword Args:
        *:                  All keywords update the lattice properties

    Returns:
        lattice (Lattice):          New :py:class:`.Lattice` object

    See Also:
        :py:meth:`.Lattice.load` for a generic lattice-loading method.
    """
    filename = Path(filename)
    with filename.open() as jsonfile:
        data = json.load(jsonfile)
        if from_at or ("atjson" in data):
            return _load_at(data, in_file=str(filename), **kwargs)
        elif from_xsuite or ("__class__" in data):
            return XsLine.from_dict(data, use=use).to_at(**kwargs)
        else:
            msg = "Cannot guess the file origin, try 'from_at or 'from_xsuite' keywords"
            raise TypeError(msg)


register_format(
    ".json",
    load_json,
    save_json,
    descr="JSON representation of a python AT Lattice. See :py:func:`.load_json`.",
)
