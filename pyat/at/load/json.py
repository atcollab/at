"""
Handling of JSON files
"""

from __future__ import annotations

__all__ = ["save_json", "load_json"]

from os.path import abspath
import json
from typing import Optional, Any

import numpy as np

from .allfiles import register_format
from .utils import element_to_dict, element_from_dict, save_filter
from ..lattice import Element, Lattice, Particle


class _AtEncoder(json.JSONEncoder):
    """JSON encoder for specific AT types"""

    def default(self, obj):
        if isinstance(obj, Element):
            return element_to_dict(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, Particle):
            return obj.to_dict()
        else:
            return super().default(obj)


def save_json(ring: Lattice, filename: Optional[str] = None) -> None:
    """Save a :py:class:`.Lattice` as a JSON file

    Parameters:
        ring:           Lattice description
        filename:       Name of the JSON file. Default: outputs on
          :py:obj:`sys.stdout`

    See Also:
        :py:func:`.save_lattice` for a generic lattice-saving function.
        :py:meth:`.Lattice.save` for a generic lattice-saving method.
    """
    data = dict(elements=list(save_filter(ring)), properties=ring.attrs)
    if filename is None:
        print(json.dumps(data, cls=_AtEncoder, indent=2))
    else:
        with open(filename, "wt") as jsonfile:
            json.dump(data, jsonfile, cls=_AtEncoder, indent=2)


def load_json(filename: str, **kwargs) -> Lattice:
    """Create a :py:class:`.Lattice`  from a JSON file

    Parameters:
        filename:           Name of a JSON file

    Keyword Args:
        *:                  All keywords update the lattice properties

    Returns:
        lattice (Lattice):          New :py:class:`.Lattice` object

    See Also:
        :py:meth:`.Lattice.load` for a generic lattice-loading method.
    """

    def json_generator(params: dict[str, Any], fn):

        with open(params.setdefault("json_file", fn), "rt") as jsonfile:
            data = json.load(jsonfile)
        elements = data["elements"]
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

    return Lattice(abspath(filename), iterator=json_generator, **kwargs)


register_format(
    ".json",
    load_json,
    save_json,
    descr="JSON representation of a python AT Lattice",
)
