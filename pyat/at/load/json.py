"""
Handling of JSON files
"""

from __future__ import annotations

__all__ = ["save_json", "load_json"]

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
            return obj.definition
        elif isinstance(obj, Particle):
            return obj.to_dict()
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
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
    data = dict(
        elements=[element_to_dict(el) for el in save_filter(ring)],
        parameters=ring.attrs,
    )
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

    def json_generator(params: dict[str, Any], elem_list):
        particle_dict = params.pop("particle", {})
        params["particle"] = Particle(**particle_dict)
        for idx, elem_dict in enumerate(elem_list):
            yield element_from_dict(elem_dict, index=idx, check=False)

    with open(filename, "rt") as jsonfile:
        data = json.load(jsonfile)

    try:
        elements = data["elements"]
        parameters = data["parameters"]
    except KeyError:
        raise TypeError("Not a Lattice")
    if not (isinstance(elements, list) and isinstance(parameters, dict)):
        raise TypeError("Not a lattice")

    parameters.update(kwargs)
    return Lattice(elements, iterator=json_generator, **parameters)


register_format(
    ".json",
    load_json,
    save_json,
    descr="JSON representation of a python AT Lattice",
)
