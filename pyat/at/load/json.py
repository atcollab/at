from __future__ import annotations

import json
from typing import Optional, Any

import numpy as np

from . import register_format
from ..lattice import Element, Lattice, Particle, get_class_map

_CLASS_MAP = get_class_map()


def elemstr(self):
    attrs = dict(self.items())
    return "\n".join(
        [self.__class__.__name__ + ":"]
        + [f"{k:>14}: {attrs.pop(k)!s}" for k in ["FamName", "Length", "PassMethod"]]
        + [f"{k:>14}: {v!s}" for k, v in attrs.items()]
    )


class _AtEncoder(json.JSONEncoder):
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
    if filename is None:
        json.dumps(("Lattice", ring, ring.attrs), cls=_AtEncoder)
    else:
        with open(filename, "wt") as jsonfile:
            json.dump(("Lattice", ring, ring.attrs), jsonfile, cls=_AtEncoder)


def load_json(filename: str, **kwargs) -> Lattice:
    """Create a :py:class:`.Lattice`  from a JSON file

    Parameters:
        filename:           Name of a JSON file

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
        :py:meth:`.Lattice.load` for a generic lattice-loading method.
    """

    def json_generator(params: dict[str, Any], elem_list):
        particle_dict = params.pop("particle", {})
        params["particle"] = Particle(**particle_dict)
        for clname, args, keys in elem_list:
            cls = _CLASS_MAP[clname]
            yield cls(*args, **keys)

    with open(filename, "rt") as jsonfile:
        data = json.load(jsonfile)

    try:
        code, elems, prms = data
    except ValueError:
        raise TypeError("Not a Lattice")
    if not (
        isinstance(code, str)
        and isinstance(elems, list)
        and isinstance(prms, dict)
        and (code == "Lattice")
    ):
        raise TypeError("Not a lattice")

    prms.update(kwargs)
    return Lattice(elems, iterator=json_generator, **prms)


register_format(
    ".json",
    load_json,
    save_json,
    descr="JSON representation of a python AT Lattice",
)
