"""This file contains common grid definitions and methods."""

import numpy as np

from at.lattice import AtError

_pdict = {"x": 0, "xp": 1, "y": 2, "yp": 3, "dp": 4, "ct": 5}


def get_plane_index(planes) -> np.array:
    """
    Convert plane to particle coordinate index.

    Args:
        planes: an iterable tuple or list containing the plane names.
          'x'. 'xp', 'y', 'yp', 'dp', 'ct'.

    Returns:
        AT index of the plane.

    Raises:
        Error when the list contains any element not above.
    """
    planesi = np.array([], dtype=np.int32)
    msg = "Allowed values for plane are " + ", ".join(key for key in _pdict)
    for p in np.atleast_1d(planes):
        if isinstance(p, str):
            try:
                planesi = np.append(planesi, _pdict[p])
            except KeyError as keyerr:
                raise AtError(msg) from keyerr
        else:
            raise AtError(msg) from None
    return planesi
