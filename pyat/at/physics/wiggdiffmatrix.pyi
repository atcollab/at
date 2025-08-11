"""Stub file for the 'diffmatrix' extension"""

from __future__ import annotations

import numpy as np

from . import Orbit
from ..lattice import Element

def FDW(element: Element, orbit: Orbit, energy: float | None = None) -> np.ndarray: ...
