"""Stub file for the 'diffmatrix' extension"""

import numpy as np
from typing import Optional
from at.lattice import Element
from . import Orbit

def find_mpole_raddiff_matrix(element: Element, orbit: Orbit,
                              energy: Optional[float] = None) -> np.ndarray: ...