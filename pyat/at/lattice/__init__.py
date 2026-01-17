"""Helper functions for working with AT lattices.

A :py:class:`.Lattice` in pyAT is a sequence of :py:class:`.Element` objects.
These functions are useful for building ad manipulating these sequences.
"""
import sys
import numpy as np
from .exceptions import *
from .axisdef import *
from .options import *
from .particle_object import Particle
from .elements import *
from .utils import *
from .variables import *
from .geometry import *
from .transformation import *
from .lattice_object import *
# from .lattice_variables import *
from .cavity_access import *
from .deprecated import *
# Define the module "lattice.constants" for backward compatibility
# noinspection PyUnresolvedReferences
from .. import constants
sys.modules['at.lattice.constants'] = sys.modules['at.constants']

# Type definitions
Orbit = np.ndarray
