"""
Helper functions for working with AT lattices.

A :py:class:`.Lattice` in pyAT is a sequence of :py:class:`.Element` objects.
These functions are useful for building ad manipulating these sequences.
"""
import sys
from .options import DConstant, random
from .particle_object import Particle
from .elements import *
from .utils import *
from .lattice_object import *
from .cavity_access import *
from .variable_elements import *
from .deprecated import *
# Define the module "lattice.constants" for backward compatibility
# noinspection PyUnresolvedReferences
from .. import constants
sys.modules['at.lattice.constants'] = sys.modules['at.constants']
