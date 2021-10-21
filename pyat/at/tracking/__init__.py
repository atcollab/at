"""
Tracking functions
"""

# noinspection PyUnresolvedReferences
from .atpass import _atpass, _elempass, isopenmp
from .patpass import patpass
from .track import *
from .particles import *
# noinspection PyProtectedMember
from ..lattice.options import _Dst


# noinspection PyUnusedLocal
def _omp(self):
    """True is AT is compiled with OpenMP"""
    return isopenmp()


_Dst.openmp = property(_omp)
