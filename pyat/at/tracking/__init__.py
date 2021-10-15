"""
Tracking functions
"""

# noinspection PyUnresolvedReferences,PyProtectedMember
from .atpass import _atpass, _elempass, isopenmp
from .patpass import patpass
from .track import *
from .particles import *
# noinspection PyProtectedMember
from ..lattice.options import _Dst


# noinspection PyUnusedLocal
def _omp(self):
    return isopenmp()


_Dst.omp = property(_omp, doc="True is AT is compiled with OpenMP")
