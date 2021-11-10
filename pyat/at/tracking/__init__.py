"""
Tracking functions
"""

# noinspection PyUnresolvedReferences
from .atpass import atpass, elempass, isopenmp, ismpi
from .patpass import patpass
from .track import *
from .particles import *
# noinspection PyProtectedMember
from ..lattice.options import _Dst


# noinspection PyUnusedLocal
def _omp(self):
    """True is AT is compiled with OpenMP"""
    return isopenmp()


# noinspection PyUnusedLocal
def _mpi(self):
    """True is AT is compiled with MPI"""
    return ismpi()


_Dst.openmp = property(_omp)
_Dst.mpi = property(_mpi)




