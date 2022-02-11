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

if ismpi():
    # if AT is compiled with mpicc this is required
    # noinspection PyUnresolvedReferences
    from mpi4py import MPI

_Dst.openmp = property(lambda self: isopenmp())
_Dst.mpi = property(lambda self: ismpi())
