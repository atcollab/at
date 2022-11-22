"""
Tracking functions
"""
# noinspection PyUnresolvedReferences
from ..cconfig import ismpi
# noinspection PyUnresolvedReferences
from .atpass import reset_rng, common_rng, thread_rng
from .patpass import patpass
from .track import *
from .particles import *

# initialise the C random generators
if ismpi():
    from mpi4py import MPI
    _comm = MPI.COMM_WORLD
    _rank = _comm.Get_rank()
else:
    _rank = 0
reset_rng(_rank)
