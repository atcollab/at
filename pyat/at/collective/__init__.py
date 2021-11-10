# noinspection PyProtectedMember
from ..lattice.options import _Dst

# when atimplib was compiled with mpicc
# this is required

if(_Dst.mpi):
    try:
        from mpi4py import MPI
    except ModuleNotFoundError:
        print('mpi4py could not be imported')
