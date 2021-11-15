# noinspection PyProtectedMember
from ..lattice.options import DConstant

# when atimplib was compiled with mpicc
# this is required

if(DConstant.mpi):
    try:
        from mpi4py import MPI
    except ModuleNotFoundError:
        print('mpi4py could not be imported')
