# noinspection PyProtectedMember
from ..lattice.options import DConstant

# when atimplib was compiled with mpicc
# this is required

if(DConstant.mpi):
    from mpi4py import MPI
