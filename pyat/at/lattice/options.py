"""Global set of constants"""
from ..cconfig import ismpi, isopenmp
from numpy.random import Generator, PCG64, SeedSequence
import os
import multiprocessing

if ismpi():
    from mpi4py import MPI
    _comm = MPI.COMM_WORLD
    _MPI_sz = _comm.Get_size()
    _MPI_rk = _comm.Get_rank()
else:
    _MPI_sz = 1
    _MPI_rk = 0


def _newgen(seed=12345):
    ss = SeedSequence(seed)
    seeds = ss.spawn(_MPI_sz + 1)
    return Generator(PCG64(seeds[0])), Generator(PCG64(seeds[_MPI_rk + 1]))


class _Dst(object):
    """Set of constants for AT numerical analysis"""
    # optimized values:
    # see https://github.com/atcollab/at/pull/199#issuecomment-745420443
    XYStep = 3.e-8           # Coordinate step for differentiation
    DPStep = 3.e-6           # Momentum step for dispersion and chromaticity
    OrbConvergence = 1.e-12  # Convergence criterion for orbit
    OrbMaxIter = 20          # Max. number of iterations for orbit
    omp_num_threads = int(os.environ.get('OMP_NUM_THREADS', '0'))
    patpass_poolsize = multiprocessing.cpu_count()
    patpass_startmethod = None
    _rank = _MPI_rk         # MPI rank

    def __setattr__(self, name, value):
        _ = getattr(self, name)     # make sure attribute exists
        object.__setattr__(self, name, value)

    def reset(self, name):
        delattr(self, name)

    @property
    def mpi(self):
        return ismpi()

    @property
    def openmp(self):
        return isopenmp()

    @property
    def rank(self):
        return self._rank


class _Random(object):
    """Random generators for AT"""
    def __init__(self):
        self.reset()

    @classmethod
    def reset(cls, seed=12345):
        cls.common, cls.thread = _newgen(seed)


DConstant = _Dst()
"""
Default values for AT algorithms

Attributes:
    XYStep:              Coordinate step for differentiation
    DPStep:              Momentum step for dispersion and chromaticity
    OrbConvergence:      Convergence criterion for orbit
    OrbMaxIter:          Max. number of iterations for orbit
    omp_num_threads:     Default number of OpenMP threads
    patpass_poolsize:    Default size of multiprocessing pool
    patpass_startmethod: Default start method for the multiprocessing
    mpi:                 :py:obj:`True` if MPI is active
    openmp:              :py:obj:`True` if OpenMP is active

Methods:
    reset(attrname):    Reset the attribute to its default value
"""

random = _Random()
"""
Random generators for AT

Attributes:
    common(:py:class:`~numpy.random.Generator`): random generator common to
        all threads
    thread(:py:class:`~numpy.random.Generator`): random generator specific to
        each thread

Methods:
    reset(seed=None):   Reset the generators.

Examples:
    >>> at.random.common.random()
    0.8699988509120198

    Generates a random number in [0, 1) which is identical in all threads

    >>> at.random.thread.normal(size=5)
    array([-0.02326275, -0.89806482, -1.69086929, -0.74399398,  0.70156743])

    Generates an array of Gaussian random values which is specific to
    each thread
"""
