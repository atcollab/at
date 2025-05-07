"""Global set of constants."""

from __future__ import annotations

__all__ = ["test_mode", "set_test_mode", "DConstant", "random"]

import multiprocessing
import os
from contextlib import contextmanager
from contextvars import ContextVar

from numpy.random import Generator, PCG64, SeedSequence

from ..cconfig import ismpi, isopenmp, iscuda, isopencl

if ismpi():
    from mpi4py import MPI

    _comm = MPI.COMM_WORLD
    _MPI_sz = _comm.Get_size()
    _MPI_rk = _comm.Get_rank()
else:
    _MPI_sz = 1
    _MPI_rk = 0

_test_mode: ContextVar[int] = ContextVar("test_mode", default=0)


@contextmanager
def set_test_mode(mode: int = 1) -> None:
    r"""Context manager to set a \"test mode\" integer value.

    Args:
        mode:   a non-zero value which will be returned by the :py:func:`test_mode`
          function.

    Example:
        >>> assert not test_mode(), "test mode should not be set"
        >>> with set_test_mode():
        ...     assert test_mode(), "test mode should be set"

        Test mode contexts can be nested:

        >>> with set_test_mode(1):
        ...     assert test_mode() == 1, "outer context"
        ...     with set_test_mode(2):
        ...         assert test_mode() == 2, "inner context"
        ...     assert test_mode() == 1, "back to outer context"
    """
    tok = _test_mode.set(mode)
    try:
        yield None
    finally:
        _test_mode.reset(tok)


def test_mode() -> int:
    """Return the test mode value.

    Returns:
        test_mode:  0 outside a :py:func:`set_test_mode` context. Within such context,
          return its non-zero value
    """
    return _test_mode.get()


def _newgen(seed=12345):
    ss = SeedSequence(seed)
    seeds = ss.spawn(_MPI_sz + 1)
    return Generator(PCG64(seeds[0])), Generator(PCG64(seeds[_MPI_rk + 1]))


class _Dst:
    """Set of constants for AT numerical analysis."""

    # optimized values:
    # see https://github.com/atcollab/at/pull/199#issuecomment-745420443
    XYStep: float = 3.0e-8  # Coordinate step for differentiation
    DPStep: float = 3.0e-6  # Momentum step for dispersion and chromaticity
    OrbConvergence: float = 1.0e-12  # Convergence criterion for orbit
    OrbMaxIter: int = 20  # Max. number of iterations for orbit
    TStol: float = 1.0e-9  # Tolerance for synchronous phase search
    omp_num_threads: int = int(os.environ.get("OMP_NUM_THREADS", "0"))
    patpass_poolsize: int = multiprocessing.cpu_count()
    patpass_startmethod: str | None = None
    _rank: int = _MPI_rk  # MPI rank

    def __setattr__(self, name: str, value):
        _ = getattr(self, name)  # make sure attribute exists
        object.__setattr__(self, name, value)

    def reset(self, name: str) -> None:
        delattr(self, name)

    @property
    def mpi(self) -> bool:
        return ismpi()

    @property
    def openmp(self) -> bool:
        return isopenmp()

    @property
    def cuda(self) -> bool:
        return iscuda()

    @property
    def opencl(self) -> bool:
        return isopencl()

    @property
    def rank(self) -> int:
        return self._rank


class _Random:
    """Random generators for AT."""

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
    TStol:               Tolerance for synchronous phase search
    omp_num_threads:     Default number of OpenMP threads
    patpass_poolsize:    Default size of multiprocessing pool
    patpass_startmethod: Default start method for the multiprocessing
    mpi:                 :py:obj:`True` if MPI is active
    openmp:              :py:obj:`True` if OpenMP is active
    cuda:                :py:obj:`True` if CUDA is active
    opencl:              :py:obj:`True` if OpenCL is active

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
