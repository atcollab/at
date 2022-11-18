"""Global set of constants"""
import os
import multiprocessing


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

    def __setattr__(self, name, value):
        _ = getattr(self, name)     # make sure attribute exists
        object.__setattr__(self, name, value)

    def reset(self, name):
        delattr(self, name)


DConstant = _Dst()
"""
Default values for AT algorithms

Attributes:
    XYStep:             Coordinate step for differentiation
    DPStep:             Momentum step for dispersion and chromaticity
    OrbConvergence:     Convergence criterion for orbit
    OrbMaxIter:         Max. number of iterations for orbit
    omp_num_threads:    Default number of OpenMP threads
"""
