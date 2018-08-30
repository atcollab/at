"""
Simple parallelisation of atpass() using multiprocessing.
"""
import multiprocessing
# noinspection PyUnresolvedReferences
from .atpass import atpass
from .lattice import uint32_refpts
import numpy


def _atpass_one(args):
    ring, rin, turns, refpts = args
    return atpass(ring, rin, turns, refpts)


def _patpass(ring, rin, nturns, refpts, pool_size):
    pool = multiprocessing.Pool(pool_size)
    args = [(ring, rin[:, i], nturns, refpts) for i in range(rin.shape[1])]
    results = pool.map(_atpass_one, args)
    return numpy.concatenate(results, axis=1)


def patpass(ring, rin, nturns, refpts=None, reuse=True, pool_size=None):
    """
    Simple parallel implementation of atpass().  If more than one particle
    is supplied, use multiprocessing to run each particle in a separate
    process.
    """
    if not reuse:
        raise ValueError('patpass does not support altering lattices')
    if refpts is None:
        refpts = len(ring)
    refs = uint32_refpts(refpts, len(ring))
    if pool_size is None:
        pool_size = multiprocessing.cpu_count()
    return _patpass(ring, rin, nturns, refs, pool_size)
