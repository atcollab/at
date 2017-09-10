"""
Simple parallelisation of atpass() using multiprocessing.
"""
import multiprocessing
from at import atpass
import numpy


def _atpass_one(args):
    ring, rin, turns = args
    ans = atpass(ring, rin, turns)
    return ans


def _patpass(ring, rin, nturns, pool_size):
    rout = numpy.zeros((rin.shape[0], rin.shape[1] * nturns))
    pool = multiprocessing.Pool(pool_size)
    args = [(ring, rin[:, i], nturns) for i in range(rin.shape[1])]
    results = pool.map(_atpass_one, args)
    for i, res in enumerate(results):
        # Fold the results back into the same shape as atpass.
        rout[:, i::rin.shape[1]] = res
    return rout


def patpass(ring, rin, nturns, reuse=True, refpts=None, pool_size=None):
    """
    Simple parallel implementation of atpass().  If more than one particle
    is supplied, use multiprocessing to run each particle in a separate
    process.
    """
    if not reuse:
        raise ValueError('patpass does not support altering lattices')
    if refpts is not None:
        raise ValueError('patpass does not yet support refpts')
    if pool_size is None:
        pool_size = multiprocessing.cpu_count()
    return _patpass(ring, rin, nturns, pool_size)
