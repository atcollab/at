"""
Simple parallelisation of atpass() using multiprocessing.
"""
import multiprocessing
from at.tracking import atpass
from at.lattice import uint32_refpts
import numpy


__all__ = ['patpass']


def _atpass_one(args):
    ring, rin, turns, refpts = args
    return atpass(ring, rin, turns, refpts)


def _patpass(ring, r_in, nturns, refpts, pool_size):
    pool = multiprocessing.Pool(pool_size)
    args = [(ring, r_in[:, i], nturns, refpts) for i in range(r_in.shape[1])]
    results = pool.map(_atpass_one, args)
    return numpy.concatenate(results, axis=1)


def patpass(ring, r_in, nturns, refpts=None, reuse=True, pool_size=None):
    """
    Simple parallel implementation of atpass().  If more than one particle
    is supplied, use multiprocessing to run each particle in a separate
    process.

    INPUT:
        ring            lattice description
        r_in:           6xN array: input coordinates of N particles
        nturns:         number of passes through the lattice line
        refpts          elements at which data is returned. It can be:
                        1) an integer in the range [-len(ring), len(ring)-1]
                           selecting the element according to python indexing
                           rules. As a special case, len(ring) is allowed and
                           refers to the end of the last element,
                        2) an ordered list of such integers without duplicates,
                        3) a numpy array of booleans of maximum length
                           len(ring)+1, where selected elements are True.
                        Defaults to None, meaning no refpts, equivelent to
                        passing an empty array for calculation purposes.
    """
    if not reuse:
        raise ValueError('patpass does not support altering lattices')
    if refpts is None:
        refpts = len(ring)
    refs = uint32_refpts(refpts, len(ring))
    if pool_size is None:
        pool_size = multiprocessing.cpu_count()
    return _patpass(ring, r_in, nturns, refs, pool_size)
