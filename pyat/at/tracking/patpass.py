"""
Simple parallelisation of atpass() using multiprocessing.
"""
from warnings import warn
import multiprocessing
from at.tracking import atpass
from sys import platform
import numpy


__all__ = ['patpass']


globring = None


def _atpass_one(args):
    return atpass(globring, *args)


def _atpass(ring, r_in, nturns, refs, pool_size=None, globvar=True):
    if platform.startswith('linux') and globvar:
        global globring
        globring = ring
        args = [(r_in[:, i], nturns, refs) for i in range(r_in.shape[1])]
        with multiprocessing.Pool(pool_size) as pool:
            results = pool.map(_atpass_one, args)
        globring = None
    else:
        args = [(ring, r_in[:, i], nturns, refs) for i in range(r_in.shape[1])]
        with multiprocessing.Pool(pool_size) as pool:
            results = pool.starmap(atpass, args)
    return numpy.concatenate(results, axis=1)


def patpass(ring, r_in, nturns=1, refpts=None, pool_size=None, globvar=True):
    """
    Simple parallel implementation of atpass().  If more than one particle
    is supplied, use multiprocessing to run each particle in a separate
    process. In case a single particle is provided or the ring contains
    ImpedanceTablePass element, atpass is returned

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
        pool_size       number of processes, if None the min(npart,nproc) is used
        globvar         For linux machines speed-up is achieved by defining a global
                        ring variable, this can be disabled using globvar=False
    """
    if refpts is None:
        refpts = len(ring)
    refs = ring.uint32_refpts(refpts)
    pm_ok = [e.PassMethod=='ImpedanceTablePass' for e in ring]
    if len(numpy.atleast_1d(r_in[0]))>1 or any(pm_ok):
        if pool_size is None:
            pool_size = min(len(r_in[0]),multiprocessing.cpu_count())
        return _atpass(ring, r_in, nturns, refs, pool_size=pool_size, globvar=globvar)
    else:
        return atpass(ring, r_in, nturns, refs)
