"""
Simple parallelisation of atpass() using multiprocessing.
"""
from functools import partial
import multiprocessing
from at.tracking import atpass
from at.lattice import uint32_refpts
from sys import platform
from warnings import warn
from at.lattice import AtWarning, elements
import numpy


__all__ = ['patpass']


globring = None


def format_results(results, r_in, losses):
    rin = [r['rin'] for r in results]
    r_in = numpy.vstack(rin).T
    if losses:
        rout = [r['results'][0] for r in results]
        rout = numpy.concatenate(rout, axis=1)
        lin = [r['results'][1] for r in results]
        lout = {}
        for k in lin[0].keys():
            lout[k] = numpy.hstack([l[k] for l in lin])
        return rout, lout
    else:
        rout = [r['results'] for r in results]
        rout = numpy.concatenate(rout, axis=1)
        return rout


def _atpass_one(ring, rin, **kwargs):
    if ring is None:
        result = atpass(globring, rin, **kwargs)
    else:
        result = atpass(ring, rin, **kwargs)
    return {'rin': rin, 'results': result}


def _atpass(ring, r_in, pool_size, globvar, **kwargs):
    if platform.startswith('linux') and globvar:
        global globring
        globring = ring
        args = [(None, r_in[:, i]) for i in range(r_in.shape[1])]
        with multiprocessing.Pool(pool_size) as pool:
            results = pool.starmap(partial(_atpass_one, **kwargs), args)
        globring = None
    else:
        args = [(ring, r_in[:, i]) for i in range(r_in.shape[1])]
        with multiprocessing.Pool(pool_size) as pool:
            results = pool.starmap(partial(_atpass_one, **kwargs), args)
    losses = kwargs.pop('losses', False)
    return format_results(results, r_in, losses)


def patpass(ring, r_in, nturns=1, refpts=None, losses=False, pool_size=None,
            globvar=True, **kwargs):
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
        losses          Activate loss maps
        pool_size       number of processes, if None the min(npart,nproc) is used
        globvar         For linux machines speed-up is achieved by defining a global
                        ring variable, this can be disabled using globvar=False

     OUTPUT:
        (6, N, R, T) array containing output coordinates of N particles
        at R reference points for T turns.
        If losses ==True: {islost,turn,elem,coord} dictionnary containing
        flag for particles lost (True -> particle lost), turn, element and
        coordinates at which the particle is lost. Set to zero for particles
        that survived
    """
    if refpts is None:
        refpts = len(ring)
    refpts = uint32_refpts(refpts, len(ring))
    pm_ok = [e.PassMethod in elements._collective for e in ring]
    if len(numpy.atleast_1d(r_in[0])) > 1 and not any(pm_ok):
        if pool_size is None:
            pool_size = min(len(r_in[0]), multiprocessing.cpu_count())
        return _atpass(ring, r_in, pool_size, globvar, nturns=nturns,
                       refpts=refpts, losses=losses)
    else:
        if any(pm_ok):
            warn(AtWarning('Collective PassMethod found: use single process'))
        return atpass(ring, r_in, nturns=nturns, refpts=refpts, losses=losses)
