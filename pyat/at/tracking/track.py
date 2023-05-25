import numpy
from .atpass import atpass as _atpass, elempass as _elempass
from .utils import fortran_align, has_collective, format_results
from .utils import initialize_args
from ..lattice import Lattice, Element, Particle, Refpts, End
from ..lattice import get_uint32_index
from ..lattice import AtWarning, DConstant, random
from typing import List, Iterable, Optional, Union
from functools import partial
import multiprocessing
from warnings import warn
from .atpass import reset_rng


__all__ = ['track_function']

_imax = numpy.iinfo(int).max
_globring: Optional[List[Element]] = None


def _atpass_fork(seed, rank, rin, **kwargs):
    """Single forked job"""
    reset_rng(rank, seed=seed)
    result = _atpass(_globring, rin, **kwargs)
    return rin, result


def _atpass_spawn(ring, seed, rank, rin, **kwargs):
    """Single spawned job"""
    reset_rng(rank, seed=seed)
    result = _atpass(ring, rin, **kwargs)
    return rin, result


def _pass(ring, r_in, pool_size, start_method, **kwargs):
    ctx = multiprocessing.get_context(start_method)
    # Split input in as many slices as processes
    args = enumerate(numpy.array_split(r_in, pool_size, axis=1))
    # Generate a new starting point for C RNGs
    seed = random.common.integers(0, high=_imax, dtype=int)
    global _globring
    _globring = ring
    if ctx.get_start_method() == 'fork':
        passfunc = partial(_atpass_fork, seed, **kwargs)
    else:
        passfunc = partial(_atpass_spawn, ring, seed, **kwargs)
    # Start the parallel jobs
    with ctx.Pool(pool_size) as pool:
        results = pool.starmap(passfunc, args)
    _globring = None
    # Gather the results
    losses = kwargs.pop('losses', False)
    return format_results(results, r_in, losses)

@fortran_align
def _element_pass(element: Element, r_in, **kwargs):
    return _elempass(element, r_in, **kwargs)


@fortran_align
def _lattice_pass(lattice: Iterable[Element], r_in, nturns: int = 1,
                  refpts: Refpts = End, **kwargs):
    if not isinstance(lattice, list):
        lattice = list(lattice)
    refs = get_uint32_index(lattice, refpts)
    kwargs['reuse'] = kwargs.pop('keep_lattice', False)
    return _atpass(lattice, r_in, nturns, refpts=refs, **kwargs)


@fortran_align
def _plattice_pass(lattice: Iterable[Element], r_in, nturns: int = 1,
                   refpts: Refpts = End, pool_size: int = None,
                   start_method: str = None, **kwargs):
    if not isinstance(lattice, list):
        lattice = list(lattice)
    refpts = get_uint32_index(lattice, refpts)
    any_collective = has_collective(lattice)
    rshape = r_in.shape
    if len(rshape) >= 2 and rshape[1] > 1 and not any_collective:
        if pool_size is None:
            pool_size = min(len(r_in[0]), multiprocessing.cpu_count(),
                            DConstant.patpass_poolsize)
        return _pass(lattice, r_in, pool_size, start_method, nturns=nturns,
                     refpts=refpts, **kwargs)
    else:
        if any_collective:
            warn(AtWarning('Collective PassMethod found: use single process'))
        else:
            warn(AtWarning('no parallel computation for a single particle'))
        return _atpass(lattice, r_in, nturns=nturns, refpts=refpts, **kwargs)


@initialize_args
def track_function(lattice: Union[Element, Iterable[Element]], r_in, nturns: int = 1,
                   refpts: Refpts = End, use_mp=False, **kwargs):

    pool_size = kwargs.pop('pool_size', None)
    start_method = kwargs.pop('start_method', None)
    losses = kwargs.pop('losses', False)

    trackdata = {}
    trackparam = {'rin': r_in.copy()}
    elem_kw = ['energy', 'particle']
    saved_kw = elem_kw + ['turn']
    try:
        npart = numpy.shape(r_in)[1]
    except IndexError:
        npart = 1
    [trackparam.update((kw, kwargs.get(kw))) for kw in kwargs if kw in saved_kw]
    trackparam.update({'refpts': lattice.get_uint32_index(refpts),
                       'nturns': nturns,
                       'npart': npart})

    ldtype = [('islost', numpy.bool_),
              ('turn', numpy.uint32),
              ('elem', numpy.uint32),
              ('coord', numpy.float64, (6, )),
              ]
    loss_map = numpy.recarray((npart,), ldtype)

    if isinstance(lattice, Element):
        kwargs = {k: v for k, v in kwargs.items() if k in elem_kw}
        rout = _element_pass(lattice, r_in, **kwargs)
    elif use_mp:
        rout = _plattice_pass(lattice, r_in, nturns=nturns, refpts=refpts,
                              pool_size=pool_size, start_method=start_method,
                              lossses=losses, **kwargs)
    else:
        rout = _lattice_pass(lattice, r_in, nturns=nturns, refpts=refpts,
                             losses=losses, **kwargs)
    if losses:
        rout, lm = rout
        for k, v in lm.items():
            loss_map[k] = v

    trackparam.update({'rout': r_in})
    trackdata.update({'loss_map': loss_map})
    return rout, trackparam, trackdata


Lattice.track = track_function
