"""
Simple parallelisation of atpass() using multiprocessing.
"""
from functools import partial
import multiprocessing
# noinspection PyProtectedMember
from ..lattice.utils import _uint32_refs
from ..lattice import AtWarning, Element, DConstant, random
from ..lattice import Refpts, End
from warnings import warn
from .atpass import reset_rng, atpass as _atpass
from .track import fortran_align
from typing import Iterable, Optional, List
import numpy as np

__all__ = ['patpass']

_imax = np.iinfo(int).max

_globring: Optional[List[Element]] = None


def format_results(results, r_in, losses):
    lin, lout = zip(*results)
    # Update r_in with values at the end of tracking
    np.concatenate(lin, out=r_in, axis=1)
    if losses:
        lout, ldic = zip(*lout)
        keys = ldic[0].keys()
        dicout = dict(((k, np.hstack([li[k] for li in ldic])) for k in keys))
        return np.concatenate(lout, axis=1), dicout
    else:
        return np.concatenate(lout, axis=1)


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
    args = enumerate(np.array_split(r_in, pool_size, axis=1))
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
def patpass(lattice: Iterable[Element], r_in, nturns: int = 1,
            refpts: Refpts = End, pool_size: int = None,
            start_method: str = None, **kwargs):
    """
    Simple parallel implementation of :py:func:`.lattice_pass`.
    If more than one particle is supplied, use multiprocessing. For a
    single particle or if the lattice contains :py:class:`.Collective`
    elements, :py:func:`.atpass` is used.

    :py:func:`patpass` tracks particles through each element of a lattice
    calling the element-specific tracking function specified in the Element's
    *PassMethod* field.

    Parameters:
        lattice:                list of elements
        r_in:                   (6, N) array: input coordinates of N particles.
          *r_in* is modified in-place and reports the coordinates at
          the end of the element. For the best efficiency, *r_in*
          should be given as F_CONTIGUOUS numpy array.
        nturns:                 number of turns to be tracked
        refpts:                 Selects the location of coordinates output.
          See ":ref:`Selecting elements in a lattice <refpts>`"
        pool_size:              number of processes. If None,
          ``min(npart,nproc)`` is used
        start_method:           python multiprocessing start method.
          :py:obj:`None` uses the python default that is considered safe.
          Available values: ``'fork'``, ``'spawn'``, ``'forkserver'``.
          Default for linux is ``'fork'``, default for macOS and  Windows is
          ``'spawn'``. ``'fork'`` may be used on macOS to speed up the
          calculation or to solve Runtime Errors, however it is considered
          unsafe.

    Keyword arguments:
        keep_lattice (bool):    Use elements persisted from a previous
          call. If :py:obj:`True`, assume that the lattice has not changed
          since the previous call.
        keep_counter (bool):    Keep the turn number from the previous
          call.
        turn (int):             Starting turn number. Ignored if
          *keep_counter* is :py:obj:`True`. The turn number is necessary to
          compute the absolute path length used in RFCavityPass.
        losses (bool):          Boolean to activate loss maps output
        omp_num_threads (int):  Number of OpenMP threads
          (default: automatic)

    The following keyword arguments overload the Lattice values

    Keyword arguments:
        particle (Particle):    circulating particle.
          Default: *lattice.particle* if existing,
          otherwise *Particle('relativistic')*
        energy (float):         lattice energy. Default 0.

    If *energy* is not available, relativistic tracking if forced,
    *rest_energy* is ignored.

    Returns:
        r_out: (6, N, R, T) array containing output coordinates of N particles
          at R reference points for T turns.
        loss_map: If *losses* is :py:obj:`True`: dictionary with the
          following key:

          ==============    ===================================================
          **islost**        (npart,) bool array indicating lost particles
          **turn**          (npart,) int array indicating the turn at
                            which the particle is lost
          **element**       ((npart,) int array indicating the element at
                            which the particle is lost
          **coord**         (6, npart) float array giving the coordinates at
                            which the particle is lost (zero for surviving
                            particles)
          ==============    ===================================================

    .. note::

       * For multiparticle tracking with large number of turn the size of
         *r_out* may increase excessively. To avoid memory issues
         :pycode:`lattice_pass(lattice, r_in, refpts=[])` can be used.
         An empty list is returned and the tracking results of the last turn
         are stored in *r_in*.
       * By default, :py:func:`patpass` will use all the available CPUs.
         To change the number of cores used in ALL functions using
         :py:func:`patpass` (:py:mod:`~at.acceptance.acceptance` module for
         example) it is possible to set ``at.DConstant.patpass_poolsize``
         to the desired value.

    """
    def collective(rg) -> bool:
        """True if any element involves collective effects"""
        for elem in rg:
            if elem.is_collective:
                return True
        return False

    if not isinstance(lattice, list):
        lattice = list(lattice)
    refpts = _uint32_refs(lattice, refpts)
    bunch_currents = getattr(lattice, 'bunch_currents', np.zeros(1))
    bunch_spos = getattr(lattice, 'bunch_spos', np.zeros(1))
    kwargs.update(bunch_currents=bunch_currents, bunch_spos=bunch_spos)
    kwargs['reuse'] = kwargs.pop('keep_lattice', False)
    any_collective = collective(lattice)
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
