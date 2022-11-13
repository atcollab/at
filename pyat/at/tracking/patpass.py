"""
Simple parallelisation of atpass() using multiprocessing.
"""
from functools import partial
import multiprocessing
from at.lattice import uint32_refpts
from warnings import warn
# noinspection PyUnresolvedReferences
from .atpass import atpass as _atpass
from .track import fortran_align
from at.lattice import AtWarning, DConstant
import numpy as np


__all__ = ['patpass']

_atpassf = fortran_align(_atpass)

globring = None


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


def _atpass_one(ring, rin, **kwargs):
    kwargs['id'] = multiprocessing.current_process().ident
    if ring is None:
        result = _atpass(globring, rin, **kwargs)
    else:
        result = _atpass(ring, rin, **kwargs)
    return rin, result


@fortran_align
def _pass(ring, r_in, pool_size, start_method, **kwargs):
    ctx = multiprocessing.get_context(start_method)
    args = np.array_split(r_in, pool_size, axis=1)
    if ctx.get_start_method() == 'fork':
        global globring
        globring = ring
        with ctx.Pool(pool_size) as pool:
            results = pool.map(partial(_atpass_one, None, **kwargs), args)
        globring = None
    else:
        with ctx.Pool(pool_size) as pool:
            results = pool.map(partial(_atpass_one, ring, **kwargs), args)
    losses = kwargs.pop('losses', False)
    return format_results(results, r_in, losses)


# noinspection PyIncorrectDocstring
def patpass(ring, r_in, nturns=1, refpts=None, pool_size=None,
            start_method=None, **kwargs):
    """
    patpass(lattice, r_in, nturns=1, refpts=None, keep_lattice=False,
    keep_counter=False, turn=0, losses=False, omp_num_threads=None,
    pool_size=None, start_method=None)

    Simple parallel implementation of atpass().  If more than one particle
    is supplied, use multiprocessing to run each particle in a separate
    process. In case a single particle is provided or the ring contains
    ImpedanceTablePass element, atpass() is returned

    ``patpass`` tracks particles through each element of a lattice
    calling the element-specific tracking function specified in the Element's
    ``PassMethod`` field.

    Parameters:
        lattice (Iterable[Element]): list of elements
        r_in:                   6 x n_particles Fortran-ordered numpy array.
          On return, rin contains the final coordinates of the particles
        nturns (int):           number of turns to be tracked
        refpts (Uint32_refs):   numpy array of indices of elements where
          output is desired:

          * 0 means entrance of the first element
          * len(line) means end of the last element

    Keyword arguments:
        keep_lattice (Optional[bool]):  use elements persisted from a previous
          call. If True, assume that the lattice has not changed since
          the previous call.
        keep_counter (Optional[bool]):  Keep the turn number from the previous
          call.
        turn (Optional[int]):           Starting turn number. Ignored if
          keep_counter is True. The turn number is necessary to compute the
          absolute path length used in RFCavityPass.
        losses (Optional[bool]):        Boolean to activate loss maps output
        pool_size (Optional[int]):      number of processes. If None,
          ``min(npart,nproc)`` is used
        start_method (Optional[str]):   This parameter allows to change the
          python multiprocessing start method, default=None uses the python
          defaults that is considered safe. Available parameters:
          '``fork'``, ``'spawn'``, ``'forkserver'``. Default for linux is
          ``'fork'``, default for macOS and  Windows is ``'spawn'``. ``'fork'``
          may be used for macOS to speed up the calculation or to solve
          Runtime Errors, however it is considered unsafe.
        omp_num_threads (Optional[int]): number of OpenMP threads
          (default: automatic)

    The following keyword arguments overload the Lattice values

    Keyword arguments:
        particle (Optional[Particle]):  circulating particle.
          Default: ``lattice.particle`` if existing,
          otherwise ``Particle('relativistic')``
        energy (Optiona[float]):        lattice energy. Default 0.

    If ``energy`` is not available, relativistic tracking if forced,
    ``rest_energy`` is ignored.

    Returns:
        r_out: (6, N, R, T) array containing output coordinates of N particles
          at R reference points for T turns.

          If losses ==True: {islost,turn,elem,coord} dictionary containing
          flag for particles lost (True -> particle lost), turn, element and
          coordinates at which the particle is lost. Set to zero for particles
          that survived

    .. note::

       * For multiparticle tracking with large number of turn the size of
         ``r_out`` may increase excessively. To avoid memory issues
         ``patpass(lattice, r_in, refpts=[])`` can be used. An empty list
         is returned and the tracking results of the last turn are stored in
         ``r_in``.
       * By default, ``patpass`` will use all the available CPUs, to change
         the number of cores used in ALL functions using ``patpass``
         (``acceptance`` module for example) it is possible to set
         ``at.DConstant.patpass_poolsize`` to the desired value

    """
    def collective(rg) -> bool:
        """True if any element involves collective effects"""
        for elem in rg:
            if elem.is_collective:
                return True
        return False

    if not isinstance(ring, list):
        ring = list(ring)
    if refpts is None:
        refpts = len(ring)
    refpts = uint32_refpts(refpts, len(ring))
    bunch_currents = getattr(ring, 'bunch_currents', np.zeros(1))
    bunch_spos = getattr(ring, 'bunch_spos', np.zeros(1))
    kwargs.update(bunch_currents=bunch_currents, bunch_spos=bunch_spos)
    any_collective = collective(ring)
    rshape = r_in.shape
    if len(rshape) >= 2 and rshape[1] > 1 and not any_collective:
        if pool_size is None:
            pool_size = min(len(r_in[0]), multiprocessing.cpu_count(),
                            DConstant.patpass_poolsize)
        return _pass(ring, r_in, pool_size, start_method, nturns=nturns,
                     refpts=refpts, **kwargs)
    else:
        if any_collective:
            warn(AtWarning('Collective PassMethod found: use single process'))
        else:
            warn(AtWarning('no parallel computation for a single particle'))
        return _atpassf(ring, r_in, nturns=nturns, refpts=refpts, **kwargs)
