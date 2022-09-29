"""
Simple parallelisation of atpass() using multiprocessing.
"""
from functools import partial
import multiprocessing
from at.lattice import uint32_refpts
from warnings import warn
# noinspection PyUnresolvedReferences
from .atpass import atpass as _atpass
from at.lattice import AtWarning, DConstant
import numpy


__all__ = ['patpass']


globring = None


def format_results(results, r_in, losses):
    rin = [r['rin'] for r in results]
    r_in[:] = numpy.vstack(rin).T[:]
    if losses:
        rout = [r['results'][0] for r in results]
        rout = numpy.concatenate(rout, axis=1)
        lin = [r['results'][1] for r in results]
        lout = {}
        for k in lin[0].keys():
            lout[k] = numpy.hstack([li[k] for li in lin])
        return rout, lout
    else:
        rout = [r['results'] for r in results]
        rout = numpy.concatenate(rout, axis=1)
        return rout


def _atpass_one(ring, rin, **kwargs):
    if ring is None:
        result = _atpass(globring, rin, **kwargs)
    else:
        result = _atpass(ring, rin, **kwargs)
    return {'rin': rin, 'results': result}


def _pass(ring, r_in, pool_size, start_method, **kwargs):
    ctx = multiprocessing.get_context(start_method)
    if ctx.get_start_method() == 'fork':
        global globring
        globring = ring
        args = [(None, r_in[:, i]) for i in range(r_in.shape[1])]
        with ctx.Pool(pool_size) as pool:
            results = pool.starmap(partial(_atpass_one, **kwargs), args)
        globring = None
    else:
        args = [(ring, r_in[:, i]) for i in range(r_in.shape[1])]
        with ctx.Pool(pool_size) as pool:
            results = pool.starmap(partial(_atpass_one, **kwargs), args)
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
          may used for macOS to speed up the calculation or to solve
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
    bunch_currents = getattr(ring, 'bunch_currents', numpy.zeros(1))
    bunch_spos = getattr(ring, 'bunch_spos', numpy.zeros(1))
    kwargs.update(bunch_currents=bunch_currents, bunch_spos=bunch_spos)
    any_collective = collective(ring)
    if len(numpy.atleast_1d(r_in[0])) > 1 and not any_collective:
        if pool_size is None:
            pool_size = min(len(r_in[0]), multiprocessing.cpu_count(),
                            DConstant.patpass_poolsize)
        return _pass(ring, r_in, pool_size, start_method, nturns=nturns,
                     refpts=refpts, **kwargs)
    else:
        if any_collective:
            warn(AtWarning('Collective PassMethod found: use single process'))
        if r_in.flags.f_contiguous:
            return _atpass(ring, r_in, nturns=nturns,
                           refpts=refpts, **kwargs)
        else:
            r_fin = numpy.asfortranarray(r_in)
            r_out = _atpass(ring, r_fin, nturns=nturns,
                            refpts=refpts, **kwargs)
            r_in[:] = r_fin[:]
            return r_out
