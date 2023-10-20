from .atpass import atpass as _atpass, elempass as _elempass
from .utils import fortran_align
from .track import internal_lpass, internal_epass, internal_plpass
from .utils import initialize_lpass
from ..lattice import Lattice, Element, Particle, Refpts, End
from typing import Iterable, Optional
from warnings import warn


__all__ = ['element_pass', 'lattice_pass', 'patpass', 'atpass', 'elempass']


# noinspection PyIncorrectDocstring
def atpass(*args, **kwargs):
    """
    atpass(line, r_in, nturns, refpts=[], reuse=False, omp_num_threads=0)

    Track input particles *r_in* along line for *nturns* turns.

    Parameters:
        line (Sequence[Element]): list of elements
        r_in:                   6 x n_particles Fortran-ordered numpy array.
          On return, rin contains the final coordinates of the particles
        nturns (int):           number of turns to be tracked

    Keyword arguments:
        refpts (Uint32_refs):   numpy array of indices of elements where
          output is desired:

          * 0 means entrance of the first element
          * len(line) means end of the last element

        energy:                 nominal energy [eV]
        rest_energy:            rest_energy of the particle [eV]
        charge:                 particle charge [elementary charge]
        reuse (bool):           if True, use previously cached description
          of the lattice.
        omp_num_threads (int):  number of OpenMP threads
          (default 0: automatic)
        losses (bool):          if True, process losses

    Returns:
        r_out:  6 x n_particles x n_refpts x n_turns Fortran-ordered
          numpy array of particle coordinates

    :meta private:
    """
    warn(UserWarning("The public interface for tracking is 'lattice_pass'"))
    return _atpass(*args, **kwargs)


# noinspection PyIncorrectDocstring
def elempass(*args, **kwargs):
    """elempass(element, r_in)

    Track input particles *r_in* through a single element.

    Parameters:
        element (Element):  AT element
        rin:                6 x n_particles Fortran-ordered numpy array.
          On return, rin contains the final coordinates of the particles

    Keyword arguments:
        energy:             nominal energy [eV]
        rest_energy:        rest_energy of the particle [eV]
        charge:             particle charge [elementary charge]

    :meta private:
    """
    warn(UserWarning("The public interface for tracking is 'element_pass'"))
    return _elempass(*args, **kwargs)


def lattice_pass(lattice: Iterable[Element], r_in, nturns: int = 1,
                 refpts: Refpts = End, **kwargs):
    """
    :py:func:`lattice_pass` tracks particles through each element of a lattice
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

        particle (Optional[Particle]):  circulating particle.
          Default: :code:`lattice.particle` if existing,
          otherwise :code:`Particle('relativistic')`
        energy (Optiona[float]):        lattice energy. Default 0.
        unfold_beam (bool): Internal beam folding activate, this
            assumes the input particles are in bucket 0, works only
            if all bucket see the same RF Voltage.
            Default: :py:obj:`True`

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

       * :pycode:`lattice_pass(lattice, r_in, refpts=len(line))` is the same as
         :pycode:`lattice_pass(lattice, r_in)` since the reference point
         len(line) is the exit of the last element.
       * :pycode:`lattice_pass(lattice, r_in, refpts=0)` is a copy of *r_in*
         since the reference point 0 is the entrance of the first element.
       * To resume an interrupted tracking (for instance to get intermediate
         results), one must use one of the *turn* or *keep_counter*
         keywords to ensure the continuity of the turn number.
       * For multiparticle tracking with large number of turn the size of
         *r_out* may increase excessively. To avoid memory issues
         :pycode:`lattice_pass(lattice, r_in, refpts=None)` can be used.
         An empty list is returned and the tracking results of the last turn
         are stored in *r_in*.
       * To model buckets with different RF voltage :pycode:`unfold_beam=False`
         has to be used. The beam can be unfolded using the function
         :py:func:`.unfold_beam`. This function takes into account
         the true voltage in each bucket and distributes the particles in the
         bunches defined by :code:`ring.fillpattern` using a 6D orbit search.
    """
    lattice = initialize_lpass(lattice, nturns, kwargs)
    return internal_lpass(lattice, r_in, nturns=nturns, refpts=refpts,
                          no_varelem=False, **kwargs)


@fortran_align
def element_pass(element: Element, r_in, **kwargs):
    """Tracks particles through a single element.

    Parameters:
        element:                AT element
        r_in:                   (6, N) array: input coordinates of N particles.
          *r_in* is modified in-place and reports the coordinates at
          the end of the element. For the best efficiency, *r_in*
          should be given as F_CONTIGUOUS numpy array.

    Keyword arguments:
        particle (Particle):    circulating particle.
          Default: :code:`lattice.particle` if existing,
          otherwise :code:`Particle('relativistic')`
        energy (float):         lattice energy. Default 0.

    If *energy* is not available, relativistic tracking if forced,
    *rest_energy* is ignored.

    Returns:
        r_out:              (6, N) array containing output the coordinates of
          the particles at the exit of the element.
    """
    return internal_epass(element, r_in, **kwargs)


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
    kwargs['pool_size'] = pool_size
    kwargs['start_method'] = start_method
    lattice = initialize_lpass(lattice, nturns, kwargs)
    return internal_plpass(lattice, r_in, nturns=nturns,
                           refpts=refpts, **kwargs)


Lattice.lattice_pass = lattice_pass
