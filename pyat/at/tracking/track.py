import numpy
import functools
from warnings import warn
# noinspection PyUnresolvedReferences
from .atpass import atpass as _atpass, elempass as _elempass
from ..lattice import Element, Particle, Refpts, uint32_refpts
from ..lattice import elements, get_elements
from typing import List, Iterable


__all__ = ['fortran_align', 'lattice_pass', 'element_pass', 'atpass',
           'elempass']

DIMENSION_ERROR = 'Input to lattice_pass() must be a 6xN array.'


def _set_beam_monitors(ring: List[Element], nbunch: int, nturns: int):
    monitors = get_elements(ring, elements.BeamMoments)
    for m in monitors:
        m.set_buffers(nturns, nbunch)
    return len(monitors) == 0


def fortran_align(func):
    # noinspection PyShadowingNames
    """decorator to ensure that *r_in* is Fortran-aligned

    :py:func:`fortran_align` ensures that the 2nd argument (usually *r_in*) of
    the decorated function is Fortran-aligned before calling the function

    Example:

        >>> @fortran_align
        ... def element_pass(element: Element, r_in, **kwargs):
        ... ...

        Ensure than *r_in* is Fortran-aligned
    """
    @functools.wraps(func)
    def wrapper(lattice, r_in, *args, **kwargs):
        assert r_in.shape[0] == 6 and r_in.ndim in (1, 2), DIMENSION_ERROR
        if r_in.flags.f_contiguous:
            return func(lattice, r_in, *args, **kwargs)
        else:
            r_fin = numpy.asfortranarray(r_in)
            r_out = func(lattice, r_fin, *args, **kwargs)
            r_in[:] = r_fin[:]
            return r_out

    return wrapper


@fortran_align
def lattice_pass(lattice: Iterable[Element], r_in, nturns: int = 1,
                 refpts: Refpts = None, **kwargs):
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
        refpts:                 numpy array of indices of elements where
          output is desired:

          * len(line) means end of the last element (default)
          * 0 means entrance of the first element

    Keyword arguments:
        keep_lattice (bool):    Use elements persisted from a previous
          call. If True, assume that the lattice has not changed since
          the previous call.
        keep_counter (bool):    Keep the turn number from the previous
          call.
        turn (int):             Starting turn number. Ignored if
          keep_counter is True. The turn number is necessary to compute the
          absolute path length used in RFCavityPass.
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

       * ``lattice_pass(lattice, r_in, refpts=len(line))`` is the same as
         ``lattice_pass(lattice, r_in)`` since the reference point len(line) is
         the exit of the last element.
       * ``lattice_pass(lattice, r_in, refpts=0)`` is a copy of *r_in* since
         the reference point 0 is the entrance of the first element.
       * To resume an interrupted tracking (for instance to get intermediate
         results), one must use one of the *turn* or *keep_counter*
         keywords to ensure the continuity of the turn number.
       * For multiparticle tracking with large number of turn the size of
         *r_out* may increase excessively. To avoid memory issues
         ``lattice_pass(lattice, r_in, refpts=[])`` can be used. An empty list
         is returned and the tracking results of the last turn are stored in
         *r_in*.

    """
    if not isinstance(lattice, list):
        lattice = list(lattice)
    if refpts is None:
        refpts = len(lattice)
    refs = uint32_refpts(refpts, len(lattice))
    # define properties if lattice is not a Lattice object
    nbunch = getattr(lattice, 'nbunch', 1)
    bunch_currents = getattr(lattice, 'bunch_currents', numpy.zeros(1))
    bunch_spos = getattr(lattice, 'bunch_spos', numpy.zeros(1))
    kwargs.update(bunch_currents=bunch_currents, bunch_spos=bunch_spos)
    no_bm = _set_beam_monitors(lattice, nbunch, nturns)
    kwargs['reuse'] = kwargs.pop('keep_lattice', False) and no_bm
    # atpass returns 6xNxRxT array
    # * N is number of particles;
    # * R is number of refpts
    # * T is the number of turns
    return _atpass(lattice, r_in, nturns, refpts=refs, **kwargs)


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
          Default: *lattice.particle* if existing,
          otherwise *Particle('relativistic')*
        energy (float):         lattice energy. Default 0.

    If *energy* is not available, relativistic tracking if forced,
    *rest_energy* is ignored.

    Returns:
        r_out:              (6, N) array containing output the coordinates of
          the particles at the exit of the element.
    """
    return _elempass(element, r_in, **kwargs)


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
