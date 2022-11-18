import numpy
from warnings import warn
# noinspection PyUnresolvedReferences
from .atpass import atpass as _atpass, elempass as _elempass
from ..lattice import DConstant, uint32_refpts, get_elements
from ..lattice import BeamMonitor


__all__ = ['lattice_pass', 'element_pass', 'atpass', 'elempass']

DIMENSION_ERROR = 'Input to lattice_pass() must be a 6xN array.'


def _set_beam_monitors(ring, nturns):
    monitors = get_elements(ring, BeamMonitor)
    for m in monitors:
        m.set_buffers(nturns)
    return len(monitors)==0


# noinspection PyIncorrectDocstring
def lattice_pass(lattice, r_in, nturns=1, refpts=None, keep_lattice=False,
                 omp_num_threads=None, **kwargs):
    """
    lattice_pass(lattice, r_in, nturns=1, refpts=None, keep_lattice=False, keep_counter=False, turn=0, losses=False, omp_num_threads=None)

    lattice_pass tracks particles through each element of a lattice
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

          If losses is True: {islost,turn,elem,coord} dictionary containing
          flag for particles lost (True -> particle lost), turn, element and
          coordinates at which the particle is lost. Set to zero for particles
          that survived

    .. note::

       * ``lattice_pass(lattice, r_in, refpts=len(line))`` is the same as
         ``lattice_pass(lattice, r_in)`` since the reference point len(line) is
         the exit of the last element.
       * ``lattice_pass(lattice, r_in, refpts=0)`` is a copy of ``r_in`` since
         the reference point 0 is the entrance of the first element.
       * To resume an interrupted tracking (for instance to get intermediate
         results), one must use one of the ``turn`` or ``keep_counter``
         keywords to ensure the continuity of the turn number.
       * For multiparticle tracking with large number of turn the size of
         ``r_out`` may increase excessively. To avoid memory issues
         ``lattice_pass(lattice, r_in, refpts=[])`` can be used. An empty list
         is returned and the tracking results of the last turn are stored in
         ``r_in``.

    """
    assert r_in.shape[0] == 6 and r_in.ndim in (1, 2), DIMENSION_ERROR
    if not isinstance(lattice, list):
        lattice = list(lattice)
    if refpts is None:
        refpts = len(lattice)
    if omp_num_threads is None:
        omp_num_threads = DConstant.omp_num_threads
    refs = uint32_refpts(refpts, len(lattice))
    no_bm = _set_beam_monitors(lattice, nturns)
    keep_lattice = keep_lattice and no_bm
    bunch_currents = getattr(lattice, 'bunch_currents', numpy.zeros(1))
    bunch_spos = getattr(lattice, 'bunch_spos', numpy.zeros(1))
    kwargs.update({'bunch_currents': bunch_currents,
                   'bunch_spos': bunch_spos})
    # atpass returns 6xAxBxC array where n = x*y*z;
    # * A is number of particles;
    # * B is number of refpts
    # * C is the number of turns
    if r_in.flags.f_contiguous:
        return _atpass(lattice, r_in, nturns, refpts=refs,
                       reuse=keep_lattice,
                       omp_num_threads=omp_num_threads,
                       **kwargs)
    else:
        r_fin = numpy.asfortranarray(r_in)
        r_out = _atpass(lattice, r_fin, nturns, refpts=refs,
                        reuse=keep_lattice,
                        omp_num_threads=omp_num_threads,
                        **kwargs)
        r_in[:] = r_fin[:]
        return r_out


def element_pass(element, r_in, **kwargs):
    """Tracks particles through a single element.

    Parameters:
        element (Element):  AT element
        r_in:               (6, N) array: input coordinates of N particles.
          r_in is modified in-place and reports the coordinates at
          the end of the element. For the best efficiency, r_in
          should be given as F_CONTIGUOUS numpy array.

    Keyword arguments:
        particle (Optional[Particle]):  circulating particle.
          Default: Particle('relativistic')
        energy (Optiona[float]):        lattice energy

    If ``energy`` is not available, relativistic tracking if forced,
    ``rest_energy`` is ignored.

    Returns:
        r_out:              (6, N) array containing output the coordinates of
          the particles at the exit of the element.
    """
    r_in = numpy.asfortranarray(r_in)
    return _elempass(element, r_in, **kwargs)


# noinspection PyIncorrectDocstring
def atpass(*args, **kwargs):
    """
    atpass(line, r in, nturns, refpts=[], reuse=False, omp_num_threads=0)

    Track input particles ``rin`` along line for ``nturns`` turns.
    Record 6D phase space at elements corresponding to refpts for each turn.

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
        reuse (Optional[bool]): if True, use previously cached description
          of the lattice.
        omp_num_threads (Optional[int]): number of OpenMP threads
          (default 0: automatic)
        losses (Optional[bool]):if True, process losses

    Returns:
        r_out:  6 x n_particles x n_refpts x n_turns Fortran-ordered
          numpy array of particle coordinates

    :meta private:
    """
    warn(UserWarning("The public interface for tracking is 'lattice_pass'"))
    return _atpass(*args, **kwargs)


# noinspection PyIncorrectDocstring
def elempass(*args, **kwargs):
    """elempass(element, rin)

    Tracks input particles rin through a single element.

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
