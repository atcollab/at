import numpy
from warnings import warn
# noinspection PyUnresolvedReferences
from .atpass import atpass as _atpass, elempass as _elempass
from ..lattice import Particle, DConstant, uint32_refpts


__all__ = ['lattice_pass', 'element_pass', 'atpass', 'elempass']

DIMENSION_ERROR = 'Input to lattice_pass() must be a 6xN array.'


def lattice_pass(lattice, r_in, nturns=1, refpts=None, keep_lattice=False,
                 omp_num_threads=None, **kwargs):
    """lattice_pass tracks particles through each element of a lattice
    calling the element-specific tracking function specified in the
    lattice[i].PassMethod field.

    Note:

     * lattice_pass(lattice, r_in, refpts=len(line)) is the same as
       lattice_pass(lattice, r_in) since the reference point len(line) is the
       exit of the last element
     * lattice_pass(lattice, r_in, refpts=0) is a copy of r_in since the
       reference point 0 is the entrance of the first element

    PARAMETERS
        lattice:    iterable of AT elements
        r_in:       (6, N) array: input coordinates of N particles.
                    r_in is modified in-place and reports the coordinates at
                    the end of the tracking. For the the best efficiency, r_in
                    should be given as F_CONTIGUOUS numpy array.
    KEYWORDS
        nturns=1:   number of passes through the lattice line
        refpts      elements at which data is returned. It can be:
                    1) an integer in the range [-len(ring), len(ring)-1]
                       selecting the element according to python indexing
                       rules. As a special case, len(ring) is allowed and
                       refers to the end of the last element,
                    2) an ordered list of such integers without duplicates,
                    3) a numpy array of booleans of maximum length
                       len(ring)+1, where selected elements are True.
                    Default: end of lattice
        keep_lattice: use elements persisted from a previous call to at.atpass.
                    If True, assume that the lattice has not changed since
                    that previous call.
        losses:     Boolean to activate loss maps output, default is False
    The following keyword overloads a value from lattice:
        particle:   circulating particle. Default: lattice.particle if
                    existing, otherwise Particle('relativistic')
    The following keywords overload values from lattice of from particle
    keyword
        energy      lattice energy
        rest_energy rest energy of the circulating particle [eV]
        charge      charge of the circulating particle [elementary charge]

    If 'energy' is not available, relativistic tracking if forced, rest_energy
    is ignored.

    OUTPUT
        (6, N, R, T) array containing output coordinates of N particles
        at R reference points for T turns.
        If losses ==True: {islost,turn,elem,coord} dictionary containing
        flag for particles lost (True -> particle lost), turn, element and
        coordinates at which the particle is lost. Set to zero for particles
        that survived
    """
    assert r_in.shape[0] == 6 and r_in.ndim in (1, 2), DIMENSION_ERROR
    if not isinstance(lattice, list):
        lattice = list(lattice)
    if refpts is None:
        refpts = len(lattice)
    if omp_num_threads is None:
        omp_num_threads = DConstant.omp_num_threads
    refs = uint32_refpts(refpts, len(lattice))
    particle = kwargs.pop('particle', getattr(lattice, 'particle', Particle()))
    try:
        # try to get 'energy' from the lattice
        kwargs.setdefault('energy', getattr(lattice, 'energy'))
    except AttributeError:
        pass
    if 'energy' in kwargs:
        # energy available, use the particle properties
        kwargs.setdefault('rest_energy', particle.rest_energy)
        kwargs.setdefault('charge', particle.charge)
    else:
        # energy no available, force relativistic tracking
        kwargs['rest_energy'] = 0.0
        kwargs['charge'] = -1.0
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
    """element_pass tracks particles through a single element.

    PARAMETERS
        element:    AT element
        r_in:       (6, N) array: input coordinates of N particles.
                    r_in is modified in-place and reports the coordinates at
                    the end of the tracking. For the the best efficiency, r_in
                    should be given as F_CONTIGUOUS numpy array.
    KEYWORDS
        particle:   circulating particle. Default: Particle('relativistic')
        energy      lattice energy
    The following keywords overload the values from the particle keyword:
        rest_energy rest energy of the circulating particle [eV]
        charge      charge of the circulating particle [elementary charge]

    If 'energy' is not available, relativistic tracking if forced, rest_energy
    is ignored.

    OUTPUT
        (6, N) array containing output the coordinates of the particles at the
        exit of the element.
    """
    particle = kwargs.pop('particle', Particle())
    if 'energy' in kwargs:
        # energy available: use the particle properties
        kwargs.setdefault('rest_energy', particle.rest_energy)
        kwargs.setdefault('charge', particle.charge)
    else:
        # energy not available: force relativistic tracking
        kwargs['rest_energy'] = 0.0
        kwargs['charge'] = -1.0
    r_in = numpy.asfortranarray(r_in)
    return _elempass(element, r_in, **kwargs)


def atpass(*args, **kwargs):
    warn(UserWarning("The public interface for tracking is 'lattice_pass'"))
    return _atpass(*args, **kwargs)


def elempass(*args, **kwargs):
    warn(UserWarning("The public interface for tracking is 'element_pass'"))
    return _elempass(*args, **kwargs)
