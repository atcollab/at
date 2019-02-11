"""
Closed orbit related functions
"""

import numpy
import scipy.constants as constants
from ..lattice import AtWarning, AtError, get_s_pos, RFCavity, uint32_refpts
from ..tracking import lattice_pass
import warnings

__all__ = ['find_orbit4', 'find_sync_orbit', 'find_orbit6']

STEP_SIZE = 1e-6
MAX_ITERATIONS = 20
CONVERGENCE = 1e-12


def find_orbit4(ring, dp=0.0, refpts=None, guess=None, **kwargs):
    """findorbit4 finds the closed orbit in the 4-d transverse phase
    space by numerically solving for a fixed point of the one turn
    map M calculated with lattice_pass.

        (X, PX, Y, PY, dP2, CT2 ) = M (X, PX, Y, PY, dP1, CT1)

    under the CONSTANT MOMENTUM constraint, dP2 = dP1 = dP and
    there is NO constraint on the 6-th coordinate CT

    IMPORTANT!!! findorbit4 imposes a constraint on dP and relaxes
    the constraint on the revolution frequency. A physical storage
    ring does exactly the opposite: the momentum deviation of a
    particle on the closed orbit settles at the value
    such that the revolution is synchronous with the RF cavity

                HarmNumber*Frev = Frf

    To impose this artificial constraint in find_orbit4, PassMethod
    used for any element SHOULD NOT
    1.  change the longitudinal momentum dP (cavities , magnets with radiation)
    2.  have any time dependence (localized impedance, fast kickers etc)

    cod = find_orbit4(ring, dP)
        cod:    6x1 vector - fixed point at the entrance of the 1-st element of the ring (x,px,y,py)

    cod ,orbit = find_orbit4(ring, dP, refpts)
        cod:    6x1 vector - fixed point at the entrance of the 1-st element of the ring (x,px,y,py)
        orbit:  6xNrefs array - orbit at each location specified in refpts

    ... = find_orbit4(ring,...,guess=initial_orbit)
        sets the initial search to initial_orbit

    PARAMETERS
        ring            lattice description
        dp              momentum deviation. Defaults to 0
        refpts          elements at which data is returned. It can be:
                        1) an integer in the range [-len(ring), len(ring)-1]
                           selecting the element according to python indexing
                           rules. As a special case, len(ring) is allowed and
                           refers to the end of the last element,
                        2) an ordered list of such integers without duplicates,
                        3) a numpy array of booleans of maximum length
                           len(ring)+1, where selected elements are True.
                        Defaults to None, if refpts is None an empty array is
                        returned for orbit.

    See also find_sync_orbit, find_orbit6.
    """
    # We seek
    #  - f(x) = x
    #  - g(x) = f(x) - x = 0
    #  - g'(x) = f'(x) - 1
    # Use a Newton-Raphson-type algorithm:
    #  - r_n+1 = r_n - g(r_n) / g'(r_n)
    #  - r_n+1 = r_n - (f(r_n) - r_n) / (f'(r_n) - 1)
    #
    # (f(r_n) - r_n) / (f'(r_n) - 1) can be seen as x = b/a where we use least
    #     squares fitting to determine x when ax = b
    # f(r_n) - r_n is denoted b
    # f'(r_n) is the 4x4 jacobian, denoted j4

    convergence = kwargs.pop('convergence', CONVERGENCE)
    max_iterations = kwargs.pop('max_iterations', MAX_ITERATIONS)
    step_size = kwargs.pop('step_size', STEP_SIZE)
    if guess is None:
        ref_in = numpy.zeros((6,), order='F')
        ref_in[4] = dp
    else:
        ref_in = guess

    delta_matrix = numpy.zeros((6, 5), order='F')
    for i in range(4):
        delta_matrix[i, i] = step_size
    id4 = numpy.asfortranarray(numpy.identity(4))
    change = 1
    itercount = 0
    keeplattice = False
    while (change > convergence) and itercount < max_iterations:
        in_mat = ref_in.reshape((6, 1)) + delta_matrix
        _ = lattice_pass(ring, in_mat, refpts=[], keep_lattice=keeplattice)
        # the reference particle after one turn
        ref_out = in_mat[:, 4]
        # 4x4 jacobian matrix from numerical differentiation:
        # f(x+d) - f(x) / d
        j4 = (in_mat[:4, :4] - in_mat[:4, 4:]) / step_size
        a = j4 - id4  # f'(r_n) - 1
        b = ref_out[:4] - ref_in[:4]
        b_over_a, _, _, _ = numpy.linalg.lstsq(a, b, rcond=-1)
        r_next = ref_in - numpy.append(b_over_a, numpy.zeros((2,)))
        # determine if we are close enough
        change = numpy.linalg.norm(r_next - ref_in)
        itercount += 1
        ref_in = r_next
        keeplattice = True

    if itercount == max_iterations:
        warnings.warn(AtWarning('Maximum number of iterations reached. Possible non-convergence'))

    uint32refs = uint32_refpts(refpts, len(ring))
    all_points = numpy.empty((0, 6), dtype=float) if (len(uint32refs) == 0) else numpy.squeeze(
        lattice_pass(ring, ref_in.copy(order='K'), refpts=uint32refs, keep_lattice=keeplattice), axis=(1, 3)).T

    return ref_in, all_points


def find_sync_orbit(ring, dct=0.0, refpts=None, guess=None, **kwargs):
    """find_sync_orbit finds the closed orbit, synchronous with the RF cavity
    and momentum deviation dP (first 5 components of the phase space vector)
    % by numerically solving  for a fixed point
    % of the one turn map M calculated with lattice_pass

        (X, PX, Y, PY, dP2, CT2 ) = M (X, PX, Y, PY, dP1, CT1)

    under constraints CT2 - CT1 =  dCT = C(1/Frev - 1/Frev0) and dP2 = dP1 , where
    Frev0 = Frf0/HarmNumber is the design revolution frequency
    Frev  = (Frf0 + dFrf)/HarmNumber is the imposed revolution frequency

    IMPORTANT!!!  find_sync_orbit imposes a constraint (CT2 - CT1) and
    dP2 = dP1 but no constraint on the value of dP1, dP2
    The algorithm assumes time-independent fixed-momentum ring
    to reduce the dimensionality of the problem.

    To impose this artificial constraint in find_sync_orbit
    PassMethod used for any element SHOULD NOT
    1.  change the longitudinal momentum dP (cavities , magnets with radiation)
    2.  have any time dependence (localized impedance, fast kickers etc).

    cod = find_sync_orbit(ring, dct)
        cod:    6x1 vector - fixed point at the entrance of the 1-st element of the ring (x,px,y,py,delta)

    cod ,orbit = find_sync_orbit(ring, dct, refpts)
        cod:    6x1 vector - fixed point at the entrance of the 1-st element of the ring (x,px,y,py,delta)
        orbit:  6xNrefs array - orbit at each location specified in refpts

    ... = find_sync_orbit(ring,...,guess=initial_orbit)
        sets the initial search to initial_orbit

    PARAMETERS
        ring            lattice description
        dct              ? Defaults to 0
        refpts          elements at which data is returned. It can be:
                        1) an integer in the range [-len(ring), len(ring)-1]
                           selecting the element according to python indexing
                           rules. As a special case, len(ring) is allowed and
                           refers to the end of the last element,
                        2) an ordered list of such integers without duplicates,
                        3) a numpy array of booleans of maximum length
                           len(ring)+1, where selected elements are True.
                        Defaults to None, if refpts is None an empty array is
                        returned for orbit.

    See also find_orbit4, find_orbit6.
    """
    convergence = kwargs.pop('convergence', CONVERGENCE)
    max_iterations = kwargs.pop('max_iterations', MAX_ITERATIONS)
    step_size = kwargs.pop('step_size', STEP_SIZE)
    ref_in = numpy.zeros((6,), order='F') if guess is None else guess

    delta_matrix = numpy.zeros((6, 6), order='F')
    for i in range(5):
        delta_matrix[i, i] = step_size
    theta5 = numpy.zeros((5,), order='F')
    theta5[4] = dct
    id5 = numpy.zeros((5, 5), order='F')
    for i in range(4):
        id5[i, i] = 1.0
    idx = numpy.array([0, 1, 2, 3, 5])
    change = 1
    itercount = 0
    keeplattice = False
    while (change > convergence) and itercount < max_iterations:
        in_mat = ref_in.reshape((6, 1)) + delta_matrix
        _ = lattice_pass(ring, in_mat, refpts=[], keep_lattice=keeplattice)
        # the reference particle after one turn
        ref_out = in_mat[:, -1]
        # 5x5 jacobian matrix from numerical differentiation:
        # f(x+d) - f(x) / d
        j5 = (in_mat[idx, :5] - in_mat[idx, 5:]) / step_size
        a = j5 - id5  # f'(r_n) - 1
        b = ref_out[idx] - numpy.append(ref_in[:4], 0.0) - theta5
        b_over_a, _, _, _ = numpy.linalg.lstsq(a, b, rcond=-1)
        r_next = ref_in - numpy.append(b_over_a, 0.0)
        # determine if we are close enough
        change = numpy.linalg.norm(r_next - ref_in)
        itercount += 1
        ref_in = r_next
        keeplattice = True

    if itercount == max_iterations:
        warnings.warn(AtWarning('Maximum number of iterations reached. Possible non-convergence'))

    uint32refs = uint32_refpts(refpts, len(ring))
    all_points = numpy.empty((0, 6), dtype=float) if (len(uint32refs) == 0) else numpy.squeeze(
        lattice_pass(ring, ref_in.copy(order='K'), refpts=uint32refs, keep_lattice=keeplattice), axis=(1, 3)).T

    return ref_in, all_points


def find_orbit6(ring, refpts=None, guess=None, **kwargs):
    """find_orbit6 finds the closed orbit in the full 6-D phase space
    by numerically solving  for a fixed point of the one turn
    map M calculated with lattice_pass

    (X, PX, Y, PY, DP, CT2 ) = M (X, PX, Y, PY, DP, CT1)

    with constraint  CT2 - CT1 = C*HarmNumber(1/Frf - 1/Frf0)

    IMPORTANT!!! find_orbit6 is a realistic simulation
    1.  The Frf frequency in the RF cavities (may be different from Frf0)
        imposes the synchronous condition
        CT2 - CT1 = C*HarmNumber(1/Frf - 1/Frf0)
    2.  The algorithm numerically calculates
        6-by-6 Jacobian matrix J6. In order for (J-E) matrix
        to be non-singular it is NECESSARY to use a realistic
        PassMethod for cavities with non-zero momentum kick
        (such as ThinCavityPass).
    3.  find_orbit6 can find orbits with radiation.
        In order for the solution to exist the cavity must supply
        adequate energy compensation.
        In the simplest case of a single cavity, it must have
        'Voltage' field set so that Voltage > Erad - energy loss per turn
    4.  find_orbit6 starts the search from [0,,0, 0, 0, 0, 0], unless
        the third argument is specified: find_orbit6(ring,...,guess=initial_orbit)
        There exist a family of solutions that correspond to different RF buckets
        They differ in the 6-th coordinate by C*Nb/Frf. Nb = 1 .. HarmNum-1
    5.  The value of the 6-th coordinate found at the cavity gives
        the equilibrium RF phase. If there is no radiation the phase is 0;

    cod = find_orbit6(ring)
        cod:    6x1 vector - fixed point at the entrance of the 1-st element of the RING (x,px,y,py,delta,ct)

    cod, orbit = find_orbit6(ring, refpts)
        cod:    6x1 vector - fixed point at the entrance of the 1-st element of the RING (x,px,y,py,delta)
        orbit:  6xNrefs array - orbit at each location specified in refpts

    ... = find_orbit6(ring,...,guess=initial_orbit)
        sets the initial search to initial_orbit

    PARAMETERS
        ring            lattice description
        refpts          elements at which data is returned. It can be:
                        1) an integer in the range [-len(ring), len(ring)-1]
                           selecting the element according to python indexing
                           rules. As a special case, len(ring) is allowed and
                           refers to the end of the last element,
                        2) an ordered list of such integers without duplicates,
                        3) a numpy array of booleans of maximum length
                           len(ring)+1, where selected elements are True.
                        Defaults to None, if refpts is None an empty array is
                        returned for orbit.

    See also find_orbit4, find_sync_orbit.
    """
    convergence = kwargs.pop('convergence', CONVERGENCE)
    max_iterations = kwargs.pop('max_iterations', MAX_ITERATIONS)
    step_size = kwargs.pop('step_size', STEP_SIZE)
    ref_in = numpy.zeros((6,), order='F') if guess is None else guess

    # Get evolution period
    l0 = get_s_pos(ring, len(ring))
    cavities = [elem for elem in ring if isinstance(elem, RFCavity)]
    if len(cavities) == 0:
        raise AtError('No cavity found in the lattice.')

    f_rf = cavities[0].Frequency
    harm_number = cavities[0].HarmNumber

    theta = numpy.zeros((6,))
    theta[5] = constants.speed_of_light * harm_number / f_rf - l0

    delta_matrix = numpy.zeros((6, 7), order='F')
    for i in range(6):
        delta_matrix[i, i] = step_size

    id6 = numpy.asfortranarray(numpy.identity(6))
    change = 1
    itercount = 0
    keeplattice = False
    while (change > convergence) and itercount < max_iterations:
        in_mat = ref_in.reshape((6, 1)) + delta_matrix
        _ = lattice_pass(ring, in_mat, refpts=[], keep_lattice=keeplattice)
        # the reference particle after one turn
        ref_out = in_mat[:, 6]
        # 6x6 jacobian matrix from numerical differentiation:
        # f(x+d) - f(x) / d
        j6 = (in_mat[:, :6] - in_mat[:, 6:]) / step_size
        a = j6 - id6  # f'(r_n) - 1
        b = ref_out[:] - ref_in[:] - theta
        b_over_a, _, _, _ = numpy.linalg.lstsq(a, b, rcond=-1)
        r_next = ref_in - b_over_a
        # determine if we are close enough
        change = numpy.linalg.norm(r_next - ref_in)
        itercount += 1
        ref_in = r_next
        keeplattice = True

    if itercount == max_iterations:
        warnings.warn(AtWarning('Maximum number of iterations reached. Possible non-convergence'))

    uint32refs = uint32_refpts(refpts, len(ring))
    all_points = numpy.empty((0, 6), dtype=float) if (len(uint32refs) == 0) else numpy.squeeze(
        lattice_pass(ring, ref_in.copy(order='K'), refpts=uint32refs, keep_lattice=keeplattice), axis=(1, 3)).T

    return ref_in, all_points
