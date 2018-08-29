import numpy
import scipy.constants as constants
import at
import math
import warnings

__all__ = ['find_orbit4', 'find_sync_orbit', 'find_orbit6', 'find_m44', 'find_m66', 'get_twiss']

XYDEFSTEP = 6.055454452393343e-006  # Optimal delta?
DPSTEP = 6.055454452393343e-006  # Optimal delta?
DDP = 1e-8
STEP_SIZE = 1e-6
MAX_ITERATIONS = 20
CONVERGENCE = 1e-12

# dtype for structured array containing Twiss parameters
TWISS_DTYPE = [('idx', numpy.uint32),
               ('s_pos', numpy.float64),
               ('closed_orbit', numpy.float64, (6,)),
               ('dispersion', numpy.float64, (4,)),
               ('alpha', numpy.float64, (2,)),
               ('beta', numpy.float64, (2,)),
               ('mu', numpy.float64, (2,)),
               ('m44', numpy.float64, (4, 4))]

# Prepare symplectic identity matrix
_s0 = numpy.zeros((2, 2), order='F')
_s2 = numpy.array([[0, 1], [-1, 0]], order='F')
_s4 = numpy.concatenate((numpy.concatenate((_s2, _s0), axis=1), numpy.concatenate((_s0, _s2), axis=1)), axis=0)
# prepare Identity matrix


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
        _ = at.lattice_pass(ring, in_mat, refpts=[], keep_lattice=keeplattice)
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
        warnings.warn(at.AtWarning('Maximum number of iterations reached. Possible non-convergence'))

    if refpts is None:
        output = ref_in
    else:
        # We know that there is one particle and one turn, so select the
        # (6, nrefs) output.
        all_points = at.lattice_pass(ring,
                                     ref_in.copy(order='K'),
                                     refpts=refpts,
                                     keep_lattice=keeplattice)[:,0,:,0]
        output = (ref_in, all_points)
    return output


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
        _ = at.lattice_pass(ring, in_mat, refpts=[], keep_lattice=keeplattice)
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
        warnings.warn(at.AtWarning('Maximum number of iterations reached. Possible non-convergence'))

    if refpts is None:
        output = ref_in
    else:
        all_points = numpy.squeeze(at.lattice_pass(ring, ref_in.copy(order='K'), refpts=refpts,
                                                   keep_lattice=keeplattice))
        output = (ref_in, all_points)
    return output


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

    See also find_orbit4, find_sync_orbit.
    """
    convergence = kwargs.pop('convergence', CONVERGENCE)
    max_iterations = kwargs.pop('max_iterations', MAX_ITERATIONS)
    step_size = kwargs.pop('step_size', STEP_SIZE)
    ref_in = numpy.zeros((6,), order='F') if guess is None else guess

    # Get evolution period
    l0 = at.get_s_pos(ring, len(ring))
    cavities = list(filter(at.checktype(at.RFCavity), ring))
    if len(cavities) == 0:
        raise at.AtError('No cavity found in the lattice.')

    fRF = cavities[0].Frequency
    harm_number = cavities[0].HarmNumber

    theta = numpy.zeros((6,))
    theta[5] = constants.speed_of_light * harm_number / fRF - l0

    delta_matrix = numpy.zeros((6, 7), order='F')
    for i in range(6):
        delta_matrix[i, i] = step_size

    id6 = numpy.asfortranarray(numpy.identity(6))
    change = 1
    itercount = 0
    keeplattice = False
    while (change > convergence) and itercount < max_iterations:
        in_mat = ref_in.reshape((6, 1)) + delta_matrix
        _ = at.lattice_pass(ring, in_mat, refpts=[], keep_lattice=keeplattice)
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
        warnings.warn(at.AtWarning('Maximum number of iterations reached. Possible non-convergence'))

    if refpts is None:
        output = ref_in
    else:
        all_points = numpy.squeeze(at.lattice_pass(ring, ref_in.copy(order='K'), refpts=refpts,
                                                   keep_lattice=keeplattice))
        output = (ref_in, all_points)
    return output


def find_m44(ring, dp=0.0, refpts=None, orbit=None, output_orbit=False, **kwargs):
    """find_m44 numerically finds the 4x4 transfer matrix of an accelerator lattice
    for a particle with relative momentum deviation DP

    IMPORTANT!!! find_m44 assumes constant momentum deviation.
    PassMethod used for any element in the lattice SHOULD NOT
    1.  change the longitudinal momentum dP
        (cavities , magnets with radiation, ...)
    2.  have any time dependence (localized impedance, fast kickers, ...)

    m44 = find_m44(lattice, dp=0.0)
        returns a full one-turn matrix at the entrance of the first element
        !!! With this syntax find_m44 assumes that the lattice
        is a ring and first finds the closed orbit

    m44, t = find_m44(lattice, dp, refpts)
        returns 4x4 transfer matrices between the entrance of the first element and each element indexed by refpts.
        t is a 4x4xNrefs array

    m44, t = find_m44(lattice, dp, refpts, full=True)
        Same as above except that matrixes returned in t are full 1-turn matrices at the entrance of each
        element indexed by refpts.

    ... = find_m44(lattice, ..., orbit=closed_orbit)
        Does not search for the closed orbit. Instead closed_orbit,a vector of initial conditions is used.
        This syntax is useful to specify the entrance orbit if lattice is not a ring or to avoid recomputing the
        closed orbit if it is already known.

    m44, t, orbit = find_m44(lattice, ..., output_orbit=True)
        Returns in addition the closed orbit at the entrance of the 1st element

    See also find_m66, find_orbit4
    """
    def mrotate(m):
        m = numpy.squeeze(m)
        return numpy.linalg.multi_dot([m, m44, _s4.T, m.T, _s4])

    XYStep = kwargs.pop('XYStep', XYDEFSTEP)
    full = kwargs.pop('full', False)
    keeplattice = False
    if orbit is None:
        orbit = find_orbit4(ring, dp)
        keeplattice = True
    # Construct matrix of plus and minus deltas
    dg = numpy.asfortranarray(0.5 * numpy.diag([XYStep] * 6)[:, :4])
    dmat=numpy.concatenate((dg, -dg, numpy.zeros((6, 1))), axis=1)
    # Add the deltas to multiple copies of the closed orbit
    in_mat = orbit.reshape(6, 1) + dmat

    refs = () if refpts is None else refpts
    out_mat = numpy.squeeze(at.lattice_pass(ring, in_mat, refpts=refs, keep_lattice=keeplattice), axis=3)
    # out_mat: 8 particles at n refpts for one turn
    # (x + d) - (x - d) / d
    m44 = (in_mat[:4, :4] - in_mat[:4, 4:-1]) / XYStep

    if refpts is not None:
        mstack = (out_mat[:4, :4, :] - out_mat[:4, 4:-1, :]) / XYStep
        if full:
            mstack = numpy.stack(map(mrotate, numpy.split(mstack, mstack.shape[2], axis=2)), axis=2)
        if output_orbit:
            return m44, mstack, out_mat[:, -1, :]
        else:
            return m44, mstack
    else:
        return m44


def find_m66(ring, refpts=None, orbit=None, output_orbit=False, **kwargs):
    """find_m66 numerically finds the 6x6 transfer matrix of an accelerator lattice
    by differentiation of lattice_pass near the closed orbit.
    FINDM66 uses find_orbit6 to search for the closed orbit in 6-D
    In order for this to work the ring MUST have a CAVITY element

    m66 = find_m66(lattice)
        returns the full one-turn 6-by-6 matrix at the entrance of the first element

    m66, t = find_m66(lattice, refpts)
        returns 6x6 transfer matrices between the entrance of the first element and each element indexed by refpts.
        t is 6x6xNrefs array.

    ... = find_m66(lattice, ..., orbit=closed_orbit) - Does not search for closed orbit.
        Does not search for the closed orbit. Instead closed_orbit,a vector of initial conditions is used.
        This syntax is useful to specify the entrance orbit if lattice is not a ring or to avoid recomputing the
        closed orbit if it is already known.

    m66, t, orbit = find_m66(lattice, ..., output_orbit=True)
        Returns in addition the closed orbit at the entrance of the 1st element

    See also find_m44, find_orbit6

    """
    XYStep = kwargs.pop('XYStep', XYDEFSTEP)
    DPStep = kwargs.pop('DPStep', DPSTEP)
    keeplattice = False
    if orbit is None:
        orbit = find_orbit6(ring)
        keeplattice = True

    # Construct matrix of plus and minus deltas
    scaling = numpy.array([XYStep, XYStep, XYStep, XYStep, DPStep, DPStep])
    dg = numpy.asfortranarray(0.5*numpy.diag(scaling))
    dmat = numpy.concatenate((dg, -dg, numpy.zeros((6,1))), axis=1)

    in_mat = orbit.reshape(6, 1) + dmat

    refs = () if refpts is None else refpts
    out_mat = numpy.squeeze(at.lattice_pass(ring, in_mat, refpts=refs, keep_lattice=keeplattice), axis=3)
    # out_mat: 12 particles at n refpts for one turn
    # (x + d) - (x - d) / d
    m66 = (in_mat[:, :6] - in_mat[:, 6:-1]) / scaling.reshape((1, 6))

    if refpts is not None:
        mstack = (out_mat[:, :6, :] - out_mat[:, 6:-1, :]) / XYStep
        if output_orbit:
            return m66, mstack, out_mat[:, -1, :]
        else:
            return m66, mstack
    else:
        return m66


def betatron_phase_unwrap(m):
    """
    Unwrap negative jumps in betatron.
    """
    dp = numpy.diff(m)
    jumps = numpy.append([0], dp) < 0
    return m + numpy.cumsum(jumps) * numpy.pi


def get_twiss(ring, dp=0.0, refpts=None, get_chrom=False, ddp=DDP):
    """
    Determine Twiss parameters by first finding the transfer matrix.
    """

    def twiss22(mat, ms):
        """
        Calculate Twiss parameters from the standard 2x2 transfer matrix
        (i.e. x or y).
        """
        sin_mu_end = (numpy.sign(mat[0, 1]) *
                      math.sqrt(-mat[0, 1] * mat[1, 0] -
                                (mat[0, 0] - mat[1, 1]) ** 2 / 4))
        alpha0 = (mat[0, 0] - mat[1, 1]) / 2.0 / sin_mu_end
        beta0 = mat[0, 1] / sin_mu_end
        beta = ((ms[0, 0, :] * beta0 - ms[0, 1, :] * alpha0) **
                2 + ms[0, 1, :] ** 2) / beta0
        alpha = -((ms[0, 0, :] * beta0 - ms[0, 1, :] * alpha0) *
                  (ms[1, 0, :] * beta0 - ms[1, 1, :] * alpha0) +
                  ms[0, 1, :] * ms[1, 1, :]) / beta0
        mu = numpy.arctan(ms[0, 1, :] / (ms[0, 0, :] * beta0 - ms[0, 1, :] * alpha0))
        mu = betatron_phase_unwrap(mu)
        return alpha, beta, mu

    chrom = None

    refpts = at.uint32_refpts(refpts, len(ring))
    nrefs = refpts.size
    if refpts[-1] != len(ring):
        refpts = numpy.append(refpts, [len(ring)])

    orbit4, orbit = find_orbit4(ring, dp, refpts)
    m44, mstack = find_m44(ring, dp, refpts, orbit=orbit4)

    ax, bx, mx = twiss22(m44[:2, :2], mstack[:2, :2, :])
    ay, by, my = twiss22(m44[2:, 2:], mstack[2:, 2:, :])

    tune = numpy.array((mx[-1], my[-1])) / (2 * numpy.pi)
    twiss = numpy.zeros(nrefs, dtype=TWISS_DTYPE)
    twiss['idx'] = refpts[:nrefs]
    twiss['s_pos'] = at.get_s_pos(ring, refpts[:nrefs])
    twiss['closed_orbit'] = numpy.rollaxis(orbit, -1)[:nrefs]
    twiss['m44'] = numpy.rollaxis(mstack, -1)[:nrefs]
    twiss['alpha'] = numpy.rollaxis(numpy.vstack((ax, ay)), -1)[:nrefs]
    twiss['beta'] = numpy.rollaxis(numpy.vstack((bx, by)), -1)[:nrefs]
    twiss['mu'] = numpy.rollaxis(numpy.vstack((mx, my)), -1)[:nrefs]
    twiss['dispersion'] = numpy.NaN
    # Calculate chromaticity by calling this function again at a slightly
    # different momentum.
    if get_chrom:
        twissb, tuneb, _ = get_twiss(ring, dp + ddp, refpts[:nrefs])
        chrom = (tuneb - tune) / ddp
        twiss['dispersion'] = (twissb['closed_orbit'] - twiss['closed_orbit'])[:,:4] / ddp

    return twiss, tune, chrom
