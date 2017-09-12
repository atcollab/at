from __future__ import print_function
import numpy
import at
from at import lattice
import math

EPS = 1e-10
XYDEFSTEP = 6.055454452393343e-006  # Optimal delta?
DDP = 1e-8


# dtype for structured array containing Twiss parameters
TWISS_DTYPE = [('idx', numpy.uint32),
               ('s_pos', numpy.float64),
               ('closed_orbit', numpy.float64, (4,)),
               ('dispersion', numpy.float64, (4,)),
               ('alpha', numpy.float64, (2,)),
               ('beta', numpy.float64, (2,)),
               ('mu', numpy.float64, (2,)),
               ('m44', numpy.float64, (4, 4))]


def find_orbit4(ring, dp=0.0, refpts=None):
    """findorbit4 finds the closed orbit in the 4-d transverse phase
    space by numerically solving for a fixed point of the one turn
    map M calculated with linepass

        (X, PX, Y, PY, dP2, CT2 ) = M (X, PX, Y, PY, dP1, CT1)

    under the CONSTANT MOMENTUM constraint, dP2 = dP1 = dP and
    there is NO constraint on the 6-th coordinate CT

    IMPORTANT!!! findorbit4 imposes a constraint on dP and relaxes
    the constraint on the revolution frequency. A physical storage
    ring does exactly the opposite: the momentum deviation of a
    particle on the closed orbit settles at the value
    such that the revolution is synchronous with the RF cavity

                HarmNumber*Frev = Frf

    To impose this artificial constraint in findorbit4, PassMethod
    used for any element SHOULD NOT
    1. change the longitudinal momentum dP (cavities , magnets with radiation)
    2. have any time dependence (localized impedance, fast kickers etc)

    findorbit4(RING, dP) is 4x1 vector - fixed point at the
    entrance of the 1-st element of the RING (x,px,y,py)

    findorbit4(RING, dP, refpts) is 4-by-Length(refpts)
    array of column vectors - fixed points (x,px,y,py)
    at the entrance of each element indexed refpts array.
    refpts is an array of increasing indexes that select elements
    from the range 0 to length(RING).
    See further explanation of refpts in the 'help' for FINDSPOS

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
    STEP_SIZE = 1e-6
    MAX_ITERATIONS = 20
    CONVERGENCE = 1e-12
    r_in = numpy.zeros((6,), order='F')
    r_in[4] = dp
    delta_matrix = numpy.zeros((6, 5), order='F')
    for i in range(4):
        delta_matrix[i, i] += STEP_SIZE
    change = 1
    itercount = 0
    keeplattice = False
    while (change > CONVERGENCE) and itercount < MAX_ITERATIONS:
        in_mat = r_in.reshape((6, 1)) + delta_matrix
        out_mat = at.linepass(ring, in_mat, keep_lattice=keeplattice)
        # the reference particle after one turn
        ref_out = out_mat[:4, 4]
        # 4x4 jacobian matrix from numerical differentiation:
        # f(x+d) - f(x) / d
        j4 = (out_mat[:4, :4] - ref_out.reshape((4, 1))) / STEP_SIZE
        a = j4 - numpy.identity(4)  # f'(r_n) - 1
        b = ref_out - r_in[:4]
        b_over_a, _, _, _ = numpy.linalg.lstsq(a, b)
        r_next = r_in - numpy.append(b_over_a, numpy.zeros((2,)))
        # determine if we are close enough
        change = numpy.linalg.norm(r_next - r_in)
        itercount += 1
        r_in = r_next
        keeplattice = True

    all_points = at.linepass(ring, r_in, refpts, keep_lattice=keeplattice)
    return r_in[:4], all_points[:4, :]


def find_m44(ring, dp=0.0, refpts=None, orbit4=None, XYStep=XYDEFSTEP):
    """
    Determine the transfer matrix for ring, by first finding the closed orbit.
    """
    refpts = lattice.uint32_refpts(refpts, len(ring))
    nrefs = refpts.size
    if refpts[-1] != len(ring):
        refpts = numpy.append(refpts, [len(ring)])
    keeplattice = False
    if orbit4 is None:
        orbit4, _ = find_orbit4(ring, dp)
        keeplattice = True
    # Append zeros to closed 4-orbit.
    orbit6 = numpy.append(orbit4, (dp, 0.0)).reshape(6, 1)
    # Construct matrix of plus and minus deltas
    dmat = numpy.zeros((6, 8), order='F')
    for i in range(4):
        dmat[i, i] = 0.5 * XYStep
        dmat[i, i + 4] = -0.5 * XYStep
    # Add the deltas to multiple copies of the closed orbit
    in_mat = orbit6 + dmat

    out_mat = at.linepass(ring, in_mat, refpts, keep_lattice=keeplattice)
    tmat3 = numpy.reshape(out_mat[:4, :], (4, 8, -1), order='F')
    # (x + d) - (x - d) / d
    mstack = (tmat3[:, :4, :] - tmat3[:, 4:8, :]) / XYStep
    m44 = mstack[:, :, -1]
    return m44, mstack[:, :, :nrefs]


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

    refpts = lattice.uint32_refpts(refpts, len(ring))
    nrefs = refpts.size
    if refpts[-1] != len(ring):
        refpts = numpy.append(refpts, [len(ring)])

    orbit4, orbit = find_orbit4(ring, dp, refpts)
    m44, mstack = find_m44(ring, dp, refpts, orbit4=orbit4)

    ax, bx, mx = twiss22(m44[:2, :2], mstack[:2, :2, :])
    ay, by, my = twiss22(m44[2:, 2:], mstack[2:, 2:, :])

    tune = numpy.array((mx[-1], my[-1])) / (2 * numpy.pi)
    twiss = numpy.zeros(nrefs, dtype=TWISS_DTYPE)
    twiss['idx'] = refpts[:nrefs]
    twiss['s_pos'] = lattice.get_s_pos(ring, refpts[:nrefs])
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
        twiss['dispersion'] = (twissb['closed_orbit'] - twiss['closed_orbit']) / ddp

    return twiss, tune, chrom


def m66(ring):
    epsmat = EPS * numpy.identity(6)
    mm = at.atpass(ring, epsmat, 1) / EPS
    return numpy.transpose(mm)
