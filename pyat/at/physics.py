import numpy
import at
from at import lattice
import math
import collections


EPS = 1e-10
DP = 1e-8
DDP = 1e-8


Twiss = collections.namedtuple('Twiss', ('refpts', 's_pos', 'closed_orbit',
                                         'm44', 'alpha', 'beta', 'mu', 'tune',
                                         'dispersion', 'chrom'))


def find_orbit4(ring, dp=DP, refpts=None):
    """
    Determine the closed orbit assuming constant momentum deviation.

    The closed orbit does not change in a pass of the ring with the
    specified momentum difference.  We seek
     - f(x) = x
     - g(x) = f(x) - x = 0
     - g'(x) = f'(x) - 1
    Use a Newton-Raphson-type algorithm:
     - r_n+1 = r_n - g(r_n) / g'(r_n)
     - r_n+1 = r_n - (f(r_n) - r_n) / (f'(r_n) - 1)

    (f(r_n) - r_n) / (f'(r_n) - 1) can be seen as x = b/a where we use least
        squares fitting to determine x when ax = b
    f(r_n) - r_n is denoted b
    f'(r_n) is the 4x4 jacobian, denoted j4
    """
    STEP_SIZE = 1e-6
    MAX_ITERATIONS = 20
    CONVERGENCE = 1e-12
    refpts = lattice.get_refpts(refpts, len(ring))
    r_in = numpy.zeros((1, 6))
    r_in[0, 4] = dp
    delta_matrix = numpy.zeros((5, 6))
    for i in range(4):
        delta_matrix[i, i] += STEP_SIZE
    change = 1
    itercount = 0
    while (change > CONVERGENCE) and itercount < MAX_ITERATIONS:
        in_mat = r_in + delta_matrix
        out_mat = at.atpass(ring, in_mat, 1)
        # the reference particle after one turn
        ref_out = out_mat[4, :4]
        # 4x4 jacobian matrix from numerical differentiation:
        # f(x+d) - f(x) / d
        j4 = (out_mat[:4, :4] - ref_out) / STEP_SIZE
        a = j4 - numpy.identity(4)  # f'(r_n) - 1
        b = ref_out - r_in[:, :4]
        # transpose the arrays to (4,4) (4,1)
        b_over_a, _, _, _ = numpy.linalg.lstsq(a.T, b.T)
        r_out = r_in - numpy.append(b_over_a, numpy.zeros((1, 2)))
        # determine if we are close enough
        change = numpy.linalg.norm(r_out - r_in)
        itercount += 1
        r_in = r_out

    all_points = at.atpass(ring, r_in, 1, refpts)
    return all_points[:, :4]


def find_m44(ring, dp=DP, refpts=None):
    """
    Determine the transfer matrix for ring, by first finding the closed orbit.
    """
    last_included = False if refpts is None else len(ring) in refpts
    refpts = lattice.get_refpts(refpts, len(ring), append_last=True)
    # Optimal delta?
    d = 6.055454452393343e-006
    orbit4 = find_orbit4(ring, dp)
    # Append zeros to closed 4-orbit.
    bottom = numpy.array([dp, 0]).reshape(1, 2)
    orbit6 = numpy.append(orbit4, bottom, axis=1)
    # Construct matrix of plus and minus deltas
    dmat = numpy.append(numpy.identity(4)*d, -numpy.identity(4)*d, axis=0)
    dmat = numpy.append(dmat, numpy.zeros((8, 2)), axis=1)
    # Add the deltas to multiple copies of the closed orbit
    in_mat = (numpy.zeros((8, 6)) + orbit6) + dmat

    out_mat = at.atpass(ring, in_mat, 1, refpts)
    out_mat = out_mat[:, :4].reshape(-1, 8, 4)
    # (x + d) - (x - d) / 2*d
    m44 = (out_mat[-1, :4, :4] - out_mat[-1, 4:8, :4]) / (2 * d)
    out_mat = out_mat if last_included else out_mat[:-1, :, :]
    mstack = (out_mat[:, :4, :4] - out_mat[:, 4:8, :4]) / (2 * d)
    return m44, mstack


def betatron_phase_unwrap(m):
    """
    Unwrap negative jumps in betatron.
    """
    dp = numpy.diff(m)
    jumps = numpy.append([0], dp) < 0
    return m + numpy.cumsum(jumps) * numpy.pi


def get_twiss(ring, dp=DP, refpts=None, get_chrom=False, ddp=DDP):
    """
    Determine Twiss parameters by first finding the transfer matrix.
    """

    def twiss22(mat, ms):
        """
        Calculate Twiss parameters from the standard 2x2 transfer matrix
        (i.e. x or y).
        """
        sin_mu_end = (numpy.sign(mat[1, 0]) *
                      math.sqrt(-mat[1, 0] * mat[0, 1] -
                                (mat[0, 0] - mat[1, 1]) ** 2 / 4))
        alpha_end = (mat[0, 0] - mat[1, 1]) / 2 / sin_mu_end
        beta_end = mat[1, 0] / sin_mu_end
        beta = ((ms[:, 0, 0] * beta_end - ms[:, 1, 0] * alpha_end) **
                2 + ms[:, 1, 0] ** 2) / beta_end
        alpha = -((ms[:, 0, 0] * beta_end - ms[:, 1, 0] * alpha_end) *
                  (ms[:, 0, 1] * beta_end - ms[:, 1, 1] * alpha_end) +
                   ms[:, 1, 0] * ms[:, 1, 1]) / beta_end
        mu = numpy.arctan(ms[:, 1, 0] /
                          (ms[:, 0, 0] * beta_end - ms[:, 1, 0] * alpha_end))
        mu = betatron_phase_unwrap(mu)
        return alpha, beta, mu

    chrom = None
    dispersion = None

    refpts = lattice.get_refpts(refpts, len(ring))
    # We're calling find_orbit4 twice here.  We should do better.
    orbit = find_orbit4(ring, dp, refpts)
    m44, mstack = find_m44(ring, dp, refpts)
    s_pos = lattice.get_s_pos(ring, refpts)

    ax, bx, mx = twiss22(m44[:2, :2], mstack[:, :2, :2])
    ay, by, my = twiss22(m44[2:, 2:], mstack[:, 2:, 2:])

    tune = numpy.array((mx[-1], my[-1])) / (2 * numpy.pi)
    beta = numpy.append(bx, by).reshape(2, -1)
    alpha = numpy.append(ax, ay).reshape(2, -1)
    mu = numpy.append(mx, my).reshape(2, -1)
    # Calculate chromaticity by calling this function again at a slightly
    # different momentum.
    if get_chrom:
        dtwiss = get_twiss(ring, dp + ddp, refpts)
        dispersion = (dtwiss.closed_orbit - orbit) / ddp
        chrom = (dtwiss.tune - tune) / ddp

    twiss = Twiss(refpts, s_pos, orbit, m44, alpha,
                  beta, mu, tune, dispersion, chrom)
    return twiss


def m66(ring):
    epsmat = EPS * numpy.identity(6)
    mm = at.atpass(ring, epsmat, 1) / EPS
    return numpy.transpose(mm)
