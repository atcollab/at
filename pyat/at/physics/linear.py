import numpy
from numpy.linalg import multi_dot as md
from math import sqrt, acos
from .matrix import find_m44, _s2
from .orbit import find_orbit4
from ..lattice import uint32_refpts, get_s_pos
from ..track import lattice_pass

__all__ = ['get_twiss', 'linopt']

DDP = 1e-8

# dtype for structured array containing Twiss parameters
TWISS_DTYPE = [('idx', numpy.uint32),
               ('s_pos', numpy.float64),
               ('closed_orbit', numpy.float64, (6,)),
               ('dispersion', numpy.float64, (4,)),
               ('alpha', numpy.float64, (2,)),
               ('beta', numpy.float64, (2,)),
               ('mu', numpy.float64, (2,)),
               ('m44', numpy.float64, (4, 4))]

# dtype for structured array containing linopt parameters
LINDATA_DTYPE = TWISS_DTYPE + [('A', numpy.float64, (2, 2)),
                               ('B', numpy.float64, (2, 2)),
                               ('C', numpy.float64, (2, 2)),
                               ('gamma', numpy.float64)]

def betatron_phase_unwrap(m):
    """
    Unwrap negative jumps in betatron.
    """
    dp = numpy.diff(m)
    jumps = numpy.append([0], dp) < 0
    return m + numpy.cumsum(jumps) * 2.0 * numpy.pi


def twiss22(ms, alpha0, beta0):
    """
    Calculate Twiss parameters from the standard 2x2 transfer matrix
    (i.e. x or y).
    """
    bbb = ms[0, 1, :]
    aaa = ms[0, 0, :] * beta0 - bbb * alpha0
    beta = (aaa * aaa + bbb * bbb) / beta0
    alpha = -(aaa * (ms[1, 0, :] * beta0 - ms[1, 1, :] * alpha0) + bbb * ms[1, 1, :]) / beta0
    mu = numpy.arctan2(bbb, aaa)
    mu = betatron_phase_unwrap(mu)
    return alpha, beta, mu


def ab(m22):
    diff = (m22[0, 0] - m22[1, 1]) / 2.0
    sinmu = numpy.sign(m22[0, 1]) * sqrt(-m22[0, 1] * m22[1, 0] - diff * diff)
    alpha = diff / sinmu
    beta = m22[0, 1] / sinmu
    return alpha, beta


def _tune(m44):
    cos_mu = 0.5 * numpy.trace(m44)
    tune = acos(cos_mu)/2.0/numpy.pi
    if m44[0, 1] < 0.0:
        tune = 1.0 - tune
    return tune


def get_twiss(ring, dp=0.0, refpts=None, get_chrom=False, orbit=None, keep_lattice=False, ddp=DDP):
    """Determine Twiss parameters by first finding the transfer matrix.

    The Twiss structured array has nrefs elements, so:
     * twiss['idx'].shape is (nrefs,)
     * twiss['closed_orbit'].shape is (nrefs, 6).

    Returns:
        twiss - structured array
        tune - numpy array of shape (2,)
        chrom - numpy array of shape (2,)
    """
    refs = () if refpts is None else refpts

    if orbit is None:
        orbit = find_orbit4(ring, dp)
        keep_lattice = True

    orbs = numpy.squeeze(lattice_pass(ring, orbit.copy(order='K'), refpts=refs,
                                      keep_lattice=keep_lattice))
    m44, mstack = find_m44(ring, dp, refs, orbit=orbit, keep_lattice=True)
    nrefs = mstack.shape[2]

    alpha_x, beta_x, mu_x = twiss22(mstack[:2, :2, :], *ab(m44[:2, :2]))
    alpha_z, beta_z, mu_z = twiss22(mstack[2:, 2:, :], *ab(m44[2:, 2:]))

    tune = numpy.array([_tune(m44[:2, :2]), _tune(m44[2:, 2:])])
    twiss = numpy.zeros(nrefs, dtype=TWISS_DTYPE)
    twiss['idx'] = uint32_refpts(refs, len(ring))
    # Use rollaxis to get the arrays in the correct shape for the twiss
    # structured array - that is, with nrefs as the first dimension.
    twiss['s_pos'] = get_s_pos(ring, refs[:nrefs])
    twiss['closed_orbit'] = numpy.rollaxis(orbs, -1)
    twiss['m44'] = numpy.rollaxis(mstack, -1)
    twiss['alpha'] = numpy.rollaxis(numpy.vstack((alpha_x, alpha_z)), -1)
    twiss['beta'] = numpy.rollaxis(numpy.vstack((beta_x, beta_z)), -1)
    twiss['mu'] = numpy.rollaxis(numpy.vstack((mu_x, mu_z)), -1)
    # Calculate chromaticity by calling this function again at a slightly
    # different momentum.
    if get_chrom:
        twissb, tuneb, _ = get_twiss(ring, dp + ddp, refs[:nrefs])
        chrom = (tuneb - tune) / ddp
        twiss['dispersion'] = (twissb['closed_orbit'] - twiss['closed_orbit'])[:, :4] / ddp
    else:
        chrom = None
        twiss['dispersion'] = numpy.NaN

    return twiss, tune, chrom


def linopt(ring, dp=0.0, refpts=None, get_chrom=False, orbit=None, keep_lattice=False, ddp=DDP):
    def analyze(r44):
        t44=r44.reshape((4, 4))
        mm = t44[:2, :2]
        nn = t44[2:, 2:]
        m = t44[:2, 2:]
        n = t44[2:, :2]
        gamma = sqrt(numpy.linalg.det(numpy.dot(n, C) + numpy.dot(G, nn)))
        msa = (md((G, mm)) - md((m, _s2, C.T, _s2.T))) / gamma
        msb = (numpy.dot(n, C) + numpy.dot(G, nn)) / gamma
        cc = md(((numpy.dot(mm, C) + numpy.dot(G, m)), _s2, msb.T, _s2.T))
        aa = md((msa, A, _s2, msa.T, _s2.T))
        bb = md((msb, B, _s2, msb.T, _s2.T))
        return msa, msb, gamma, cc, aa, bb

    refs = () if refpts is None else refpts
    if orbit is None:
        orbit = find_orbit4(ring, dp)
        keep_lattice = True
    orbs = numpy.squeeze(lattice_pass(ring, orbit.copy(order='K'), refpts=refpts,
                                      keep_lattice=keep_lattice))
    m44, mstack = find_m44(ring, dp, refs, orbit=orbit, keep_lattice=True)
    nrefs = mstack.shape[2]

    # Calculate A, B, C, gamma at the first element
    M = m44[:2, :2]
    N = m44[2:, 2:]
    m = m44[:2, 2:]
    n = m44[2:, :2]

    H = m + md((_s2, n.T, _s2.T))
    t = numpy.trace(M - N)
    t2 = t * t
    t2h = t2 + 4.0 * numpy.linalg.det(H)

    g = sqrt(1.0 + sqrt(t2 / t2h)) / sqrt(2.0)
    G = numpy.diag((g, g))
    C = -H * numpy.sign(t) / (g * sqrt(t2h))
    A = md((G, G, M)) - numpy.dot(G, (md((m, _s2, C.T, _s2.T)) + md((C, n)))) + md((C, N, _s2, C.T, _s2.T))
    B = md((G, G, N)) + numpy.dot(G, (md((_s2, C.T, _s2.T, m)) + md((n, C)))) + md((_s2, C.T, _s2.T, M, C))
    MSA, MSB, gamma, CL, AL, BL = zip(*map(analyze, numpy.split(mstack, mstack.shape[2], axis=2)))
    alpha_a, beta_a, mu_a = twiss22(numpy.stack(MSA, axis=2), *ab(A))
    alpha_b, beta_b, mu_b = twiss22(numpy.stack(MSB, axis=2), *ab(B))

    tune = numpy.array([_tune(A), _tune(B)])
    lindata = numpy.zeros(nrefs, dtype=LINDATA_DTYPE)
    lindata['idx'] = uint32_refpts(refs, len(ring))
    # Use rollaxis to get the arrays in the correct shape for the lindata
    # structured array - that is, with nrefs as the first dimension.
    lindata['s_pos'] = get_s_pos(ring, refs)
    lindata['closed_orbit'] = numpy.rollaxis(orbs, -1)
    lindata['m44'] = numpy.rollaxis(mstack, -1)
    lindata['alpha'] = numpy.rollaxis(numpy.vstack((alpha_a, alpha_b)), -1)
    lindata['beta'] = numpy.rollaxis(numpy.vstack((beta_a, beta_b)), -1)
    lindata['mu'] = numpy.rollaxis(numpy.vstack((mu_a, mu_b)), -1)
    lindata['A'] = AL
    lindata['B'] = BL
    lindata['C'] = CL
    lindata['gamma'] = gamma
    lindata['dispersion'] = numpy.NaN
    if get_chrom:
        l_up, tune_up, _ = linopt(ring, dp+0.5*ddp, refs, keep_lattice=True)
        l_down, tune_down, _ = linopt(ring, dp-0.5*ddp, refs, keep_lattice=True)
        chrom = (tune_up-tune_down)/ddp
        lindata['dispersion'] = (l_up['closed_orbit'] - l_down['closed_orbit'])[:, :4] / ddp
    else:
        chrom = None
        lindata['dispersion'] = numpy.NaN

    return lindata, tune, chrom
