import numpy
from numpy.linalg import multi_dot as md
from math import sqrt, atan2, pi
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


def _twiss22(ms, alpha0, beta0):
    """
    Calculate Twiss parameters from the standard 2x2 transfer matrix
    (i.e. x or y).
    """
    bbb = ms[0, 1, :]
    aaa = ms[0, 0, :] * beta0 - bbb * alpha0
    beta = (aaa * aaa + bbb * bbb) / beta0
    alpha = -(aaa * (ms[1, 0, :] * beta0 - ms[1, 1, :] * alpha0) + bbb * ms[1, 1, :]) / beta0
    mu = numpy.arctan2(bbb, aaa)
    # Unwrap negative jumps in betatron phase advance
    dmu = numpy.diff(mu)
    jumps = numpy.append([0], dmu) < 0
    mu += numpy.cumsum(jumps) * 2.0 * numpy.pi
    return alpha, beta, mu


def _closure(m22):
    diff = (m22[0, 0] - m22[1, 1]) / 2.0
    sinmu = numpy.sign(m22[0, 1]) * sqrt(-m22[0, 1] * m22[1, 0] - diff * diff)
    cosmu = 0.5 * numpy.trace(m22)
    alpha = diff / sinmu
    beta = m22[0, 1] / sinmu
    tune = (atan2(sinmu, cosmu) / 2.0 / pi) % 1
    return alpha, beta, tune


def get_twiss(ring, dp=0.0, refpts=None, get_chrom=False, orbit=None, keep_lattice=False, ddp=DDP):
    """
    Perform linear analysis of the NON-COUPLED lattices

    PARAMETERS
        ring            lattice description
        dp              momentum deviation. Defaults to 0
        refpts          elements at which data is returned. It can be
                        1) an integer (0 indicating the first element)
                        2) a list of integers
                        3) a numpy array of booleans as long as ring where
                           selected elements are true
                        Defaults to None

    KEYWORDS
        orbit           avoids looking for the colsed orbit if is already known ((6,) array)
        get_chrom       compute dispersion and chromaticities. Needs computing the optics
                        at 2 different momentum deviations around the central one.
                        Defaults to False
        keep_lattice    Assume no lattice change since the previous tracking.
                        Defaults to False

    OUTPUT
        twiss           linear optics data
        tune            [tune_h, tune_v], fractional part of the linear tunes
        chrom           [ksi_h , ksi_v], vector of chromaticities ksi = d(nu)/(dP/P).
                        Only computed if 'get_chrom' is True

        twiss is a structured array with fields:
        idx             element index in the ring                           (nrefs,)
        s_pos           longitudinal position [m]                           (nrefs,)
        closed_orbit    closed orbit vector with                            (nrefs, 6)
        dispersion      dispersion vector.                                  (nrefs, 4)
                        Only computed if 'get_chrom' is True                (nrefs, 4)
        m44             4x4 transfer matrix M from the beginning of ring    (nrefs, 4, 4)
                        to the entrance of the element [2]
        mu              [mux, muy], A and B betatron phase                  (nrefs, 2)
        beta            [betax, betay] vector                               (nrefs, 2)
        alpha           [alphax, alphay] vector                             (nrefs, 2)
        All values are given at the entrance of each element specified in refpts.

    See also linopts
    """
    uintrefs = uint32_refpts([] if refpts is None else refpts, len(ring))

    if orbit is None:
        orbit = find_orbit4(ring, dp)
        keep_lattice = True

    orbs = numpy.squeeze(lattice_pass(ring, orbit.copy(order='K'), refpts=uintrefs,
                                      keep_lattice=keep_lattice))
    m44, mstack = find_m44(ring, dp, uintrefs, orbit=orbit, keep_lattice=True)
    nrefs = uintrefs.size

    # Get initial twiss parameters
    a0_x, b0_x, tune_x = _closure(m44[:2, :2])
    a0_y, b0_y, tune_y = _closure(m44[2:, 2:])
    tune = numpy.array([tune_x, tune_y])

    # Calculate chromaticity by calling this function again at a slightly
    # different momentum.
    if get_chrom:
        l_up, tune_up, _ = get_twiss(ring, dp + 0.5 * ddp, uintrefs, keep_lattice=True)
        l_down, tune_down, _ = get_twiss(ring, dp - 0.5 * ddp, uintrefs, keep_lattice=True)
        chrom = (tune_up - tune_down) / ddp
        dispersion = (l_up['closed_orbit'] - l_down['closed_orbit'])[:, :4] / ddp
        d0 = numpy.NaN
    else:
        chrom = None
        dispersion = numpy.NaN
        d0 = numpy.NaN

    data0 = numpy.array(
        (0, 0.0, orbit, d0, numpy.array([a0_x, a0_y]), numpy.array([b0_x, b0_y]), 0.0, m44), dtype=TWISS_DTYPE)

    # Propagate to reference points
    if nrefs > 0:
        alpha_x, beta_x, mu_x = _twiss22(mstack[:2, :2, :], a0_x, b0_x)
        alpha_z, beta_z, mu_z = _twiss22(mstack[2:, 2:, :], a0_y, b0_y)

        twiss = numpy.zeros(nrefs, dtype=TWISS_DTYPE)
        twiss['idx'] = uintrefs
        # Use rollaxis to get the arrays in the correct shape for the twiss
        # structured array - that is, with nrefs as the first dimension.
        twiss['s_pos'] = get_s_pos(ring, uintrefs[:nrefs])
        twiss['closed_orbit'] = numpy.rollaxis(orbs, -1)
        twiss['m44'] = numpy.rollaxis(mstack, -1)
        twiss['alpha'] = numpy.stack((alpha_x, alpha_z), axis=1)
        twiss['beta'] = numpy.stack((beta_x, beta_z), axis=1)
        twiss['mu'] = numpy.stack((mu_x, mu_z), axis=1)
        twiss['dispersion'] = dispersion
    else:
        twiss = numpy.array([], dtype=TWISS_DTYPE)

    return twiss, tune, chrom


def linopt(ring, dp=0.0, refpts=None, get_chrom=False, orbit=None, keep_lattice=False, ddp=DDP):
    """
    Perform linear analysis of the COUPLED lattices

    lindata, tune, chrom = linopt(ring, dp, refpts)

    PARAMETERS
        ring            lattice description
        dp              momentum deviation. Defaults to 0
        refpts          Optional: elements at which data is returned. It can be
                        1) an integer (0 indicating the first element)
                        2) a list of integers
                        3) a numpy array of booleans as long as ring where
                           selected elements are true
                        Defaults to None

    KEYWORDS
        orbit           avoids looking for the colsed orbit if is already known ((6,) array)
        get_chrom       compute dispersion and chromaticities. Needs computing the optics
                        at 2 different momentum deviations around the central one.
                        Defaults to False
        keep_lattice    Assume no lattice change since the previous tracking.
                        Defaults to False

    OUTPUT
        lindata0        linear optics data at the entrance of the ring
        tune            [tune_A, tune_B], linear tunes for the two normal modes of linear motion [1]
        chrom           [ksi_A , ksi_B], vector of chromaticities ksi = d(nu)/(dP/P).
                        Only computed if 'get_chrom' is True
        lindata         Only returned if refpts is not None:
                        linear optics at the points refered to by refpts

        lindata is a structured array with fields:
        idx             element index in the ring                           (nrefs,)
        s_pos           longitudinal position [m]                           (nrefs,)
        closed_orbit    closed orbit vector with                            (nrefs, 6)
        dispersion      dispersion vector.                                  (nrefs, 4)
                        Only computed if 'get_chrom' is True                (nrefs, 4)
        m44             4x4 transfer matrix M from the beginning of ring    (nrefs, 4, 4)
                        to the entrance of the element [2]
        A               (2, 2) matrix A in [3]                              (nrefs, 2, 2)
        B               (2, 2) matrix B in [3]                              (nrefs, 2, 2)
        C               (2, 2) matrix C in [3]                              (nrefs, 2, 2)
        gamma           gamma parameter of the transformation to eigenmodes (nrefs,)
        mu              [mux, muy], A and B betatron phase                  (nrefs, 2)
        beta            [betax, betay] vector                               (nrefs, 2)
        alpha           [alphax, alphay] vector                             (nrefs, 2)
        All values are given at the entrance of each element specified in refpts.

    REFERENCES
        [1] D.Edwars,L.Teng IEEE Trans.Nucl.Sci. NS-20, No.3, p.885-888, 1973
        [2] E.Courant, H.Snyder
        [3] D.Sagan, D.Rubin Phys.Rev.Spec.Top.-Accelerators and beams, vol.2 (1999)

    See also get_twiss

    """

    def analyze(r44):
        t44 = r44.reshape((4, 4))
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

    uintrefs = uint32_refpts([] if refpts is None else refpts, len(ring))

    if orbit is None:
        orbit = find_orbit4(ring, dp)
        keep_lattice = True
    orbs = numpy.squeeze(lattice_pass(ring, orbit.copy(order='K'), refpts=refpts,
                                      keep_lattice=keep_lattice))
    m44, mstack = find_m44(ring, dp, uintrefs, orbit=orbit, keep_lattice=True)
    nrefs = uintrefs.size

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

    # Get initial twiss parameters
    a0_a, b0_a, tune_a = _closure(A)
    a0_b, b0_b, tune_b = _closure(B)
    tune = numpy.array([tune_a, tune_b])

    if get_chrom:
        d0_up, tune_up, _, l_up = linopt(ring, dp + 0.5 * ddp, uintrefs, keep_lattice=True)
        d0_down, tune_down, _, l_down = linopt(ring, dp - 0.5 * ddp, uintrefs, keep_lattice=True)
        chrom = (tune_up - tune_down) / ddp
        dispersion = (l_up['closed_orbit'] - l_down['closed_orbit'])[:, :4] / ddp
        d0 = (d0_up['closed_orbit'] - d0_down['closed_orbit'])[:4] / ddp
    else:
        chrom = None
        dispersion = numpy.NaN
        d0 = numpy.NaN

    data0 = numpy.array(
        (0, 0.0, orbit, d0, numpy.array([a0_a, a0_b]), numpy.array([b0_a, b0_b]), 0.0, m44, A, B, C, g),
        dtype=LINDATA_DTYPE)

    # Propagate to reference points
    if nrefs > 0:
        MSA, MSB, gamma, CL, AL, BL = zip(*map(analyze, numpy.split(mstack, mstack.shape[2], axis=2)))
        alpha_a, beta_a, mu_a = _twiss22(numpy.stack(MSA, axis=2), a0_a, b0_a)
        alpha_b, beta_b, mu_b = _twiss22(numpy.stack(MSB, axis=2), a0_b, b0_b)

        lindata = numpy.zeros(nrefs, dtype=LINDATA_DTYPE)
        lindata['idx'] = uintrefs
        # Use rollaxis to get the arrays in the correct shape for the lindata
        # structured array - that is, with nrefs as the first dimension.
        lindata['s_pos'] = get_s_pos(ring, uintrefs)
        lindata['closed_orbit'] = numpy.rollaxis(orbs, -1)
        lindata['m44'] = numpy.rollaxis(mstack, -1)
        lindata['alpha'] = numpy.stack((alpha_a, alpha_b), axis=1)
        lindata['beta'] = numpy.stack((beta_a, beta_b), axis=1)
        lindata['mu'] = numpy.stack((mu_a, mu_b), axis=1)
        lindata['A'] = AL
        lindata['B'] = BL
        lindata['C'] = CL
        lindata['gamma'] = gamma
        lindata['dispersion'] = dispersion
    else:
        lindata = numpy.array([], dtype=LINDATA_DTYPE)

    if refpts is None:
        return data0, tune, chrom
    else:
        return data0, tune, chrom, lindata
