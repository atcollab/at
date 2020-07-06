"""
Coupled or non-coupled 4x4 linear motion
"""
import numpy
from math import sqrt, atan2, pi, factorial
from at.lattice import Lattice, check_radiation, uint32_refpts, get_s_pos, \
    bool_refpts
from at.tracking import lattice_pass
from at.physics import find_orbit4, find_m44, jmat
from .harmonic_analysis import get_tunes_harmonic

__all__ = ['get_twiss', 'linopt', 'avlinopt', 'get_mcf', 'get_tune',
           'get_chrom']

DDP = 1e-8

_jmt = jmat(1)

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
    bbb = ms[:, 0, 1]
    aaa = ms[:, 0, 0] * beta0 - bbb * alpha0
    beta = (aaa * aaa + bbb * bbb) / beta0
    alpha = -(aaa * (ms[:, 1, 0] * beta0 - ms[:, 1, 1] * alpha0) +
              bbb * ms[:, 1, 1]) / beta0
    mu = numpy.arctan2(bbb, aaa)
    # Unwrap negative jumps in betatron phase advance
    dmu = numpy.diff(numpy.append([0], mu))
    jumps = dmu < 0
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


@check_radiation(False)
def get_twiss(ring, dp=0.0, refpts=None, get_chrom=False, orbit=None,
              keep_lattice=False, ddp=DDP):
    """
    Perform linear analysis of the NON-COUPLED lattices

    twiss0, tune, chrom, twiss = get_twiss(ring, dp[, refpts])

    PARAMETERS
        ring            lattice description.
        dp=0.0          momentum deviation.
        refpts=None     elements at which data is returned. It can be:
                        1) an integer in the range [-len(ring), len(ring)-1]
                           selecting the element according to python indexing
                           rules. As a special case, len(ring) is allowed and
                           refers to the end of the last element,
                        2) an ordered list of such integers without duplicates,
                        3) a numpy array of booleans of maximum length
                           len(ring)+1, where selected elements are True.

    KEYWORDS
        orbit           avoids looking for the closed orbit if is already known
                        ((6,) array)
        get_chrom=False compute dispersion and chromaticities. Needs computing
                        the optics at 2 different momentum deviations around
                        the central one.
        keep_lattice    Assume no lattice change since the previous tracking.
                        Defaults to False
        ddp=1.0E-8      momentum deviation used for computation of
                        chromaticities and dispersion

    OUTPUT
        twiss0          linear optics data at the entrance/end of the ring
        tune            [tune_h, tune_v], fractional part of the linear tunes
        chrom           [ksi_h , ksi_v], chromaticities ksi = d(nu)/(dP/P).
                        Only computed if 'get_chrom' is True
        twiss           linear optics at the points refered to by refpts, if
                        refpts is None an empty twiss structure is returned.

        twiss is a record array with fields:
        idx             element index in the ring
        s_pos           longitudinal position [m]
        closed_orbit    (6,) closed orbit vector
        dispersion      (4,) dispersion vector
                        Only computed if 'get_chrom' is True
        m44             (4, 4) transfer matrix M from the beginning of ring
                        to the entrance of the element
        mu              (2,) betatron phase (modulo 2*pi)
        beta            [betax, betay] vector
        alpha           [alphax, alphay] vector
        All values given at the entrance of each element specified in refpts.

        Field values can be obtained with either
        twiss['idx']    or
        twiss.idx


    See also linopt
    """
    uintrefs = uint32_refpts(refpts, len(ring))

    if orbit is None:
        orbit, _ = find_orbit4(ring, dp, keep_lattice=keep_lattice)
        keep_lattice = True

    orbs = numpy.squeeze(
        lattice_pass(ring, orbit.copy(order='K'), refpts=uintrefs,
                     keep_lattice=keep_lattice),
        axis=(1, 3)).T
    m44, mstack = find_m44(ring, dp, uintrefs, orbit=orbit, keep_lattice=True)
    nrefs = uintrefs.size

    # Get initial twiss parameters
    a0_x, b0_x, tune_x = _closure(m44[:2, :2])
    a0_y, b0_y, tune_y = _closure(m44[2:, 2:])
    tune = numpy.array([tune_x, tune_y])

    # Calculate chromaticity by calling this function again at a slightly
    # different momentum.
    if get_chrom:
        d0_up, tune_up, _, l_up = get_twiss(ring, dp + 0.5 * ddp,
                                            uintrefs, keep_lattice=True)
        d0_down, tune_down, _, l_down = get_twiss(ring, dp - 0.5 * ddp,
                                                  uintrefs, keep_lattice=True)
        chrom = (tune_up - tune_down) / ddp
        dispersion = (l_up['closed_orbit'] -
                      l_down['closed_orbit'])[:, :4] / ddp
        disp0 = (d0_up['closed_orbit'] - d0_down['closed_orbit'])[:4] / ddp
    else:
        chrom = numpy.array([numpy.NaN, numpy.NaN])
        dispersion = numpy.NaN
        disp0 = numpy.NaN

    twiss0 = numpy.rec.fromarrays(
        (len(ring), get_s_pos(ring, len(ring))[0], orbit, disp0,
         numpy.array([a0_x, a0_y]), numpy.array([b0_x, b0_y]),
         2.0 * pi * tune, m44), dtype=TWISS_DTYPE)

    # Propagate to reference points
    twiss = numpy.rec.array(numpy.zeros(nrefs, dtype=TWISS_DTYPE))
    if nrefs > 0:
        alpha_x, beta_x, mu_x = _twiss22(mstack[:, :2, :2], a0_x, b0_x)
        alpha_z, beta_z, mu_z = _twiss22(mstack[:, 2:, 2:], a0_y, b0_y)

        twiss['idx'] = uintrefs
        twiss['s_pos'] = get_s_pos(ring, uintrefs[:nrefs])
        twiss['closed_orbit'] = orbs
        twiss['m44'] = mstack
        twiss['alpha'] = numpy.stack((alpha_x, alpha_z), axis=1)
        twiss['beta'] = numpy.stack((beta_x, beta_z), axis=1)
        twiss['mu'] = numpy.stack((mu_x, mu_z), axis=1)
        twiss['dispersion'] = dispersion

    return twiss0, tune, chrom, twiss


# noinspection PyPep8Naming
@check_radiation(False)
def linopt(ring, dp=0.0, refpts=None, get_chrom=False, orbit=None,
           keep_lattice=False, ddp=DDP, coupled=True):
    """
    Perform linear analysis of a lattice

    lindata0, tune, chrom, lindata = linopt(ring, dp[, refpts])

    PARAMETERS
        ring            lattice description.
        dp=0.0          momentum deviation.
        refpts=None     elements at which data is returned. It can be:
                        1) an integer in the range [-len(ring), len(ring)-1]
                           selecting the element according to python indexing
                           rules. As a special case, len(ring) is allowed and
                           refers to the end of the last element,
                        2) an ordered list of such integers without duplicates,
                        3) a numpy array of booleans of maximum length
                           len(ring)+1, where selected elements are True.

    KEYWORDS
        orbit           avoids looking for the closed orbit if is already known
                        ((6,) array)
        get_chrom=False compute dispersion and chromaticities. Needs computing
                        the optics at 2 different momentum deviations around
                        the central one.
        keep_lattice    Assume no lattice change since the previous tracking.
                        Defaults to False
        ddp=1.0E-8      momentum deviation used for computation of
                        chromaticities and dispersion
        coupled=True    if False, simplify the calculations by assuming
                        no H/V coupling

    OUTPUT
        lindata0        linear optics data at the entrance/end of the ring
        tune            [tune_A, tune_B], linear tunes for the two normal modes
                        of linear motion [1]
        chrom           [ksi_A , ksi_B], chromaticities ksi = d(nu)/(dP/P).
                        Only computed if 'get_chrom' is True
        lindata         linear optics at the points refered to by refpts, if
                        refpts is None an empty lindata structure is returned.

        lindata is a record array with fields:
        idx             element index in the ring
        s_pos           longitudinal position [m]
        closed_orbit    (6,) closed orbit vector
        dispersion      (4,) dispersion vector
                        Only computed if 'get_chrom' is True
        m44             (4, 4) transfer matrix M from the beginning of ring
                        to the entrance of the element [2]
        A               (2, 2) matrix A in [3]
        B               (2, 2) matrix B in [3]
        C               (2, 2) matrix C in [3]
        gamma           gamma parameter of the transformation to eigenmodes
        mu              [mux, muy], betatron phase (modulo 2*pi)
        beta            [betax, betay] vector
        alpha           [alphax, alphay] vector
        All values given at the entrance of each element specified in refpts.

        Field values can be obtained with either
        lindata['idx']    or
        lindata.idx

    REFERENCES
        [1] D.Edwars,L.Teng IEEE Trans.Nucl.Sci. NS-20, No.3, p.885-888, 1973
        [2] E.Courant, H.Snyder
        [3] D.Sagan, D.Rubin Phys.Rev.Spec.Top.-Accelerators and beams,
            vol.2 (1999)

    See also get_twiss

    """

    # noinspection PyShadowingNames
    def analyze(r44):
        t44 = r44.reshape((4, 4))
        mm = t44[:2, :2]
        nn = t44[2:, 2:]
        m = t44[:2, 2:]
        n = t44[2:, :2]
        gamma = sqrt(numpy.linalg.det(numpy.dot(n, C) + numpy.dot(G, nn)))
        msa = (G.dot(mm) - m.dot(_jmt.dot(C.T.dot(_jmt.T)))) / gamma
        msb = (numpy.dot(n, C) + numpy.dot(G, nn)) / gamma
        cc = (numpy.dot(mm, C) + numpy.dot(G, m)).dot(
            _jmt.dot(msb.T.dot(_jmt.T)))
        return msa, msb, gamma, cc

    uintrefs = uint32_refpts([] if refpts is None else refpts, len(ring))

    if orbit is None:
        orbit, _ = find_orbit4(ring, dp, keep_lattice=keep_lattice)
        keep_lattice = True
    orbs = numpy.squeeze(
        lattice_pass(ring, orbit.copy(order='K'), refpts=uintrefs,
                     keep_lattice=keep_lattice),
        axis=(1, 3)).T
    m44, mstack = find_m44(ring, dp, uintrefs, orbit=orbit, keep_lattice=True)
    nrefs = uintrefs.size

    M = m44[:2, :2]
    N = m44[2:, 2:]
    m = m44[:2, 2:]
    n = m44[2:, :2]

    if coupled:
        # Calculate A, B, C, gamma at the first element
        H = m + _jmt.dot(n.T.dot(_jmt.T))
        t = numpy.trace(M - N)
        t2 = t * t
        t2h = t2 + 4.0 * numpy.linalg.det(H)

        g = sqrt(1.0 + sqrt(t2 / t2h)) / sqrt(2.0)
        G = numpy.diag((g, g))
        C = -H * numpy.sign(t) / (g * sqrt(t2h))
        A = G.dot(G.dot(M)) - numpy.dot(G, (
            m.dot(_jmt.dot(C.T.dot(_jmt.T))) + C.dot(n))) + C.dot(
            N.dot(_jmt.dot(C.T.dot(_jmt.T))))
        B = G.dot(G.dot(N)) + numpy.dot(G, (
            _jmt.dot(C.T.dot(_jmt.T.dot(m))) + n.dot(C))) + _jmt.dot(
            C.T.dot(_jmt.T.dot(M.dot(C))))
    else:
        A = M
        B = N
        C = numpy.zeros((2, 2))
        g = 1.0

    # Get initial twiss parameters
    a0_a, b0_a, tune_a = _closure(A)
    a0_b, b0_b, tune_b = _closure(B)
    tune = numpy.array([tune_a, tune_b])

    if get_chrom:
        d0_up, tune_up, _, l_up = linopt(ring, dp + 0.5 * ddp, uintrefs,
                                         keep_lattice=True,
                                         coupled=coupled)
        d0_down, tune_down, _, l_down = linopt(ring, dp - 0.5 * ddp, uintrefs,
                                               keep_lattice=True,
                                               coupled=coupled)
        chrom = (tune_up - tune_down) / ddp
        dispersion = (l_up['closed_orbit'] -
                      l_down['closed_orbit'])[:, :4] / ddp
        disp0 = (d0_up['closed_orbit'] - d0_down['closed_orbit'])[:4] / ddp
    else:
        chrom = numpy.array([numpy.NaN, numpy.NaN])
        dispersion = numpy.NaN
        disp0 = numpy.NaN

    lindata0 = numpy.rec.fromarrays(
        (len(ring), get_s_pos(ring, len(ring))[0], orbit, disp0,
         numpy.array([a0_a, a0_b]), numpy.array([b0_a, b0_b]),
         2.0 * pi * tune, m44, A, B, C, g),
        dtype=LINDATA_DTYPE)

    # Propagate to reference points
    lindata = numpy.rec.array(numpy.zeros(nrefs, dtype=LINDATA_DTYPE))
    if nrefs > 0:
        if coupled:
            MSA, MSB, gamma, CL = zip(*[analyze(ms44) for ms44 in mstack])
            msa = numpy.stack(MSA, axis=0)
            msb = numpy.stack(MSB, axis=0)
        else:
            msa = mstack[:, :2, :2]
            msb = mstack[:, 2:, 2:]
            gamma = 1.0
            CL = numpy.zeros((1, 2, 2))

        alpha_a, beta_a, mu_a = _twiss22(msa, a0_a, b0_a)
        alpha_b, beta_b, mu_b = _twiss22(msb, a0_b, b0_b)

        lindata['idx'] = uintrefs
        lindata['s_pos'] = get_s_pos(ring, uintrefs)
        lindata['closed_orbit'] = orbs
        lindata['m44'] = mstack
        lindata['alpha'] = numpy.stack((alpha_a, alpha_b), axis=1)
        lindata['beta'] = numpy.stack((beta_a, beta_b), axis=1)
        lindata['mu'] = numpy.stack((mu_a, mu_b), axis=1)
        lindata['A'] = [ms.dot(A.dot(_jmt.dot(ms.T.dot(_jmt.T))))
                        for ms in msa]
        lindata['B'] = [ms.dot(B.dot(_jmt.dot(ms.T.dot(_jmt.T))))
                        for ms in msb]
        lindata['C'] = CL
        lindata['gamma'] = gamma
        lindata['dispersion'] = dispersion

    return lindata0, tune, chrom, lindata


# noinspection PyPep8Naming
@check_radiation(False)
def avlinopt(ring, dp=0.0, refpts=None, **kwargs):
    """
    Perform linear analysis of a lattice and returns average beta, dispersion
    and phase advance

    lindata,avebeta,avemu,avedisp,tune,chrom = avlinopt(ring, dp[, refpts])

    PARAMETERS
        ring            lattice description.
        dp=0.0          momentum deviation.
        refpts=None     elements at which data is returned. It can be:
                        1) an integer in the range [-len(ring), len(ring)-1]
                           selecting the element according to python indexing
                           rules. As a special case, len(ring) is allowed and
                           refers to the end of the last element,
                        2) an ordered list of such integers without duplicates,
                        3) a numpy array of booleans of maximum length
                           len(ring)+1, where selected elements are True.

    KEYWORDS
        orbit           avoids looking for the closed orbit if is already known
                        ((6,) array)
        keep_lattice    Assume no lattice change since the previous tracking.
                        Defaults to False
        ddp=1.0E-8      momentum deviation used for computation of
                        chromaticities and dispersion
        coupled=True    if False, simplify the calculations by assuming
                        no H/V coupling

    OUTPUT
        lindata         linear optics at the points refered to by refpts, if
                        refpts is None an empty lindata structure is returned.
                        See linopt for details
        avebeta         Average beta functions [betax,betay] at refpts
        avemu           Average phase advances [mux,muy] at refpts
        avedisp         Average dispersion [Dx,Dx',Dy,Dy',muy] at refpts
        avespos         Average s position at refpts
        tune            [tune_A, tune_B], linear tunes for the two normal modes
                        of linear motion [1]
        chrom           [ksi_A , ksi_B], chromaticities ksi = d(nu)/(dP/P).
                        Only computed if 'get_chrom' is True


    See also get_twiss,linopt

    """
    def get_strength(elem):
        try:
            k = elem.PolynomB[1]
        except (AttributeError, IndexError):
            k = 0.0
        return k

    def betadrift(beta0, beta1, alpha0, lg):
        gamma0 = (alpha0 * alpha0 + 1) / beta0
        return 0.5 * (beta0 + beta1) - gamma0 * lg * lg / 6

    def betafoc(beta1, alpha0, alpha1, k2, lg):
        gamma1 = (alpha1 * alpha1 + 1) / beta1
        return 0.5 * ((gamma1 + k2 * beta1) * lg + alpha1 - alpha0) / k2 / lg

    def dispfoc(dispp0, dispp1, k2, lg):
        return (dispp0 - dispp1) / k2 / lg

    boolrefs = bool_refpts([] if refpts is None else refpts, len(ring))
    length = numpy.array([el.Length for el in ring[boolrefs]])
    strength = numpy.array([get_strength(el) for el in ring[boolrefs]])
    longelem = bool_refpts([], len(ring))
    longelem[boolrefs] = (length != 0)

    shorti_refpts = (~longelem) & boolrefs
    longi_refpts = longelem & boolrefs
    longf_refpts = numpy.roll(longi_refpts, 1)

    all_refs = shorti_refpts | longi_refpts | longf_refpts
    _, tune, chrom, d_all = linopt(ring, dp=dp, refpts=all_refs,
                                   get_chrom=True, **kwargs)
    lindata = d_all[boolrefs[all_refs]]

    avebeta = lindata.beta.copy()
    avemu = lindata.mu.copy()
    avedisp = lindata.dispersion.copy()
    aves = lindata.s_pos.copy()

    di = d_all[longi_refpts[all_refs]]
    df = d_all[longf_refpts[all_refs]]

    long = (length != 0.0)
    kfoc = (strength != 0.0)
    foc = long & kfoc
    nofoc = long & (~kfoc)
    K2 = numpy.stack((strength[foc], -strength[foc]), axis=1)
    fff = foc[long]
    length = length.reshape((-1, 1))

    avemu[long] = 0.5 * (di.mu + df.mu)
    aves[long] = 0.5 * (df.s_pos + di.s_pos)
    avebeta[nofoc] = \
        betadrift(di.beta[~fff], df.beta[~fff], di.alpha[~fff], length[nofoc])
    avebeta[foc] = \
        betafoc(df.beta[fff], di.alpha[fff], df.alpha[fff], K2, length[foc])
    avedisp[numpy.ix_(long, [1, 3])] = \
        (df.dispersion[:, [0, 2]] - di.dispersion[:, [0, 2]]) / length[long]
    idx = numpy.ix_(~fff, [0, 2])
    avedisp[numpy.ix_(nofoc, [0, 2])] = (di.dispersion[idx] +
                                         df.dispersion[idx]) * 0.5
    idx = numpy.ix_(fff, [1, 3])
    avedisp[numpy.ix_(foc, [0, 2])] = \
        dispfoc(di.dispersion[idx], df.dispersion[idx], K2, length[foc])
    return lindata, avebeta, avemu, avedisp, aves, tune, chrom


@check_radiation(False)
def get_mcf(ring, dp=0.0, ddp=DDP, keep_lattice=False):
    """Compute momentum compaction factor

    PARAMETERS
        ring            lattice description
        dp              momentum deviation. Defaults to 0

    KEYWORDS
        keep_lattice    Assume no lattice change since the previous tracking.
                        Defaults to False
        ddp=1.0E-8      momentum deviation used for differentiation
    """
    fp_a, _ = find_orbit4(ring, dp=dp - 0.5 * ddp, keep_lattice=keep_lattice)
    fp_b, _ = find_orbit4(ring, dp=dp + 0.5 * ddp, keep_lattice=True)
    fp = numpy.stack((fp_a, fp_b),
                     axis=0).T  # generate a Fortran contiguous array
    b = numpy.squeeze(lattice_pass(ring, fp, keep_lattice=True), axis=(2, 3))
    ring_length = get_s_pos(ring, len(ring))
    return (b[5, 1] - b[5, 0]) / ddp / ring_length[0]


@check_radiation(False)
def get_tune(ring, method='linopt', dp=0, **kwargs):
    """
    gets the tune using several available methods

    method can be 'linopt', 'fft' or 'laskar'
    linopt: returns the tune from the linopt function
    fft: tracks a single particle (one for x and one for y)
    and computes the tune from the fft
    laskar: tracks a single particle (one for x and one for y)
    and computes the harmonic components


    INPUT
    for linopt:
        no input needed

    for harmonic:
        nturns: number of turn
        amplitude: amplitude of oscillation
        method: laskar or fft
        num_harmonics: number of harmonic components to compute
        (before mask applied, default=20)
        fmin/fmax: determine the boundaries within which the tune is
        located [default 0->1]
        hann: flag to turn on hanning window [default-> False]
        remove_dc: Removes the mean offset of oscillation data

    OUTPUT
        tunes = np.array([Qx,Qy])
    """

    def gen_centroid(ring, ampl, nturns, dp, remove_dc):
        orbit, _ = find_orbit4(ring,dp)
        p0 = numpy.array([orbit,]*2).T
        p0[0, 0] += ampl
        p0[2, 1] += ampl
        p1 = lattice_pass(ring, p0, refpts=len(ring), nturns=nturns)
        cent_x = p1[0, 0, 0, :]
        cent_xp = p1[1, 0, 0, :]
        cent_y = p1[2, 1, 0, :]
        cent_yp = p1[3, 1, 0, :]
        if remove_dc:
            cent_x -= numpy.mean(cent_x)
            cent_y -= numpy.mean(cent_y)
            cent_xp -= numpy.mean(cent_xp)
            cent_yp -= numpy.mean(cent_yp)
        return cent_x - 1j*cent_xp, cent_y - 1j*cent_yp

    if method == 'linopt':
        _, tunes, _, _ = linopt(ring, dp=dp)
    else:
        nturns = kwargs.pop('nturns', 512)
        ampl = kwargs.pop('ampl', 1.0e-6)
        remove_dc = kwargs.pop('remove_dc', True)
        cent_x, cent_y = gen_centroid(ring, ampl, nturns, dp, remove_dc)
        cents = numpy.vstack((cent_x, cent_y))
        tunes = get_tunes_harmonic(cents, method,**kwargs)
    return tunes


@check_radiation(False)
def get_chrom(ring, method='linopt', dp=0, **kwargs):
    """
    gets the chromaticty using several available methods

    method can be 'linopt', 'fft' or 'laskar'
    linopt: returns the chromaticity from the linopt function
    fft: tracks a single particle (one for x and one for y)
    and computes the tune from the fft
    harmonic: tracks a single particle (one for x and one for y)
    and computes the harmonic components

    see get_tune for kwargs inputs

    OUTPUT
        chromaticities = np.array([Q'x,Q'y])
    """

    if method=='fft':
        print('Warning fft method not accurate to get the '+ 
              'chromaticity')

    ddp = kwargs.pop('ddp', DDP)
    tune_up = get_tune(ring, method=method, dp=dp + 0.5 * ddp, **kwargs)
    tune_down = get_tune(ring, method=method, dp=dp - 0.5 * ddp, **kwargs)
    chrom = (tune_up - tune_down) / ddp
    return numpy.array(chrom)


Lattice.linopt = linopt
Lattice.avlinopt = avlinopt
Lattice.get_mcf = get_mcf
Lattice.avlinopt = avlinopt
Lattice.get_tune = get_tune
Lattice.get_chrom = get_chrom
