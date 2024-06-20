"""
Coupled or non-coupled 4x4 linear motion
"""
from __future__ import annotations
import numpy
from math import sqrt, pi, sin, cos, atan2
from collections.abc import Callable
import warnings
from scipy.linalg import solve
from ..constants import clight
from ..lattice import DConstant, Refpts, get_bool_index, get_uint32_index
from ..lattice import AtWarning, Lattice, Orbit, check_6d, get_s_pos
from ..lattice import AtError
from ..lattice import frequency_control
from ..tracking import internal_lpass
from .orbit import find_orbit4, find_orbit6
from .matrix import find_m44, find_m66
from .amat import a_matrix, jmat, jmatswap
from .harmonic_analysis import get_tunes_harmonic

__all__ = ['linopt', 'linopt2', 'linopt4', 'linopt6', 'avlinopt',
           'get_optics', 'get_tune', 'get_chrom']

_S = jmat(1)
_S2 = numpy.array([[0, 1], [-1, 0]], dtype=numpy.float64)

# dtype for structured array containing linopt parameters
_DATA2_DTYPE = [('alpha', numpy.float64, (2,)),
                ('beta', numpy.float64, (2,)),
                ('mu', numpy.float64, (2,))
                ]

_DATA4_DTYPE = [('alpha', numpy.float64, (2,)),
                ('beta', numpy.float64, (2,)),
                ('mu', numpy.float64, (2,)),
                ('gamma', numpy.float64),
                ('A', numpy.float64, (2, 2)),
                ('B', numpy.float64, (2, 2)),
                ('C', numpy.float64, (2, 2))
                ]

_DATA6_DTYPE = [('alpha', numpy.float64, (2,)),
                ('beta', numpy.float64, (2,)),
                ('mu', numpy.float64, (3,)),
                ('R', numpy.float64, (3, 6, 6)),
                ('A', numpy.float64, (6, 6)),
                ('dispersion', numpy.float64, (4,)),
                ]

_DATAX_DTYPE = [('alpha', numpy.float64, (2,)),
                ('beta', numpy.float64, (2,)),
                ('mu', numpy.float64, (2,)),
                ('R', numpy.float64, (2, 4, 4)),
                ('A', numpy.float64, (4, 4)),
                ]

_W24_DTYPE = [('W', numpy.float64, (2,)),
              ('Wp', numpy.float64, (2,)),
              ('dalpha', numpy.float64, (2,)),
              ('dbeta', numpy.float64, (2,)),
              ('dmu', numpy.float64, (2,)),
              ('ddispersion', numpy.float64, (4,)),
              ]


_W6_DTYPE = [('W', numpy.float64, (2,)),
             ('Wp', numpy.float64, (2,)),
             ('dalpha', numpy.float64, (2,)),
             ('dbeta', numpy.float64, (2,)),
             ('dmu', numpy.float64, (3,)),
             ('dR', numpy.float64, (3, 6, 6)),
             ('ddispersion', numpy.float64, (4,)),
             ]

_WX_DTYPE = [('W', numpy.float64, (2,)),
             ('Wp', numpy.float64, (2,)),
             ('dalpha', numpy.float64, (2,)),
             ('dbeta', numpy.float64, (2,)),
             ('dmu', numpy.float64, (2,)),
             ('dR', numpy.float64, (2, 4, 4)),
             ('ddispersion', numpy.float64, (4,)),
             ]


_IDX_DTYPE = [('idx', numpy.uint32)]


def _twiss22(t12, alpha0, beta0):
    """Propagate Twiss parameters"""
    bbb = t12[:, 0, 1]
    aaa = t12[:, 0, 0] * beta0 - bbb * alpha0
    beta = (aaa * aaa + bbb * bbb) / beta0
    alpha = -(aaa * (t12[:, 1, 0] * beta0 - t12[:, 1, 1] * alpha0) +
              bbb * t12[:, 1, 1]) / beta0
    mu = numpy.arctan2(bbb, aaa)
    # Unwrap negative jumps in betatron phase advance
    # dmu = numpy.diff(numpy.append([0], mu))
    # jumps = dmu < -1.0e-3
    # mu += numpy.cumsum(jumps) * 2.0 * numpy.pi
    return alpha, beta, mu


def _closure(m22):
    diff = (m22[0, 0] - m22[1, 1]) / 2.0
    try:
        sinmu = numpy.sign(m22[0, 1]) * sqrt(-m22[0, 1]*m22[1, 0] - diff*diff)
        cosmu = 0.5 * numpy.trace(m22)
        alpha = diff / sinmu
        beta = m22[0, 1] / sinmu
        return alpha, beta, cosmu + sinmu*1j
    except ValueError:          # Unstable motion
        return numpy.nan, numpy.nan, numpy.nan


# noinspection PyShadowingNames,PyPep8Naming
def _tunes(ring, **kwargs):
    """"""
    if ring.is_6d:
        nd = 3
        mt, _ = find_m66(ring, **kwargs)
    else:
        nd = 2
        mt, _ = find_m44(ring, **kwargs)
    try:
        _, vps = a_matrix(mt)
        tunes = numpy.mod(numpy.angle(vps) / 2.0 / pi, 1.0)
    except AtError:
        warnings.warn(AtWarning('Unstable ring'))
        tunes = numpy.empty(nd)
        tunes[:] = numpy.nan
    return tunes


def _analyze2(mt, ms):
    """Analysis of a 2D 1-turn transfer matrix"""
    A = mt[:2, :2]
    B = mt[2:, 2:]
    alp0_a, bet0_a, vp_a = _closure(A)
    alp0_b, bet0_b, vp_b = _closure(B)
    vps = numpy.array([vp_a, vp_b])
    el0 = (numpy.array([alp0_a, alp0_b]),
           numpy.array([bet0_a, bet0_b]),
           0.0)
    alpha_a, beta_a, mu_a = _twiss22(ms[:, :2, :2], alp0_a, bet0_a)
    alpha_b, beta_b, mu_b = _twiss22(ms[:, 2:, 2:], alp0_b, bet0_b)
    els = (numpy.stack((alpha_a, alpha_b), axis=1),
           numpy.stack((beta_a, beta_b), axis=1),
           numpy.stack((mu_a, mu_b), axis=1))
    return vps, _DATA2_DTYPE, el0, els, _W24_DTYPE


def _analyze4(mt, ms):
    """Analysis of a 4D 1-turn transfer matrix according to Sagan, Rubin"""
    def propagate(t12):
        M = t12[:2, :2]
        N = t12[2:, 2:]
        m = t12[:2, 2:]
        n = t12[2:, :2]
        ff = n @ C + g * N
        gamma = sqrt(numpy.linalg.det(ff))
        e12 = (g * M - m @ _S @ C.T @ _S.T) / gamma
        f12 = ff / gamma
        a12 = e12 @ A @ _S @ e12.T @ _S.T
        b12 = f12 @ B @ _S @ f12.T @ _S.T
        c12 = M @ C + g * m @ _S @ f12.T @ _S.T
        return e12, f12, gamma, a12, b12, c12

    M = mt[:2, :2]
    N = mt[2:, 2:]
    m = mt[:2, 2:]
    n = mt[2:, :2]
    H = m + _S @ n.T @ _S.T
    detH = numpy.linalg.det(H)
    if detH == 0.0:
        g = 1.0
        C = -H
        A = M
        B = N
    else:
        t = numpy.trace(M - N)
        t2 = t * t
        t2h = t2 + 4.0 * detH
        g2 = (1.0 + sqrt(t2 / t2h)) / 2
        g = sqrt(g2)
        C = -H * numpy.sign(t) / (g * sqrt(t2h))
        A = g2 * M - g * (m @ _S @ C.T @ _S.T + C @ n) + \
            C @ N @ _S @ C.T @ _S.T
        B = g2 * N + g * (_S @ C.T @ _S.T @ m + n @ C) + \
            _S @ C.T @ _S.T @ M @ C
    alp0_a, bet0_a, vp_a = _closure(A)
    alp0_b, bet0_b, vp_b = _closure(B)
    vps = numpy.array([vp_a, vp_b])
    el0 = (numpy.array([alp0_a, alp0_b]),
           numpy.array([bet0_a, bet0_b]),
           0.0, g, A, B, C)
    if ms.shape[0] > 0:
        e, f, g, ai, bi, ci = zip(*[propagate(mi) for mi in ms])
        alp_a, bet_a, mu_a = _twiss22(numpy.array(e), alp0_a, bet0_a)
        alp_b, bet_b, mu_b = _twiss22(numpy.array(f), alp0_b, bet0_b)
        els = (numpy.stack((alp_a, alp_b), axis=1),
               numpy.stack((bet_a, bet_b), axis=1),
               numpy.stack((mu_a, mu_b), axis=1), numpy.array(g),
               numpy.stack(ai, axis=0), numpy.stack(bi, axis=0),
               numpy.stack(ci, axis=0))
    else:
        els = (numpy.empty((0, 2)), numpy.empty((0, 2)),
               numpy.empty((0, 2)), numpy.empty((0,)),
               numpy.empty((0, 2, 2)), numpy.empty((0, 2, 2)),
               numpy.empty((0, 2, 2)))
    return vps, _DATA4_DTYPE, el0, els, _W24_DTYPE


def _analyze6(mt, ms):
    """Analysis of a 2D, 4D, 6D 1-turn transfer matrix
    according to Wolski"""
    def get_phase(a22):
        """Return the phase for A standardization"""
        return atan2(a22[0, 1], a22[0, 0])

    def standardize(aa, slcs):
        """Apply rotation to put A in std form"""
        def rot2(slc):
            rot = -get_phase(aa[slc, slc])
            cs = cos(rot)
            sn = sin(rot)
            return aa[:, slc] @ numpy.array([[cs, sn], [-sn, cs]])

        return numpy.concatenate([rot2(slc) for slc in slcs], axis=1)

    def r_matrices(ai):
        # Rk = A * S * Ik * inv(A) * S.T
        def mul2(slc):
            return ai[:, slc] @ tt[slc, slc]

        ais = numpy.concatenate([mul2(slc) for slc in slices], axis=1)
        invai = solve(ai, ss.T)
        ri = numpy.array(
            [ais[:, sl] @ invai[sl, :] for sl in slices])
        mui = numpy.array([get_phase(ai[sl, sl]) for sl in slices])
        return mui, ri, ai

    def propagate4(ri, phi, ai):
        betai = numpy.stack((ri[..., 0, 0, 0], ri[..., 1, 2, 2]), axis=-1)
        alphai = -numpy.stack((ri[..., 0, 1, 0], ri[..., 1, 3, 2]), axis=-1)
        return alphai, betai, phi, ri, ai

    def propagate6(ri, phi, ai):
        alphai, betai, phi, ri, ai = propagate4(ri, phi, ai)
        dispersion = ri[..., 2, :4, 4] / ri[..., 2, 4, 4, numpy.newaxis]
        return alphai, betai, phi, ri, ai, dispersion

    nv = mt.shape[0]
    dms = nv // 2
    slices = [slice(2 * i, 2 * (i + 1)) for i in range(dms)]
    ss = jmat(dms)
    tt = jmatswap(dms)
    if dms >= 3:
        propagate = propagate6
        dtype = _DATA6_DTYPE
        wtype = _W6_DTYPE
    else:
        propagate = propagate4
        dtype = _DATAX_DTYPE
        wtype = _WX_DTYPE
    a0, vps = a_matrix(mt)

    astd = standardize(a0, slices)
    phi0, r0, _ = r_matrices(astd)
    el0 = propagate(r0, phi0, astd)
    if ms.shape[0] > 0:
        ps, rs, aas = zip(*[r_matrices(mi @ astd) for mi in ms])
        els = propagate(numpy.array(rs), numpy.array(ps), numpy.array(aas))
    elif dms >= 3:
        els = (numpy.empty((0, dms)), numpy.empty((0, dms)),
               numpy.empty((0, dms)),
               numpy.empty((0, dms, nv, nv)), numpy.empty((0, nv, nv)),
               numpy.empty((0, 4)))
    else:
        els = (numpy.empty((0, dms)), numpy.empty((0, dms)),
               numpy.empty((0, dms)),
               numpy.empty((0, dms, nv, nv)), numpy.empty((0, nv, nv)))
    return vps, dtype, el0, els, wtype


# noinspection PyShadowingNames,PyPep8Naming
def _linopt(ring: Lattice, analyze, refpts=None, dp=None, dct=None, df=None,
            orbit=None, twiss_in=None, get_chrom=False, get_w=False,
            keep_lattice=False, mname='M', add0=(), adds=(), cavpts=None,
            **kwargs):
    """"""
    def build_sigma(orbit, dp=None):
        """Build the initial distribution at entrance of the transfer line"""
        try:
            d0 = twiss_in['dispersion']
        except (ValueError, KeyError):  # record arrays throw ValueError !
            d0 = numpy.zeros((4,))

        try:
            rmat = twiss_in['R']
        except (ValueError, KeyError):  # record arrays throw ValueError !
            rmat = None

        try:
            alphas = twiss_in['alpha']
            betas = twiss_in['beta']
        except (ValueError, KeyError):  # record arrays throw ValueError !
            alphas = numpy.zeros((2,))
            betas = numpy.ones((2,))

        dorbit = numpy.hstack((d0, 1.0, 0.0))
        if orbit is None:
            try:
                orbit = twiss_in['closed_orbit']
            except (ValueError, KeyError):  # record arrays throw ValueError !
                orbit = numpy.zeros((6,))
            if dp is not None:
                orbit = orbit + dorbit * dp

        if dp is not None:
            try:
                drmat = twiss_in['dR']
            except (ValueError, KeyError):  # record arrays throw ValueError !
                drmat = None

            try:
                dd0 = twiss_in['ddispersion']
                dalpha = twiss_in['dalpha']
                dbeta = twiss_in['dbeta']
            except (ValueError, KeyError):  # record arrays throw ValueError !
                msg = ("'get_w' option for a line requires 'twiss_in' calculated "
                       "with 'get_w' activated")
                raise AtError(msg)

            orbit = orbit + numpy.hstack((dd0, 1.0, 0.0)) * dp * dp
            dorbit = numpy.hstack((d0+dd0*dp, 1.0, 0.0))

            if (rmat is not None) and (drmat is not None):
                rmat = rmat + drmat * dp
            else:
                alphas = alphas + dalpha * dp
                betas = betas + dbeta * dp

        if rmat is not None:
            # For some reason, "emittances" must be different...
            sigm = rmat[0, ...]+10.0*rmat[1, ...]
            if rmat.shape[0] >= 3:
                sigm = sigm+0.1*rmat[2, ...]
                dorbit = rmat[2, :, 4] / rmat[2, 4, 4, numpy.newaxis]
        else:
            slices = [slice(2 * i, 2 * (i + 1)) for i in range(2)]
            ab = numpy.stack((alphas, betas), axis=1)
            sigm = numpy.zeros((4, 4))
            for slc, (alpha, beta) in zip(slices, ab):
                gamma = (1.0+alpha*alpha)/beta
                sigm[slc, slc] = numpy.array([[beta, -alpha], [-alpha, gamma]])

        return orbit, sigm, dorbit

    def chrom_w(ringup, ringdn, orbitup, orbitdn,
                refpts=None, **kwargs):
        """Compute the chromaticity and W-functions"""
        # noinspection PyShadowingNames
        def off_momentum(rng, orb0, dp=None, **kwargs):
            if twiss_in is None:
                mt, ms = get_matrix(rng, refpts=refpts, orbit=orb0, **kwargs)
                mxx = mt
                dpup = orb0[4] + 0.5 * dp_step
                dpdn = orb0[4] - 0.5 * dp_step
                o0up = None
                o0dn = None
            else:
                orbit, sigma, dorbit = build_sigma(orb0, dp=dp)
                mt, ms = get_matrix(ring, refpts=refpts, orbit=orbit, **kwargs)
                mxx = sigma @ jmat(sigma.shape[0] // 2)
                o0up = orbit + dorbit * 0.5 * dp_step
                o0dn = orbit - dorbit * 0.5 * dp_step
                dpup = None
                dpdn = None
            vps, _, el0, els, wtype = analyze(mxx, ms)
            tunes = _tunes(rng, orbit=orb0)
            o0up, oup = get_orbit(ring, refpts=refpts, guess=orb0, dp=dpup,
                                  orbit=o0up, **kwargs)
            o0dn, odn = get_orbit(ring, refpts=refpts, guess=orb0, dp=dpdn,
                                  orbit=o0dn, **kwargs)
            d0 = (o0up - o0dn)[:4] / dp_step
            ds = numpy.array([(up - dn)[:4] / dp_step for up, dn in zip(oup, odn)])
            return tunes, el0, els, d0, ds, wtype

        def wget(ddp, elup, eldn, has_r):
            """Compute the chromatic amplitude function"""
            *data_up, = elup  # Extract alpha and beta
            *data_dn, = eldn
            alpha_up, beta_up, mu_up = data_up[:3]
            alpha_dn, beta_dn, mu_dn = data_dn[:3]
            db = numpy.array(beta_up - beta_dn) / ddp
            mb = (beta_up + beta_dn) / 2
            da = numpy.array(alpha_up - alpha_dn) / ddp
            ma = (alpha_up + alpha_dn) / 2
            wa = da - ma / mb * db
            wb = db / mb
            ww = numpy.sqrt(wa ** 2 + wb ** 2)
            wp = numpy.arctan2(wa, wb)
            dmu = numpy.array(mu_up - mu_dn) / ddp
            data_out = (ww, wp, da, db, dmu)
            if has_r:
                r_up = data_up[3]
                r_dn = data_dn[3]
                data_out += (numpy.array(r_up - r_dn) / ddp, )
            return data_out

        deltap = orbitup[4] - orbitdn[4]
        *data_up, = off_momentum(ringup, orbitup, dp=0.5*deltap,
                                 **kwargs)
        *data_dn, = off_momentum(ringdn, orbitdn, dp=-0.5*deltap,
                                 **kwargs)
        tunesup, el0up, elsup, d0up, dsup, wtype = data_up
        tunesdn, el0dn, elsdn, d0dn, dsdn, _ = data_dn
        has_r = len(wtype) == 7
        # in 6D, dp comes out of find_orbit6
        chrom = (tunesup-tunesdn) / deltap
        dd0 = numpy.array(d0up - d0dn) / deltap
        dds = numpy.array(dsup - dsdn) / deltap
        data0 = wget(deltap, el0up, el0dn, has_r)
        datas = wget(deltap, elsup, elsdn, has_r)
        data0 = data0 + (dd0,)
        datas = datas + (dds,)
        return chrom, data0, datas

    def unwrap(mu):
        """Remove the phase jumps"""
        dmu = numpy.diff(numpy.concatenate((numpy.zeros((1, dms)),
                                            mu)), axis=0)
        jumps = dmu < -1.e-3
        mu += numpy.cumsum(jumps, axis=0) * 2.0 * numpy.pi

    dp_step = kwargs.get('DPStep', DConstant.DPStep)
    addtype = kwargs.pop('addtype', [])

    if ring.is_6d:
        get_matrix = find_m66
        get_orbit = find_orbit6
    else:
        get_matrix = find_m44
        get_orbit = find_orbit4

    o0up = None
    o0dn = None
    if twiss_in is None:   # Ring
        if orbit is None:
            orbit, _ = get_orbit(ring, dp=dp, dct=dct, df=df,
                                 keep_lattice=keep_lattice, **kwargs)
            keep_lattice = True
        # Get 1-turn transfer matrix
        mt, ms = get_matrix(ring, refpts=refpts, orbit=orbit, **kwargs)
        mxx = mt
    else:                       # Transfer line
        orbit, sigma, dorbit = build_sigma(orbit)
        # Get 1-turn transfer matrix
        mt, ms = get_matrix(ring, refpts=refpts, orbit=orbit, **kwargs)
        mxx = sigma @ jmat(sigma.shape[0] // 2)
        o0up = orbit + dorbit * 0.5 * dp_step
        o0dn = orbit - dorbit * 0.5 * dp_step

    # Perform analysis
    vps, dtype, el0, els, wtype = analyze(mxx, ms)
    tunes = _tunes(ring, orbit=orbit)

    if (get_chrom or get_w) and mt.shape == (6, 6):
        f0 = ring.get_rf_frequency(cavpts=cavpts)
        df = dp_step * ring.disable_6d(copy=True).slip_factor * f0
        rgup = ring.set_rf_frequency(f0 + 0.5 * df, cavpts=cavpts, copy=True)
        rgdn = ring.set_rf_frequency(f0 - 0.5 * df, cavpts=cavpts, copy=True)
        if o0up is None:
            o0up, _ = get_orbit(rgup, guess=orbit, orbit=o0up, **kwargs)
            o0dn, _ = get_orbit(rgdn, guess=orbit, orbit=o0dn, **kwargs)
    else:
        rgup = ring
        rgdn = ring

    # Propagate the closed orbit
    orb0, orbs = get_orbit(ring, refpts=refpts, orbit=orbit,
                           keep_lattice=keep_lattice)
    spos = ring.get_s_pos(refpts)

    nrefs = orbs.shape[0]
    dms = vps.size
    if dms >= 3:            # 6D processing
        dtype = dtype + [('closed_orbit', numpy.float64, (6,)),
                         ('M', numpy.float64, (2*dms, 2*dms)),
                         ('s_pos', numpy.float64)]
        data0 = (orb0, numpy.identity(2*dms), 0.0)
        datas = (orbs, ms, spos)
        length = ring.get_s_pos(len(ring))[0]
        damping_rates = -numpy.log(numpy.absolute(vps))
        damping_times = length / clight / damping_rates
    else:               # 4D processing
        kwargs['keep_lattice'] = True
        dpup = orb0[4] + 0.5 * dp_step
        dpdn = orb0[4] - 0.5 * dp_step
        o0up, oup = get_orbit(ring, refpts=refpts, guess=orb0, dp=dpup,
                              orbit=o0up, **kwargs)
        o0dn, odn = get_orbit(ring, refpts=refpts, guess=orb0, dp=dpdn,
                              orbit=o0dn, **kwargs)
        d0 = (o0up - o0dn)[:4] / dp_step
        ds = numpy.array([(up - dn)[:4] / dp_step for up, dn in zip(oup, odn)])
        dtype = dtype + [('dispersion', numpy.float64, (4,)),
                         ('closed_orbit', numpy.float64, (6,)),
                         (mname, numpy.float64, (2*dms, 2*dms)),
                         ('s_pos', numpy.float64)]
        data0 = (d0, orb0, mt, get_s_pos(ring, len(ring))[0])
        datas = (ds, orbs, ms, spos)
        damping_times = numpy.nan

    if get_w:
        dtype = dtype + wtype
        chrom, ddata0, ddatas = chrom_w(rgup, rgdn, o0up, o0dn,
                                        refpts, **kwargs)
        data0 = data0 + ddata0
        datas = datas + ddatas
    elif get_chrom:
        tunesup = _tunes(rgup, orbit=o0up)
        tunesdn = _tunes(rgdn, orbit=o0dn)
        deltap = o0up[4] - o0dn[4]
        chrom = (tunesup - tunesdn) / deltap
    else:
        chrom = numpy.nan

    beamdata = numpy.array((tunes, chrom, damping_times),
                           dtype=[('tune', numpy.float64, (dms,)),
                                  ('chromaticity', numpy.float64, (dms,)),
                                  ('damping_time', numpy.float64, (dms,))
                                  ]).view(numpy.recarray)

    dtype = dtype + addtype
    elemdata0 = numpy.array(el0+data0+add0, dtype=dtype).view(numpy.recarray)
    elemdata = numpy.recarray((nrefs,), dtype=dtype)
    if nrefs > 0:
        for name, value in zip(numpy.dtype(dtype).names, els+datas+adds):
            elemdata[name] = value
        unwrap(elemdata.mu)
    return elemdata0, beamdata, elemdata


@check_6d(False)
def linopt2(ring: Lattice, *args, **kwargs):
    r"""Linear analysis of an uncoupled lattice

    Parameters:
        ring:   Lattice description.

    Keyword Args:
        refpts (Refpts):        Elements at which data is returned.
          It can be:

          1. an integer in the range [-len(ring), len(ring)-1]
             selecting the element according to python indexing rules.
             As a special case, len(ring) is allowed and refers to the end
             of the last element,
          2. an ordered list of such integers without duplicates,
          3. a numpy array of booleans of maximum length len(ring)+1,
             where selected elements are :py:obj:`True`.
        dp (float):             Momentum deviation. Defaults to :py:obj:`None`
        dct (float):            Path lengthening. Defaults to :py:obj:`None`
        df (float):             Deviation of RF frequency. Defaults to
          :py:obj:`None`
        orbit (Orbit):          Avoids looking for the closed orbit if it is
          already known ((6,) array)
        get_chrom (bool):       Compute chromaticities. Needs computing
          the tune at 2 different momentum deviations around the central one.
        get_w (bool):           Computes chromatic amplitude functions
          (W, WP) [4]_, and derivatives of the dispersion and twiss parameters
          versus dp. Needs to compute the optics at 2 different momentum
          deviations around the central one.
        keep_lattice (bool):    Assume no lattice change since the
          previous tracking. Defaults to :py:obj:`False`
        XYStep (float):         Step size.
          Default: :py:data:`DConstant.XYStep <.DConstant>`
        DPStep (float):         Momentum step size.
          Default: :py:data:`DConstant.DPStep <.DConstant>`
        twiss_in:               Initial conditions for transfer line optics.
          Record array as output by :py:func:`.linopt6`, or dictionary. Keys:

          R or alpha, beta
            mandatory (2,) arrays
          closed_orbit
            Optional (6,) array, default 0
          dispersion
            Optional (6,) array, default 0

          If present, the attribute **R** will be used, otherwise the
          attributes **alpha** and **beta** will be used. All other attributes
          are ignored.

    Returns:
        elemdata0:      Linear optics data at the entrance of the ring
        ringdata:       Lattice properties
        elemdata:       Linear optics at the points refered to by *refpts*,
          if refpts is :py:obj:`None` an empty lindata structure is returned.

    **elemdata** is a record array with fields:

    ================    ===================================================
    **s_pos**           longitudinal position [m]
    **M**               (4, 4) transfer matrix M from the beginning of ring
                        to the entrance of the element [2]_
    **closed_orbit**    (6,) closed orbit vector
    **dispersion**      (4,) dispersion vector
    **beta**            :math:`\left[ \beta_x,\beta_y \right]` vector
    **alpha**           :math:`\left[ \alpha_x,\alpha_y \right]` vector
    **mu**              :math:`\left[ \mu_x,\mu_y \right]`, betatron phase
                        (modulo :math:`2\pi`)
    **W**               :math:`\left[ W_x,W_y \right]` only if *get_w*
                        is :py:obj:`True`: chromatic amplitude function
    **Wp**              :math:`\left[ Wp_x,Wp_y \right]` only if *get_w*
                        is :py:obj:`True`: chromatic phase function
    **dalpha**          (2,) alpha derivative vector
                        (:math:`\Delta \alpha/ \delta_p`)
    **dbeta**           (2,) beta derivative vector
                        (:math:`\Delta \beta/ \delta_p`)
    **dmu**             (2,) mu derivative vector
                        (:math:`\Delta \mu/ \delta_p`)
    **ddispersion**     (4,) dispersion derivative vector
                        (:math:`\Delta D/ \delta_p`)
    ================    ===================================================

    All values given at the entrance of each element specified in refpts.
    Field values can be obtained with either *lindata['idx']* or *lindata.idx*

    **ringdata** is a record array with fields:

    =================   ======
    **tune**            Fractional tunes
    **chromaticity**    Chromaticities, only computed if *get_chrom* is
                        :py:obj:`True`
    =================   ======

    References:
        **[1]** D.Edwards,L.Teng IEEE Trans.Nucl.Sci. NS-20, No.3, p.885-888
        , 1973

        .. [2] E.Courant, H.Snyder

        **[3]** D.Sagan, D.Rubin Phys.Rev.Spec.Top.-Accelerators and beams,
        vol.2 (1999)

        .. [4] Brian W. Montague Report LEP Note 165, CERN, 1979
    """
    return _linopt(ring, _analyze2, *args, **kwargs)


@check_6d(False)
def linopt4(ring: Lattice, *args, **kwargs):
    r"""Linear analysis of a H/V coupled lattice

    4D-analysis of coupled motion following Sagan/Rubin

    Parameters:
        ring:   Lattice description.

    Keyword Args:
        refpts (Refpts):        Elements at which data is returned.
          It can be:

          1. an integer in the range [-len(ring), len(ring)-1]
             selecting the element according to python indexing rules.
             As a special case, len(ring) is allowed and refers to the end
             of the last element,
          2. an ordered list of such integers without duplicates,
          3. a numpy array of booleans of maximum length len(ring)+1,
             where selected elements are :py:obj:`True`.
        dp (float):             Momentum deviation. Defaults to :py:obj:`None`
        dct (float):            Path lengthening. Defaults to :py:obj:`None`
        df (float):             Deviation of RF frequency. Defaults to
          :py:obj:`None`
        orbit (Orbit):          Avoids looking for the closed orbit if it is
          already known ((6,) array)
        get_chrom (bool):       Compute chromaticities. Needs computing
          the tune at 2 different momentum deviations around the central one.
        get_w (bool):           Computes chromatic amplitude functions
          (W) [8]_. Needs to compute the optics at 2 different momentum
          deviations around the central one.
        keep_lattice (bool):    Assume no lattice change since the
          previous tracking. Defaults to :py:obj:`False`
        XYStep (float):         Step size.
          Default: :py:data:`DConstant.XYStep <.DConstant>`
        DPStep (float):         Momentum step size.
          Default: :py:data:`DConstant.DPStep <.DConstant>`
        twiss_in:               Initial conditions for transfer line optics.
          Record array as output by :py:func:`.linopt6`, or dictionary. Keys:

          R or alpha, beta
            mandatory (2,) arrays
          closed_orbit
            Optional (6,) array, default 0
          dispersion
            Optional (6,) array, default 0

          If present, the attribute **R** will be used, otherwise the
          attributes **alpha** and **beta** will be used. All other attributes
          are ignored.

    Returns:
        elemdata0:      Linear optics data at the entrance of the ring
        ringdata:       Lattice properties
        elemdata:       Linear optics at the points refered to by *refpts*,
          if refpts is :py:obj:`None` an empty lindata structure is returned.

    **elemdata** is a record array with fields:

    ================    ===================================================
    **s_pos**           longitudinal position [m]
    **M**               (4, 4) transfer matrix M from the beginning of ring
                        to the entrance of the element [6]_
    **closed_orbit**    (6,) closed orbit vector
    **dispersion**      (4,) dispersion vector
    **beta**            :math:`\left[ \beta_x,\beta_y \right]` vector
    **alpha**           :math:`\left[ \alpha_x,\alpha_y \right]` vector
    **mu**              :math:`\left[ \mu_x,\mu_y \right]`, betatron phase
                        (modulo :math:`2\pi`)
    **gamma**           gamma parameter of the transformation to
                        eigenmodes [7]_
    **W**               :math:`\left[ W_x,W_y \right]` only if *get_w*
                        is :py:obj:`True`: chromatic amplitude function
    **Wp**              :math:`\left[ Wp_x,Wp_y \right]` only if *get_w*
                        is :py:obj:`True`: chromatic phase function
    **dalpha**          (2,) alpha derivative vector
                        (:math:`\Delta \alpha/ \delta_p`)
    **dbeta**           (2,) beta derivative vector
                        (:math:`\Delta \beta/ \delta_p`)
    **dmu**             (2,) mu derivative vector
                        (:math:`\Delta \mu/ \delta_p`)
    **ddispersion**     (4,) dispersion derivative vector
                        (:math:`\Delta D/ \delta_p`)
    ================    ===================================================

    All values given at the entrance of each element specified in refpts.
    Field values can be obtained with either
    *lindata['idx']* or *lindata.idx*

    **ringdata** is a record array with fields:

    =================   ======
    **tune**            Fractional tunes
    **chromaticity**    Chromaticities, only computed if *get_chrom* is
                        :py:obj:`True`
    =================   ======

    References:
        **[5]** D.Edwards,L.Teng IEEE Trans.Nucl.Sci. NS-20, No.3, p.885-888
        , 1973

        .. [6] E.Courant, H.Snyder

        .. [7] D.Sagan, D.Rubin Phys.Rev.Spec.Top.-Accelerators and beams,
           vol.2 (1999)

        .. [8] Brian W. Montague Report LEP Note 165, CERN, 1979
    """
    return _linopt(ring, _analyze4, *args, **kwargs)


@frequency_control
def linopt6(ring: Lattice, *args, **kwargs):
    r"""Linear analysis of a fully coupled lattice using normal modes

    For circular machines, :py:func:`linopt6` analyses

    * the 4x4 1-turn transfer matrix if *ring* is 4D, or
    * the 6x6 1-turn transfer matrix if  *ring* is 6D.

    For a transfer line, The "twiss_in" intput must contain either:

    *  a field **R**, as provided by ATLINOPT6, or
    * the fields **beta** and **alpha**, as provided by linopt and linopt6

    Parameters:
        ring:   Lattice description.

    Keyword Args:
        refpts (Refpts):        Elements at which data is returned.
          It can be:

          1. an integer in the range [-len(ring), len(ring)-1]
             selecting the element according to python indexing rules.
             As a special case, len(ring) is allowed and refers to the end
             of the last element,
          2. an ordered list of such integers without duplicates,
          3. a numpy array of booleans of maximum length len(ring)+1,
             where selected elements are :py:obj:`True`.
        dp (float):             Momentum deviation. Defaults to :py:obj:`None`
        dct (float):            Path lengthening. Defaults to :py:obj:`None`
        df (float):             Deviation of RF frequency. Defaults to
          :py:obj:`None`
        orbit (Orbit):          Avoids looking for the closed orbit if it is
          already known ((6,) array)
        get_chrom (bool):       Compute chromaticities. Needs computing
          the tune at 2 different momentum deviations around the central one.
        get_w (bool):           Computes chromatic amplitude functions
          (W) [11]_. Needs to compute the optics at 2 different momentum
          deviations around the central one.
        keep_lattice (bool):    Assume no lattice change since the
          previous tracking. Defaults to :py:obj:`False`
        XYStep (float):         Step size.
          Default: :py:data:`DConstant.XYStep <.DConstant>`
        DPStep (float):         Momentum step size.
          Default: :py:data:`DConstant.DPStep <.DConstant>`
        twiss_in:               Initial conditions for transfer line optics.
          Record array as output by :py:func:`.linopt6`, or dictionary. Keys:

          R or alpha, beta
            mandatory (2,) arrays
          closed_orbit
            Optional (6,) array, default 0
          dispersion
            Optional (6,) array, default 0

          If present, the attribute **R** will be used, otherwise the
          attributes **alpha** and **beta** will be used. All other attributes
          are ignored.
        cavpts (Refpts):        Cavity location for off-momentum tuning

    Returns:
        elemdata0:      Linear optics data at the entrance of the ring
        ringdata:       Lattice properties
        elemdata:       Linear optics at the points refered to by *refpts*,
          if refpts is :py:obj:`None` an empty lindata structure is returned.

    **elemdata** is a record array with fields:

    ================    ===================================================
    **s_pos**           longitudinal position [m]
    **M**               (6, 6) transfer matrix M from the beginning of ring
                        to the entrance of the element
    **closed_orbit**    (6,) closed orbit vector
    **dispersion**      (4,) dispersion vector
    **A**               (6,6) A-matrix
    **R**               (3, 6, 6) R-matrices
    **beta**            :math:`\left[ \beta_x,\beta_y \right]` vector
    **alpha**           :math:`\left[ \alpha_x,\alpha_y \right]` vector
    **mu**              :math:`\left[ \mu_x,\mu_y \right]`, betatron phase
                        (modulo :math:`2\pi`)
    **W**               :math:`\left[ W_x,W_y \right]` only if *get_w*
                        is :py:obj:`True`: chromatic amplitude function
    **Wp**              :math:`\left[ Wp_x,Wp_y \right]` only if *get_w*
                        is :py:obj:`True`: chromatic phase function
    **dalpha**          (2,) alpha derivative vector
                        (:math:`\Delta \alpha/ \delta_p`)
    **dbeta**           (2,) beta derivative vector
                        (:math:`\Delta \beta/ \delta_p`)
    **dmu**             (2,) mu derivative vector
                        (:math:`\Delta \mu/ \delta_p`)
    **ddispersion**     (4,) dispersion derivative vector
                        (:math:`\Delta D/ \delta_p`)
    **dR**              (3, 6, 6) R derivative vector
                        (:math:`\Delta R/ \delta_p`)
    ================    ===================================================

    All values given at the entrance of each element specified in refpts.
    Field values can be obtained with either
    *lindata['idx']* or *lindata.idx*

    **ringdata** is a record array with fields:

    =================   ======
    **tune**            Fractional tunes
    **chromaticity**    Chromaticities, only computed if *get_chrom* is
                        :py:obj:`True`
    **damping_time**    Damping times [s] (only if radiation is ON)
    =================   ======

    References:
        **[9]** Etienne Forest, Phys. Rev. E 58, 2481 – Published 1 August 1998

        **[10]** Andrzej Wolski, Phys. Rev. ST Accel. Beams 9, 024001 –
        Published 3 February 2006

        .. [11] Brian W. Montague Report LEP Note 165, CERN, 1979
    """
    return _linopt(ring, _analyze6, *args, **kwargs)


def linopt_auto(ring: Lattice, *args, **kwargs):
    """
    This is a convenience function to automatically switch to the faster
    :py:func:`linopt2` in case the *coupled* keyword argument is
    :py:obj:`False` **and** ring.is_6d is :py:obj:`False`.
    Otherwise the default :py:func:`linopt6` is used

    Parameters: Same as :py:func:`.linopt2` or :py:func:`.linopt6`

    Keyword Args;
        coupled (bool):     If set to :py:obj:`False`, H/V coupling
          will be ignored to simplify the calculation (needs ring.is_6d
          :py:obj:`False`)


    Returns:
        elemdata0:      Linear optics data at the entrance of the ring
        ringdata:       Lattice properties
        elemdata:       Linear optics at the points refered to by *refpts*,
          if refpts is :py:obj:`None` an empty lindata structure is returned.

    Warning:
        The output varies depending whether :py:func:`.linopt2` or
        :py:func:`.linopt6` is called. To be used with care!
    """
    if not (kwargs.pop('coupled', True) or ring.is_6d):
        return linopt2(ring, *args, **kwargs)
    else:
        return linopt6(ring, *args, **kwargs)


def get_optics(ring: Lattice, refpts: Refpts = None,
               dp: float = None,
               method: Callable = linopt6,
               **kwargs):
    """Linear analysis of a fully coupled lattice

    Parameters:
        ring:                   Lattice description.
        refpts:                 Elements at which data is returned.
          It can be:

          1. an integer in the range [-len(ring), len(ring)-1]
             selecting the element according to python indexing rules.
             As a special case, len(ring) is allowed and refers to the end
             of the last element,
          2. an ordered list of such integers without duplicates,
          3. a numpy array of booleans of maximum length len(ring)+1,
             where selected elements are :py:obj:`True`.
        dp:                     Momentum deviation.
        method (Callable):      Method for linear optics:

          :py:obj:`~.linear.linopt2`: no longitudinal motion, no H/V coupling,

          :py:obj:`~.linear.linopt4`: no longitudinal motion, Sagan/Rubin
          4D-analysis of coupled motion,

          :py:obj:`~.linear.linopt6` (default): with or without longitudinal
          motion, normal mode analysis

    Keyword Args:
        dct (float):            Path lengthening. Defaults to :py:obj:`None`
        df (float):             Deviation of RF frequency. Defaults to
          :py:obj:`None`
        orbit (Orbit):          Avoids looking for the closed orbit if it is
          already known ((6,) array)
        get_chrom (bool):       Compute chromaticities. Needs computing
          the tune at 2 different momentum deviations around the central one.
        get_w (bool):           Computes chromatic amplitude functions
          (W) [4]_. Needs to compute the optics at 2 different momentum
          deviations around the central one.
        keep_lattice (bool):    Assume no lattice change since the
          previous tracking. Defaults to :py:obj:`False`
        XYStep (float):         Step size.
          Default: :py:data:`DConstant.XYStep <.DConstant>`
        DPStep (float):         Momentum step size.
          Default: :py:data:`DConstant.DPStep <.DConstant>`
        twiss_in:               Initial conditions for transfer line optics.
          Record array as output by :py:func:`.linopt6`, or dictionary. Keys:

          R or alpha, beta
            mandatory (2,) arrays
          closed_orbit
            Optional (6,) array, default 0
          dispersion
            Optional (6,) array, default 0

          If present, the attribute **R** will be used, otherwise the
          attributes **alpha** and **beta** will be used. All other attributes
          are ignored.

    Returns:
        elemdata0:      Linear optics data at the entrance of the ring
        ringdata:       Lattice properties
        elemdata:       Linear optics at the points refered to by *refpts*,
          if refpts is :py:obj:`None` an empty lindata structure is returned.

    Warning:
        The format of output record arrays depends on the selected method.
        See :py:func:`linopt2`, :py:func:`linopt4`, :py:func:`linopt6`.
    """
    return method(ring, refpts=refpts, dp=dp, **kwargs)


# noinspection PyPep8Naming
@check_6d(False)
def linopt(ring: Lattice, dp: float = 0.0, refpts: Refpts = None,
           get_chrom: bool = False, **kwargs):
    """Linear analysis of a H/V coupled lattice (deprecated)

    Parameters:
        ring:           lattice description.
        dp:             momentum deviation.
        refpts:         elements at which data is returned. It can be:
                        1) an integer in the range [-len(ring), len(ring)-1]
                           selecting the element according to python indexing
                           rules. As a special case, len(ring) is allowed and
                           refers to the end of the last element,
                        2) an ordered list of such integers without duplicates,
                        3) a numpy array of booleans of maximum length
                           len(ring)+1, where selected elements are
                           :py:obj:`True`.
    Keyword Args:
        orbit           avoids looking for the closed orbit if is already known
                        ((6,) array)
        get_chrom=False compute chromaticities. Needs computing the tune at
                        2 different momentum deviations around the central one.
        get_w=False     computes chromatic amplitude functions (W) [4].
                        Needs to compute the optics at 2 different momentum
                        deviations around the central one.
        keep_lattice    Assume no lattice change since the previous tracking.
                        Defaults to :py:obj:`False`
        XYStep=1.0e-8   transverse step for numerical computation
        DPStep=1.0E-6   momentum deviation used for computation of
                        chromaticities and dispersion
        coupled=True    if :py:obj:`False`, simplify the calculations by
                        assuming no H/V coupling
        twiss_in=None   Initial conditions for transfer line optics. Record
                        array as output by linopt, or dictionary. Keys:
                        'alpha' and 'beta'  (mandatory)
                        'closed_orbit',     (default 0)
                        'dispersion'        (default 0)
                        All other attributes are ignored.
    OUTPUT
        lindata0        linear optics data at the entrance of the ring
        tune            [tune_A, tune_B], linear tunes for the two normal modes
                        of linear motion [1]
        chrom           [ksi_A , ksi_B], chromaticities ksi = d(nu)/(dP/P).
                        Only computed if 'get_chrom' is :py:obj:`True`
        lindata         linear optics at the points refered to by refpts, if
                        refpts is None an empty lindata structure is returned.

        lindata is a record array with fields:
        idx             element index in the ring
        s_pos           longitudinal position [m]
        m44             (4, 4) transfer matrix M from the beginning of ring
                        to the entrance of the element [2]
        closed_orbit    (6,) closed orbit vector
        dispersion      (4,) dispersion vector
        beta            [betax, betay] vector
        alpha           [alphax, alphay] vector
        mu              [mux, muy], betatron phase (modulo 2*pi)
        W               (2,) chromatic amplitude function (only if get_w==True)
        All values given at the entrance of each element specified in refpts.
        In case coupled == :py:obj:`True` additional outputs are available:
        gamma           gamma parameter of the transformation to eigenmodes
        A               (2, 2) matrix A in [3]
        B               (2, 2) matrix B in [3]
        C               (2, 2) matrix C in [3]
        Field values can be obtained with either
        lindata['idx']    or
        lindata.idx
    REFERENCES
        [1] D.Edwards,L.Teng IEEE Trans.Nucl.Sci. NS-20, No.3, p.885-888, 1973
        [2] E.Courant, H.Snyder
        [3] D.Sagan, D.Rubin Phys.Rev.Spec.Top.-Accelerators and beams,
            vol.2 (1999)
        [4] Brian W. Montague Report LEP Note 165, CERN, 1979

    :meta private:
    """
    analyze = _analyze4 if kwargs.pop('coupled', True) else _analyze2
    eld0, bd, eld = _linopt(ring, analyze, refpts, dp=dp, get_chrom=get_chrom,
                            add0=(0,), adds=(get_uint32_index(ring, refpts),),
                            addtype=[('idx', numpy.uint32)],
                            mname='m44', **kwargs)
    return eld0, bd.tune, bd.chromaticity, eld


# noinspection PyPep8Naming
@check_6d(False)
def avlinopt(ring, dp, refpts, **kwargs):
    r"""Linear analysis of a lattice with average values

    :py:func:`avlinopt` returns average beta, mu, dispersion over the lattice
    elements.

    Parameters:
        ring:       Lattice description.
        dp:         Momentum deviation.
        refpts:     Elements at which data is returned.
          It can be:

          1. an integer in the range [-len(ring), len(ring)-1]
             selecting the element according to python indexing rules.
             As a special case, len(ring) is allowed and refers to the end
             of the last element,
          2. an ordered list of such integers without duplicates,
          3. a numpy array of booleans of maximum length len(ring)+1,
             where selected elements are :py:obj:`True`.

    Keyword Args:
        dct (float):            Path lengthening. Defaults to :py:obj:`None`
        df (float):             Deviation of RF frequency. Defaults to
          :py:obj:`None`
        orbit (Orbit):          Avoids looking for the closed orbit if it is
          already known ((6,) array)
        get_chrom (bool):       Compute chromaticities. Needs computing
          the tune at 2 different momentum deviations around the central one.
        get_w (bool):           Computes chromatic amplitude functions
          (W) [4]_. Needs to compute the optics at 2 different momentum
          deviations around the central one.
        keep_lattice (bool):    Assume no lattice change since the
          previous tracking. Defaults to :py:obj:`False`
        XYStep (float):         Step size.
          Default: :py:data:`DConstant.XYStep <.DConstant>`
        DPStep (float):         Momentum step size.
          Default: :py:data:`DConstant.DPStep <.DConstant>`
        twiss_in:               Initial conditions for transfer line optics.
          Record array as output by :py:func:`.linopt6`, or dictionary. Keys:

          R or alpha, beta
            mandatory (2,) arrays
          closed_orbit
            Optional (6,) array, default 0
          dispersion
            Optional (6,) array, default 0

          If present, the attribute **R** will be used, otherwise the
          attributes **alpha** and **beta** will be used. All other attributes
          are ignored.

    Returns:
        elemdata:   Linear optics at the points refered to by *refpts*,
          if refpts is :py:obj:`None` an empty lindata structure is returned.
        avebeta:    Average beta functions
          [:math:`\hat{\beta_x},\hat{\beta_y}`] at *refpts*
        avemu:      Average phase advances
          [:math:`\hat{\mu_x},\hat{\mu_y}`] at *refpts*
        avedisp:    Average dispersion [:math:`\hat{\eta_x}, \hat{\eta'_x},
          \hat{\eta_y}, \hat{\eta'_y}`] at *refpts*
        avespos:    Average s position at *refpts*
        tune:       [:math:`\nu_1,\nu_2`], linear tunes for the two normal
          modes of linear motion [1]
        chrom:      [:math:`\xi_1,\xi_2`], chromaticities

    See also:
        :py:func:`linopt4`, :py:func:`get_optics`
    """
    def get_strength(elem):
        try:
            k = elem.PolynomB[1]
        except (AttributeError, IndexError):
            k = numpy.finfo(numpy.float64).eps
        return k

    def get_sext_strength(elem):
        try:
            k = elem.PolynomB[2]
        except (AttributeError, IndexError):
            k = 0.0
        return k

    def get_roll(elem):
        try:
            k = elem.R2[0][0]
        except (AttributeError, IndexError):
            k = 1.0
        return k

    def get_dx(elem):
        try:
            k = (elem.T2[0]-elem.T1[0])/2.0
        except (AttributeError, IndexError):
            k = 0.0
        return k

    def get_bendingangle(elem):
        try:
            k = elem.BendingAngle
        except (AttributeError, IndexError):
            k = 0.0
        return k

    def get_e1(elem):
        try:
            k = elem.EntranceAngle
        except (AttributeError, IndexError):
            k = 0.0
        return k

    def get_fint(elem):
        try:
            k = elem.FringeInt1
        except (AttributeError, IndexError):
            k = 0.0
        return k

    def get_gap(elem):
        try:
            k = elem.FullGap
        except (AttributeError, IndexError):
            k = 0.0
        return k

    def sini(x, L):
        r = x.copy()
        r[x > 0] = numpy.sin(numpy.sqrt(x[x > 0])*L[x > 0]
                             )/numpy.sqrt(x[x > 0])
        r[x < 0] = numpy.sinh(numpy.sqrt(-x[x < 0])*L[x < 0]
                              )/numpy.sqrt(-x[x < 0])
        return r

    def cosi(x, L):
        r = x.copy()
        r[x > 0] = numpy.cos(numpy.sqrt(x[x > 0])*L[x > 0])
        r[x < 0] = numpy.cosh(numpy.sqrt(-x[x < 0])*L[x < 0])
        return r

    def betadrift(beta0, alpha0, lg):
        gamma0 = (alpha0 * alpha0 + 1.0) / beta0
        return beta0-alpha0*lg+gamma0*lg*lg/3

    def betafoc(beta0, alpha0, kk, lg):
        gamma0 = (alpha0 * alpha0 + 1.0) / beta0
        return ((beta0+gamma0/kk)*lg +
                (beta0-gamma0/kk)*sini(kk, 2.0*lg)/2.0 +
                (cosi(kk, 2.0*lg)-1.0)*alpha0/kk)/lg/2.0

    def dispfoc(disp0, ir, k2, lg):
        avedisp = disp0.copy()
        avedisp[:, 0::2] = (disp0[:, 0::2]*(sini(k2, lg)) +
                            disp0[:, 1::2]*(1.0-cosi(k2, lg))/k2 +
                            ir*(lg-sini(k2, lg))/k2)/lg
        avedisp[:, 1::2] = (disp0[:, 0::2]*(sini(k2, lg)) -
                            disp0[:, 0::2]*(1.0-cosi(k2, lg)) +
                            ir*(lg-cosi(k2, lg))/k2)/lg
        return avedisp

    # selected list
    boolrefs = get_bool_index(ring, refpts)
    ring_selected = ring[refpts]
    L = numpy.array([el.Length for el in ring_selected])
    K = numpy.array([get_strength(el) for el in ring_selected])
    sext_strength = numpy.array([get_sext_strength(el)
                                 for el in ring_selected])
    roll = numpy.array([get_roll(el) for el in ring_selected])
    ba = numpy.array([get_bendingangle(el) for el in ring_selected])
    e1 = numpy.array([get_e1(el) for el in ring_selected])
    Fint = numpy.array([get_fint(el) for el in ring_selected])
    gap = numpy.array([get_gap(el) for el in ring_selected])
    dx = numpy.array([get_dx(el) for el in ring_selected])
    irho = ba.copy()
    d_csi = ba.copy()

    # whole ring list
    longelem = get_bool_index(ring, None)
    longelem[boolrefs] = (L != 0)

    shorti_refpts = (~longelem) & boolrefs
    longi_refpts = longelem & boolrefs
    longf_refpts = numpy.roll(longi_refpts, 1)
    all_refs = shorti_refpts | longi_refpts | longf_refpts
    _, bd, d_all = linopt4(ring, refpts=all_refs, dp=dp,
                           get_chrom=True, **kwargs)
    lindata = d_all[boolrefs[all_refs]]

    avebeta = lindata.beta.copy()
    avemu = lindata.mu.copy()
    avedisp = lindata.dispersion.copy()
    aves = lindata.s_pos.copy()
    ClosedOrbit = lindata.closed_orbit.copy()

    di = d_all[longi_refpts[all_refs]]
    df = d_all[longf_refpts[all_refs]]

    b_long = (L != 0.0)
    b_foc = (K != 0.0)
    b_foc_long = b_long & b_foc
    b_drf = b_long & (~b_foc)

    dx0 = ClosedOrbit[:, 0]
    dx0[b_long] = (di.closed_orbit[:, 0]+df.closed_orbit[:, 0])/2
    K = K*roll + 2*sext_strength*(dx0-dx)

    K2 = numpy.stack((K[b_foc_long], -K[b_foc_long]), axis=1)
    fff = b_foc_long[b_long]
    irho[b_long] = ba[b_long]/L[b_long]
    d_csi[b_long] = ba[b_long]*(gap[b_long]*Fint[b_long] *
                                (1.0 + numpy.sin(e1[b_long])**2) /
                                numpy.cos(e1[b_long])/L[b_long])
    L2 = numpy.stack((L[b_foc_long], L[b_foc_long]), axis=1)
    irho = irho.reshape((-1, 1))
    L = L.reshape((-1, 1))

    avemu[b_long] = 0.5 * (di.mu + df.mu)
    aves[b_long] = 0.5 * (df.s_pos + di.s_pos)
    avebeta[b_drf] = betadrift(di.beta[~fff], di.alpha[~fff], L[b_drf])
    avebeta[b_foc_long] = betafoc(di.beta[fff], di.alpha[fff], K2, L2)
    avedisp[numpy.ix_(b_long, [1, 3])] = (df.dispersion[:, [0, 2]] -
                                          di.dispersion[:, [0, 2]]) / L[b_long]
    idx = numpy.ix_(~fff, [0, 2])
    avedisp[numpy.ix_(b_drf, [0, 2])] = (di.dispersion[idx] +
                                         df.dispersion[idx]) * 0.5
    avedisp[b_foc_long, :] = dispfoc(di.dispersion[fff, :],
                                     irho[b_foc_long], K2, L2)

    return lindata, avebeta, avemu, avedisp, aves, bd.tune, bd.chromaticity


@frequency_control
def get_tune(ring: Lattice, *, method: str = 'linopt',
             dp: float = None, dct: float = None, df: float = None,
             orbit: Orbit = None, **kwargs):
    r"""Computes the tunes using several available methods

    Parameters:
        ring:                   Lattice description
        method:                 ``'linopt'`` returns the tunes from the
          :py:func:`.linopt6` function,

          ``'fft'`` tracks a single particle and computes the tunes with fft,

          ``'interp_fft'`` tracks a single particle and computes the tunes with
          interpolated FFT.
        dp (float):             Momentum deviation.
        dct (float):            Path lengthening.
        df (float):             Deviation of RF frequency.
        orbit (Orbit):          Avoids looking for the closed orbit if it is
          already known ((6,) array)

    for the ``'fft'`` and ``'interp_fft'`` methods only:

    Keyword Args:
        nturns (int):           Number of turns. Default: 512
        amplitude (float):      Amplitude of oscillation.
          Default: 1.E-6
        remove_dc (bool):       Remove the mean of oscillation data.
          Default: :py:obj:`True`
        num_harmonics (int):    Number of harmonic components to
          compute (before mask applied, default: 20)
        fmin (float):           Lower tune bound. Default: 0
        fmax (float):           Upper tune bound. Default: 1
        hann (bool):            Turn on Hanning window.
          Default: :py:obj:`False`. Work only for ``'fft'``
        get_integer(bool):   Turn on integer tune (slower)

    Returns:
        tunes (ndarray):        array([:math:`\nu_x,\nu_y`])
    """
    # noinspection PyShadowingNames
    def gen_centroid(ring, ampl, nturns, remove_dc, ld):
        nv = ld.A.shape[0]
        p0 = ld.closed_orbit.copy()
        p0[0] += ampl
        p0[2] += ampl
        if nv >= 6:
            p0[4] += ampl
        p1 = numpy.squeeze(internal_lpass(ring, p0, nturns, len(ring)))
        if remove_dc:
            p1 -= numpy.mean(p1, axis=1, keepdims=True)
        p2 = solve(ld.A, p1[:nv, :])
        return numpy.conjugate(p2.T.view(dtype=complex).T)
    get_integer = kwargs.pop('get_integer', False)
    if get_integer:
        assert method == 'linopt', \
           'Integer tune only accessible with method=linopt'
    if method == 'linopt':
        if get_integer:
            _, _, c = get_optics(ring, refpts=range(len(ring)+1),
                                 dp=dp, dct=dct, df=df, orbit=orbit)
            tunes = c.mu[-1]/(2*numpy.pi)
        else:
            tunes = _tunes(ring, dp=dp, dct=dct, df=df, orbit=orbit)
    else:
        nturns = kwargs.pop('nturns', 512)
        ampl = kwargs.pop('ampl', 1.0e-6)
        remove_dc = kwargs.pop('remove_dc', True)
        ld, _, _ = linopt6(ring, dp=dp, dct=dct, df=df, orbit=orbit)
        cents = gen_centroid(ring, ampl, nturns, remove_dc, ld)
        tunes = get_tunes_harmonic(cents, method=method, **kwargs)
    return tunes


@frequency_control
def get_chrom(ring: Lattice, *, method: str = 'linopt',
              dp: float = None, dct: float = None, df: float = None,
              cavpts: Refpts = None, **kwargs):
    r"""Computes the chromaticities using several available methods

    Parameters:
        ring:               Lattice description.
        method:             ``'linopt'`` returns the tunes from the
          :py:func:`linopt6` function,

          ``'fft'`` tracks a single particle and computes the tunes with
          :py:func:`~scipy.fftpack.fft`,

          ``'interp_fft'`` tracks a single particle and computes the tunes with
          interpolated FFT.
        dp (float):         Momentum deviation.
        dct (float):        Path lengthening.
        df (float):         Deviation of RF frequency.
        cavpts:     If :py:obj:`None`, look for ring.cavpts, or
          otherwise take all cavities.

    Keyword Args:
        DPStep (float):     Momentum step for differentiation
          Default: :py:data:`DConstant.DPStep <.DConstant>`

    for the ``'fft'`` and ``'interp_fft'`` methods only:

    Keyword Args:
        nturns (int):       Number of turns. Default: 512
        amplitude (float):  Amplitude of oscillation.
          Default: 1.E-6
        remove_dc (bool):   Remove the mean of oscillation data.
          Default: :py:obj:`True`
        num_harmonics (int):Number of harmonic components to
          compute (before mask applied, default: 20)
        fmin (float):       Lower tune bound. Default: 0
        fmax (float):       Upper tune bound. Default: 1
        hann (bool):        Turn on Hanning window.
          Default: :py:obj:`False`, Work only for ``'fft'``

    Returns:
        chromaticities (ndarray):   array([:math:`\xi_x,\xi_y`])
    """

    dp_step = kwargs.pop('DPStep', DConstant.DPStep)
    if method == 'fft':
        print('Warning fft method not accurate to get the ' +
              'chromaticity')

    if ring.is_6d:
        f0 = ring.get_rf_frequency(cavpts=cavpts)
        df = dp_step * ring.disable_6d(copy=True).slip_factor * f0
        rgup = ring.set_rf_frequency(f0 + 0.5 * df, cavpts=cavpts, copy=True)
        o0up, _ = find_orbit6(rgup, **kwargs)
        tune_up = get_tune(rgup,  method=method, orbit=o0up, **kwargs)
        rgdn = ring.set_rf_frequency(f0 - 0.5 * df, cavpts=cavpts, copy=True)
        o0dn, _ = find_orbit6(rgdn, **kwargs)
        tune_down = get_tune(rgdn,  method=method, orbit=o0dn, **kwargs)
        dp_step = o0up[4] - o0dn[4]
    else:
        if dct is not None or df is not None:
            dp = find_orbit4(ring, dct=dct, df=df)[0][4]
        elif dp is None:
            dp = 0.0
        tune_up = get_tune(ring, method=method, dp=dp + 0.5*dp_step, **kwargs)
        tune_down = get_tune(ring, method=method,
                             dp=dp - 0.5*dp_step, **kwargs)

    return (tune_up - tune_down) / dp_step


Lattice.linopt = linopt
Lattice.linopt2 = linopt2
Lattice.linopt4 = linopt4
Lattice.linopt6 = linopt6
Lattice.get_optics = get_optics
Lattice.avlinopt = avlinopt
Lattice.get_tune = get_tune
Lattice.get_chrom = get_chrom
