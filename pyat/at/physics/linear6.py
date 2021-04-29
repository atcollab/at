from math import sin, cos, atan2, pi
from itertools import repeat
from warnings import warn
import numpy
from scipy.constants import c as clight
from scipy.linalg import solve, block_diag
from at.lattice import get_rf_frequency, set_rf_frequency, DConstant, get_s_pos
from at.lattice import AtWarning
from at.physics import a_matrix, jmat, get_mcf, find_m44, find_m66
from at.physics import find_orbit4, find_sync_orbit, find_orbit6
from at.tracking import lattice_pass

_S2 = numpy.array([[0, 1], [-1, 0]], dtype=numpy.float64)
_Tn = [_S2, _S2, _S2.T]

W_DTYPE6 = [('R', numpy.float64, (3, 6, 6)),
            ('A', numpy.float64, (6, 6)),
            ('alpha', numpy.float64, (2,)),
            ('beta', numpy.float64, (2,)),
            ('mu', numpy.float64, (3,)),
            ('dispersion', numpy.float64, (4,)),
            ('M', numpy.float64, (6, 6)),
            ('s_pos', numpy.float64),
            ('closed_orbit', numpy.float64, (6,)),
            ]

W_DTYPE4 = [('R', numpy.float64, (2, 4, 4)),
            ('A', numpy.float64, (4, 4)),
            ('alpha', numpy.float64, (2,)),
            ('beta', numpy.float64, (2,)),
            ('mu', numpy.float64, (2,)),
            ('dispersion', numpy.float64, (4,)),
            ('M', numpy.float64, (4, 4)),
            ('s_pos', numpy.float64),
            ('closed_orbit', numpy.float64, (6,)),
            ]

W_DTYPEW = [('W', numpy.float64, (2,))]

__all__ = ['linopt6']


def _r_analysis(a0, mstack):
    """Compute phase advance and R-matrices in 2D, 4D or 6D
    R0, R = _bk_analysis(A0)

    PARAMETERS
    A0      A-matrix at the beginning of the lattice
    mstack  transfer matrix from the entrance of the lattice to the selected
            element

    OUTPUT
    R0      R-matrices at the entrance of the lattice
    R       Iterator over the results at the selected elements
    """
    def get_phase(a22):
        """Return the phase for A standardization"""
        return atan2(a22[0, 1], a22[0, 0])

    def standardize(aa, slcs):
        """Apply rotation to put A in std form"""
        rotmat = numpy.zeros((nv, nv))
        for slc in slcs:
            rot = -get_phase(aa[slc, slc])
            cs = cos(rot)
            sn = sin(rot)
            rotmat[slc, slc] = numpy.array([[cs, sn], [-sn, cs]])
        return numpy.dot(aa, rotmat)

    def r_matrices(ai):
        # Rk = A * S * Ik * inv(A) * S.T
        ais = ai.dot(tt)            # Longitudinal swap
        invai = solve(ai, ss.T)
        ri = numpy.array([numpy.dot(ais[:, sl], invai[sl, :]) for sl in slices])
        return ri

    def propagate(ai):
        ri = r_matrices(ai)
        phi = numpy.array([get_phase(ai[sl, sl]) for sl in slices])
        return ri, ai, phi

    nv = a0.shape[0]
    dms = nv // 2
    slices = [slice(2*i, 2*(i+1)) for i in range(dms)]
    ss = jmat(dms)
    tt = block_diag(*_Tn[:dms])     # Used instead of ss for longitudinal swap

    astd = standardize(a0, slices)
    inival = (r_matrices(astd), astd, numpy.zeros((dms,)))
    return inival, (propagate(mi.dot(astd)) for mi in mstack)


def _find_orbit(ring, dp=None, refpts=None, orbit=None, ct=None, **kwargs):
    """"""
    if ring.radiation:
        if dp is not None:
            warn(AtWarning('In 6D, "dp" and "ct" are ignored'))
        if orbit is None:
            orbit, _ = find_orbit6(ring, **kwargs)
    else:
        if dp is None:
            dp = 0.0
        if orbit is None:
            if ct is not None:
                orbit, _ = find_sync_orbit(ring, ct, **kwargs)
            else:
                orbit, _ = find_orbit4(ring, dp, **kwargs)

    if refpts is None:
        orbs = []
    else:
        orbs = numpy.squeeze(
            lattice_pass(ring, orbit.copy(order='K'), refpts=refpts,
                         keep_lattice=True), axis=(1, 3)).T
    return orbit, orbs


def linopt6(ring, dp=None, refpts=None, orbit=None, cavpts=None, twiss_in=None,
            get_chrom=False, get_w=False, keep_lattice=False, **kwargs):
    """Perform linear analysis of a fully coupled lattice
    elemdata0, beamdata, elemdata = linopt6(lattice)

    For circular machines, linopt6 analyses
    the 4x4 1-turn transfer matrix if radiation is OFF, or
    the 6x6 1-turn transfer matrix if radiation is ON.

    For a transfer line, The "twiss_in" intput must contain either:
     - a field 'R', as provided by ATLINOPT6, or
      - the fields 'beta' and 'alpha', as provided by linopt and linopt6

    PARAMETERS
        ring            lattice description.
        refpts=None     elements at which data is returned.

    KEYWORDS
        orbit           avoids looking for the closed orbit if is already known
                        ((6,) array)
        keep_lattice    Assume no lattice change since the previous tracking.
                        Defaults to False
        XYStep=1.0e-8   transverse step for numerical computation
        DPStep=1.0E-6   momentum deviation used for computation of
                        the closed orbit
        twiss_in=None   Initial conditions for transfer line optics. Record
                        array, as output by linopt or linopt6.
                        If present, the attribute 'R' will be used, otherwise
                        The attributes 'alpha' and 'beta' will be used. All
                        other attributes are ignored.
        cavpts=None     Cavity location for off-momentum tuning

    OUTPUT
        elemdata0       linear optics data at the entrance/end of the ring
        beamdata        lattice properties
        elemdata        linear optics at the points refered to by refpts, if
                        refpts is None an empty elemdata structure is returned.

        elemdata is a record array with fields:
        R               R-matrices (3, 6, 6)
        A               A-matrix (6, 6)
        M               Transfer matrix from the entrance of the line (6, 6)
        beta            [betax, betay] vector
        alpha           [alphax, alphay] vector
        dispersion      (4,) dispersion vector
        mu              [mux, muy], betatron phases
        s_pos           longitudinal position [m]
        closed_orbit    (6,) closed orbit vector

        All values given at the entrance of each element specified in refpts.
        Field values can be obtained with either
        elemdata['beta']    or
        elemdata.beta

        beamdata is a record with fields:
        tunes           Fractional tunes
        damping_times   Damping times [s]

    REFERENCES
        [1] Etienne Forest, Phys. Rev. E 58, 2481 – Published 1 August 1998
        [2] Andrzej Wolski, Phys. Rev. ST Accel. Beams 9, 024001 –
            Published 3 February 2006
    """
    def get_alphabeta(r):
        beta = numpy.array([r[0, 0, 0], r[1, 2, 2]])
        alpha = numpy.array([r[0, 1, 0], r[1, 3, 2]])
        return alpha, beta

    # noinspection PyShadowingNames
    def _lin6(ring, dp, orbit, refpts=None, mxx=None, **kwargs):
        """"""
        if ring.radiation:
            mt, ms = find_m66(ring, refpts=refpts, orbit=orbit, **kwargs)
        else:
            mt, ms = find_m44(ring, dp, refpts=refpts, orbit=orbit, **kwargs)

        a0, vps = a_matrix(mt if mxx is None else mxx)
        val0, vals = _r_analysis(a0, ms)
        return ms, vps, val0, vals

    def build_sigma(twin):
        """Build the initial distribution at entrance of the transfer line"""
        try:
            sigm = numpy.sum(twin.R, axis=0)
        except AttributeError:
            slices = [slice(2 * i, 2 * (i + 1)) for i in range(4)]
            ab = numpy.stack((twin.alpha, twin.beta), axis=1)
            sigm = numpy.zeros((4, 4))
            for slc, (alpha, beta) in zip(slices, ab):
                gamma = (1.0+alpha*alpha)/beta
                sigm[slc, slc] = numpy.array([[beta, -alpha], [-alpha, gamma]])
        return sigm

    def output6(r123, a, mu, *args):
        """Extract output parameters from Bk matrices"""
        alpha, beta = get_alphabeta(r123)
        dispersion = r123[2, :4, 4] / r123[2, 4, 4]
        return (r123, a, alpha, beta, mu, dispersion) + args

    def output4(r12, a, *args):
        """Extract output parameters from Bk matrices"""
        alpha, beta = get_alphabeta(r12)
        return (r12, a, alpha, beta) + args

    def get_disp(ringup, ringdn, dpup, dpdn, refpts=None, matpts=None,
                 keep_lattice=False, **kwargs):

        def off_momentum(rng, dp):
            orb0, orbs = _find_orbit(ring, dp, refpts=refpts,
                                     keep_lattice=keep_lattice, **kwargs)
            dp = orb0[4]      # in 6D, dp comes out of find_orbit6
            _, vps, el0, els = _lin6(rng, dp, orb0, refpts=matpts,
                                     keep_lattice=True, **kwargs)
            tunes = numpy.mod(numpy.angle(vps) / 2.0 / pi, 1.0)
            return dp, tunes, orb0, orbs, el0, els

        def chromfunc(ddp, elup, eldn):
            aup, bup = get_alphabeta(elup[0])
            adn, bdn = get_alphabeta(eldn[0])
            db = (bup - bdn) / ddp
            mb = (bup + bdn) / 2
            da = (aup - adn) / ddp
            ma = (aup + adn) / 2
            w = numpy.sqrt((da - ma / mb * db) ** 2 + (db / mb) ** 2)
            return w

        dpup, tunesup, o0up, orbup, el0up, elup = off_momentum(ringup, dpup)
        dpdn, tunesdn, o0dn, orbdn, el0dn, eldn = off_momentum(ringdn, dpdn)
        deltap = dpup-dpdn
        chrom = (tunesup-tunesdn) / deltap
        disp0 = (o0up-o0dn)[:4]/deltap
        disp = ((oup-odn)[:4]/deltap for oup, odn in zip(orbup, orbdn))
        w0 = chromfunc(deltap, el0up, el0dn)
        w = (chromfunc(deltap, *v) for v in zip(elup, eldn))
        return chrom, disp0, disp, w0, w

    def unwrap(mu):
        """Remove the phase jumps"""
        dmu = numpy.diff(numpy.concatenate((numpy.zeros((1, dms)), mu)), axis=0)
        jumps = dmu < 0
        mu += numpy.cumsum(jumps, axis=0) * 2.0 * numpy.pi

    dp_step = kwargs.get('DPStep', DConstant.DPStep)

    if twiss_in is None:            # Circular machine
        mxx = None
    else:                           # Transfer line
        if orbit is None:
            orbit = numpy.zeros((6,))
        sigma = build_sigma(twiss_in)
        mxx = sigma.dot(jmat(sigma.shape[0] // 2))

    orb0, orbs = _find_orbit(ring, dp, refpts=refpts, orbit=orbit, **kwargs)
    dp = orb0[4]
    ms, vps, el0, els = _lin6(ring, dp, orb0, refpts=refpts, mxx=mxx,
                              keep_lattice=keep_lattice, **kwargs)

    dms = vps.size
    length = ring.get_s_pos(len(ring))[0]
    tunes = numpy.mod(numpy.angle(vps) / 2.0 / pi, 1.0)
    damping_rates = -numpy.log(numpy.absolute(vps))
    damping_times = length / clight / damping_rates

    spos = get_s_pos(ring, refpts)
    if dms >= 3:            # 6D processsing
        output = output6
        dtype = W_DTYPE6
        data0 = [*el0, numpy.identity(2*dms), 0.0, orb0]
        datas = [els, iter(ms), iter(spos), iter(orbs)]
        if get_chrom or get_w:
            f0 = get_rf_frequency(ring, cavpts=cavpts)
            df = -dp_step * get_mcf(ring.radiation_off(copy=True)) * f0
            rgup = set_rf_frequency(ring, f0 + 0.5*df, cavpts=cavpts, copy=True)
            rgdn = set_rf_frequency(ring, f0 - 0.5*df, cavpts=cavpts, copy=True)
            if get_w:
                dtype = W_DTYPE6 + W_DTYPEW
                chrom, _, _, w0, ws = get_disp(rgup, rgdn, None, None,
                                               matpts=refpts, **kwargs)
                data0 = data0.append(w0)
                datas.append(ws)
            else:
                chrom, _, _, _, _ = get_disp(rgup, rgdn, None, None)
        else:
            chrom = numpy.NaN
    else:               # 4D processsing
        output = output4
        dtype = W_DTYPE4
        data0 = [*el0]
        datas = [els]
        if get_w:
            dtype = W_DTYPE4 + W_DTYPEW
            chrom, d0, ds, w0, ws = get_disp(ring, ring, dp + 0.5*dp_step,
                                             dp - 0.5*dp_step, refpts=refpts,
                                             matpts=refpts,
                                             keep_lattice=True, **kwargs)
            data0 += [d0, numpy.identity(2*dms), 0.0, orb0, w0]
            datas += [ds, iter(ms), iter(spos), iter(orbs), ws]
        elif get_chrom:
            chrom, d0, ds, w0, ws = get_disp(ring, ring, dp + 0.5*dp_step,
                                             dp - 0.5*dp_step, refpts=refpts,
                                             keep_lattice=True, **kwargs)
            data0 += [d0, numpy.identity(6, 6), 0.0, orb0]
            datas += [ds, iter(ms), iter(spos), iter(orbs)]
        else:
            chrom = numpy.NaN
            data0 += [numpy.NaN, numpy.identity(6, 6), 0.0, orb0]
            datas += [repeat(numpy.NaN), iter(ms), iter(spos), iter(orbs)]

    elemdata0 = numpy.array(output(*data0), dtype=dtype).view(numpy.recarray)
    elemdata = numpy.fromiter((output(*el, *v) for el, *v in zip(*datas)),
                              dtype, count=ring.refcount(refpts)
                              ).view(numpy.recarray)

    beamdata = numpy.array((tunes, chrom, damping_times),
                           dtype=[('tunes', numpy.float64, (dms,)),
                                  ('chromaticities', numpy.float64, (dms,)),
                                  ('damping_times', numpy.float64, (dms,))
                                  ]).view(numpy.recarray)

    unwrap(elemdata.mu)
    return elemdata0, beamdata, elemdata
