from math import sin, cos, atan2, pi
import numpy
from scipy.constants import c as clight
from scipy.linalg import solve
from at.physics import a_matrix, jmat, find_m44, find_m66

W_DTYPE6 = [('A', numpy.float64, (6, 6)),
            ('R', numpy.float64, (3, 6, 6)),
            ('M', numpy.float64, (6, 6)),
            ('alpha', numpy.float64, (2,)),
            ('beta', numpy.float64, (2,)),
            ('dispersion', numpy.float64, (4,)),
            ('mu', numpy.float64, (3,)),
            ]

W_DTYPE4 = [('A', numpy.float64, (4, 4)),
            ('R', numpy.float64, (2, 4, 4)),
            ('M', numpy.float64, (4, 4)),
            ('alpha', numpy.float64, (2,)),
            ('beta', numpy.float64, (2,)),
            ('mu', numpy.float64, (2,)),
            ]

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

    def propagate(ai, slcs):
        """Propagate the phase and mode matrices
        Rk = A * S * Ik * inv(A) * S.T
        """
        ais = ai.dot(s)
        invai = solve(ai, s.T)
        ri = numpy.array([numpy.dot(ais[:, s], invai[s, :]) for s in slcs])
        phi = numpy.array([get_phase(ai[slc, slc]) for slc in slcs])
        return ai, ri, phi

    nv = a0.shape[0]
    slices = [slice(2*i, 2*(i+1)) for i in range(nv // 2)]
    s = jmat(nv // 2)

    astd = standardize(a0, slices)
    _, r0, _ = propagate(astd, slices)
    return r0, (propagate(mi.dot(astd), slices) for mi in mstack)


def linopt6(ring, dp=None, refpts=None, orbit=None, twiss_in=None, **kwargs):
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
    OUTPUT
        elemdata0       linear optics data at the entrance/end of the ring
        beamdata        lattice properties
        elemdata        linear optics at the points refered to by refpts, if
                        refpts is None an empty elemdata structure is returned.

        elemdata is a record array with fields:
        R               R-matrices (3, 6, 6)
        M               Transfer matrix from the entrance of the line (6, 6)
        beta            [betax, betay] vector
        alpha           [alphax, alphay] vector
        dispersion      (4,) dispersion vector
        mu              [mux, muy], betatron phases
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
    def build_1turn_map(is6d):
        """Build (6,6) or (4,4) one-turn map depending of flag is6d"""
        if is6d:
            t, ts = find_m66(ring, refpts=refpts, orbit=orbit, **kwargs)
        else:
            t, ts = find_m44(ring, ddp=dp, refpts=refpts, orbit=orbit, **kwargs)
        return t, ts

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

    def output6(ms, mu, r123, a):
        """Extract output parameters from Bk matrices"""
        beta = numpy.array([r123[0, 0, 0], r123[1, 2, 2]])
        alpha = numpy.array([r123[0, 1, 0], r123[1, 3, 2]])
        dispersion = r123[2, :4, 4] / r123[2, 4, 4]
        return a, r123, ms, alpha, beta, dispersion, mu

    def output4(ms, mu, r12, a):
        """Extract output parameters from Bk matrices"""
        beta = numpy.array([r12[0, 0, 0], r12[1, 2, 2]])
        alpha = numpy.array([r12[0, 1, 0], r12[1, 3, 2]])
        return a, r12, ms, alpha, beta, mu

    def unwrap(mu):
        """Remove the phase jumps"""
        dmu = numpy.diff(numpy.concatenate((numpy.zeros((1, dms)), mu)), axis=0)
        jumps = dmu < 0
        mu += numpy.cumsum(jumps, axis=0) * 2.0 * numpy.pi

    if twiss_in is None:        # Circular machine
        mxx, mstack = build_1turn_map(ring.radiation)
        dms = mxx.shape[0] // 2
    else:                       # Transfer line
        if orbit is None:
            orbit = numpy.zeros((6,))
        sigma = build_sigma(twiss_in)
        dms = sigma.shape[0] // 2
        _, mstack = build_1turn_map(dms >= 3)
        mxx = sigma.dot(jmat(dms))

    if dms >= 3:
        dtype = W_DTYPE6
        output = output6
    else:
        dtype = W_DTYPE4
        output = output4

    a0, vps = a_matrix(mxx)
    length = ring.get_s_pos(len(ring))[0]
    tunes = numpy.mod(numpy.angle(vps) / 2.0 / pi, 1.0)
    damping_rates = -numpy.log(numpy.absolute(vps))
    damping_times = length / clight / damping_rates

    r0, phi_rr = _r_analysis(a0, mstack)
    elemdata0 = numpy.array(output(numpy.identity(2*dms), 0.0, r0, a0),
                            dtype=dtype).view(numpy.recarray)
    elemdata = numpy.array([output(ms, phi, ri, ai) for ms, (ai, ri, phi)
                            in zip(mstack, phi_rr)],
                           dtype=dtype).view(numpy.recarray)
    beamdata = numpy.array((tunes, damping_times),
                           dtype=[('tunes', numpy.float64, (dms,)),
                                  ('damping_times', numpy.float64, (dms,))
                                  ]).view(numpy.recarray)
    unwrap(elemdata.mu)
    return elemdata0, beamdata, elemdata
