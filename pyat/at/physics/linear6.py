from math import sin, cos, atan2, pi
import numpy
from at.lattice import check_radiation
from at.physics import a_matrix, find_m44, find_m66

W_DTYPE6 = [('R1', numpy.float64, (6, 6)),
            ('R2', numpy.float64, (6, 6)),
            ('R3', numpy.float64, (6, 6)),
            ('alpha', numpy.float64, (2,)),
            ('beta', numpy.float64, (2,)),
            ('dispersion', numpy.float64, (4,)),
            ('mu', numpy.float64, (3,)),
            ]

W_DTYPE4 = [('R1', numpy.float64, (4, 4)),
            ('R2', numpy.float64, (4, 4)),
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
        """Propagate the phase and mode matrices"""
        ri = [numpy.dot(ai[:, s], ai.T[s, :]) for s in slcs]
        phi = numpy.array([get_phase(ai[slc, slc]) for slc in slcs])
        return phi, ri

    nv = a0.shape[0]
    slices = [slice(2*i, 2*(i+1)) for i in range(nv // 2)]

    astd = standardize(a0, slices)
    _, R0 = propagate(astd, slices)
    return R0, (propagate(mi.dot(astd), slices) for mi in mstack)


def linopt6(ring, dp=None, refpts=None, twiss_in=None, **kwargs):
    """Perform linear analysis of a lattice
    elemdata0, beamdata, elemdata = linopt6(lattice)

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
        twiss_in=None   Initial twiss to compute transfer line optics of the
                        type lindata, the initial orbit in twiss_in is ignored,
                        only the beta and alpha are required other quatities
                        set to 0 if absent
    OUTPUT
        elemdata0       linear optics data at the entrance/end of the ring
        beamdata        lattice properties
        elemdata        linear optics at the points refered to by refpts, if
                        refpts is None an empty elemdata structure is returned.

        elemdata is a record array with fields:
        R1              R-matrix for mode 1 (~horizontal)
        R2              R-matrix for mode 2 (~vertical)
        beta            [betax, betay] vector
        alpha           [alphax, alphay] vector
        dispersion      (4,) dispersion vector
        mu              [mux, muy], betatron phases
        All values given at the entrance of each element specified in refpts.
        Field values can be obtained with either
        elemdata['beta']    or
        elemdata.beta

        beamdata is a record with fields:
        tunes           fractional tunes
        damping_rates   damping rates

    REFERENCES
        [1] Etienne Forest, Phys. Rev. E 58, 2481 – Published 1 August 1998
        [2] Andrzej Wolski, Phys. Rev. ST Accel. Beams 9, 024001 –
            Published 3 February 2006
    """
    def output6(mu, r1, r2, r3):
        """Extract output parameters from Bk matrices"""
        beta = numpy.array([r1[0, 0], r2[2, 2]])
        alpha = numpy.array([r1[1, 0], r2[3, 2]])
        dispersion = r3[:4, 4] / r3[4, 4]
        return r1, r2, r3, alpha, beta, dispersion, mu

    def output4(mu, r1, r2):
        """Extract output parameters from Bk matrices"""
        beta = numpy.array([r1[0, 0], r2[2, 2]])
        alpha = numpy.array([r1[1, 0], r2[3, 2]])
        return r1, r2, alpha, beta, mu

    def unwrap(mu):
        """Remove the phase jumps"""
        dmu = numpy.diff(numpy.concatenate((numpy.zeros((1, dms)), mu)), axis=0)
        jumps = dmu < 0
        mu += numpy.cumsum(jumps, axis=0) * 2.0 * numpy.pi

    if ring.radiation:
        mxx, mstack = find_m66(ring, refpts=refpts, **kwargs)
        dtype = W_DTYPE6
        output = output6
    else:
        mxx, mstack = find_m44(ring, ddp=dp, refpts=refpts, **kwargs)
        dtype = W_DTYPE4
        output = output4

    dms = mxx.shape[0] // 2
    a0, vps = a_matrix(mxx)
    tunes = numpy.mod(numpy.angle(vps) / 2.0 / pi, 1.0)
    damping_rates = -numpy.log(numpy.absolute(vps))

    r0, phi_rr = _r_analysis(a0, mstack)
    elemdata0 = numpy.array(output(0.0, *r0),
                            dtype=dtype).view(numpy.recarray)
    elemdata = numpy.array([output(phi, *ri) for phi, ri in phi_rr],
                           dtype=dtype).view(numpy.recarray)
    beamdata = numpy.array((tunes, damping_rates),
                           dtype=[('tunes', numpy.float64, (dms,)),
                                  ('damping_rates', numpy.float64, (dms,))
                                  ]).view(numpy.recarray)
    unwrap(elemdata.mu)
    return elemdata0, beamdata, elemdata
