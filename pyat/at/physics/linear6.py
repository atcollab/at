from math import sin, cos, atan2
import numpy
from at.lattice import check_radiation
from at.physics import amat, find_m66

W_DTYPE = [('B1', numpy.float64, (6, 6)),
           ('B2', numpy.float64, (6, 6)),
           ('alpha', numpy.float64, (2,)),
           ('beta', numpy.float64, (2,)),
           ('dispersion', numpy.float64, (4,)),
           ('mu', numpy.float64, (2,)),
           ]

__all__ = ['bk_analysis', 'linopt6']


def bk_analysis(m0, mstack):
    """Compute phase advance and BK matrices in 2D, 4D or 6D
    """
    def make_Tk(rng):
        """Build the Tk selection matrices"""
        tk = numpy.zeros((nv, nv))
        tk[rng, rng] = numpy.identity(2)
        return tk

    def get_phase(a22):
        """Return the phase for A standardization"""
        return atan2(a22[0, 1], a22[0, 0])

    def standardize(aa, slcs):
        """Apply rotation to set A in std form"""
        rotmat = numpy.zeros((nv, nv))
        for slc in slcs:
            r1 = -get_phase(aa[slc, slc])
            cs = cos(r1)
            sn = sin(r1)
            rotmat[slc, slc] = numpy.array([[cs, sn], [-sn, cs]])
        return numpy.dot(aa, rotmat)

    def propagate(ai, tks, slcs):
        """Propagate the A matrices"""
        bk = [ai.dot(tk.dot(ai.T)) for tk in tks]
        phi = numpy.array([get_phase(ai[slc, slc]) for slc in slcs])
        return phi, bk

    nv = m0.shape[0]
    slices = [slice(2*i, 2*(i+1)) for i in range(nv // 2)]
    Tk = [make_Tk(slc) for slc in slices]

    astd = standardize(amat(m0), slices)
    return (propagate(mi.dot(astd), Tk, slices) for mi in mstack)


@check_radiation(True)
def linopt6(ring, refpts=None, twiss_in=None, **kwargs):
    """"""
    def output(mu, b1, b2, b3):
        """Extract output parameters from Bk matrices"""
        beta = numpy.array([b1[0, 0], b2[2, 2]])
        alpha = numpy.array([b1[1, 0], b2[3, 2]])
        dispersion = b3[:4, 4]/b3[4, 4]
        return b1, b2, alpha, beta, dispersion, mu[:2]

    def unwrap(mu):
        """Remove the phase jumps"""
        dmu = numpy.diff(numpy.concatenate((numpy.zeros((1, 2)), mu)), axis=0)
        jumps = dmu < 0
        mu += numpy.cumsum(jumps, axis=0) * 2.0 * numpy.pi

    mxx, mstack = find_m66(ring, refpts=refpts, **kwargs)
    phi_bk = bk_analysis(mxx, mstack)
    reca = numpy.array([output(phi, *bk) for phi, bk in phi_bk],
                       dtype=W_DTYPE).view(numpy.recarray)
    unwrap(reca.mu)
    return reca
