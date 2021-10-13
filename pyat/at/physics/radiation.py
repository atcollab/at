"""
Radiation and equilibrium emittances
"""
from math import sin, cos, tan, sqrt, sinh, cosh, pi
import numpy
from scipy.linalg import inv, det, solve_sylvester
from at.lattice.constants import clight, e_mass
from at.lattice import Lattice, check_radiation
from at.lattice import Element, Dipole, Wiggler, DConstant, Multipole
from at.lattice import get_refpts, get_value_refpts
from at.lattice import uint32_refpts, set_value_refpts
from at.tracking import lattice_pass
from at.physics import find_orbit6, find_m66, find_elem_m66
from at.physics import find_mpole_raddiff_matrix, get_tunes_damp

__all__ = ['ohmi_envelope', 'get_radiation_integrals', 'quantdiffmat',
           'gen_quantdiff_elem', 'tapering']

_NSTEP = 60  # nb slices in a wiggler period

_submat = [slice(0, 2), slice(2, 4), slice(6, 3, -1)]

# dtype for structured array containing optical parameters
ENVELOPE_DTYPE = [('r66', numpy.float64, (6, 6)),
                  ('r44', numpy.float64, (4, 4)),
                  ('m66', numpy.float64, (6, 6)),
                  ('orbit6', numpy.float64, (6,)),
                  ('emitXY', numpy.float64, (2,)),
                  ('emitXYZ', numpy.float64, (3,))]


def _cumulb(it):
    """accumulate diffusion matrices"""
    cumul = numpy.zeros((6, 6))
    yield cumul
    for el, orbin, b in it:
        m = find_elem_m66(el, orbin)
        cumul = m.dot(cumul).dot(m.T) + b
        yield cumul


def _dmatr(ring, orbit=None, keep_lattice=False):
    """
    compute the cumulative diffusion and orbit
    matrices over the ring
    """
    nelems = len(ring)
    energy = ring.energy
    allrefs = uint32_refpts(range(nelems + 1), nelems)

    if orbit is None:
        orbit, _ = find_orbit6(ring, keep_lattice=keep_lattice)
        keep_lattice = True

    orbs = numpy.squeeze(
        lattice_pass(ring, orbit.copy(order='K'), refpts=allrefs,
                     keep_lattice=keep_lattice), axis=(1, 3)).T
    b0 = numpy.zeros((6, 6))
    bb = [find_mpole_raddiff_matrix(elem, elemorb, energy)
          if elem.PassMethod.endswith('RadPass') else b0
          for elem, elemorb in zip(ring, orbs)]
    bbcum = numpy.stack(list(_cumulb(zip(ring, orbs, bb))), axis=0)
    return bbcum, orbs


def _lmat(dmat):
    """
    lmat is Cholesky decomp of dmat unless diffusion is 0 in
    vertical.  Then do chol on 4x4 hor-long matrix and put 0's
    in vertical
    """
    lmat = numpy.zeros((6, 6))
    try:
        lmat = numpy.linalg.cholesky(dmat)
    except numpy.linalg.LinAlgError:
        nz = numpy.where(dmat != 0)
        cmat = numpy.reshape(dmat[nz], (4, 4))
        cmat = numpy.linalg.cholesky(cmat)
        lmat[nz] = numpy.reshape(cmat, (16,))
    return lmat


@check_radiation(True)
def ohmi_envelope(ring, refpts=None, orbit=None, keep_lattice=False):
    """
    Calculate the equilibrium beam envelope in a
    circular accelerator using Ohmi's beam envelope formalism [1]

    emit0, beamdata, emit = ohmi_envelope(ring[, refpts])

    PARAMETERS
        ring            Lattice object.
        refpts=None     elements at which data is returned. It can be:
                        1) an integer in the range [-len(ring), len(ring)-1]
                           selecting the element according to python indexing
                           rules. As a special case, len(ring) is allowed and
                           refers to the end of the last element,
                        2) an ordered list of such integers without duplicates,
                        3) a numpy array of booleans of maximum length
                           len(ring)+1, where selected elements are True.

    KEYWORDS
        orbit=None          Avoids looking for the closed orbit if it is
                            already known ((6,) array)
        keep_lattice=False  Assume no lattice change since the previous
                            tracking

    OUTPUT
        emit0               emittance data at the start/end of the ring
        beamdata            beam parameters at the start of the ring
        emit                emittance data at the points refered to by refpts,
                            if refpts is None an empty structure is returned.

        emit is a record array with fields:
        r66                 (6, 6) equilibrium envelope matrix R
        r44                 (4, 4) betatron emittance matrix (dpp = 0)
        m66                 (6, 6) transfer matrix from the start of the ring
        orbit6              (6,) closed orbit
        emitXY              (2,) betatron emittance projected on xxp and yyp
        emitXYZ             (3,) 6x6 emittance projected on xxp, yyp, ldp

        beamdata is a record array with fields:
        tunes               tunes of the 3 normal modes
        damping_rates       damping rates of the 3 normal modes
        mode_matrices       R-matrices of the 3 normal modes
        mode_emittances     equilibrium emittances of the 3 normal modes

        Field values can be obtained with either
        emit['r66']    or
        emit.r66

    REFERENCES
        [1] K.Ohmi et al. Phys.Rev.E. Vol.49. (1994)
    """

    def process(r66):
        # projections on xx', zz', ldp
        emit3sq = numpy.array([det(r66[s, s]) for s in _submat])
        # Prevent from unrealistic negative values of the determinant
        emit3 = numpy.sqrt(numpy.maximum(emit3sq, 0.0))
        # Emittance cut for dpp=0
        if emit3[0] < 1.E-13:  # No equilibrium emittance
            r44 = numpy.nan * numpy.ones((4, 4))
        elif emit3[1] < 1.E-13:  # Uncoupled machine
            minv = inv(r66[[0, 1, 4, 5], :][:, [0, 1, 4, 5]])
            r44 = numpy.zeros((4, 4))
            r44[:2, :2] = inv(minv[:2, :2])
        else:  # Coupled machine
            minv = inv(r66)
            r44 = inv(minv[:4, :4])
        # betatron emittances (dpp=0)
        emit2sq = numpy.array(
            [det(r44[s, s], check_finite=False) for s in _submat[:2]])
        # Prevent from unrealistic negative values of the determinant
        emit2 = numpy.sqrt(numpy.maximum(emit2sq, 0.0))
        return r44, emit2, emit3

    def propag(m, cumb, orbit6):
        """Propagate the beam matrix to refpts"""
        sigmatrix = m.dot(rr).dot(m.T) + cumb
        m44, emit2, emit3 = process(sigmatrix)
        return sigmatrix, m44, m, orbit6, emit2, emit3

    nelems = len(ring)
    uint32refs = uint32_refpts(refpts, nelems)
    bbcum, orbs = _dmatr(ring, orbit=orbit, keep_lattice=keep_lattice)
    mring, ms = find_m66(ring, uint32refs, orbit=orbs[0], keep_lattice=True)
    # ------------------------------------------------------------------------
    # Equation for the moment matrix R is
    #         R = MRING*R*MRING' + BCUM;
    # We rewrite it in the form of Lyapunov-Sylvester equation to use scipy's
    # solve_sylvester function
    #            A*R + R*B = Q
    # where
    #               A =  inv(MRING)
    #               B = -MRING'
    #               Q = inv(MRING)*BCUM
    # ------------------------------------------------------------------------
    aa = inv(mring)
    bb = -mring.T
    qq = numpy.dot(aa, bbcum[-1])
    rr = solve_sylvester(aa, bb, qq)
    rr = 0.5 * (rr + rr.T)
    rr4, emitxy, emitxyz = process(rr)
    r66data = get_tunes_damp(mring, rr)

    data0 = numpy.rec.fromarrays(
        (rr, rr4, mring, orbs[0], emitxy, emitxyz),
        dtype=ENVELOPE_DTYPE)
    if uint32refs.shape == (0,):
        data = numpy.recarray((0,), dtype=ENVELOPE_DTYPE)
    else:
        data = numpy.rec.fromrecords(
            list(map(propag, ms, bbcum[uint32refs], orbs[uint32refs, :])),
            dtype=ENVELOPE_DTYPE)

    return data0, r66data, data


def get_radiation_integrals(ring, dp=None, twiss=None, **kwargs):
    """
    Compute the 5 radiation integrals for uncoupled lattices.

    PARAMETERS
        ring            lattice description.

    KEYWORDS
        twiss=None      linear optics at all points (from linopt). If None,
                        it will be computed.
        dp=0.0          Ignored if radiation is ON. Momentum deviation.
        dct=None        Ignored if radiation is ON. Path lengthening.
                        If specified, dp is ignored and the off-momentum is
                        deduced from the path lengthening.
        method=linopt6  Method used for the analysis of the transfer matrix.
                        See get_optics.
                        linopt6: default
                        linopt2: faster if no longitudinal motion and
                                 no H/V coupling,

    OUTPUT
        i1, i2, i3, i4, i5
    """

    def dipole_radiation(elem, vini, vend):
        """Analytically compute the radiation integrals in dipoles"""
        beta0 = vini.beta[0]
        alpha0 = vini.alpha[0]
        eta0 = vini.dispersion[0]
        etap0 = vini.dispersion[1]

        ll = elem.Length
        theta = elem.BendingAngle
        rho = ll / theta
        rho2 = rho * rho
        k2 = elem.K + 1.0 / rho2
        eps1 = tan(elem.EntranceAngle) / rho
        eps2 = tan(elem.ExitAngle) / rho

        eta3 = vend.dispersion[0]
        alpha1 = alpha0 - beta0 * eps1
        gamma1 = (1.0 + alpha1 * alpha1) / beta0
        etap1 = etap0 + eta0 * eps1
        etap2 = vend.dispersion[1] - eta3 * eps2

        h0 = gamma1*eta0*eta0 + 2.0*alpha1*eta0*etap1 + beta0*etap1*etap1

        if k2 != 0.0:
            if k2 > 0.0:  # Focusing
                kl = ll * sqrt(k2)
                ss = sin(kl) / kl
                cc = cos(kl)
            else:  # Defocusing
                kl = ll * sqrt(-k2)
                ss = sinh(kl) / kl
                cc = cosh(kl)
            eta_ave = (theta - (etap2 - etap1)) / k2 / ll
            bb = 2.0 * (alpha1 * eta0 + beta0 * etap1) * rho
            aa = -2.0 * (alpha1 * etap1 + gamma1 * eta0) * rho
            h_ave = h0 + (aa * (1.0 - ss) + bb * (1.0 - cc) / ll
                          + gamma1 * (3.0 - 4.0 * ss + ss * cc) / 2.0 / k2
                          - alpha1 * (1.0 - cc) ** 2 / k2 / ll
                          + beta0 * (1.0 - ss * cc) / 2.0
                          ) / k2 / rho2
        else:
            eta_ave = 0.5 * (eta0 + eta3) - ll * ll / 12.0 / rho
            hp0 = 2.0 * (alpha1 * eta0 + beta0 * etap1) / rho
            h2p0 = 2.0 * (-alpha1 * etap1 + beta0 / rho - gamma1 * eta0) / rho
            h_ave = h0 + hp0 * ll / 2.0 + h2p0 * ll * ll / 6.0 \
                    - alpha1 * ll ** 3 / 4.0 / rho2 \
                    + gamma1 * ll ** 4 / 20.0 / rho2

        di1 = eta_ave * ll / rho
        di2 = ll / rho2
        di3 = ll / abs(rho) / rho2
        di4 = eta_ave * ll * (2.0 * elem.K + 1.0 / rho2) / rho \
              - (eta0 * eps1 + eta3 * eps2) / rho
        di5 = h_ave * ll / abs(rho) / rho2
        return numpy.array([di1, di2, di3, di4, di5])

    def wiggler_radiation(elem, dini):
        """Compute the radiation integrals in wigglers with the following
        approximations:

        - The wiggler is aligned with the closed orbit
        - The self-induced dispersion is neglected in I4 and I5, but is is used
          as a lower limit for the I5 contribution
        - I1, I2 are integrated analytically
        - I3 is integrated analytically for a single harmonic, numerically
          otherwise
        """

        def b_on_axis(wiggler, s):
            """On-axis wiggler field"""

            def harm(coef, h, phi):
                return -Bmax * coef * numpy.cos(h*kws + phi)

            kw = 2 * pi / wiggler.Lw
            Bmax = wiggler.Bmax
            kws = kw * s
            zz = [numpy.zeros(kws.shape)]
            vh = zz + [harm(pb[1], pb[4], pb[5]) for pb in wiggler.By.T]
            vv = zz + [-harm(pb[1], pb[4], pb[5]) for pb in wiggler.Bx.T]
            bys = numpy.sum(numpy.stack(vh), axis=0)
            bxs = numpy.sum(numpy.stack(vv), axis=0)
            return bxs, bys

        le = elem.Length
        alphax0 = dini.alpha[0]
        betax0 = dini.beta[0]
        gammax0 = (alphax0 * alphax0 + 1) / betax0
        eta0 = dini.dispersion[0]
        etap0 = dini.dispersion[1]
        H0 = gammax0*eta0*eta0 + 2*alphax0*eta0*etap0 + betax0*etap0*etap0
        avebetax = betax0 + alphax0*le + gammax0*le*le/3

        kw = 2 * pi / elem.Lw
        rhoinv = elem.Bmax / Brho
        coefh = elem.By[1, :] * rhoinv
        coefv = elem.Bx[1, :] * rhoinv
        coef2 = numpy.concatenate((coefh, coefv))
        if len(coef2) == 1:
            di3 = le * coef2[0] ** 3 * 4 / 3 / pi
        else:
            bx, bz = b_on_axis(elem, numpy.linspace(0, elem.Lw, _NSTEP + 1))
            rinv = numpy.sqrt(bx*bx + bz*bz) / Brho
            di3 = numpy.trapz(rinv ** 3) * le / _NSTEP
        di2 = le * (numpy.sum(coefh * coefh) + numpy.sum(coefv * coefv)) / 2
        di1 = -di2 / kw / kw
        di4 = 0
        if len(coefh) > 0:
            d5lim = 4 * avebetax * le * coefh[0] ** 5 / 15 / pi / kw / kw
        else:
            d5lim = 0
        di5 = max(H0 * di3, d5lim)
        return numpy.array([di1, di2, di3, di4, di5])

    Brho = sqrt(ring.energy**2 - e_mass**2) / clight
    integrals = numpy.zeros((5,))

    if twiss is None:
        _, _, twiss = ring.get_optics(refpts=range(len(ring) + 1), dp=dp,
                                      get_chrom=True, **kwargs)
    elif len(twiss) != len(ring) + 1:
        raise ValueError('length of Twiss data should be {0}'
                         .format(len(ring) + 1))
    for (el, vini, vend) in zip(ring, twiss[:-1], twiss[1:]):
        if isinstance(el, Dipole) and el.BendingAngle != 0.0:
            integrals += dipole_radiation(el, vini, vend)
        elif isinstance(el, Wiggler) and el.PassMethod != 'DriftPass':
            integrals += wiggler_radiation(el, vini)
    return tuple(integrals)


@check_radiation(True)
def quantdiffmat(ring, orbit=None):
    """
    This function computes the diffusion matrix of the whole ring

    PARAMETERS
        ring            lattice description.
        orbit=None      initial orbit

    OUTPUT
        diffusion matrix (6,6)
    """
    bbcum, _ = _dmatr(ring, orbit=orbit)
    diffmat = [(bbc + bbc.T) / 2 for bbc in bbcum]
    return numpy.round(diffmat[-1], 24)


@check_radiation(True)
def gen_quantdiff_elem(ring, orbit=None):
    """
    Generates a quantum diffusion element
    """
    dmat = quantdiffmat(ring, orbit=orbit)
    lmat = numpy.asfortranarray(_lmat(dmat))
    diff_elem = Element('Diffusion', Lmatp=lmat, PassMethod='QuantDiffPass')
    return diff_elem


@check_radiation(True)
def tapering(ring, multipoles=True, niter=1, **kwargs):
    """
    Scales magnet strength with local energy to cancel the closed orbit
    and optics errors due to synchrotron radiations. PolynomB is used for
    dipoles such that the machine geometry is maintained. This is the ideal
    tapering scheme where magnets and multipoles components (PolynomB and
    PolynomA) are scaled individually.
    !!! WARNING: This method works only for lattices without errors and
    corrections: if not all corrections and field errors will also be
    scaled !!!
    tapering(ring) or ring.tapering()
    PARAMETERS
        ring            lattice description.

    KEYWORDS
        multipoles=True scale all multipoles
        niter=1         number of iteration
        XYStep=1.0e-8   transverse step for numerical computation
        DPStep=1.0E-6   momentum deviation used for computation of orbit6
    """

    xy_step = kwargs.pop('XYStep', DConstant.XYStep)
    dp_step = kwargs.pop('DPStep', DConstant.DPStep)
    dipoles = get_refpts(ring, Dipole)
    b0 = get_value_refpts(ring, dipoles, 'BendingAngle')
    k0 = get_value_refpts(ring, dipoles, 'PolynomB', index=0)
    ld = get_value_refpts(ring, dipoles, 'Length')

    for i in range(niter):
        _, o6 = find_orbit6(ring, refpts=range(len(ring)+1),
                            XYStep=xy_step, DPStep=dp_step)
        dpps = (o6[dipoles, 4] + o6[dipoles+1, 4]) / 2
        set_value_refpts(ring, dipoles, 'PolynomB', b0/ld*dpps+k0*(1+dpps),
                         index=0)

    if multipoles:
        mults = get_refpts(ring, Multipole)
        k0 = get_value_refpts(ring, dipoles, 'PolynomB', index=0)
        _, o6 = find_orbit6(ring, refpts=range(len(ring)+1),
                            XYStep=xy_step, DPStep=dp_step)
        dpps = (o6[mults, 4] + o6[mults+1, 4]) / 2
        for dpp, el in zip(dpps, ring[mults]):
            el.PolynomB *= 1+dpp
            el.PolynomA *= 1+dpp
        set_value_refpts(ring, dipoles, 'PolynomB', k0, index=0)


Lattice.ohmi_envelope = ohmi_envelope
Lattice.get_radiation_integrals = get_radiation_integrals
Lattice.tapering = tapering
