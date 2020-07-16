"""
Radiation and equilibrium emittances
"""
from math import sin, cos, tan, sqrt, sinh, cosh, pi
import numpy
from scipy.linalg import inv, det, solve_sylvester
from at.lattice import Lattice, check_radiation, uint32_refpts
from at.lattice import Element, elements
from at.tracking import lattice_pass
from at.physics import find_orbit6, find_m66, find_elem_m66, get_tunes_damp
from at.physics import Cgamma, linopt, find_mpole_raddiff_matrix

__all__ = ['ohmi_envelope', 'get_radiation_integrals', 'quantdiffmat',
           'gen_quantdiff_elem']

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
    bb = [find_mpole_raddiff_matrix(elem, orbs[i], energy)
          if elem.PassMethod.endswith('RadPass') else b0
          for i, elem in enumerate(ring)]

    bbcum = numpy.stack(list(_cumulb(zip(ring, orbs, bb))), axis=0)
    return bbcum, orbs


def _lmat(dmat):
    '''
    lmat is Cholesky decomp of dmat unless diffusion is 0 in
    vertical.  Then do chol on 4x4 hor-long matrix and put 0's
    in vertical
    '''
    lmat = numpy.zeros((6, 6))
    try:
        lmat = numpy.linalg.cholesky(dmat)
    except numpy.linalg.LinAlgError:
        nz = numpy.where(dmat != 0)
        cmat = numpy.reshape(dmat[nz], (4, 4))
        cmat = numpy.linalg.cholesky(cmat)
        lmat[nz] = numpy.reshape(cmat, (16, ))
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
                            already known                           (6,) array)
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


@check_radiation(False)
def get_radiation_integrals(ring, dp=0.0, twiss=None):
    """
    Compute the 5 radiation integrals for uncoupled lattices. No RF cavity or
    radiating element is allowed.

    PARAMETERS
        ring            lattice description.
        dp=0.0          momentum deviation

    KEYWORDS
        twiss=None      linear optics at all points (from linopt). If None,
                        it will be computed.

    OUTPUT
        i1, i2, i3, i4, i5
    """
    i1 = 0.0
    i2 = 0.0
    i3 = 0.0
    i4 = 0.0
    i5 = 0.0

    if twiss is None:
        _, _, _, twiss = linopt(ring, dp, range(len(ring) + 1),
                                get_chrom=True, coupled=False)
    elif len(twiss) != len(ring)+1:
        raise ValueError('length of Twiss data should be {0}'
                         .format(len(ring)+1))
    for (elem, vini, vend) in zip(ring, twiss[:-1], twiss[1:]):
        if isinstance(elem, elements.Dipole) and elem.BendingAngle != 0.0:
            beta0 = vini.beta[0]
            alpha0 = vini.alpha[0]
            gamma0 = (1.0 + alpha0 * alpha0) / beta0
            eta0 = vini.dispersion[0]
            etap0 = vini.dispersion[1]

            ll = elem.Length
            theta = elem.BendingAngle
            rho = ll / theta
            rho2 = rho * rho
            k2 = elem.K + 1.0/rho2
            eps1 = tan(elem.EntranceAngle) / rho
            eps2 = tan(elem.ExitAngle) / rho

            eta3 = vend.dispersion[0]
            alpha1 = alpha0 - beta0*eps1
            etap1 = etap0 + eta0*eps1
            etap2 = vend.dispersion[1] - eta3*eps2

            h0 = gamma0*eta0*eta0 + 2.0*alpha1*eta0*etap1 + beta0*etap1*etap1

            if k2 != 0.0:
                if k2 > 0.0:        # Focusing
                    kl = ll * sqrt(k2)
                    ss = sin(kl) / kl
                    cc = cos(kl)
                else:               # Defocusing
                    kl = ll * sqrt(-k2)
                    ss = sinh(kl) / kl
                    cc = cosh(kl)
                eta_ave = (theta - (etap2 - etap1)) / k2 / ll
                bb = 2.0 * (alpha1*eta0 + beta0*etap1) * rho
                aa = -2.0 * (alpha1*etap1 + gamma0*eta0) * rho
                h_ave = h0 + (aa * (1.0-ss) + bb * (1.0-cc) / ll
                              + gamma0 * (3.0-4.0*ss+ss*cc) / 2.0 / k2
                              - alpha1 * (1.0-cc)**2 / k2 / ll
                              + beta0 * (1.0-ss*cc) / 2.0
                              ) / k2 / rho2
            else:
                eta_ave = 0.5 * (eta0 + eta3) - ll*ll / 12.0 / rho
                hp0 = 2.0 * (alpha1 * eta0 + beta0 * etap1) / rho
                h2p0 = 2.0 * (-alpha1*etap1 + beta0/rho - gamma0*eta0) / rho
                h_ave = h0 + hp0*ll/2.0 + h2p0*ll*ll/6.0 \
                    - alpha1*ll**3/4.0/rho2 \
                    + gamma0*ll**4/20.0/rho2

            i1 += eta_ave * ll / rho
            i2 += ll / rho2
            i3 += ll / abs(rho) / rho2
            i4 += eta_ave * ll * (2.0*elem.K+1.0/rho2) / rho \
                - (eta0*eps1 + eta3*eps2)/rho
            i5 += h_ave * ll / abs(rho) / rho2

    return i1, i2, i3, i4, i5


def get_energy_loss(ring):
    """Energy loss per turn [eV]

    Losses = Cgamma / 2pi * EGeV^4 * i2
    """
    lenthe = numpy.array(
        [(elem.Length, elem.BendingAngle) for elem in ring if
         isinstance(elem, elements.Dipole)])
    lendp = lenthe[:, 0]
    theta = lenthe[:, 1]

    i2 = ring.periodicity * (numpy.sum(theta * theta / lendp))
    e_loss = Cgamma / 2.0 / pi * ring.energy**4 * i2
    return e_loss


@check_radiation(True)
def quantdiffmat(ring, orbit=None):
    '''
    This function computes the diffusion matrix of the whole ring

    PARAMETERS
        ring            lattice description.
        orbit=None      initial orbit

    OUTPUT
        diffusion matrix (6,6)
    '''
    bbcum, _ = _dmatr(ring, orbit=orbit)
    diffmat = [(bbc + bbc.T)/2 for bbc in bbcum]
    return numpy.round(diffmat[-1], 24)


@check_radiation(True)
def gen_quantdiff_elem(ring, orbit=None):
    '''
    Generates a quantum diffusion element
    '''
    dmat = quantdiffmat(ring, orbit=orbit)
    lmat = numpy.asfortranarray(_lmat(dmat))
    diff_elem = Element('Diffusion',
                        Lmatp=lmat,
                        PassMethod='QuantDiffPass')
    return diff_elem


Lattice.ohmi_envelope = ohmi_envelope
Lattice.get_radiation_integrals = get_radiation_integrals
Lattice.energy_loss = property(get_energy_loss)
