"""
Radiation and equilibrium emittances
"""
import numpy
from numpy.linalg import multi_dot as md
from scipy.linalg import inv, det, solve_sylvester
from ..lattice import uint32_refpts
from ..tracking import lattice_pass
from ..physics import find_orbit6, find_m66, find_elem_m66, get_tunes_damp
# noinspection PyUnresolvedReferences
from .diffmatrix import find_mpole_raddiff_matrix

__all__ = ['ohmi_envelope']

_submat = [slice(0, 2), slice(2, 4), slice(6, 3, -1)]

# dtype for structured array containing optical parameters
ENVELOPE_DTYPE = [('r66', numpy.float64, (6, 6)),
                  ('r44', numpy.float64, (4, 4)),
                  ('m66', numpy.float64, (6, 6)),
                  ('orbit6', numpy.float64, (6,)),
                  ('emitXY', numpy.float64, (2,)),
                  ('emitXYZ', numpy.float64, (3,))]


def ohmi_envelope(ring, refpts=None, orbit=None, keep_lattice=False):
    """
    Calculate the equilibrium beam envelope in a
    circular accelerator using Ohmi's beam envelope formalism [1]

    emit0, mode_emit, damping_rates, tunes, emit = ohmi_envelope(ring[, refpts])

    PARAMETERS
        ring                at.Lattice object
        refpts              elements at which data is returned. It can be
                            1) an integer (0 indicating the first element)
                            2) a list of integers
                            3) a numpy array of booleans as long as ring where selected elements are true
                            Defaults to None

    KEYWORDS
        orbit=None          Avoids looking for the colsed orbit if is already known ((6,) array)
        keep_lattice=False  Assume no lattice change since the previous tracking.
                            Defaults to False

    OUTPUT
        emit0               emittance data at the start/end of the ring
        beamdata            beam parameters at the start of the ring
        emit                emittance data at the points refered to by refpts

        emit is a structured array with fields:
        R66                 (6, 6) equilibrium envelope matrix R
        R44                 (4, 4) betatron emittance matrix (dpp = 0)
        T66                 (6, 6) transfer matrix from the start of the ring
        orbit6              (6,) closed orbit
        emitXY              betatron emittance projected on xxp and yyp
        emitXYZ             6x6 emittance projected on xxp, yyp, ldp

        beamdata is a named tuple with attributes:
        tunes               tunes of the 3 normal modes
        damping_rates       damping rates of the 3 normal modes
        mode_matrices       R-matrices of the 3 normal modes
        mode_emittances     equilibrium emittances of the 3 normal modes

    REFERENCES
        [1] K.Ohmi et al. Phys.Rev.E. Vol.49. (1994)
    """

    def cumulb(it):
        """accumulate diffusion matrices"""
        cumul = numpy.zeros((6, 6))
        yield cumul
        for elem, orbin, b in it:
            m = find_elem_m66(elem, orbin)
            cumul = md((m, cumul, m.T)) + b
            yield cumul

    def process(r66):
        # projections on xx', zz', ldp
        emit3 = numpy.sqrt(numpy.array([det(r66[s, s]) for s in _submat]))
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
        emit2 = numpy.sqrt(numpy.array([det(r44[s, s], check_finite=False) for s in _submat[:2]]))
        return r44, emit2, emit3

    def propag(m, cumb, orbit6):
        """Propagate the beam matrix to refpts"""
        sigmatrix = md((m, rr, m.T)) + cumb
        m44, emit2, emit3 = process(sigmatrix)
        return sigmatrix, m44, m, orbit6, emit2, emit3

    nelems = len(ring)
    uint32refs = uint32_refpts(refpts, nelems)
    allrefs = uint32_refpts(range(nelems+1), nelems)

    if orbit is None:
        orbit, _ = find_orbit6(ring, keep_lattice=keep_lattice)
        keep_lattice = True

    orbs = numpy.squeeze(lattice_pass(ring, orbit.copy(order='K'), refpts=allrefs, keep_lattice=keep_lattice),
                         axis=(1, 3)).T
    mring, ms = find_m66(ring, uint32refs, orbit=orbit, keep_lattice=True)
    b0 = numpy.zeros((6, 6))
    bb = [find_mpole_raddiff_matrix(elem, orbit, ring.energy) if elem.PassMethod.endswith('RadPass') else b0
          for elem in ring]
    bbcum = numpy.stack(cumulb(zip(ring, orbs, bb)), axis=0)
    # ------------------------------------------------------------------------
    # Equation for the moment matrix R is
    #         R = MRING*R*MRING' + BCUM;
    # We rewrite it in the form of Lyapunov-Sylvester equation to use scipy's solve_sylvester function
    #            A*R + R*B = Q
    # where
    #               A =  inv(MRING)
    #               B = -MRING'
    #               Q = inv(MRING)*BCUM
    # -----------------------------------------------------------------------
    aa = inv(mring)
    bb = -mring.T
    qq = numpy.dot(aa, bbcum[-1])
    rr = solve_sylvester(aa, bb, qq)
    rr = 0.5 * (rr + rr.T)
    rr4, emitxy, emitxyz = process(rr)
    r66data = get_tunes_damp(mring, rr)

    data0 = numpy.array((rr, rr4, mring, orbit, emitxy, emitxyz), dtype=ENVELOPE_DTYPE)
    data = numpy.array(list(map(propag, ms, bbcum[uint32refs], orbs[uint32refs, :])), dtype=ENVELOPE_DTYPE)

    return data0, r66data, data
