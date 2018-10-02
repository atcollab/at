"""
Radiation and equilibrium emittances
"""
import numpy
from numpy.linalg import multi_dot as md
from scipy.linalg import inv, eig, solve_sylvester
from ..lattice import uint32_refpts
from ..tracking import lattice_pass
from .orbit import find_orbit6
from .matrix import find_m66, find_elem_m66
# noinspection PyUnresolvedReferences
from .diffmatrix import find_mpole_raddiff_matrix

__all__ = ['ohmi_envelope']

# dtype for structured array containing Twiss parameters
ENVELOPE_DTYPE = [('R', numpy.float64, (6, 6)),
                  ('tilt', numpy.float64),
                  ('sigma', numpy.float64, (2,))]


def ohmi_envelope(ring, radindex, refpts=None, orbit=None, keep_lattice=False):
    """
    Calculate the equilibrium beam envelope in a
    circular accelerator using Ohmi's beam envelope formalism [1]

    envelope, rmsdp, rmsbl = ohmi_envelope(ring, ,radelemindex, refpts)

    PARAMETERS
        ring            lattice description
        radelemindex    elements producing radiation (same format as refpts)
        refpts          elements at which data is returned. It can be
                        1) an integer (0 indicating the first element)
                        2) a list of integers
                        3) a numpy array of booleans as long as ring where
                           selected elements are true
                        Defaults to ring entrance ([0])

    OUTPUT
        envelope        envelope description at the entrance of each element specified in refpts.
        rmsdp           RMS momentum spread
        rmsbl           RMS bunch length [m]

        envelope is a structured array with fields:
        sigma           [sigma_a, sigma_b] - RMS size [m] along
                        the principal axis of a tilted ellips
                        Assuming normal distribution exp(-(z^2)/(2*sigma_z))
        tilt            Tilt angle of the XY ellipse [rad]
                        Positive Tilt corresponds to Corkscrew (right)
                        rotation of XY plane around s-axis
        R               (6, 6) equilibrium envelope matrix R

    REFERENCES
        [1] K.Ohmi et al. Phys.Rev.E. Vol.49. (1994)
    """
    def cumulb(it):
        """accumulate diffusion matrices"""
        bcum = numpy.zeros((6, 6))
        yield bcum
        for elem, orbin, b in it:
            m = find_elem_m66(elem, orbin)
            bcum = md((m, bcum, m.T)) + b
            yield bcum

    def propag(m, cumb):
        """Propagate the beam matrix to refpts"""
        sigmatrix = md((m, rr, m.T)) + cumb
        # noinspection PyTupleAssignmentBalance
        dr, u = eig(sigmatrix[0:3:2, 0:3:2])
        tilt = numpy.arcsin(0.5 * (u[1, 0] - u[0, 1]))
        sigma = numpy.sqrt(dr)
        return sigmatrix, tilt, sigma

    nelems = len(ring)
    uint32refs = uint32_refpts([0] if refpts is None else refpts, nelems)
    allrefs = uint32_refpts(range(nelems), nelems)

    if orbit is None:
        orbit = find_orbit6(ring, keep_lattice=keep_lattice)
        keep_lattice = True

    orb = numpy.rollaxis(numpy.squeeze(lattice_pass(ring, orbit.copy(order='K'), refpts=allrefs,
                                                    keep_lattice=keep_lattice)), -1)
    mring, ms = find_m66(ring, uint32refs, orbit=orbit, keep_lattice=True)
    bb = [numpy.zeros((6, 6))] * nelems
    for idx in radindex:
        bb[idx] = find_mpole_raddiff_matrix(ring[idx], numpy.squeeze(orb[idx, :]), 6.0e9)
    batbeg = numpy.stack(cumulb(zip(ring, orb, bb)), axis=0)
    # ------------------------------------------------------------------------
    # Equation for the moment matrix RR is
    #         RR = MRING*RR*MRING' + BCUM;
    # We rewrite it in the form of Lyapunov-Sylvester equation to use scipy's solve_sylvester function
    #            AA*RR + RR*BB = QQ
    # where
    #               AA =  inv(MRING)
    #               BB = -MRING'
    #               QQ = inv(MRING)*BCUM
    # -----------------------------------------------------------------------
    aa = inv(mring)
    bb = -mring.T
    qq = numpy.dot(aa, batbeg[-1])
    rr = solve_sylvester(aa, bb, qq)
    rmsdp = numpy.sqrt(rr[4, 4])
    rmsbl = numpy.sqrt(rr[5, 5])
    lindata = numpy.array(list(map(propag, numpy.rollaxis(ms, -1), batbeg[uint32refs])), dtype=ENVELOPE_DTYPE)
    return lindata, rmsdp, rmsbl
