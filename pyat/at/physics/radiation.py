"""
Radiation and equilibrium emittances
"""

from __future__ import annotations

__all__ = [
    "ohmi_envelope",
    "get_radiation_integrals",
    "quantdiffmat",
    "gen_quantdiff_elem",
    "tapering",
]

from math import sin, cos, tan, sqrt, sinh, cosh, pi

import numpy as np
from scipy.linalg import inv, det, solve_sylvester

from ..lattice import Dipole, Wiggler, EnergyLoss, DConstant, test_mode
from ..lattice import Lattice, Element, check_radiation, Refpts, All
from ..lattice import Quadrupole, Multipole, QuantumDiffusion
from ..lattice import Collective, SimpleQuantDiff
from ..lattice import frequency_control, set_value_refpts
from ..constants import Cgamma
from . import ELossMethod
from . import find_mpole_raddiff_matrix, FDW, get_tunes_damp
from . import find_orbit6, find_m66, find_elem_m66, Orbit
from ..tracking import internal_lpass, diffusion_matrix

_new_methods = {
    "BndMPoleSymplectic4RadPass",
    "StrMPoleSymplectic4RadPass",
    "ExactMultipoleRadPass",
    "GWigSymplecticRadPass",
    "EnergyLossRadPass",
}

_NSTEP = 60  # nb slices in a wiggler period

_submat = [slice(0, 2), slice(2, 4), slice(6, 3, -1)]

# dtype for structured array containing optical parameters
ENVELOPE_DTYPE = [
    ("r66", np.float64, (6, 6)),
    ("r44", np.float64, (4, 4)),
    ("m66", np.float64, (6, 6)),
    ("orbit6", np.float64, (6,)),
    ("emitXY", np.float64, (2,)),
    ("emitXYZ", np.float64, (3,)),
]

_b0 = np.zeros((6, 6), dtype=np.float64)


def _dmatr(ring: Lattice, orbit: Orbit = None, keep_lattice: bool = False):
    """
    compute the cumulative diffusion and orbit
    matrices over the ring
    """

    def _cumulb(it):
        """accumulate diffusion matrices"""
        cumul = np.zeros((6, 6))
        yield cumul
        for el, orbin, b in it:
            m = find_elem_m66(el, orbin, energy=energy, particle=ring.particle)
            cumul = m.dot(cumul).dot(m.T) + b
            yield cumul

    def substitute(elem):
        if elem.PassMethod not in _new_methods:
            elem = elem.copy()
            elem.PassMethod = "StrMPoleSymplectic4RadPass"
        return elem

    def elem_diffusion(elem: Element, elemorb):
        if elem.PassMethod.endswith("RadPass"):
            if not test_mode():
                return diffusion_matrix(substitute(elem), elemorb, energy=energy)
            elif hasattr(elem, "Bmax"):
                return FDW(elem, orbit, energy)
            else:
                return find_mpole_raddiff_matrix(elem, elemorb, energy)
        else:
            return _b0

    energy = ring.energy

    if orbit is None:
        orbit, _ = find_orbit6(ring, keep_lattice=keep_lattice)
        keep_lattice = True

    orbs = np.squeeze(
        internal_lpass(
            ring, orbit.copy(order="K"), refpts=All, keep_lattice=keep_lattice
        ),
        axis=(1, 3),
    ).T

    bb = [elem_diffusion(elem, orb) for elem, orb in zip(ring, orbs)]

    bbcum = np.stack(list(_cumulb(zip(ring, orbs, bb))), axis=0)
    return bbcum, orbs


def _lmat(dmat):
    """
    lmat is Cholesky decomp of dmat unless diffusion is 0 in
    vertical.  Then do chol on 4x4 hor-long matrix and put 0's
    in vertical
    """
    lmat = np.zeros((6, 6))
    try:
        lmat = np.linalg.cholesky(dmat)
    except np.linalg.LinAlgError:
        nz = np.where(dmat != 0)
        cmat = np.reshape(dmat[nz], (4, 4))
        cmat = np.linalg.cholesky(cmat)
        lmat[nz] = np.reshape(cmat, (16,))
    return lmat


@check_radiation(True)
def ohmi_envelope(
    ring: Lattice,
    refpts: Refpts = None,
    orbit: Orbit = None,
    keep_lattice: bool = False,
):
    """Calculates the equilibrium beam envelope

    Computation based on Ohmi's beam envelope formalism [1]_

    Parameters:
        ring:           Lattice description. Radiation must be ON
        refpts:         Observation points
        orbit:          Avoids looking for initial the closed orbit if it is
          already known ((6,) array).
        keep_lattice:   Assume no lattice change since the previous tracking.
          Default: False

    Returns:
        emit0 (np.recarray):     Emittance data at the start/end of the ring
        beamdata (np.recarray):  Beam parameters at the start of the ring
        emit (np.recarray):      Emittance data at the points selected to by
          ``refpts``

    **emit** is a :py:obj:`record array <numpy.recarray>` with the following
    fields:

    ================    ===================================================
    **r66**             (6, 6) equilibrium envelope matrix R
    **r44**             (4, 4) betatron emittance matrix (dpp = 0)
    **m66**             (6, 6) transfer matrix from the start of the ring
    **orbit6**          (6,) closed orbit
    **emitXY**          (2,) betatron emittance projected on xxp and yyp
    **emitXYZ**         (3,) 6x6 emittance projected on xxp, yyp, ldp
    ================    ===================================================

    Values given at the entrance of each element specified in ``refpts``.

    Field values can be obtained with either
    ``emit['r66']`` or ``emit.r66``

    **beamdata** is a :py:obj:`record array <numpy.recarray>` with the
    following fields:

    ====================  ===================================================
    **tunes**             tunes of the 3 normal modes
    **damping_rates**     damping rates of the 3 normal modes
    **mode_matrices**     R-matrices of the 3 normal modes
    **mode_emittances**   equilibrium emittances of the 3 normal modes
    ====================  ===================================================

    References:
        .. [1] K.Ohmi et al. Phys.Rev.E. Vol.49. (1994)
    """

    def process(r66):
        # projections on xx', zz', ldp
        emit3sq = np.array([det(r66[s, s]) for s in _submat])
        # Prevent from unrealistic negative values of the determinant
        emit3 = np.sqrt(np.maximum(emit3sq, 0.0))
        # Emittance cut for dpp=0
        if emit3[0] < 1.0e-13:  # No equilibrium emittance
            r44 = np.nan * np.ones((4, 4))
        elif emit3[1] < 1.0e-13:  # Uncoupled machine
            minv = inv(r66[[0, 1, 4, 5], :][:, [0, 1, 4, 5]])
            r44 = np.zeros((4, 4))
            r44[:2, :2] = inv(minv[:2, :2])
        else:  # Coupled machine
            minv = inv(r66)
            r44 = inv(minv[:4, :4])
        # betatron emittances (dpp=0)
        emit2sq = np.array([det(r44[s, s], check_finite=False) for s in _submat[:2]])
        # Prevent from unrealistic negative values of the determinant
        emit2 = np.sqrt(np.maximum(emit2sq, 0.0))
        return r44, emit2, emit3

    def propag(m, cumb, orbit6):
        """Propagate the beam matrix to refpts"""
        sigmatrix = m.dot(rr).dot(m.T) + cumb
        m44, emit2, emit3 = process(sigmatrix)
        return sigmatrix, m44, m, orbit6, emit2, emit3

    rtmp = ring.disable_6d(QuantumDiffusion, Collective, SimpleQuantDiff, copy=True)
    uint32refs = rtmp.get_uint32_index(refpts)
    bbcum, orbs = _dmatr(rtmp, orbit=orbit, keep_lattice=keep_lattice)
    mring, ms = find_m66(rtmp, uint32refs, orbit=orbs[0], keep_lattice=True)
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
    qq = np.dot(aa, bbcum[-1])
    rr = solve_sylvester(aa, bb, qq)
    rr = 0.5 * (rr + rr.T)
    rr4, emitxy, emitxyz = process(rr)
    r66data = get_tunes_damp(mring, rr)

    data0 = np.rec.fromarrays(
        (rr, rr4, mring, orbs[0], emitxy, emitxyz), dtype=ENVELOPE_DTYPE
    )
    if uint32refs.shape == (0,):
        data = np.recarray((0,), dtype=ENVELOPE_DTYPE)
    else:
        data = np.rec.fromrecords(
            list(map(propag, ms, bbcum[uint32refs], orbs[uint32refs, :])),
            dtype=ENVELOPE_DTYPE,
        )

    return data0, r66data, data


@frequency_control
def get_radiation_integrals(
    ring, dp: float = None, twiss=None, **kwargs
) -> tuple[float, float, float, float, float]:
    r"""Computes the 5 radiation integrals for uncoupled lattices.

    Parameters:
        ring:   Lattice description
        twiss:  Linear optics at all points (from :py:func:`.linopt6`).
          If ``None``, it will be computed.

    Keyword Args:
        dp (float):             Momentum deviation. Defaults to :py:obj:`None`
        dct (float):            Path lengthening. Defaults to :py:obj:`None`
        df (float):             Deviation of RF frequency. Defaults to
        method (Callable):  Method for linear optics:

          :py:obj:`~.linear.linopt2`: no longitudinal motion, no H/V coupling,

          :py:obj:`~.linear.linopt6` (default): with or without longitudinal
          motion, normal mode analysis

    Returns:
        i1 (float): Radiation integrals - :math:`I_1 \quad [m]`
        i2 (float): :math:`I_2 \quad [m^{-1}]`
        i3 (float): :math:`I_3 \quad [m^{-2}]`
        i4 (float): :math:`I_4 \quad [m^{-1}]`
        i5 (float): :math:`I_5 \quad [m^{-1}]`
    """

    def element_radiation(elem: Dipole | Quadrupole, vini, vend):
        """Analytically compute the radiation integrals in dipoles"""
        beta0 = vini.beta[0]
        alpha0 = vini.alpha[0]
        eta0 = vini.dispersion[0]
        etap0 = vini.dispersion[1]
        theta = getattr(elem, "BendingAngle", None)
        if theta is None:
            xpi = vini.closed_orbit[1] / (1.0 + vini.closed_orbit[4])
            xpo = vend.closed_orbit[1] / (1.0 + vend.closed_orbit[4])
            theta = xpi - xpo
        if abs(theta) < 1.0e-7:
            return np.zeros(5)
        angin = getattr(elem, "EntranceAngle", 0.0)
        angout = getattr(elem, "ExitAngle", 0.0)

        ll = elem.Length
        rho = ll / theta
        rho2 = rho * rho
        k2 = elem.K + 1.0 / rho2
        eps1 = tan(angin) / rho
        eps2 = tan(angout) / rho
        eta3 = vend.dispersion[0]
        alpha1 = alpha0 - beta0 * eps1
        gamma1 = (1.0 + alpha1 * alpha1) / beta0
        etap1 = etap0 + eta0 * eps1
        etap2 = vend.dispersion[1] - eta3 * eps2

        h0 = gamma1 * eta0 * eta0 + 2.0 * alpha1 * eta0 * etap1 + beta0 * etap1 * etap1

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
            h_ave = (
                h0
                + (
                    aa * (1.0 - ss)
                    + bb * (1.0 - cc) / ll
                    + gamma1 * (3.0 - 4.0 * ss + ss * cc) / 2.0 / k2
                    - alpha1 * (1.0 - cc) ** 2 / k2 / ll
                    + beta0 * (1.0 - ss * cc) / 2.0
                )
                / k2
                / rho2
            )
        else:
            eta_ave = 0.5 * (eta0 + eta3) - ll * ll / 12.0 / rho
            hp0 = 2.0 * (alpha1 * eta0 + beta0 * etap1) / rho
            h2p0 = 2.0 * (-alpha1 * etap1 + beta0 / rho - gamma1 * eta0) / rho
            h_ave = (
                h0
                + hp0 * ll / 2.0
                + h2p0 * ll * ll / 6.0
                - alpha1 * ll**3 / 4.0 / rho2
                + gamma1 * ll**4 / 20.0 / rho2
            )

        di1 = eta_ave * ll / rho
        di2 = ll / rho2
        di3 = ll / abs(rho) / rho2
        di4 = (
            eta_ave * ll * (2.0 * elem.K + 1.0 / rho2) / rho
            - (eta0 * eps1 + eta3 * eps2) / rho
        )
        di5 = h_ave * ll / abs(rho) / rho2
        return np.array([di1, di2, di3, di4, di5])

    def wiggler_radiation(elem: Wiggler, dini):
        """Compute the radiation integrals in wigglers with the following
        approximations:

        - The wiggler is aligned with the closed orbit
        - The self-induced dispersion is neglected in I4 and I5, but is is used
          as a lower limit for the I5 contribution
        - I1, I2 are integrated analytically
        - I3 is integrated analytically for a single harmonic, numerically
          otherwise
        """

        def b_on_axis(wiggler: Wiggler, s):
            """On-axis wiggler field"""

            def harm(coef, h, phi):
                return -Bmax * coef * np.cos(h * kws + phi)

            kw = 2 * pi / wiggler.Lw
            Bmax = wiggler.Bmax
            kws = kw * s
            zz = [np.zeros(kws.shape)]
            vh = zz + [harm(pb[1], pb[4], pb[5]) for pb in wiggler.By.T]
            vv = zz + [-harm(pb[1], pb[4], pb[5]) for pb in wiggler.Bx.T]
            bys = np.sum(np.stack(vh), axis=0)
            bxs = np.sum(np.stack(vv), axis=0)
            return bxs, bys

        le = elem.Length
        alphax0 = dini.alpha[0]
        betax0 = dini.beta[0]
        gammax0 = (alphax0 * alphax0 + 1) / betax0
        eta0 = dini.dispersion[0]
        etap0 = dini.dispersion[1]
        H0 = gammax0 * eta0 * eta0 + 2 * alphax0 * eta0 * etap0 + betax0 * etap0 * etap0
        avebetax = betax0 + alphax0 * le + gammax0 * le * le / 3

        kw = 2 * pi / elem.Lw
        rhoinv = elem.Bmax / Brho
        coefh = elem.By[1, :] * rhoinv
        coefv = elem.Bx[1, :] * rhoinv
        coef2 = np.concatenate((coefh, coefv))
        if len(coef2) == 1:
            di3 = le * coef2[0] ** 3 * 4 / 3 / pi
        else:
            bx, bz = b_on_axis(elem, np.linspace(0, elem.Lw, _NSTEP + 1))
            rinv = np.sqrt(bx * bx + bz * bz) / Brho
            di3 = np.trapz(rinv**3) * le / _NSTEP
        di2 = le * (np.sum(coefh * coefh) + np.sum(coefv * coefv)) / 2
        di1 = -di2 / kw / kw
        di4 = 0
        if len(coefh) > 0:
            d5lim = 4 * avebetax * le * coefh[0] ** 5 / 15 / pi / kw / kw
        else:
            d5lim = 0
        di5 = max(H0 * di3, d5lim)
        return np.array([di1, di2, di3, di4, di5])

    def eloss_radiation(elem: EnergyLoss, coef):
        # Assuming no diffusion
        di2 = elem.EnergyLoss / coef
        return np.array([0.0, di2, 0.0, 0.0, 0.0])

    Brho = ring.BRho
    elosscoef = Cgamma / 2.0 / pi * ring.energy**4
    integrals = np.zeros((5,))

    if twiss is None:
        _, _, twiss = ring.get_optics(
            refpts=range(len(ring) + 1), dp=dp, get_chrom=True, **kwargs
        )
    elif len(twiss) != len(ring) + 1:
        raise ValueError(f"length of Twiss data should be {len(ring) + 1}")
    for el, vini, vend in zip(ring, twiss[:-1], twiss[1:]):
        if isinstance(el, (Dipole, Quadrupole)):
            integrals += element_radiation(el, vini, vend)
        elif isinstance(el, Wiggler) and el.PassMethod != "DriftPass":
            integrals += wiggler_radiation(el, vini)
        elif isinstance(el, EnergyLoss):
            integrals += eloss_radiation(el, elosscoef)
    return tuple(integrals)


@check_radiation(True)
def quantdiffmat(ring: Lattice, orbit: Orbit = None) -> np.ndarray:
    """Computes the diffusion matrix of the whole ring

    Parameters:
        ring:           Lattice description. Radiation must be ON
        orbit:          Avoids looking for initial the closed orbit if it is
          already known ((6,) array).

    Returns:
        diffmat (ndarray):  Diffusion matrix (6,6)
    """
    bbcum, _ = _dmatr(ring, orbit=orbit)
    diffmat = [(bbc + bbc.T) / 2 for bbc in bbcum]
    return np.round(diffmat[-1], 24)


@check_radiation(True)
def gen_quantdiff_elem(ring: Lattice, orbit: Orbit = None) -> QuantumDiffusion:
    """Generates a quantum diffusion element

    Parameters:
        ring:           Lattice description. Radiation must be ON
        orbit:          Avoids looking for initial the closed orbit if it is
          already known ((6,) array).

    Returns:
        diffElem (QuantumDiffusion): Quantum diffusion element
    """
    dmat = quantdiffmat(ring, orbit=orbit)
    lmat = np.asfortranarray(_lmat(dmat))
    diff_elem = QuantumDiffusion("Diffusion", lmat)
    return diff_elem


@check_radiation(True)
def tapering(ring: Lattice, multipoles: bool = True, niter: int = 1, **kwargs) -> None:
    """Scales magnet strengths

    Scales magnet strengths with local energy to cancel the closed orbit
    and optics errors due to synchrotron radiations. The dipole bending angle
    is changed by changing the reference momentum. This is the ideal
    tapering scheme where magnets and multipoles components (PolynomB and
    PolynomA) are scaled individually.

    Warning:
        This method works only for lattices without errors and
        corrections: if not all corrections and field errors will also be
        scaled !!!

    Parameters:
        ring:           Lattice description
        multipoles:     Scales all multipoles
        niter:          Number of iteration

    Keyword Args:
        method (ELossMethod):   Method for energy loss computation.
          See :py:class:`.ELossMethod`
        XYStep (float):         Step size.
          Default: :py:data:`DConstant.XYStep <.DConstant>`
        DPStep (float):         Momentum step size.
          Default: :py:data:`DConstant.DPStep <.DConstant>`
    """

    xy_step = kwargs.pop("XYStep", DConstant.XYStep)
    dp_step = kwargs.pop("DPStep", DConstant.DPStep)
    method = kwargs.pop("method", ELossMethod.TRACKING)
    dipin = ring.get_bool_index(Dipole)
    dipout = np.roll(dipin, 1)
    multin = ring.get_bool_index(Multipole) & ~dipin
    multout = np.roll(multin, 1)

    for _i in range(niter):
        _, o6 = find_orbit6(
            ring,
            refpts=range(len(ring) + 1),
            XYStep=xy_step,
            DPStep=dp_step,
            method=method,
        )
        dpps = (o6[dipin, 4] + o6[dipout, 4]) / 2.0
        set_value_refpts(ring, dipin, "FieldScaling", 1 + dpps)
        if multipoles:
            _, o6 = find_orbit6(
                ring,
                refpts=range(len(ring) + 1),
                XYStep=xy_step,
                DPStep=dp_step,
                method=method,
            )
            dppm = (o6[multin, 4] + o6[multout, 4]) / 2
            set_value_refpts(ring, multin, "FieldScaling", 1 + dppm)


Lattice.ohmi_envelope = ohmi_envelope
Lattice.get_radiation_integrals = get_radiation_integrals
Lattice.tapering = tapering
