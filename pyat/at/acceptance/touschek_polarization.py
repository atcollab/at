r"""Touschek lifetime including beam polarization, following Fu et al. [1]_.

The inverse lifetime is
:math:`1/\tau(P) = \langle C\rangle + \langle F\rangle P^2`, where P is the
polarization degree and the brackets denote length-weighted averages.
The local contributions C(s) and F(s) have units of inverse seconds.

This model assumes ultrarelativistic electrons in the laboratory frame and
nonrelativistic collisions in the center-of-mass frame. Its accuracy must be
checked when the latter approximation is not justified.

The implementation is kept separate from :mod:`.touschek` for comparison before
an eventual replacement of the existing Touschek lifetime calculation.

References:
    .. [1] H. Fu, Z. Duan, and J. Wang, "New analytical formula on the Touschek
       lifetime of the polarized electron beam in electron storage rings",
       Phys. Rev. Accel. Beams 27, 124401 (2024).
       https://doi.org/10.1103/PhysRevAccelBeams.27.124401
"""

from numbers import Real

import numpy as np
from scipy import integrate
from scipy.special import i0

from .touschek import get_beam_sizes, _init_ma_rp
from ..constants import qe, clight, _e_radius
from ..lattice import Lattice

__all__ = ["get_lifetime_polarized"]


def _int_c(k, um, B1, B2):
    r"""Integrand of C(s), from Fu et al., Eq. (27).

    With :math:`u_m = \beta^2\delta_m^2`, the lower bound is
    :math:`t_m = u_m/(1-u_m)`. The substitution :math:`t = \tan^2 k`
    maps :math:`[t_m, \infty)` to :math:`[\arctan\sqrt{t_m}, \pi/2)`,
    with Jacobian :math:`dt/dk = 2\sqrt{t}(1+t)`. This is the same finite-interval
    transformation used after Eq. (42) in Piwinski's
    https://arxiv.org/abs/physics/9903034.

    To avoid overflow for large Bessel arguments, at :math:`x = B_2 t \geq 500`
    use :math:`I_0(x) \sim \exp(x)/\sqrt{2\pi x}` and combine its exponential
    with :math:`\exp(-B_1 t)` before evaluation.
    """
    t = np.tan(k) ** 2
    fact = np.sqrt((1 + t) ** 5 / t**3) * (
        t / (um * (1 + t)) - 1 + 0.5 * np.log(um * (1 + t) / t)
    )

    if B2 * t < 500:
        intc = fact * np.exp(-B1 * t) * i0(B2 * t) * 2 * np.sqrt(t) * (1 + t)
    else:
        intc = (
            fact
            * np.exp((B2 - B1) * t)
            / np.sqrt(2 * np.pi * B2 * t)
            * 2
            * np.sqrt(t)
            * (1 + t)
        )

    return intc


def _int_f(k, um, B1, B2):
    r"""Integrand of F(s), from Fu et al., Eq. (28).

    The substitution :math:`t = \tan^2 k`, bounds, and Jacobian are as in
    :func:`_int_c`. For :math:`x = B_2 t \geq 500`, substitute
    :math:`I_0(x) \sim \exp(x)/\sqrt{2\pi x}` and combine exponentials to
    avoid overflow.
    """
    t = np.tan(k) ** 2

    fact = 0.5 * np.sqrt((1 + t) ** 5 / t**3) * np.log(um * (1 + t) / t)

    if B2 * t < 500:
        intf = fact * np.exp(-B1 * t) * i0(B2 * t) * 2 * np.sqrt(t) * (1 + t)
    else:
        intf = (
            fact
            * np.exp((B2 - B1) * t)
            / np.sqrt(2 * np.pi * B2 * t)
            * 2
            * np.sqrt(t)
            * (1 + t)
        )

    return intf


def _get_vals(
    ring,
    rp,
    ma,
    emity,
    bunch_curr,
    emitx=None,
    sigs=None,
    sigp=None,
    zn=None,
    epsabs=1.0e-16,
    epsrel=1.0e-12,
    compute_f=False,
):
    """Calculate C(s) and, only when requested, F(s) in inverse seconds."""
    emitx, sigs, sigp = get_beam_sizes(
        ring,
        bunch_curr,
        zn=zn,
        emitx=emitx,
        sigs=sigs,
        sigp=sigp,
    )

    nc = bunch_curr / ring.revolution_frequency / qe
    beta2 = ring.beta * ring.beta
    gamma2 = ring.gamma * ring.gamma

    emit = np.array([emitx, emity])
    _, _, ld = ring.get_optics(refpts=rp)

    bxy = ld.beta
    bxy2 = bxy * bxy
    axy = ld.alpha
    dxy = ld.dispersion[:, [0, 2]]
    dxy2 = dxy * dxy
    dpxy = ld.dispersion[:, [1, 3]]

    sigb = np.sqrt(bxy * emit)
    sigb2 = sigb * sigb
    sigp2 = sigp * sigp

    dt = dxy * axy + dpxy * bxy
    dt2 = dt * dt

    sigh2 = 1 / (1 / sigp2 + np.sum((dxy2 + dt2) / sigb2, axis=1))

    bs = bxy2 / sigb2 * (1 - (sigh2 * (dt2 / sigb2).T).T)

    bg2i = 1 / (2 * beta2 * gamma2)

    B1 = bg2i * np.sum(bs, axis=1)

    B2sq = bg2i**2 * (
        np.diff(bs, axis=1).T ** 2
        + 4 * sigh2**2 * np.prod(bxy2 * dt2, axis=1) / np.prod(sigb2**2, axis=1)
    )

    B2 = np.sqrt(B2sq).ravel()

    val_c = np.zeros((2, len(rp)))
    val_f = np.zeros_like(val_c) if compute_f else None

    for i in range(2):
        dpp = ma[:, i]
        um = beta2 * dpp * dpp
        km = np.arctan(np.sqrt(um / (1 - um)))

        for ii in range(len(rp)):
            args = (um[ii], B1[ii], B2[ii])
            val_c[i, ii], *_ = integrate.quad(
                _int_c,
                km[ii],
                np.pi / 2,
                args=args,
                epsabs=epsabs,
                epsrel=epsrel,
            )

            if val_f is not None:
                val_f[i, ii], *_ = integrate.quad(
                    _int_f,
                    km[ii],
                    np.pi / 2,
                    args=args,
                    epsabs=epsabs,
                    epsrel=epsrel,
                )

    prefactor = (
        clight
        * _e_radius**2
        * np.prod(bxy, axis=1)
        * np.sqrt(sigh2)
        * nc
        / (
            8
            * np.sqrt(np.pi)
            * beta2
            * gamma2**2
            * np.prod(sigb2, axis=1)
            * sigs
            * sigp
        )
    )

    val_c *= prefactor
    if val_f is not None:
        val_f *= prefactor

    return val_c, val_f


def get_lifetime_polarized(
    ring,
    emity,
    bunch_curr,
    emitx=None,
    sigs=None,
    sigp=None,
    zn=None,
    momap=None,
    refpts=None,
    offset=None,
    polarization: float = 0.0,
    **kwargs,
):
    r"""Compute the Touschek lifetime at a given beam polarization.

    Uses the approximation of
    `Fu et al. <https://doi.org/10.1103/PhysRevAccelBeams.27.124401>`_,
    assuming ultrarelativistic electrons in the laboratory frame and
    nonrelativistic collisions in the center-of-mass frame. The latter
    assumption must be checked before applying this model at higher energies.
    This is also the model used at zero polarization, unlike the Piwinski
    formula in :func:`.get_lifetime`.

    The local inverse lifetime is :math:`C(s) + F(s)P^2`. Both contributions
    are averaged over the two momentum-loss directions, then weighted by the
    lengths of the evaluated elements:

    .. math::

        \frac{1}{\tau(P)} = \langle C\rangle + \langle F\rangle P^2,
        \qquad A_\mathrm{avg} = -\frac{\langle F\rangle}{\langle C\rangle}.

    Args:
        ring: Lattice used for the calculation.
        emity: Vertical emittance [m].
        bunch_curr: Bunch current [A].
        emitx: Horizontal emittance [m]. If None, computed from the radiation
            equilibrium emittance.
        sigs: RMS bunch length [m]. If None, computed from the lattice and zn.
        sigp: RMS relative energy spread (dimensionless). If None, computed
            from the lattice.
        zn: Full-ring longitudinal impedance Z/n [Ohm]. If None, use the
            zero-current bunch length; otherwise include bunch lengthening.
        momap: Momentum aperture, a (N, 2) array of relative momentum bounds
            corresponding to refpts. If None, calculate it by tracking.
        refpts: Reference points at which momap is supplied or calculated.
            Defaults to all element entrances. Zero-length elements are
            excluded from the lifetime average.
        offset: Initial orbit, broadcastable to (N, 6). If None, use the
            closed orbit when calculating momentum acceptance.
        polarization: Finite real scalar degree between 0 and 1, not a
            percentage. Default: 0.0. Booleans are not accepted. At zero,
            F and A_avg are not computed.

    Keyword Args:
        interpolate: Interpolate the aperture onto finite-length elements
            between the first and last selected refpts. Default: True.
        epsabs: Absolute quadrature tolerance. Default: 1.0e-16.
        epsrel: Relative quadrature tolerance. Default: 1.0e-12.

    Other keyword arguments are passed to :func:`.get_momentum_acceptance`
    when momap is not supplied, including resolution, amplitude, nturns,
    grid_mode, and use_mp.

    Returns:
        tl: Touschek lifetime at the requested polarization [s].
        ma: Momentum aperture at the returned refpts, shape (N, 2).
        refpts: Reference points used for the lifetime calculation, shape (N,).
        C: Local unpolarized inverse-lifetime contribution [1/s], shape (N,).
            Returned only for nonzero polarization.
        F: Local polarization-dependent inverse-lifetime coefficient [1/s],
            shape (N,). Returned only for nonzero polarization.
        A_avg: Dimensionless average analyzing power. Returned only for
            nonzero polarization.

    Raises:
        TypeError: If polarization is not a real scalar degree or is a boolean.
        ValueError: If polarization is nonfinite or outside [0, 1].

    Note:
        Returns ``(tl, ma, refpts)`` at zero polarization and
        ``(tl, ma, refpts, C, F, A_avg)`` otherwise. The returned C, F, and A_avg
        do not depend on the nonzero polarization degree used to obtain them.
        Once A_avg and the unpolarized lifetime are known, another degree can
        be evaluated without repeating the integrals:

        .. math::

            \tau(P) = \frac{\tau(0)}{1 - A_\mathrm{avg} P^2}.
    """
    if isinstance(polarization, (bool, np.bool_)) or not isinstance(polarization, Real):
        msg = "polarization must be a real scalar degree, not a boolean"
        raise TypeError(msg)
    if not 0.0 <= polarization <= 1.0:
        msg = "polarization must be finite and between 0 and 1"
        raise ValueError(msg)
    polarization = float(polarization)

    interpolate = kwargs.pop("interpolate", True)
    epsabs = kwargs.pop("epsabs", 1.0e-16)
    epsrel = kwargs.pop("epsrel", 1.0e-12)

    ma, rp = _init_ma_rp(
        ring,
        refpts=refpts,
        offset=offset,
        momap=momap,
        interpolate=interpolate,
        check_zero=True,
        **kwargs,
    )

    val_c, val_f = _get_vals(
        ring,
        rp,
        ma,
        emity,
        bunch_curr,
        emitx=emitx,
        sigs=sigs,
        sigp=sigp,
        zn=zn,
        epsabs=epsabs,
        epsrel=epsrel,
        compute_f=polarization > 0.0,
    )

    lengths = np.array([ring[i].Length for i in rp])
    c = np.mean(val_c, axis=0)
    c_avg = np.average(c, weights=lengths)
    if val_f is None:
        return 1 / c_avg, ma, rp

    f = np.mean(val_f, axis=0)
    f_avg = np.average(f, weights=lengths)
    tl = 1 / (c_avg + f_avg * polarization**2)
    a_avg = -f_avg / c_avg

    return tl, ma, rp, c, f, a_avg


Lattice.get_lifetime_polarized = get_lifetime_polarized
