"""
Implementation of the calculation of Touschek lifetime including polarization
Reference: Fu et al. doi/10.1103/PhysRevAccelBeams.27.124401.

The Touschek lifetime can be given as 1/t = <C> + <F>P**2,
where t is the Touschek lifetime, P is the polarization, and <..> is the average.
C and F are a function of the longitudinal position.
"""

import numpy
from scipy import integrate
from scipy.special import iv

import at

from .touschek import (
    get_beam_sizes,
    _init_ma_rp,
)
from ..constants import qe, clight, _e_radius
from ..lattice import Lattice


def int_C(k, km, um, B1, B2):
    r"""
    Integrand of C(s) function.
    The infinite integral has been turned into a finite integral with the
    transformation t = tan(k)**2.
    """
    t = numpy.tan(k) ** 2
    tm = numpy.tan(km) ** 2

    one_pl_t5 = (1 + t) * (1 + t) * (1 + t) * (1 + t) * (1 + t)
    t3 = t * t * t

    beta_ma2 = um

    fact = (
        numpy.sqrt(one_pl_t5 / t3)
        * (
            t / (beta_ma2 * (1 + t))
            - 1
            + 0.5 * numpy.log(beta_ma2 * (1 + t) / t)
        )
    )

    if B2 * t < 500:
        intc = (
            fact
            * numpy.exp(-B1 * t)
            * iv(0, B2 * t)
            * 2
            * numpy.sqrt(t)
            * (1 + t)
        )
    else:
        intc = (
            fact
            * numpy.exp(B2 * t - B1 * t)
            / numpy.sqrt(2 * numpy.pi * B2 * t)
            * 2
            * numpy.sqrt(t)
            * (1 + t)
        )

    return intc


def int_F(k, km, um, B1, B2):
    r"""
    Integrand of F(s) function.
    The infinite integral has been turned into a finite integral with the
    transformation t = tan(k)**2.
    """
    t = numpy.tan(k) ** 2
    tm = numpy.tan(km) ** 2

    one_pl_t5 = (1 + t) * (1 + t) * (1 + t) * (1 + t) * (1 + t)
    t3 = t * t * t

    beta_ma2 = um

    fact = (
        0.5
        * numpy.sqrt(one_pl_t5 / t3)
        * numpy.log(beta_ma2 * (1 + t) / t)
    )

    if B2 * t < 500:
        intf = (
            fact
            * numpy.exp(-B1 * t)
            * iv(0, B2 * t)
            * 2
            * numpy.sqrt(t)
            * (1 + t)
        )
    else:
        intf = (
            fact
            * numpy.exp(B2 * t - B1 * t)
            / numpy.sqrt(2 * numpy.pi * B2 * t)
            * 2
            * numpy.sqrt(t)
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
    polarization=True,
):
    r"""
    Calculate the integrals to obtain C(s).

    If polarization=True, also calculate F(s).
    """
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

    emit = numpy.array([emitx, emity])
    _, _, ld = ring.get_optics(refpts=rp)

    bxy = ld.beta
    bxy2 = bxy * bxy
    axy = ld.alpha
    dxy = ld.dispersion[:, [0, 2]]
    dxy2 = dxy * dxy
    dpxy = ld.dispersion[:, [1, 3]]

    sigb = numpy.sqrt(bxy * emit)
    sigb2 = sigb * sigb
    sigp2 = sigp * sigp

    dt = dxy * axy + dpxy * bxy
    dt2 = dt * dt

    sigh2 = 1 / (
        1 / sigp2
        + numpy.sum((dxy2 + dt2) / sigb2, axis=1)
    )

    bs = bxy2 / sigb2 * (
        1 - (sigh2 * (dt2 / sigb2).T).T
    )

    bg2i = 1 / (2 * beta2 * gamma2)

    B1 = bg2i * numpy.sum(bs, axis=1)

    B2sq = (
        bg2i**2
        * (
            numpy.diff(bs, axis=1).T**2
            + 4
            * sigh2**2
            * numpy.prod(bxy2 * dt2, axis=1)
            / numpy.prod(sigb2**2, axis=1)
        )
    )

    B2 = numpy.squeeze(numpy.sqrt(B2sq))

    
    val_C = numpy.zeros((2, len(rp)))

    
    if polarization:
        val_F = numpy.zeros((2, len(rp)))

    for i in range(2):

        dpp = ma[:, i]
        um = beta2 * dpp * dpp
        km = numpy.arctan(numpy.sqrt(um / (1 - um)))

        for ii in range(len(rp)):

            args = (
                km[ii],
                um[ii],
                B1[ii],
                B2[ii],
            )

            
            val_C[i, ii], *_ = integrate.quad(
                int_C,
                args[0],
                numpy.pi / 2,
                args=args,
                epsabs=epsabs,
                epsrel=epsrel,
            )

          
            if polarization:
                val_F[i, ii], *_ = integrate.quad(
                    int_F,
                    args[0],
                    numpy.pi / 2,
                    args=args,
                    epsabs=epsabs,
                    epsrel=epsrel,
                )

    prefactor = (
        clight
        * _e_radius**2
        * numpy.prod(bxy, axis=1)
        * numpy.sqrt(sigh2)
        * nc
        / (
            8
            * numpy.sqrt(numpy.pi)
            * beta2
            * gamma2**2
            * numpy.prod(sigb2, axis=1)
            * sigs
            * sigp
        )
    )

    val_C *= prefactor

    if polarization:
        val_F *= prefactor
        return val_C, val_F

    return val_C


def get_touschek_complete(
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
    polarization=False,
    **kwargs,
):
    """Touschek lifetime calculation.

    Computes the Touschek scattering C(s) function.

    If polarization=True, also computes the polarization-dependent
    F(s) function and A_avg.

    If polarization=False, F(s) and A_avg are not calculated.

    Returns:
        polarization=True:
            tl, ma, rp, C, F, A_avg

        polarization=False:
            tl, ma, rp
    """

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

    length_all = numpy.array([e.Length for e in ring[rp]])

    if polarization:
        val_C, val_F = _get_vals(
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
            polarization=True,
        )
    else:
        val_C = _get_vals(
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
            polarization=False,
        )

    tl = 1 / numpy.mean(
        [
            sum(v * length_all.T) / sum(length_all)
            for v in val_C
        ]
    )

    C = (val_C[0, :] + val_C[1, :]) / 2

    if not polarization:
        return tl, ma, rp

   
    F = (val_F[0, :] + val_F[1, :]) / 2

    lengths = numpy.array([ring[i].Length for i in rp])

    C_avg = numpy.sum(C * lengths) / numpy.sum(lengths)
    F_avg = numpy.sum(F * lengths) / numpy.sum(lengths)

    A_avg = -F_avg / C_avg

    return tl, ma, rp, C, F, A_avg


Lattice.get_touschek_complete = get_touschek_complete
