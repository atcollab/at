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
    The infinite integral has been turned into a finite integral with the transformation
    t = tan(k)**2
    """
    t = numpy.tan(k) ** 2
    tm = numpy.tan(km) ** 2
    
    one_pl_t5 = (1+t) * (1+t) * (1+t) * (1+t) * (1+t)
    t3 = t * t * t
   
    beta_ma2 = um
    
    fact =(
    numpy.sqrt(one_pl_t5 / t3 )  *
    (
    t / (beta_ma2 *(1 + t))
    - 1
    + 0.5 * numpy.log(beta_ma2 *(1 + t)/ t)
    )
    )
    
    if B2 * t < 500:
        intc = fact * numpy.exp(-B1 * t) * iv(0, B2 * t)  * 2 *  numpy.sqrt( t) * (1+t)
    else:
        intc = (
            fact
            * numpy.exp(B2 * t - B1 * t)
            / numpy.sqrt(2 * numpy.pi * B2 * t)
           * 2 *  numpy.sqrt( t) * (1+t)
        )
    return intc


def int_F(k, km, um, B1, B2):
    r"""
    Integrand of F(s) function.
    The infinite integral has been turned into a finite integral with the transformation
    t = tan(k)**2
    """
    t = numpy.tan(k) ** 2
    tm = numpy.tan(km) ** 2
    
    one_pl_t5 = (1+t) * (1+t) * (1+t) * (1+t) * (1+t)
    t3 = t * t * t
    beta_ma2 = um
    fact =(
    0.5 *    numpy.sqrt(one_pl_t5 / t3 )     
    * numpy.log(beta_ma2 *(1 + t)/ t)
     )
    
    if B2 * t < 500:
        intf = fact * numpy.exp(-B1 * t) * iv(0, B2 * t)  * 2 *  numpy.sqrt( t) * (1+t)
    else:
        intf = (
            fact
            * numpy.exp(B2 * t - B1 * t)
            / numpy.sqrt(2 * numpy.pi * B2 * t)
           * 2 *  numpy.sqrt( t) * (1+t)
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
):
    r"""
    Calculate the integrals to obtain C(s) and F(s)
    """
    emitx, sigs, sigp = get_beam_sizes(
        ring, bunch_curr, zn=zn, emitx=emitx, sigs=sigs, sigp=sigp
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

    sigh2 = 1 / (1 / sigp2 + numpy.sum((dxy2 + dt2) / sigb2, axis=1))

    bs = bxy2 / sigb2 * (1 - (sigh2 * (dt2 / sigb2).T).T)

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
    val_F = numpy.zeros((2, len(rp)))

    for i in range(2):

        dpp = ma[:, i]
        um = beta2 * dpp * dpp
        km = numpy.arctan(numpy.sqrt(um / (1 - um)))

        for ii in range(len(rp)):

            args = (km[ii], um[ii], B1[ii], B2[ii])

            val_C[i, ii], *_ = integrate.quad(
                int_C,
                args[0],
                numpy.pi / 2,
                args=args,
                epsabs=epsabs,
                epsrel=epsrel,
            )

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
    val_F *= prefactor

    return val_C, val_F
    
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
    **kwargs,
):
    """Touschek lifetime calculation

    Computes the touschek scattering C(s) function and the 
    polarization-dependent F(s) function.

    The Touschek lifetime is given as 1/t = <C> + <F>P**2,
    where t is the Touschek lifetime, P is the polarization of the beam, and <..> is the average.

    C and F functions are taken from Fu et al. doi/10.1103/PhysRevAccelBeams.27.124401.

    

    args:
        ring:            ring use for tracking
        emity:           vertical emittance [m]
        bunch_curr:      bunch current [A]

    keyword args:
        emitx=None:      horizontal emittance [m]
        sigs=None:       rms bunch length [m]
        sigp=None:       energy spread
        zn=None:         full ring :math:`Z/n` [Ohm]
        momap=None:      momentum aperture, has to be consistent with
                         ``refpts`` if provided the momentum aperture is
                         not calculated
        refpts=None:     ``refpts`` where the momentum aperture is calculated,
                         the default is to compute it for all elements in the
                         ring, ``len(refpts)>2`` is required
        resolution:      minimum distance between 2 grid points, default=1.0e-3
        amplitude:       max. amplitude for ``RADIAL`` and ``CARTESIAN`` or
                         initial step in ``RECURSIVE`` [m]
                         default = 0.1
        nturns=1024:     Number of turns for the tracking
        dp=None:         static momentum offset
        offset:         initial orbit. Default: closed orbit
        grid_mode:       ``at.GridMode.CARTESIAN/RADIAL`` track full vector
                         (default). ``at.GridMode.RECURSIVE``: recursive search
        use_mp=False:    Use python multiprocessing (``patpass``, default use
                         ``lattice_pass``). In case multi-processing is not
                         enabled ``GridMode`` is forced to
                         ``RECURSIVE`` (most efficient in single core)
        verbose=False:   Print out some inform
        epsabs, epsrel:  integral absolute and relative tolerances

    Returns:
        tl: touschek lifetime, double expressed in seconds, for the beam with 0 polarization
        ma: momentum aperture (len(refpts), 2) array
        refpts: refpts used for momentum aperture calculation
                (len(refpts), ) array
        C: C(s) scattering function. 
        F: F(s) polarization-dependency of Touschek scattering.
        A_avg: average analysing power, used to fit the polarization of the beam.

    .. note::

       * When``use_mp=True`` all the available CPUs will be used.
         This behavior can be changed by setting
         ``at.DConstant.patpass_poolsize`` to the desired value
       * When multiple ``refpts`` are provided particles are first
         projected to the beginning of the ring with tracking. Then,
         all particles are tracked up to ``nturns``. This allows to
         do most of the work in a single function call and allows for
         full parallelization.
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
    )
    tl = 1 / numpy.mean([sum(v * length_all.T) / sum(length_all) for v in val_C])
    C = (val_C[0,:] + val_C[1,:])/2
    F = (val_F[0,:] + val_F[1,:])/2
    lengths = numpy.array([ring[i].Length for i in rp])
    C_avg = numpy.sum( C * lengths ) / numpy.sum(lengths)
    F_avg = numpy.sum( F * lengths ) / numpy.sum(lengths)
    A_avg = -F_avg / C_avg
    
    return tl, ma, rp, C, F, A_avg

Lattice.get_touschek_complete = get_touschek_complete

