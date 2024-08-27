"""
Module to compute the touschek lifetime of the ring
"""

import at
import numpy
from ..lattice import Lattice, AtError, AtWarning
import warnings
from scipy.special import iv
from scipy import integrate
from scipy.optimize import fsolve
from ..constants import qe, clight, _e_radius


__all__ = ["get_bunch_length_espread", "get_lifetime", "get_scattering_rate"]


def get_bunch_length_espread(ring, zn=None, bunch_curr=None, espread=None):
    """Haissinski equation solver

    Solves the Haissinski formula and returns the bunch length and energy
    spread for given bunch current and :math:`Z/n`. If both ``zn`` and
    ``bunch_curr`` are ``None``, zero current case, otherwise both are needed
    for the calculation

    args:
        ring:             ring use for tracking

    keyword args:
        zn=None:          :math:`Z/n` for the full ring
        bunch_curr=None:  Bunch current
        espread=None:     Energy spread, if ``None`` use lattice parameter

    Returns:
        Bunch length, energy spread
    """

    def haissinski(x, cst):
        return x**3 - x - cst

    ep = at.envelope_parameters(ring.radiation_on(copy=True))
    bl0 = ep.sigma_l
    if espread is None:
        espread = ep.sigma_e
    if zn is None and bunch_curr is None:
        bl = bl0
    elif zn is None or bunch_curr is None:
        raise AtError(
            "Please provide both current and Z/n for bunch " "length calculation"
        )
    else:
        vrf = ring.get_rf_voltage()
        h = ring.harmonic_number
        etac = at.get_slip_factor(ring.radiation_off(copy=True))
        cs = numpy.cos(ep.phi_s)
        nus = ep.f_s / ring.revolution_frequency
        cst = (
            -0.5
            * numpy.sqrt(numpy.pi)
            * bunch_curr
            * zn
            / (vrf * h * cs * (abs(etac) * espread / nus) ** 3)
        )
        bl = bl0 * fsolve(haissinski, numpy.array([1.0]), args=cst)[0]
    return bl, espread


def get_beam_sizes(ring, bunch_curr, zn=None, emitx=None, sigs=None, sigp=None):
    if zn is None:
        bunch_curr = None
    if sigs is None or sigp is None:
        sigsi, sigpi = get_bunch_length_espread(
            ring, zn=zn, bunch_curr=bunch_curr, espread=sigp
        )
        if sigs is None:
            sigs = sigsi
        if sigp is None:
            sigp = sigpi
    if emitx is None:
        _, bd, _ = at.ohmi_envelope(ring.radiation_on(copy=True))
        emitx = bd.mode_emittances[0]
    return emitx, sigs, sigp


def int_piwinski(k, km, B1, B2):
    r"""
    Integrand of the piwinski formula
    In case the Bessel function has too large value
    (more than :math:`10^251`) it
    is substituted by its exponential approximation:
    :math:`I_0(x)~\frac{\exp(x)}{\sqrt{2 \pi x}}`
    """
    t = numpy.tan(k) ** 2
    tm = numpy.tan(km) ** 2
    fact = (
        (2 * t + 1) ** 2 * (t / tm / (1 + t) - 1) / t
        + t
        - numpy.sqrt(t * tm * (1 + t))
        - (2 + 1 / (2 * t)) * numpy.log(t / tm / (1 + t))
    )
    if B2 * t < 500:
        intp = fact * numpy.exp(-B1 * t) * iv(0, B2 * t) * numpy.sqrt(1 + t)
    else:
        intp = (
            fact
            * numpy.exp(B2 * t - B1 * t)
            / numpy.sqrt(2 * numpy.pi * B2 * t)
            * numpy.sqrt(1 + t)
        )
    return intp


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
    sig = numpy.sqrt(sigb2 + sigp2 * dxy2)
    sig2 = sig * sig
    dt = dxy * axy + dpxy * bxy
    dt2 = dt * dt
    sigh2 = 1 / (1 / sigp2 + numpy.sum((dxy2 + dt2) / sigb2, axis=1))

    bs = bxy2 / sigb2 * (1 - (sigh2 * (dt2 / sigb2).T).T)
    bg2i = 1 / (2 * beta2 * gamma2)
    B1 = bg2i * numpy.sum(bs, axis=1)
    B2sq = (
        bg2i
        * bg2i
        * (
            numpy.diff(bs, axis=1).T ** 2
            + 4
            * sigh2
            * sigh2
            * numpy.prod(bxy2 * dt2, axis=1)
            / numpy.prod(sigb2 * sigb2, axis=1)
        )
    )
    B2 = numpy.squeeze(numpy.sqrt(B2sq))

    val = numpy.zeros((2, len(rp)))
    for i in range(2):
        dpp = ma[:, i]
        um = beta2 * dpp * dpp
        km = numpy.arctan(numpy.sqrt(um))

        for ii in range(len(rp)):
            args = (km[ii], B1[ii], B2[ii])
            val[i, ii], *_ = integrate.quad(
                int_piwinski,
                args[0],
                numpy.pi / 2,
                args=args,
                epsabs=epsabs,
                epsrel=epsrel,
            )

        val[i] *= (
            _e_radius**2
            * clight
            * nc
            / (
                8
                * numpy.pi
                * gamma2
                * sigs
                * numpy.sqrt(
                    numpy.prod(sig2, axis=1) - sigp2 * sigp2 * numpy.prod(dxy2, axis=1)
                )
            )
            * 2
            * numpy.sqrt(numpy.pi * (B1 * B1 - B2 * B2))
        )
    return val


def _init_ma_rp(
    ring,
    refpts=None,
    offset=None,
    momap=None,
    interpolate=True,
    check_zero=True,
    **kwargs,
):
    if refpts is None:
        refpts = ring.get_uint32_index(at.All, endpoint=False)
    else:
        refpts = ring.get_uint32_index(refpts)

    if offset is not None:
        try:
            offset = numpy.broadcast_to(offset, (len(refpts), 6))
        except ValueError:
            msg = "offset and refpts have incoherent " "shapes: {0}, {1}".format(
                offset.shape, refpts.shape
            )
            raise AtError(msg)

    if momap is not None:
        assert len(momap) == len(
            refpts
        ), "Input momap and refpts have different lengths"

    if check_zero:
        mask = [ring[r].Length > 0.0 for r in refpts]
    else:
        mask = [True for _ in refpts]

    if not numpy.all(mask):
        zerolength_warning = "zero-length elements removed " "from lifetime calculation"
        warnings.warn(AtWarning(zerolength_warning))
        assert len(refpts[mask]) > 2, \
            'After removing zero-length elements: len(refpts)<2'

    refpts = refpts[mask]
    if offset is not None:
        offset = offset[mask]

    if momap is None:
        resolution = kwargs.pop("resolution", 1.0e-3)
        amplitude = kwargs.pop("amplitude", 0.1)
        kwargs.update({"refpts": refpts})
        kwargs.update({"offset": offset})
        momap, _, _ = ring.get_momentum_acceptance(resolution, amplitude, **kwargs)
    else:
        momap = momap[mask]

    if interpolate:
        refpts_all = numpy.array(
            [i for i in range(refpts[0], refpts[-1] + 1) if ring[i].Length > 0]
        )
        spos = numpy.squeeze(ring.get_s_pos(refpts))
        spos_all = numpy.squeeze(ring.get_s_pos(refpts_all))
        momp = numpy.interp(spos_all, spos, momap[:, 0])
        momn = numpy.interp(spos_all, spos, momap[:, 1])
        momap_all = numpy.vstack((momp, momn)).T
        ma, rp = momap_all, refpts_all
    else:
        ma, rp = momap, refpts
    return ma, rp


def get_lifetime(
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

    Computes the touschek lifetime using the Piwinski formula

    args:
        ring:            ring use for tracking
        emity:           verticla emittance
        bunch_curr:      bunch current

    keyword args:
        emitx=None:      horizontal emittance
        sigs=None:       rms bunch length
        sigp=None:       energy spread
        zn=None:         full ring :math:`Z/n`
        momap=None:      momentum aperture, has to be consistent with
                         ``refpts`` if provided the momentum aperture is
                         not calculated
        refpts=None:     ``refpts`` where the momentum aperture is calculated,
                         the default is to compute it for all elements in the
                         ring, ``len(refpts)>2`` is required
        resolution:      minimum distance between 2 grid points, default=1.0e-3
        amplitude:       max. amplitude for ``RADIAL`` and ``CARTESIAN`` or
                         initial step in ``RECURSIVE``
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
        tl: touschek lifetime, double expressed in seconds
        ma: momentum aperture (len(refpts), 2) array
        refpts: refpts used for momentum aperture calculation
                (len(refpts), ) array

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
        **kwargs,
    )
    length_all = numpy.array([e.Length for e in ring[rp]])
    vals = _get_vals(
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
    tl = 1 / numpy.mean([sum(v * length_all.T) / sum(length_all) for v in vals])
    return tl, ma, rp


def get_scattering_rate(
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
    """Touschek scattering rate calculation

    Computes the touschek scattering using the Piwinski formula

    args:
        ring:            ring use for tracking
        emity:           verticla emittance
        bunch_curr:      bunch current

    keyword args:
        emitx=None:      horizontal emittance
        sigs=None:       rms bunch length
        sigp=None:       energy spread
        zn=None:         full ring :math:`Z/n`
        momap=None:      momentum aperture, has to be consistent with
                         ``refpts`` if provided the momentum aperture
                         is not calculated
        refpts=None:     ``refpts`` where the momentum aperture is calculated,
                         the default is to compute it for all elements in the
                         ring, ``len(refpts)>2`` is required
        resolution:      minimum distance between 2 grid points, default=1.0e-3
        amplitude:       max. amplitude for ``RADIAL`` and ``CARTESIAN`` or
                         initial step in ``RECURSIVE``
                         default = 0.1
        nturns=1024:     Number of turns for the tracking
        dp=None:         static momentum offset
        grid_mode:       ``at.GridMode.CARTESIAN/RADIAL`` track full vector
                         (default). ``at.GridMode.RECURSIVE``: recursive search
        use_mp=False:    Use python multiprocessing (``patpass``, default use
                         ``lattice_pass``). In case multi-processing is not
                         enabled ``GridMode`` is forced to
                         ``RECURSIVE`` (most efficient in single core)
        verbose=False:   Print out some inform
        epsabs, epsrel:  integral absolute and relative tolerances

    Returns:
        scattering_rate: scattering rate double array (len(refpts), )
                         expressed in event/seconds
        ma: momentum aperture (len(refpts), 2) array
        refpts: refpts used for momentum aperture calculation
                (len(refpts), ) array
    """
    interpolate = kwargs.pop("interpolate", True)
    epsabs = kwargs.pop("epsabs", 1.0e-16)
    epsrel = kwargs.pop("epsrel", 1.0e-12)
    ma, rp = _init_ma_rp(
        ring,
        refpts=refpts,
        momap=momap,
        interpolate=interpolate,
        check_zero=False,
        offset=offset,
        **kwargs,
    )
    vals = _get_vals(
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
    scattering_rate = (
        numpy.mean(vals, axis=0) * bunch_curr / ring.revolution_frequency / qe
    )
    return scattering_rate, ma, rp


Lattice.get_bunch_length_espread = get_bunch_length_espread
Lattice.get_lifetime = get_lifetime
Lattice.get_scattering_rate = get_scattering_rate
