"""
Non-linear optics
"""
import numpy
from typing import Optional, Sequence
from scipy.special import factorial
from ..lattice import Element, Lattice
from ..tracking import internal_lpass
from .orbit import Orbit, find_orbit
from .linear import get_tune, get_chrom, linopt6
from .harmonic_analysis import get_tunes_harmonic

__all__ = [
    "detuning",
    "chromaticity",
    "gen_detuning_elem",
    "tunes_vs_amp",
    "feed_down_polynomba",
    "feed_down_from_nth_ordera",
]


def tunes_vs_amp(
    ring: Lattice,
    amp: Optional[Sequence[float]] = None,
    dim: Optional[int] = 0,
    nturns: Optional[int] = 512,
    **kwargs,
) -> numpy.ndarray:
    r"""Generates a range of tunes for varying x, y, or z amplitudes

    Parameters:
        ring:               Lattice description
        amp:                Oscillation amplitude. If amp is None the
                            tunes are calculated from :py:func:`.linopt6`
                            else tracking is used
        dim:                Coordinate index of the oscillation
        nturns:             Number of turns
        method:             ``'laskar'`` or ``'fft'``. Default: ``'laskar'``
        num_harmonics:      Number of harmonic components to compute (before
                            mask applied)
        fmin:               Lower bound for tune
        fmax:               Upper bound for tune
        hann:               Turn on Hanning window. Default: :py:obj:`False`

    Returns:
        tunes (ndarray):    Array of tunes qx, qy with shape (len(amp), 2)
    """

    def _gen_part(ring, amp, dim, orbit, ld, nturns):
        invx = numpy.array(
            [
                [1 / numpy.sqrt(ld["beta"][0]), 0],
                [ld["alpha"][0] / numpy.sqrt(ld["beta"][0]), numpy.sqrt(ld["beta"][0])],
            ]
        )

        invy = numpy.array(
            [
                [1 / numpy.sqrt(ld["beta"][1]), 0],
                [ld["alpha"][1] / numpy.sqrt(ld["beta"][1]), numpy.sqrt(ld["beta"][1])],
            ]
        )
        part = (
            numpy.array(
                [
                    orbit,
                ]
                * len(amp)
            ).T
            + 1.0e-6
        )
        part[dim, :] += amp
        part = internal_lpass(ring, part, nturns=nturns)
        sh = part.shape
        partx = numpy.reshape(part[0, :], (sh[1], sh[3]))
        partxp = numpy.reshape(part[1, :], (sh[1], sh[3]))
        party = numpy.reshape(part[2, :], (sh[1], sh[3]))
        partyp = numpy.reshape(part[3, :], (sh[1], sh[3]))
        px = numpy.array(
            [numpy.matmul(invx, [partx[i], partxp[i]]) for i in range(len(amp))]
        )
        py = numpy.array(
            [numpy.matmul(invy, [party[i], partyp[i]]) for i in range(len(amp))]
        )
        return px[:, 0, :] - 1j * px[:, 1, :], py[:, 0, :] - 1j * py[:, 1, :]

    l0, bd, _ = linopt6(ring)
    orbit = l0["closed_orbit"]
    tunes = bd["tune"]

    if amp is not None:
        partx, party = _gen_part(ring, amp, dim, orbit, l0, nturns)
        qx = get_tunes_harmonic(partx, **kwargs)
        qy = get_tunes_harmonic(party, **kwargs)
        tunes = numpy.vstack((qx, qy)).T

    return numpy.array(tunes)


def detuning(
    ring: Lattice,
    xm: Optional[float] = 0.3e-4,
    ym: Optional[float] = 0.3e-4,
    npoints: Optional[int] = 3,
    nturns: Optional[int] = 512,
    **kwargs,
):
    """Computes the tunes for a sequence of amplitudes

    This function uses :py:func:`tunes_vs_amp` to compute the tunes for
    the specified amplitudes. Then it fits this data and returns
    the detuning coefficiant dQx/dx, dQy/dx, dQx/dy, dQy/dy and the
    qx, qy arrays versus x, y arrays

    Parameters:
        ring:       Lattice description
        xm:         Maximum x amplitude
        ym:         Maximum y amplitude
        npoints:    Number of points in each plane
        nturns:     Number of turns for tracking
        method:     ``'laskar'`` or ``'fft'``. Default: ``'laskar'``
        num_harmonics:  Number of harmonic components to compute
                       (before mask applied)
        fmin:       Lower bound for tune
        fmax:       Upper bound for tune
        hann:       Turn on Hanning window. Default: :py:obj:`False`

    Returns:
        q0 (ndarray): qx, qy from horizontal and vertical amplitude scans.
                      Fit results with shape (2, 2)
        q1 (ndarray): Amplitude detuning coefficients with shape (2, 2)
                      ``[[dQx/dx, dQy/dx], [dQx/dy, dQy/dy]]``
        x (ndarray): x amplitudes (npoints, )
        q_dx (ndarray): qx, qy tunes as a function of x amplitude (npoints, 2)
        y (ndarray): y amplitudes (npoints, )
        q_dy (ndarray): qx, qy tunes as a function of y amplitude (npoints, 2)
    """
    lindata0, _, _ = linopt6(ring)
    gamma = (1 + lindata0.alpha * lindata0.alpha) / lindata0.beta

    x = numpy.linspace(-xm, xm, npoints)
    y = numpy.linspace(-ym, ym, npoints)
    x2 = x * x
    y2 = y * y

    q_dx = tunes_vs_amp(ring, amp=x, dim=0, nturns=nturns, **kwargs)
    q_dy = tunes_vs_amp(ring, amp=y, dim=2, nturns=nturns, **kwargs)

    idx = numpy.isfinite(q_dx[:, 0]) & numpy.isfinite(q_dx[:, 1])
    idy = numpy.isfinite(q_dy[:, 0]) & numpy.isfinite(q_dy[:, 1])

    fx = numpy.polyfit(x2[idx], q_dx[idx], 1)
    fy = numpy.polyfit(y2[idy], q_dy[idy], 1)

    q0 = [[fx[1, 0], fx[1, 1]], [fy[1, 0], fy[1, 1]]]
    q1 = [
        [2 * fx[0, 0] / gamma[0], 2 * fx[0, 1] / gamma[0]],
        [2 * fy[0, 0] / gamma[1], 2 * fy[0, 1] / gamma[1]],
    ]

    return numpy.array(q0), numpy.array(q1), x, q_dx, y, q_dy


def chromaticity(
    ring: Lattice,
    method: Optional[str] = "linopt",
    dpm: Optional[float] = 0.02,
    npoints: Optional[int] = 11,
    order: Optional[int] = 3,
    dp: Optional[float] = 0,
    **kwargs,
):
    r"""Computes the non-linear chromaticities

    This function computes the tunes for the specified momentum offsets.
    Then it fits this data and returns the chromaticity up to the given
    order (npoints>order)

    Parameters:
        ring:       Lattice description
        method:     ``'linopt'`` (dfault) returns the tunes from the
          :py:func:`linopt6` function,

          ``'fft'`` tracks a single particle and computes the tunes with fft,

          ``'laskar'`` tracks a single particle and computes the tunes with
          NAFF.
        dpm:        Maximum momentum deviation
        npoints:    Number of momentum deviations
        order:      Order of the fit
        dp:         Central momentum deviation

    Returns:
        fit (ndarray): Horizontal Vertical polynomial coefficients
                       from ``numpy.polyfit`` of shape (2, order+1)
        dpa: dp array of shape (npoints,)
        qz: tune array of shape (npoints, 2)
    """
    if order == 0:
        return get_chrom(ring, dp=dp, **kwargs)
    elif order > npoints - 1:
        raise ValueError("order should be smaller than npoints-1")
    else:
        dpa = numpy.linspace(-dpm, dpm, npoints)
        qz = []
        for dpi in dpa:
            qz.append(
                get_tune(ring, method=method, dp=dp + dpi, remove_dc=True, **kwargs)
            )
        fit = numpy.polyfit(dpa, qz, order)[::-1]
        fitx = fit[:, 0] / factorial(numpy.arange(order + 1))
        fity = fit[:, 1] / factorial(numpy.arange(order + 1))
        return numpy.array([fitx, fity]), dpa, numpy.array(qz)


def gen_detuning_elem(ring: Lattice, orbit: Optional[Orbit] = None) -> Element:
    """Generates a detuning element

    Parameters:
        ring:           Lattice description
        orbit:          Avoids looking for initial the closed orbit if it is
          already known ((6,) array).

    Returns:
        detuneElem:     Element reproducing the detuning of the ring with
          amplitude and momentum
    """
    if orbit is None:
        orbit, _ = find_orbit(ring)
    lindata0, _, _ = ring.linopt6(get_chrom=False, orbit=orbit)
    xsi = get_chrom(ring.radiation_off(copy=True))
    r0, r1, x, q_dx, y, q_dy = detuning(
        ring.radiation_off(copy=True), xm=1.0e-4, ym=1.0e-4, npoints=3
    )
    nonlin_elem = Element(
        "NonLinear",
        PassMethod="DeltaQPass",
        Betax=lindata0.beta[0],
        Betay=lindata0.beta[1],
        Alphax=lindata0.alpha[0],
        Alphay=lindata0.alpha[1],
        Qpx=xsi[0],
        Qpy=xsi[1],
        A1=r1[0][0],
        A2=r1[0][1],
        A3=r1[1][1],
        T1=-orbit,
        T2=orbit,
    )
    return nonlin_elem


def feed_down_polynomba(
    polb: numpy.ndarray = numpy.zeros(0),
    pola: numpy.ndarray = numpy.zeros(0),
    x0: float = 0,
    y0: float = 0,
    verbose=False,
    debug=False,
) -> dict[str, numpy.ndarray]:
    # check debug and verbose flags only once
    debugprint = print if debug else lambda *a, **k: None
    verboseprint = print if verbose else lambda *a, **k: None

    if pola is None and polb is None:
        print("At least one polynom is needed")

    # verify polynoms length
    maxorda = 0
    maxordb = 0
    if pola is not None:
        maxorda = len(pola)
    if polb is not None:
        maxordb = len(polb)
    maxord = max(maxorda, maxordb)
    debugprint(f"maxord={maxord},maxorda={maxorda},maxordb={maxordb}")

    polbpad = numpy.pad(polb, (0, maxord - maxordb), "constant", constant_values=(0, 0))
    polapad = numpy.pad(pola, (0, maxord - maxorda), "constant", constant_values=(0, 0))
    debugprint(f"polbpad={polbpad}, polapad={polapad}")

    verboseprint(f"polb={polb},pola={pola}")
    verboseprint(f"x0={x0},y0={y0}")
    polasum = numpy.zeros(maxord - 1)
    polbsum = numpy.zeros(maxord - 1)
    debugprint(f"ith=1, first order nothing to do")
    for ith in range(2, maxord + 1):
        polbaux_b, polaaux_b = feed_down_from_nth_order(
            ith, polbpad[ith - 1], x0, y0, magtype="normal"
        )
        polbaux_a, polaaux_a = feed_down_from_nth_order(
            ith, polapad[ith - 1], x0, y0, magtype="skew"
        )
        polbshort = polbaux_b + polbaux_a
        polashort = polaaux_b + polaaux_a
        polbsum = polbsum + numpy.pad(
            polbshort, (0, maxord - ith), "constant", constant_values=(0, 0)
        )
        polasum = polasum + numpy.pad(
            polashort, (0, maxord - ith), "constant", constant_values=(0, 0)
        )
        debugprint(f"ith={ith},polbsum={polbsum},polasum={polasum}")
    poldict = {"PolynomB": polbsum, "PolynomA": polasum}
    verboseprint(f"poldict={poldict}")
    return poldict


def feed_down_from_nth_order(
    nthorder: int,
    nthfieldcomp: float,
    x0: float,
    y0: float,
    magtype: str = "normal",
    debug: bool = False,
    verbose=False,
):
    # print debugging output, equivalent to extra verbose
    debugprint = print if debug else lambda *a, **k: None
    verboseprint = print if verbose else lambda *a, **k: None
    verboseprint(f"nthorder={nthorder},nthfieldcomp={nthfieldcomp},magtype={magtype}")
    verboseprint(f"x0={x0},y0={y0}")
    fakeimag = {0: 1, 1: 0, 2: -1, 3: 0}
    polbaux = numpy.zeros(nthorder - 1)
    polaaux = numpy.zeros(nthorder - 1)
    for k in range(1, nthorder):
        debugprint(f"nthorder={nthorder}, k={k}, nthfieldcomp={nthfieldcomp}")
        ichoosek = comb(nthorder - 1, k)
        debugprint(f"ichoosek={ichoosek}")
        pascalsn = comb(k, numpy.arange(k + 1))
        debugprint(f"pascalsn={pascalsn}")
        powk = numpy.arange(k + 1)
        powkflip = numpy.arange(k, -1, -1)
        debugprint(f"powk={powk}")
        debugprint(f"powkflip={powkflip}")
        recoefs = numpy.array(list(map(lambda x: fakeimag[x], numpy.mod(powk, 4))))
        imcoefs = numpy.array(list(map(lambda x: fakeimag[x], numpy.mod(powk + 3, 4))))
        debugprint(f"recoefs={recoefs}, imcoefs={imcoefs}")
        commonfactor = nthfieldcomp * ichoosek * (-1) ** k
        debugprint(f"commonfactor={commonfactor}")
        repart = recoefs * commonfactor * (pascalsn * (x0**powkflip * y0**powk))
        impart = imcoefs * commonfactor * (pascalsn * (x0**powkflip * y0**powk))
        debugprint(f"repart={repart}")
        debugprint(f"impart={impart}")
        polbaux[nthorder - k - 1] = polbaux[nthorder - k - 1] + repart.sum()
        polaaux[nthorder - k - 1] = polaaux[nthorder - k - 1] + impart.sum()
    verboseprint(f"polbaux={polbaux},polaaux={polaaux}")
    polbout, polaout = polbaux, polaaux
    # skew components swap the imaginary and real feed-down and change the sign
    # of the real part, i.e. j*j = -1
    if magtype == "skew":
        polbout, polaout = -1 * polaaux, polbaux
    return polbout, polaout
