"""
Non-linear optics
"""

from __future__ import annotations

from typing import Optional, Sequence

import numpy
from scipy.special import comb

__all__ = [
    "feeddown_polynomba",
    "feeddown_from_nth_order",
]


def feeddown_polynomba(
    polb: Sequence[float],
    pola: Sequence[float],
    xoffset: float = 0,
    yoffset: float = 0,
    verbose: bool = False,
    debug: bool = False,
) -> dict[str, numpy.ndarray]:
    """
    Return the feeddown due to a transverse offset.

    Parameters:
        polb: PolynomB.
        pola: PolynomA.
        xoffset: Default zero. Horizontal offset in meters.
        yoffset: Default zero. Vertical offset in meters.
        verbose: prints additional info.

    Returns:
        Dictionary with PolynomB and PolynomA feeddown components.

    Raises:
        ValueError: if none of the two polynoms is passed.

    Note:
        If only one polynom is passed, it is assumed to be PolynomB.
    """
    # check and verbose flags only once
    verboseprint = print if verbose else lambda *a, **k: None

    # verify polynoms length
    maxorda = len(pola)
    maxordb = len(polb)
    if maxorda == 0 and maxordb == 0:
        raise ValueError("At least one polynom is needed")
    maxord = max(maxorda, maxordb)

    polbpad = numpy.pad(polb, (0, maxord - maxordb), "constant", constant_values=(0, 0))
    polapad = numpy.pad(pola, (0, maxord - maxorda), "constant", constant_values=(0, 0))

    verboseprint(f"polb={polb},pola={pola}")
    verboseprint(f"xoffset={xoffset},yoffset={yoffset}")
    polasum = numpy.zeros(maxord - 1)
    polbsum = numpy.zeros(maxord - 1)
    for ith in range(2, maxord + 1):
        polbaux_b, polaaux_b = feeddown_from_nth_order(
            ith, polbpad[ith - 1], xoffset, yoffset, poltype="B"
        )
        polbaux_a, polaaux_a = feeddown_from_nth_order(
            ith, polapad[ith - 1], xoffset, yoffset, poltype="A"
        )
        polbshort = polbaux_b + polbaux_a
        polashort = polaaux_b + polaaux_a
        polbsum = polbsum + numpy.pad(
            polbshort, (0, maxord - ith), "constant", constant_values=(0, 0)
        )
        polasum = polasum + numpy.pad(
            polashort, (0, maxord - ith), "constant", constant_values=(0, 0)
        )
    poldict = {"PolynomB": polbsum, "PolynomA": polasum}
    verboseprint(f"poldict={poldict}")
    return poldict


def feeddown_from_nth_order(
    nthorder: int,
    nthpolcomp: float,
    xoffset: float,
    yoffset: float,
    poltype: str = "B",
    debug: bool = False,
    verbose: bool = False,
) -> tuple[numpy.ndarray, numpy.ndarray]:
    """
    Return the feeddown polynoms from an nth-order magnet component transverse offset.

    Parameters:
        nthorder: integer order of the magnet component, e.g. 1 for dipole,
            2 for quadrupoles, 3 for sextupoles, etc.
        nthpolcomp: float. nth component of the polynom, i.e. the value of
            PolynomA/B[nthorder-1].
        xoffset: float. Horizontal offset in meters.
        yoffset: float. Vertical offset in meters.
        poltype: Default 'B'. Could be 'B' or 'A', otherwise ignored.
        debug: print info on the feeddown polynom construction.
        verbose: print info on input and output parameters.

    Returns:
        Tuple of two numpy arrays with PolynomB and PolynomA from feeddown.
    """
    # print debugging output, equivalent to extra verbose
    debugprint = print if debug else lambda *a, **k: None
    verboseprint = print if verbose else lambda *a, **k: None
    verboseprint(f"nthorder={nthorder},nthpolcomp={nthpolcomp},poltype={poltype}")
    verboseprint(f"xoffset={xoffset},yoffset={yoffset}")
    fakeimag = {0: 1, 1: 0, 2: -1, 3: 0}
    polbaux = numpy.zeros(nthorder - 1)
    polaaux = numpy.zeros(nthorder - 1)
    for kterm in range(1, nthorder):
        debugprint(f"nthorder={nthorder}, kterm={kterm}, nthpolcomp={nthpolcomp}")
        ichoosek = comb(nthorder - 1, kterm)
        debugprint(f"ichoosek={ichoosek}")
        pascalsn = comb(kterm, numpy.arange(kterm + 1))
        debugprint(f"pascalsn={pascalsn}")
        powk = numpy.arange(kterm + 1)
        powkflip = numpy.arange(kterm, -1, -1)
        debugprint(f"powk={powk}")
        debugprint(f"powkflip={powkflip}")
        recoefs = numpy.array([fakeimag[numpy.mod(idx, 4)] for idx in powk])
        imcoefs = numpy.array([fakeimag[numpy.mod(idx + 3, 4)] for idx in powk])
        debugprint(f"recoefs={recoefs}, imcoefs={imcoefs}")
        commonfactor = nthpolcomp * ichoosek * (-1) ** kterm
        debugprint(f"commonfactor={commonfactor}")
        repart = (
            recoefs
            * commonfactor
            * (pascalsn * (xoffset**powkflip * yoffset**powk))
        )
        impart = (
            imcoefs
            * commonfactor
            * (pascalsn * (xoffset**powkflip * yoffset**powk))
        )
        debugprint(f"repart={repart}")
        debugprint(f"impart={impart}")
        polbaux[nthorder - kterm - 1] = polbaux[nthorder - kterm - 1] + repart.sum()
        polaaux[nthorder - kterm - 1] = polaaux[nthorder - kterm - 1] + impart.sum()
    verboseprint(f"polbaux={polbaux},polaaux={polaaux}")
    polbout, polaout = polbaux, polaaux
    # skew components swap the imaginary and real feed-down and change the sign
    # of the real part, i.e. j*j = -1
    if poltype == "A":
        polbout, polaout = -1 * polaaux, polbaux
    return polbout, polaout