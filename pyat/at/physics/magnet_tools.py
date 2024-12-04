"""Magnet tools."""

from __future__ import annotations

from typing import Sequence

import numpy
from scipy.special import comb

from ..lattice.elements import Element

__all__ = [
    "feeddown_polynomba",
    "feeddown_from_nth_order",
    "feeddown_pol_from_element",
]


def feeddown_pol_from_element(
    ele: Element, verbose: bool = False, **kwargs: dict[str, any]
) -> dict[str, numpy.ndarray]:
    """
    Return the feed down polynoms due to a transverse offset on an element.

    Parameters:
        ele: a ring element.
        verbose: prints additional info. Default: False.
        kwargs:
            offset:(x,y) horizontal and vertical offset in meters.
                Default taken from T2 or T1. Otherwise, zero.

    Returns:
        Dictionary with PolynomB and PolynomA feeddown components.

    Note: Thin lens approximation, T2=-T1.
          Bending angles are ignored.
    """
    # check verbose flags only once
    verboseprint = print if verbose else lambda *a, **k: None
    # handle cases where polynom a and b are not defined
    polb = numpy.array([])
    pola = numpy.array([])
    if hasattr(ele, "PolynomB"):
        polb = ele.PolynomB
    else:
        verboseprint(f"Element {ele.FamName} has no PolynomB")
    if hasattr(ele, "PolynomA"):
        pola = ele.PolynomA
    else:
        verboseprint(f"Element {ele.FamName} has no PolynomA")

    # check if user has set offsets
    xoffset, yoffset = kwargs.get("offset", (0, 0))
    # check if element has T1 and T2. Use one.
    if "offset" not in kwargs:
        if hasattr(ele, "T2"):
            xoffset = ele.T2[0]
            yoffset = ele.T2[2]
        elif hasattr(ele, "T1"):
            xoffset = -ele.T1[0]
            yoffset = -ele.T1[2]
        else:
            verboseprint(f"Element {ele.FamName} has no T1 or T2.")
    verboseprint(f"Using offsets xoffset={xoffset}, yoffset={yoffset}.")
    # Return the polynoms
    return feeddown_polynomba(polb=polb, pola=pola, xoffset=xoffset, yoffset=yoffset)


def feeddown_polynomba(
    polb: Sequence[float],
    pola: Sequence[float],
    xoffset: float = 0,
    yoffset: float = 0,
    verbose: bool = False,
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

    Note:
        If only one polynom is passed, it is assumed to be PolynomB.
    """
    # check verbose flags only once
    verboseprint = print if verbose else lambda *a, **k: None

    # verify polynoms length
    maxorda = len(pola)
    maxordb = len(polb)
    maxord = max(maxorda, maxordb)
    polasum = numpy.zeros(max(0, maxord - 1))
    polbsum = numpy.zeros(max(0, maxord - 1))
    if maxorda == 0 and maxordb == 0:
        verboseprint("Both polynoms are zero.")
    else:
        polbpad = numpy.pad(
            polb, (0, maxord - maxordb), "constant", constant_values=(0, 0)
        )
        polapad = numpy.pad(
            pola, (0, maxord - maxorda), "constant", constant_values=(0, 0)
        )

        verboseprint(f"polb={polb},pola={pola}")
        verboseprint(f"xoffset={xoffset},yoffset={yoffset}")
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
        verbose: print info on input and output parameters.

    Returns:
        Tuple of two numpy arrays with PolynomB and PolynomA from feeddown.
    """
    verboseprint = print if verbose else lambda *a, **k: None
    verboseprint(f"nthorder={nthorder},nthpolcomp={nthpolcomp},poltype={poltype}")
    verboseprint(f"xoffset={xoffset},yoffset={yoffset}")
    fakeimag = {0: 1, 1: 0, 2: -1, 3: 0}
    polbaux = numpy.zeros(nthorder - 1)
    polaaux = numpy.zeros(nthorder - 1)
    if nthpolcomp != 0:
        for kterm in range(1, nthorder):
            ichoosek = comb(nthorder - 1, kterm)
            pascalsn = comb(kterm, numpy.arange(kterm + 1))
            powk = numpy.arange(kterm + 1)
            powkflip = numpy.arange(kterm, -1, -1)
            recoefs = numpy.array([fakeimag[numpy.mod(idx, 4)] for idx in powk])
            imcoefs = numpy.array([fakeimag[numpy.mod(idx + 3, 4)] for idx in powk])
            commonfactor = nthpolcomp * ichoosek * (-1) ** kterm
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
            polbaux[nthorder - kterm - 1] = polbaux[nthorder - kterm - 1] + repart.sum()
            polaaux[nthorder - kterm - 1] = polaaux[nthorder - kterm - 1] + impart.sum()
    verboseprint(f"polbaux={polbaux},polaaux={polaaux}")
    polbout, polaout = polbaux, polaaux
    # skew components swap the imaginary and real feed-down and change the sign
    # of the real part, i.e. j*j = -1
    if poltype == "A":
        polbout, polaout = -1 * polaaux, polbaux
    return polbout, polaout
