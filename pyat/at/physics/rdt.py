"""Resonance Driving Terms"""
from __future__ import annotations

import numpy as np
import multiprocessing
from functools import partial
from ..lattice import Lattice, Refpts, Multipole, All
from enum import Enum
from collections.abc import Container


__all__ = ["get_rdts", "RDTType"]


class RDTType(Enum):
    """Enum class for RDT type"""

    ALL = 0  #: all available RDTs
    FOCUSING = 1  #: Normal quadrupole RDTs
    COUPLING = 2  #: Linear coupling RDTs
    CHROMATIC = 3  #: Chromatic RDTs
    #: Geometric RDTs from sextupoles
    GEOMETRIC1 = 4
    #: Geometric RDTs from octupoles
    #: optionally includes the second order contribution of sextupoles
    GEOMETRIC2 = 5
    #: Amplitude detuning coefficients
    #: optionally includes the second order contribution of sextupoles
    TUNESHIFT = 6
    #: Chromatic RDTs from octupoles
    #: optionally includes the second order contribution of sextupoles
    CHROMATIC2 = 7


def _get_polynom(elem, attr, index):
    try:
        val = getattr(elem, attr)[index]
        return val
    except (IndexError, AttributeError):
        return 0


def _compute_pf(tune, nperiods):
    """This uses the formula Sum(x^k, k=1->p) = x(x^p-1)/(x-1)"""
    pf = np.ones((9, 9), dtype=complex)
    if nperiods != 1:
        for i in range(9):
            for j in range(9):
                a1 = np.pi * 2 * (tune[0] * (i - 4) + tune[1] * (j - 4))
                a2 = a1 / nperiods
                pf[i][j] = (np.exp(1j * a1) - 1.0) / (np.exp(1j * a2) - 1.0)
    return pf


def _computedrivingterms(
    s,
    betax,
    betay,
    phix,
    phiy,
    etax,
    a2l,
    b2l,
    b3l,
    b4l,
    tune,
    rdttype,
    nperiods,
    refpt,
    second_order,
    periodic_factor,
):
    """
    Original implementation from ELEGANT
    """

    def pf(i, j):
        return periodic_factor[4 + i][4 + j]

    rdts = {}
    rdts2 = {}
    rbetax = np.sqrt(betax)
    rbetay = np.sqrt(betay)
    px = np.exp(1j * phix)
    py = np.exp(1j * phiy)

    rdts["refpts"] = refpt
    mask_a2l = np.absolute(a2l) > 1.0e-6
    mask_b2l = np.absolute(b2l) > 1.0e-6
    mask_b3l = np.absolute(b3l) > 1.0e-6
    mask_b4l = np.absolute(b4l) > 1.0e-6

    if (RDTType.FOCUSING in rdttype) or (RDTType.ALL in rdttype):
        b2lm = b2l[mask_b2l]
        betaxm = betax[mask_b2l]
        betaym = betay[mask_b2l]
        pxm = px[mask_b2l]
        pym = py[mask_b2l]
        rdts["h20000"] = 1 / 8 * pf(2, 0) * sum(b2lm * betaxm * pxm * pxm)
        rdts["h00200"] = -1 / 8 * pf(0, 2) * sum(b2lm * betaym * pym * pym)

    if (RDTType.COUPLING in rdttype) or (RDTType.ALL in rdttype):
        a2lm = a2l[mask_a2l]
        rbetaxm = rbetax[mask_a2l]
        rbetaym = rbetay[mask_a2l]
        pxm = px[mask_a2l]
        pym = py[mask_a2l]
        rdts["h10010"] = 1 / 4 * pf(1, -1) * sum(a2lm * rbetaxm * rbetaym * pxm / pym)
        rdts["h10100"] = 1 / 4 * pf(1, 1) * sum(a2lm * rbetaxm * rbetaym * pxm * pym)

    if (RDTType.CHROMATIC in rdttype) or (RDTType.ALL in rdttype):
        mask_b23l = mask_b2l | mask_b3l
        b2lm = b2l[mask_b23l]
        b3lm = b3l[mask_b23l]
        betaxm = betax[mask_b23l]
        rbetaxm = rbetax[mask_b23l]
        betaym = betay[mask_b23l]
        etaxm = etax[mask_b23l]
        pxm = px[mask_b23l]
        pym = py[mask_b23l]
        rdts["h11001"] = nperiods * sum(b3lm * betaxm * etaxm / 2 - b2lm * betaxm / 4)
        rdts["h00111"] = nperiods * sum(b2lm * betaym / 4 - b3lm * betaym * etaxm / 2)
        # fmt: off
        rdts["h20001"] = (
            1 / 2 * pf(2, 0)
            * sum((b3lm * betaxm * etaxm / 2 - b2lm * betaxm / 4) * pxm * pxm)
        )
        rdts["h00201"] = (
            1 / 2 * pf(0, 2)
            * sum((b2lm * betaym / 4 - b3lm * betaym * etaxm / 2) * pym * pym)
        )
        rdts["h10002"] = 1 / 2 * pf(1, 0 ) * sum(
            (b3lm * rbetaxm * etaxm * etaxm - b2lm * rbetaxm * etaxm)
            * pxm
        )
        # fmt: on

    if (RDTType.GEOMETRIC1 in rdttype) or (RDTType.ALL in rdttype):
        b3lm = b3l[mask_b3l]
        betaxm = betax[mask_b3l]
        betaym = betay[mask_b3l]
        rbetaxm = rbetax[mask_b3l]
        pxm = px[mask_b3l]
        pym = py[mask_b3l]
        rdts["h21000"] = 1 / 8 * pf(1, 0) * sum(b3lm * rbetaxm * betaxm * pxm)
        rdts["h30000"] = (
            1 / 24 * pf(3, 0) * sum(b3lm * rbetaxm * betaxm * pxm * pxm * pxm)
        )
        rdts["h10110"] = -1 / 4 * pf(1, 0) * sum(b3lm * rbetaxm * betaym * pxm)
        rdts["h10020"] = (
            -1 / 8 * pf(1, -2) * sum(b3lm * rbetaxm * betaym * pxm * np.conj(pym * pym))
        )
        rdts["h10200"] = (
            -1 / 8 * pf(1, 2) * sum(b3lm * rbetaxm * betaym * pxm * pym * pym)
        )

    if (RDTType.TUNESHIFT in rdttype) or (RDTType.ALL in rdttype):
        b4lm = b4l[mask_b4l]
        betaxm = betax[mask_b4l]
        betaym = betay[mask_b4l]
        rdts["dnux_dJx"] = 3 / (8 * np.pi) * nperiods * sum(b4lm * betaxm * betaxm)
        rdts["dnux_dJy"] = -3 / (4 * np.pi) * nperiods * sum(b4lm * betaxm * betaym)
        rdts["dnuy_dJy"] = 3 / (8 * np.pi) * nperiods * sum(b4lm * betaym * betaym)

    if (RDTType.GEOMETRIC2 in rdttype) or (RDTType.ALL in rdttype):
        b4lm = b4l[mask_b4l]
        betaxm = betax[mask_b4l]
        betaym = betay[mask_b4l]
        pxm = px[mask_b4l]
        pym = py[mask_b4l]
        rdts["h22000"] = 3 / 32 * nperiods * sum(b4lm * betaxm * betaxm)
        rdts["h11110"] = -3 / 8 * nperiods * sum(b4lm * betaxm * betaym)
        rdts["h00220"] = 3 / 32 * nperiods * sum(b4lm * betaym * betaym)
        rdts["h31000"] = 1 / 16 * pf(2, 0) * sum(b4lm * betaxm * betaxm * pxm * pxm)
        rdts["h40000"] = (
            1 / 64 * pf(4, 0) * sum(b4lm * betaxm * betaxm * pxm * pxm * pxm * pxm)
        )
        rdts["h20110"] = -3 / 16 * pf(2, 0) * sum(b4lm * betaxm * betaym * pxm * pxm)
        rdts["h11200"] = -3 / 16 * pf(0, 2) * sum(b4lm * betaxm * betaym * pym * pym)
        # fmt: off
        rdts["h20020"] = (
            -3 / 32 * pf(2, -2)
            * sum(b4lm * betaxm * betaym * pxm * pxm * np.conj(pym * pym))
        )
        # fmt: on
        rdts["h20200"] = (
            -3 / 32 * pf(2, 2) * sum(b4lm * betaxm * betaym * pxm * pxm * pym * pym)
        )
        rdts["h00310"] = 1 / 16 * pf(0, 2) * sum(b4lm * betaym * betaym * pym * pym)
        rdts["h00400"] = (
            1 / 64 * pf(0, 4) * sum(b4lm * betaym * betaym * pym * pym * pym * pym)
        )

    if (RDTType.CHROMATIC2 in rdttype) or (RDTType.ALL in rdttype):
        b4lm = b4l[mask_b4l]
        betaxm = betax[mask_b4l]
        betaym = betay[mask_b4l]
        etaxm = etax[mask_b4l]
        pxm = px[mask_b4l]
        pym = py[mask_b4l]
        rbetaxm = np.sqrt(betaxm)
        etaxm2 = etaxm * etaxm
        cpym = np.conj(pym)
        rdts["h21001"] = 3 / 8 * pf(1, 0) * sum(b4lm * etaxm * rbetaxm * betaxm * pxm)
        rdts["h30001"] = (
            1 / 8 * pf(3, 0) * sum(b4lm * etaxm * rbetaxm * betaxm * pxm * pxm * pxm)
        )
        # fmt: off
        rdts["h10021"] = (
            -3 / 8 * pf(1, -2)
            * sum(b4lm * etaxm * rbetaxm * betaym * pxm * cpym * cpym)
        )
        # fmt: on
        rdts["h10111"] = -3 / 4 * pf(1, 0) * sum(b4lm * etaxm * rbetaxm * betaym * pxm)
        rdts["h10201"] = (
            -3 / 8 * pf(1, 2) * sum(b4lm * etaxm * rbetaxm * betaym * pxm * pym * pym)
        )
        rdts["h11002"] = 3 / 4 * nperiods * sum(b4lm * etaxm2 * betaxm)
        rdts["h20002"] = 3 / 8 * pf(2, 0) * sum(b4lm * etaxm2 * betaxm * pxm * pxm)
        rdts["h00112"] = -3 / 4 * nperiods * sum(b4lm * etaxm2 * betaym)
        rdts["h00202"] = -3 / 8 * pf(0, 2) * sum(b4lm * etaxm2 * betaym * pym * pym)
        rdts["h10003"] = 1 / 2 * pf(1, 0) * sum(b4lm * etaxm2 * etaxm * rbetaxm * pxm)
        rdts["h00004"] = 1 / 4 * nperiods * sum(b4lm * etaxm2 * etaxm2)

    if second_order:
        assert nperiods == 1, "Second order available only for nperiods=1"

        if (RDTType.GEOMETRIC2 in rdttype) or (RDTType.ALL in rdttype):
            nelem = sum(mask_b3l)
            sm = s[mask_b3l]
            b3lm = b3l[mask_b3l]
            betaxm = betax[mask_b3l]
            betaym = betay[mask_b3l]
            rbetaxm = rbetax[mask_b3l]
            rbbetaxm = betaxm * rbetaxm
            phixm = phix[mask_b3l]
            phiym = phiy[mask_b3l]
            pxm = px[mask_b3l]
            pxm2 = pxm * pxm
            pxm3 = pxm2 * pxm
            cpxm = np.conj(pxm)
            cpxm3 = np.conj(pxm3)
            pym = py[mask_b3l]
            pym2 = pym * pym
            cpym2 = np.conj(pym2)
            cpxym2 = np.conj(pxm * pym2)
            rdts2.update(
                {
                    "h22000": 0.0,
                    "h31000": 0.0,
                    "h11110": 0.0,
                    "h11200": 0.0,
                    "h40000": 0.0,
                    "h20020": 0.0,
                    "h20110": 0.0,
                    "h20200": 0.0,
                    "h00220": 0.0,
                    "h00310": 0.0,
                    "h00400": 0.0,
                }
            )
            # fmt: off
            bsign = np.array([np.sign(sm[i] - sm) * 1j * b3lm[i] * b3lm
                              for i in range(nelem)])
            rbbx = np.array([bsign[i] * rbbetaxm[i] * rbbetaxm
                             for i in range(nelem)])
            rbxy = np.array([bsign[i] * rbetaxm[i] * rbetaxm * betaym[i]
                             for i in range(nelem)])
            ppxm = np.array([pxm[i] * pxm for i in range(nelem)])

            rdts2["h22000"] += (1.0 / 64) * np.sum([
                rbbx[i] * (pxm3[i] * cpxm3 + 3 * pxm[i] * cpxm)
                for i in range(nelem)]
            )
            rdts2["h31000"] += (1.0 / 32) * np.sum([
                rbbx[i] * pxm3[i] * cpxm
                for i in range(nelem)]
            )
            t1 = np.array([np.conj(pxm[i]) * pxm for i in range(nelem)])
            t2 = np.conj(t1)
            rdts2["h11110"] += (1.0 / 16) * np.sum([
                rbxy[i] * (betaxm * (t1[i] - t2[i]) + betaym * pym2[i]
                           * cpym2 * (t2[i] + t1[i]))
                for i in range(nelem)]
            )
            t1 = np.array([np.exp(-1j * (phixm[i] - phixm))
                           for i in range(nelem)])
            t2 = np.conj(t1)
            rdts2["h11200"] += (1.0 / 32) * np.sum([
                rbxy[i] * np.exp(1j * (2 * phiym[i]))
                * (betaxm * (t1[i] - t2[i]) + 2 * betaym * (t2[i] + t1[i]))
                for i in range(nelem)]
            )
            rdts2["h40000"] += (1.0 / 64) * np.sum([
                rbbx[i] * pxm3[i] * pxm for i in range(nelem)]
            )
            rdts2["h20020"] += (1.0 / 64) * np.sum([
                rbxy[i] * (betaxm * cpxym2[i] * pxm3 - (betaxm + 4 * betaym)
                           * ppxm[i] * cpym2[i])
                for i in range(nelem)]
            )
            rdts2["h20110"] += (1.0 / 32) * np.sum([
                rbxy[i] * (betaxm * (cpxm[i] * pxm3 - ppxm[i])
                           + 2 * betaym * ppxm[i] * pym2[i] * cpym2)
                for i in range(nelem)]
            )
            rdts2["h20200"] += (1.0 / 64) * np.sum([
                rbxy[i] * (betaxm * cpxm[i] * pxm3 * pym2[i]
                           - (betaxm - 4 * betaym)
                           * ppxm[i] * pym2[i])
                for i in range(nelem)]
            )
            rdts2["h00220"] += (1.0 / 64) * np.sum([
                rbxy[i] * betaym * (pxm[i] * pym2[i] * cpxym2 + 4 * pxm[i] * cpxm
                                    - np.conj(pxm[i] * pym2) * pxm * pym2[i])
                for i in range(nelem)]
            )
            rdts2["h00310"] += (1.0 / 32) * np.sum([
                rbxy[i] * betaym * pym2[i] * (pxm[i] * cpxm - pxm * cpxm[i])
                for i in range(nelem)]
            )
            rdts2["h00400"] += (1.0 / 64) * np.sum([
                rbxy[i] * betaym * pxm[i] * cpxm * pym2[i] * pym2
                for i in range(nelem)]
            )
            # fmt: on

        if (RDTType.CHROMATIC2 in rdttype) or (RDTType.ALL in rdttype):
            mask_b23l = mask_b2l | mask_b3l
            nelem = sum(mask_b23l)
            b2lm = b2l[mask_b23l]
            b3lm = b3l[mask_b23l]
            betxm = betax[mask_b23l]
            rbetxm = rbetax[mask_b23l]
            betym = betay[mask_b23l]
            etaxm = etax[mask_b23l]
            pxm = px[mask_b23l]
            pym = py[mask_b23l]
            betxm3o2 = betxm * rbetxm
            pxm2 = pxm * pxm
            pxm3 = pxm2 * pxm
            cpxm = np.conj(pxm)
            cpxm2 = np.conj(pxm2)
            pym2 = pym * pym
            cpym = np.conj(pym)
            cpym2 = np.conj(pym2)
            rdts2.update(
                {
                    "h21001": 0.0,
                    "h30001": 0.0,
                    "h10021": 0.0,
                    "h10111": 0.0,
                    "h10201": 0.0,
                    "h11002": 0.0,
                    "h20002": 0.0,
                    "h00112": 0.0,
                    "h00202": 0.0,
                    "h10003": 0.0,
                    "h00004": 0.0,
                }
            )
            # The variable imag_and_sign multiplies by -1 the expression in
            # ANL/APS/LS-330 March 10, 2012. Chun-xi Wang. Eq (46)
            # in order to follow AT sign convention.
            imag_and_sign = 1j * (np.tri(nelem, nelem, -1) - 1 + np.tri(nelem))
            b2b2lm = np.array([imag_and_sign[i] * b2lm[i] * b2lm for i in range(nelem)])
            b3b2lm = np.array([imag_and_sign[i] * b3lm[i] * b2lm for i in range(nelem)])
            b2b3lm = np.array([imag_and_sign[i] * b2lm[i] * b3lm for i in range(nelem)])
            b3b3lm = np.array([imag_and_sign[i] * b3lm[i] * b3lm for i in range(nelem)])
            bx3o2bx = np.array([betxm3o2[i] * betxm for i in range(nelem)])
            bxbx = np.array([betxm[i] * betxm for i in range(nelem)])
            byby = np.array([betym[i] * betym for i in range(nelem)])
            bxbx3o2etax = np.array(
                [betxm[i] * betxm3o2 * etaxm[i] for i in range(nelem)]
            )
            rbxbxby = np.array([rbetxm[i] * betxm * betym[i] for i in range(nelem)])
            rbxbyby = np.array([rbetxm[i] * betym[i] * betym for i in range(nelem)])
            bxrbxbyetax = np.array(
                [betxm[i] * rbetxm * betym * etaxm[i] for i in range(nelem)]
            )
            rbxbybyetax = np.array(
                [rbetxm * betym[i] * betym * etaxm[i] for i in range(nelem)]
            )
            rbxbx3o2etax = np.array(
                [rbetxm[i] * betxm3o2 * etaxm[i] for i in range(nelem)]
            )
            rbxrbxbyetax = np.array(
                [rbetxm[i] * rbetxm * betym * etaxm[i] for i in range(nelem)]
            )
            bxrbxetaxetax = np.array(
                [betxm[i] * rbetxm * etaxm[i] * etaxm for i in range(nelem)]
            )
            rbxbxetax = np.array([rbetxm[i] * betxm * etaxm[i] for i in range(nelem)])
            rbxrbxetaxetax = np.array(
                [rbetxm[i] * rbetxm * etaxm[i] * etaxm for i in range(nelem)]
            )

            # fmt: off
            rdts2["h21001"] += (1.0 / 32) * np.sum([
                - 1 * b3b2lm[i] * bx3o2bx[i] * (pxm[i] + pxm3[i] * cpxm2 - 2 * cpxm[i] * pxm2)
                - 2 * b3b3lm[i] * bxbx3o2etax[i] * (pxm - 2 * pxm2[i] * cpxm + cpxm2[i] * pxm3)
                for i in range(nelem)]
            )
            rdts2["h30001"] += (1.0 / 32) * np.sum([
                - 1 * b3b2lm[i] * bx3o2bx[i] * (pxm3[i] - pxm[i] * pxm2)
                - 2 * b3b3lm[i] * bxbx3o2etax[i] * (pxm3 - pxm2[i] * pxm)
                for i in range(nelem)]
            )
            rdts2["h10021"] += (1.0 / 32) * np.sum([
                + 1 * b3b2lm[i] * rbxbxby[i] * (pxm[i] * cpym2[i] - cpxm[i] * pxm2 * cpym2[i])
                + 2 * b3b2lm[i] * rbxbyby[i] * (pxm[i] * cpym2[i] - pxm[i] * cpym2)
                + 2 * b3b3lm[i] * bxrbxbyetax[i] * (pxm * cpym2 - pxm2[i] * cpxm * cpym2)
                - 4 * b3b3lm[i] * rbxbybyetax[i] * (pxm * cpym2[i] - pxm * cpym2)
                for i in range(nelem)]
            )
            rdts2["h10111"] += (1.0 / 16) * np.sum([
                + 1 * b3b2lm[i] * rbxbxby[i] * (pxm[i] - cpxm[i] * pxm2)
                + 1 * b3b2lm[i] * rbxbyby[i] * (pxm[i] * cpym2[i] * pym2 - pxm[i] * pym2[i] * cpym2)
                + 2 * b3b3lm[i] * bxrbxbyetax[i] * (pxm - pxm2[i] * cpxm)
                - 2 * b3b3lm[i] * rbxbybyetax[i] * (pxm * cpym2[i] * pym2 - pxm * pym2[i] * cpym2)
                for i in range(nelem)]
            )
            rdts2["h10201"] += (1.0 / 32) * np.sum([
                + 1 * b3b2lm[i] * rbxbxby[i] * (pxm[i] * pym2[i] - cpxm[i] * pxm2 * pym2[i])
                - 2 * b3b2lm[i] * rbxbyby[i] * (pxm[i] * pym2[i] - pxm[i] * pym2)
                + 2 * b3b3lm[i] * bxrbxbyetax[i] * (pxm * pym2 - pxm2[i] * cpxm * pym2)
                + 4 * b3b3lm[i] * rbxbybyetax[i] * (pxm * pym2[i] - pxm * pym2)
                for i in range(nelem)]
            )
            rdts2["h11002"] += (1.0 / 16) * np.sum([
                + 1 * bxbx[i] * ((b2b2lm[i] - 2 * b3b2lm[i] * etaxm[i] + 4 * b3b3lm[i] * etaxm[i] * etaxm) * pxm2[i] * cpxm2
                                 + 2 * b3b2lm[i] * etaxm[i] * cpxm2[i] * pxm2)
                + 2 * rbxbx3o2etax[i] * (b3b3lm[i] * etaxm[i] - b2b3lm[i]) * (pxm[i] * cpxm - cpxm[i] * pxm)
                for i in range(nelem)]
            )
            rdts2["h20002"] += (1.0 / 16) * np.sum([
                + bxbx[i] * ((b2b2lm[i] - 2 * b3b2lm[i] * etaxm[i] + 4 * b3b3lm[i] * etaxm[i] * etaxm) * pxm2[i]
                               + 2 * b3b2lm[i] * etaxm[i] * pxm2)
                + rbxbx3o2etax[i] * (b3b3lm[i] * etaxm[i] - b2b3lm[i]) * (pxm[i] * pxm - cpxm[i] * pxm3)
                    for i in range(nelem)]
            )
            rdts2["h00112"] += (1.0 / 16) * np.sum([
                + 1 * byby[i] * ((b2b2lm[i] - 2 * b3b2lm[i] * etaxm[i] + 4 * b3b3lm[i] * etaxm[i] * etaxm) * pym2[i] * cpym2
                                 + 2 * b3b2lm[i] * etaxm[i] * cpym2[i] * pym2)
                - 2 * rbxrbxbyetax[i] * (b3b3lm[i] * etaxm[i] - b2b3lm[i]) * (pxm[i] * cpxm - cpxm[i] * pxm)
                for i in range(nelem)]
            )
            rdts2["h00202"] += (1.0 / 16) * np.sum([
                + byby[i] * ((b2b2lm[i] - 2 * b3b2lm[i] * etaxm[i] + 4 * b3b3lm[i] * etaxm[i] * etaxm) * pym2[i]
                                + 2 * b3b2lm[i] * etaxm[i] * pym2)
                - rbxrbxbyetax[i] * (b3b3lm[i] * etaxm[i] - b2b3lm[i]) * (pxm[i] * cpxm * pym2 - cpxm[i] * pxm * pym2)
                for i in range(nelem)]
            )
            rdts2["h10003"] += (1.0 / 8) * np.sum([
                + 2 * b3b2lm[i] * bxrbxetaxetax[i] * (pxm - pxm2[i] * cpxm)
                + 1 * rbxbxetax[i] * (b2b2lm[i] - b3b2lm[i] * etaxm[i] + 2 * b3b3lm[i] * etaxm[i] * etaxm) * (pxm[i] - cpxm[i] * pxm2)
                for i in range(nelem)]
            )
            rdts2["h00004"] += (1.0 / 4) * np.sum([
                + rbxrbxetaxetax[i] * ((b2b2lm[i] - b3b2lm[i] * etaxm[i] + b3b3lm[i] * etaxm[i] * etaxm) * pxm[i] * cpxm
                                     + b3b2lm[i] * etaxm[i] * cpxm[i] * pxm)
                for i in range(nelem)]
            )
            # fmt: on

        if (RDTType.TUNESHIFT in rdttype) or (RDTType.ALL in rdttype):
            nelem = sum(mask_b3l)
            b3lm = b3l[mask_b3l]
            betaxm = betax[mask_b3l]
            betaym = betay[mask_b3l]
            phixm = phix[mask_b3l]
            phiym = phiy[mask_b3l]
            nux = tune[0]
            nuy = tune[1]
            # fmt: off
            b3f = [b3lm[i] * b3lm / np.pi for i in range(nelem)]
            sqbx = [np.sqrt(betaxm[i] * betaxm) for i in range(nelem)]
            bbx = [betaxm[i] * betaxm for i in range(nelem)]
            bby = [betaym[i] * betaym for i in range(nelem)]
            dphix = [abs(phixm[i] - phixm) for i in range(nelem)]
            dphiy = [abs(phiym[i] - phiym) for i in range(nelem)]
            rdts["dnux_dJx"] += np.sum([
                -b3f[i] / 16 * sqbx[i] * bbx[i]
                * (3 * np.cos(dphix[i] - np.pi * nux) / np.sin(np.pi * nux)
                   + np.cos(3 * dphix[i] - 3 * np.pi * nux)
                   / np.sin(3 * np.pi * nux))
                for i in range(nelem)]
            )
            rdts["dnux_dJy"] += np.sum([
                b3f[i] / 8 * sqbx[i] * betaym[i]
                * (2 * betaxm * np.cos(dphix[i] - np.pi * nux)
                   / np.sin(np.pi * nux)
                   - betaym * np.cos(dphix[i] + 2 * dphiy[i]
                                     - np.pi * (nux + 2 * nuy))
                   / np.sin(np.pi * (nux + 2 * nuy))
                   + betaym * np.cos(dphix[i] - 2 * dphiy[i]
                                     - np.pi * (nux - 2 * nuy))
                   / np.sin(np.pi * (nux - 2 * nuy)))
                for i in range(nelem)]
            )
            rdts["dnuy_dJy"] += np.sum([
                -b3f[i] / 16 * sqbx[i] * bby[i]
                * (4 * np.cos(dphix[i] - np.pi * nux)
                   / np.sin(np.pi * nux)
                   + np.cos(dphix[i] + 2 * dphiy[i]
                            - np.pi * (nux + 2 * nuy))
                   / np.sin(np.pi * (nux + 2 * nuy))
                   + np.cos(dphix[i] - 2 * dphiy[i]
                            - np.pi * (nux - 2 * nuy))
                   / np.sin(np.pi * (nux - 2 * nuy)))
                for i in range(nelem)]
            )
            # fmt: on
    return rdts, rdts2


def _get_rdtlist(
    idx_mag,
    beta,
    etax,
    phi,
    avemu,
    mu,
    smag,
    sall,
    pols,
    rdt_type,
    tune,
    nperiods,
    second_order,
    refpts,
):
    rdtlist = []
    rdtlist2 = []
    pf = _compute_pf(tune, nperiods)
    for ii in refpts:
        start_idx = sum(idx_mag < ii)
        beta_rot = np.roll(beta, -start_idx, axis=0)
        etax_rot = np.roll(etax, -start_idx)
        phi_rot = np.roll(phi, -start_idx, axis=0) - avemu[ii]
        if start_idx > 0:
            phi_rot[-start_idx:] += mu
        s_rot = np.roll(smag, -start_idx) - sall[ii]
        s_rot[-start_idx:] += sall[-1]
        pols_rot = np.roll(pols, -start_idx, axis=0)
        rdt, rdt2 = _computedrivingterms(
            s_rot,
            beta_rot[:, 0],
            beta_rot[:, 1],
            phi_rot[:, 0],
            phi_rot[:, 1],
            etax_rot,
            pols_rot[:, 0],
            pols_rot[:, 1],
            pols_rot[:, 2],
            pols_rot[:, 3],
            tune,
            rdt_type,
            nperiods,
            ii,
            second_order,
            pf,
        )
        rdtlist.append(rdt)
        rdtlist2.append(rdt2)
    return rdtlist, rdtlist2


def get_rdts(
    ring: Lattice,
    refpts: Refpts,
    rdt_type: Container[RDTType] | RDTType,
    second_order: bool = False,
    use_mp: bool = False,
    pool_size: int = None,
):
    """Get the lattice Resonance Driving Terms

    :py:func:`get_rdts` computes the ring RDTs based on the original implementation
    from ELEGANT. For consistency, pyAT keeps the sign convention of the AT MATLAB
    interface.

    Refs [#]_ and [#]_ were used to calculate the magnitude of first and second order
    rdts and Ref. [#]_ for the focusing rdt, however, different sign conventions are
    used along the references and in the original implementation. To overcome this
    issue, first and second order rdts are provided as a total and also separately
    so the user can redefine these conventions if needed.

    The resonance base of :py:obj:`~RDTType.GEOMETRIC2` rdts is explained in Eq. (54)
    of Ref [#]_.

    The periodicity property of the lattice is automatically taken into account in the
    rdt calculation, however the calculation of the second order contribution of
    sextupoles to the :py:obj:`~RDTType.GEOMETRIC2` and :py:obj:`~RDTType.TUNESHIFT`
    RDT types can only be derived for periodicity=1 (i.e. full ring provided).

    Parameters:
        ring: :py:class:`.Lattice` object
        refpts: Element refpts at which the RDTs are calculated
        rdt_type: Type of RDTs to be calculated.
          Possible RDT types are:

          * :py:obj:`RDTType.ALL`: all available RDTs
          * :py:obj:`RDTType.FOCUSING`: Normal quadrupole RDTs
          * :py:obj:`RDTType.COUPLING`: Linear coupling RDTs from skew quadrupoles
          * :py:obj:`RDTType.CHROMATIC`: Chromatic RDTs from sextupoles and normal
            quadrupoles
          * :py:obj:`RDTType.CHROMATIC2`: Chromatic RDTs from octupoles. The second
            order  contribution of sextupoles is added when *second_order* is True
          * :py:obj:`RDTType.GEOMETRIC1`: Geometric RDTs from sextupoles
          * :py:obj:`RDTType.GEOMETRIC2`: Geometric RDTs from octupoles. The second
            order contribution of sextupoles is added when *second_order* is True
          * :py:obj:`RDTType.TUNESHIFT`: Amplitude detuning coefficients. The second
            order contribution of sextupoles is added when *second_order* is True

    Keyword Args:
        second_order (bool): Compute second order terms. Computation is significantly
          longer using this method,
        use_mp (bool):       Activate parallel calculation,
        pool_size (int):     Number of processes used for parallelization.

    Returns:
        rdts (complex): rdt data at refpts,
        rdts2 (complex): contribution from sextupole second order terms
          Available only for :py:obj:`~RDTType.GEOMETRIC2` and
          :py:obj:`~RDTType.TUNESHIFT` terms,
        rdttot (complex): total rdts.

    **rdts** is a :py:class:`record array <numpy.recarray>` with fields:

    for :py:obj:`~RDTType.FOCUSING`:
        `h20000`, `h00200`

    for :py:obj:`~RDTType.COUPLING`:
        `h10010`, `h10100`

    for :py:obj:`~RDTType.CHROMATIC`:
        `h11001`, `h00111`, `h20001`, `h00201`, `h10002`

    for :py:obj:`~RDTType.GEOMETRIC1`:
        `h21000`, `h30000`, `h10110`, `h10020`, `h10200`

    for :py:obj:`~RDTType.GEOMETRIC2`:
        `h22000`, `h11110`, `h00220`, `h31000`, `h40000`, `h20110`
        `h11200`, `h20020`, `h20200`, `h00310`, `h00400`

    for :py:obj:`~RDTType.CHROMATIC2`:
        `h21001`, `h30001`, `h10021`, `h10111`, `h10201`, `h11002`
        `h20002`, `h00112`, `h00202`, `h10003`, `h00004`

    for :py:obj:`~RDTType.TUNESHIFT`:
        `dnux_dJx`, `dnux_dJy`, `dnuy_dJy`

    Example:

        >>> get_rdts(ring, reftps, [RDTType.COUPLING, RDTType.CHROMATIC])

    References:
        .. [#] J.Bengtsson, SLS Note 9 / 97, March 7, 1997, with corrections per
           W.Guo (NSLS)

        .. [#] Revised to follow C.X.Wang AOP - TN - 2009 - 020 for second - order
           terms

        .. [#] A.Franchi et al. arxiv 1711.06589, PRAB 17.074001

        .. [#] Chunxi Wang and Alex Chao. Notes on Lie algebraic analysis of achromats.
           SLAC/AP-100. Jan/1995
    """
    if not isinstance(rdt_type, Container):
        rdt_type = {rdt_type}
    nperiods = ring.periodicity
    if second_order:
        assert nperiods == 1, "Second order available only for ring.periodicity=1"

    refpts = ring.uint32_refpts(refpts)
    idx_mag = ring.get_uint32_index(Multipole)
    lo, avebeta, avemu, avedisp, *_ = ring.avlinopt(refpts=All)

    sall = ring.get_s_pos(All)
    smag = sall[idx_mag]
    beta = avebeta[idx_mag]
    etax = avedisp[idx_mag, 0]
    phi = avemu[idx_mag]
    pols = [
        [
            _get_polynom(e, "PolynomA", 1) * e.Length,
            _get_polynom(e, "PolynomB", 1) * e.Length,
            _get_polynom(e, "PolynomB", 2) * e.Length,
            _get_polynom(e, "PolynomB", 3) * e.Length,
        ]
        for e in ring[idx_mag]
    ]

    mu = lo[-1].mu
    tune = mu / 2.0 / np.pi * nperiods
    fun = partial(
        _get_rdtlist,
        idx_mag,
        beta,
        etax,
        phi,
        avemu,
        mu,
        smag,
        sall,
        pols,
        rdt_type,
        tune,
        nperiods,
        second_order,
    )
    if use_mp:
        if pool_size is None:
            pool_size = min(len(refpts), multiprocessing.cpu_count())
        ctx = multiprocessing.get_context()
        refpts = np.array_split(refpts, pool_size)
        with ctx.Pool(pool_size) as pool:
            results = pool.map(fun, refpts)
            rdtlist, rdtlist2 = zip(*results)
            rdtlist = np.concatenate(rdtlist)
            rdtlist2 = np.concatenate(rdtlist2)
    else:
        rdtlist, rdtlist2 = fun(refpts)
    rdts = {}
    rdts2 = {}
    rdttot = {}
    keylist = rdtlist[0].keys()
    for k in keylist:
        val = [rdt.get(k, None) for rdt in rdtlist]
        val2 = [rdt.get(k, None) for rdt in rdtlist2]
        if val[0] is not None:
            rdts[k] = np.array(val)
            rdttot[k] = np.array(val, dtype=complex)
        if val2[0] is not None:
            rdts2[k] = np.array(val2)
            rdttot[k] += np.array(val2, dtype=complex)
    rdts2["refpts"] = rdts["refpts"]
    rdttot["refpts"] = rdts["refpts"]
    ardts = np.rec.fromarrays(rdts.values(), names=list(rdts.keys()))
    ardts2 = np.rec.fromarrays(rdts2.values(), names=list(rdts2.keys()))
    ardttot = np.rec.fromarrays(rdttot.values(), names=list(rdttot.keys()))
    return ardts, ardts2, ardttot


Lattice.get_rdts = get_rdts
