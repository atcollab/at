from __future__ import annotations

import numpy
import multiprocessing
from functools import partial
from ..lattice import Lattice, Refpts, Multipole, All
from enum import Enum
from collections.abc import Sequence
from dataclasses import dataclass, asdict


__all__ = ["get_rdts", "RDTType"]


class RDTType(Enum):
    """Enum class for RDT type"""

    ALL = 0
    FOCUSING = 1
    COUPLING = 2
    CHROMATIC = 3
    GEOMETRIC1 = 4
    GEOMETRIC2 = 5
    TUNESHIFT = 6


@dataclass
class _RDT:
    # Location
    refpts: Sequence[int] | int = None
    # Normal quadrupoles rdts
    h20000: Sequence[complex] | complex = None
    h00200: Sequence[complex] | complex = None
    # linear coupling rdts
    h10010: Sequence[complex] | complex = None
    h10100: Sequence[complex] | complex = None
    # chromatic rdts
    h11001: Sequence[complex] | complex = None
    h00111: Sequence[complex] | complex = None
    h20001: Sequence[complex] | complex = None
    h00201: Sequence[complex] | complex = None
    h10002: Sequence[complex] | complex = None
    # sextupole geometric rdts
    h21000: Sequence[complex] | complex = None
    h30000: Sequence[complex] | complex = None
    h10110: Sequence[complex] | complex = None
    h10020: Sequence[complex] | complex = None
    h10200: Sequence[complex] | complex = None
    # octupole geometric rdts
    h22000: Sequence[complex] | complex = None
    h11110: Sequence[complex] | complex = None
    h00220: Sequence[complex] | complex = None
    h31000: Sequence[complex] | complex = None
    h40000: Sequence[complex] | complex = None
    h20110: Sequence[complex] | complex = None
    h11200: Sequence[complex] | complex = None
    h20020: Sequence[complex] | complex = None
    h20200: Sequence[complex] | complex = None
    h00310: Sequence[complex] | complex = None
    h00400: Sequence[complex] | complex = None
    # Detuning
    dnux_dJx: Sequence[float] | float = None
    dnux_dJy: Sequence[float] | float = None
    dnuy_dJy: Sequence[float] | float = None

    def __getattr__(self, item):
        return asdict(self)[item]

    def __getitem__(self, item):
        return asdict(self)[item]


def _get_polynom(elem, attr, index):
    try:
        val = getattr(elem, attr)[index]
        return val
    except (IndexError, AttributeError):
        return 0


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
):
    """
    Original implementation from ELEGANT
    Based on J.Bengtsson, SLS Note 9 / 97, March 7, 1997, with corrections per W.Guo (NSLS)
    Revised to follow C.X.Wang AOP - TN - 2009 - 020 for second - order terms
    """
    periodicfactor = numpy.ones((9, 9), dtype=complex)
    rdts = {}

    def pf(i, j):
        return periodicfactor[4 + i][4 + j]

    if nperiods != 1:
        for i in range(9):
            for j in range(9):
                a1 = numpy.pi * 2 * (tune[0] * (i - 4) + tune[1] * (j - 4))
                a2 = a1 / nperiods
                periodicfactor[i][j] = (numpy.exp(1j * a1) - 1.0) / (
                    numpy.exp(1j * a2) - 1.0
                )

    rbetax = numpy.sqrt(betax)
    rbetay = numpy.sqrt(betay)
    px = numpy.exp(1j * phix)
    py = numpy.exp(1j * phiy)

    rdts["refpts"] = refpt
    mask_a2l = numpy.absolute(a2l) > 1.0e-6
    mask_b2l = numpy.absolute(b2l) > 1.0e-6
    mask_b3l = numpy.absolute(b3l) > 1.0e-6
    mask_b4l = numpy.absolute(b4l) > 1.0e-6

    if (RDTType.FOCUSING in rdttype) or (RDTType.ALL in rdttype):
        b2lm = b2l[mask_b2l]
        betaxm = betax[mask_b2l]
        betaym = betay[mask_b2l]
        pxm = px[mask_b2l]
        pym = py[mask_b2l]
        rdts["h20000"] = sum((b2lm / 8) * betaxm * pxm * pxm * pf(2, 0))
        rdts["h00200"] = sum((b2lm / 8) * betaym * pym * pym * pf(0, 2))

    if (RDTType.COUPLING in rdttype) or (RDTType.ALL in rdttype):
        a2lm = a2l[mask_a2l]
        rbetaxm = rbetax[mask_a2l]
        rbetaym = rbetay[mask_a2l]
        pxm = px[mask_a2l]
        pym = py[mask_a2l]
        rdts["h10010"] = sum((a2lm / 4) * rbetaxm * rbetaym * pxm / pxm * pf(1, -1))
        rdts["h10100"] = sum((a2lm / 4) * rbetaxm * rbetaym * pxm / pym * pf(1, 1))

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
        rdts["h11001"] = sum((b3lm * betaxm * etaxm / 2 - b2lm * betaxm / 4) * nperiods)
        rdts["h00111"] = sum((b2lm * betaym / 4 - b3lm * betaym * etaxm / 2) * nperiods)
        rdts["h20001"] = sum(
            (b3lm * betaxm * etaxm / 2 - b2lm * betaxm / 4) / 2 * pxm * pxm * pf(2, 0)
        )
        rdts["h00201"] = sum(
            (b2lm * betaym / 4 - b3lm * betaym * etaxm / 2) / 2 * pym * pym * pf(0, 2)
        )
        # fmt: off
        rdts["h10002"] = sum(
            (b3lm * rbetaxm * etaxm * etaxm - b2lm * rbetaxm * etaxm)
            / 2 * pxm * pf(1, 0)
        )
        # fmt: on

    if (RDTType.GEOMETRIC1 in rdttype) or (RDTType.ALL in rdttype):
        b3lm = b3l[mask_b3l]
        betaxm = betax[mask_b3l]
        betaym = betay[mask_b3l]
        rbetaxm = rbetax[mask_b3l]
        pxm = px[mask_b3l]
        pym = py[mask_b3l]
        rdts["h21000"] = sum(b3lm * rbetaxm * betaxm / 8 * pxm * pf(1, 0))
        rdts["h30000"] = sum(b3lm * rbetaxm * betaxm / 24 * pxm * pxm * pxm * pf(3, 0))
        rdts["h10110"] = sum(-b3lm * rbetaxm * betaym / 4 * pxm * pf(1, 0))
        rdts["h10020"] = sum(
            -b3lm * rbetaxm * betaym / 8 * pxm * numpy.conj(pym * pym) * pf(1, -2)
        )
        rdts["h10200"] = sum(-b3lm * rbetaxm * betaym / 8 * pxm * pym * pym * pf(1, 2))

    if (RDTType.TUNESHIFT in rdttype) or (RDTType.ALL in rdttype):
        b4lm = b4l[mask_b4l]
        betaxm = betax[mask_b4l]
        betaym = betay[mask_b4l]
        rdts["dnux_dJx"] = sum(3 * b4lm * betaxm * betaxm / (8 * numpy.pi) * nperiods)
        rdts["dnux_dJy"] = sum(-3 * b4lm * betaxm * betaym / (4 * numpy.pi) * nperiods)
        rdts["dnuy_dJy"] = sum(3 * b4lm * betaym * betaym / (8 * numpy.pi) * nperiods)

    if (RDTType.GEOMETRIC2 in rdttype) or (RDTType.ALL in rdttype):
        b4lm = b4l[mask_b4l]
        betaxm = betax[mask_b4l]
        betaym = betay[mask_b4l]
        pxm = px[mask_b4l]
        pym = py[mask_b4l]
        rdts["h22000"] = sum(3 * b4lm * betaxm * betaxm / 32 * nperiods)
        rdts["h11110"] = sum(-3 * b4lm * betaxm * betaym / 8 * nperiods)
        rdts["h00220"] = sum(3 * b4lm * betaym * betaym / 32 * nperiods)
        rdts["h31000"] = sum(b4lm * betaxm * betaxm / 16 * pxm * pxm * pf(2, 0))
        rdts["h40000"] = sum(
            b4lm * betaxm * betaxm / 64 * pxm * pxm * pxm * pxm * pf(4, 0)
        )
        rdts["h20110"] = sum(-3 * b4lm * betaxm * betaym / 16 * pxm * pxm * pf(2, 0))
        rdts["h11200"] = sum(-3 * b4lm * betaxm * betaym / 16 * pym * pym * pf(0, 2))
        # fmt: off
        rdts["h20020"] = sum(
            -3 * b4lm * betaxm * betaym / 32 * pxm * pxm
            * numpy.conj(pym * pym) * pf(2, -2)
        )
        # fmt: on
        rdts["h20200"] = sum(
            -3 * b4lm * betaxm * betaym / 32 * pxm * pxm * pym * pym * pf(2, 2)
        )
        rdts["h00310"] = sum(b4lm * betaym * betaym / 16 * pym * pym * pf(0, 2))
        rdts["h00400"] = sum(
            b4lm * betaym * betaym / 64 * pym * pym * pym * pym * pf(0, 4)
        )

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
            cpxm = numpy.conj(pxm)
            cpxm3 = numpy.conj(pxm3)
            pym = py[mask_b3l]
            pym2 = pym * pym
            cpym2 = numpy.conj(pym2)
            cpxym2 = numpy.conj(pxm * pym2)
            # fmt: off
            bsign = numpy.array([numpy.sign(sm[i] - sm) * 1j * b3lm[i] * b3lm
                     for i in range(nelem)])
            rbbx = numpy.array([bsign[i] * rbbetaxm[i] * rbbetaxm
                                for i in range(nelem)])
            rbxy = numpy.array([bsign[i] * rbetaxm[i] * rbetaxm * betaym[i]
                                for i in range(nelem)])
            ppxm = numpy.array([pxm[i] * pxm for i in range(nelem)])

            rdts["h22000"] += (1.0 / 64) * numpy.sum([
                rbbx[i] * (pxm3[i] * cpxm3 + 3 * pxm[i] * cpxm)
                for i in range(nelem)]
            )
            rdts["h31000"] += (1.0 / 32) * numpy.sum([
                rbbx[i] * pxm3[i] * cpxm
                for i in range(nelem)]
            )
            t1 = numpy.array([numpy.conj(pxm[i]) * pxm for i in range(nelem)])
            t2 = numpy.conj(t1)
            rdts["h11110"] += (1.0 / 16) * numpy.sum([
                rbxy[i] * (betaxm * (t1[i] - t2[i]) + betaym * pym2[i]
                           * cpym2 * (t2[i] + t1[i]))
                for i in range(nelem)]
            )
            t1 = numpy.array([numpy.exp(-1j * (phixm[i] - phixm))
                              for i in range(nelem)])
            t2 = numpy.conj(t1)
            rdts["h11200"] += (1.0 / 32) * numpy.sum([
                rbxy[i] * numpy.exp(1j * (2 * phiym[i]))
                * (betaxm * (t1[i] - t2[i]) + 2 * betaym * (t2[i] + t1[i]))
                for i in range(nelem)]
            )
            rdts["h40000"] += (1.0 / 64) * numpy.sum([
                rbbx[i] * pxm3[i] * pxm
                for i in range(nelem)]
            )
            rdts["h20020"] += (1.0 / 64) * numpy.sum([
                rbxy[i] * (betaxm * cpxym2[i] * pxm3 - (betaxm + 4 * betaym)
                           * ppxm[i] * cpym2[i])
                for i in range(nelem)]
            )
            rdts["h20110"] += (1.0 / 32) * numpy.sum([
                rbxy[i] * (betaxm * (cpxm[i] * pxm3 - ppxm[i])
                           + 2 * betaym * ppxm[i] * pym2[i] * cpym2)
                for i in range(nelem)]
            )
            rdts["h20200"] += (1.0 / 64) * numpy.sum([
                rbxy[i] * (betaxm * cpxm[i] * pxm3 * pym2[i]
                           - (betaxm - 4 * betaym)
                           * ppxm[i] * pym2[i])
                for i in range(nelem)]
            )
            rdts["h00220"] += (1.0 / 64) * numpy.sum([
                rbxy[i] * betaym * (pxm[i] * pym2[i] * cpxym2 + 4 * pxm[i] * cpxm
                                    - numpy.conj(pxm[i] * pym2) * pxm * pym2[i])
                for i in range(nelem)]
            )
            rdts["h00310"] += (1.0 / 32) * numpy.sum([
                rbxy[i] * betaym * pym2[i] * (pxm[i] * cpxm - pxm * cpxm[i])
                for i in range(nelem)]
            )
            rdts["h00400"] += (1.0 / 64) * numpy.sum([
                rbxy[i] * betaym * pxm[i] * cpxm * pym2[i] * pym2
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
            b3f = [b3lm[i] * b3lm / numpy.pi for i in range(nelem)]
            sqbx = [numpy.sqrt(betaxm[i] * betaxm) for i in range(nelem)]
            bbx = [betaxm[i] * betaxm for i in range(nelem)]
            bby = [betaym[i] * betaym for i in range(nelem)]
            dphix = [abs(phixm[i] - phixm) for i in range(nelem)]
            dphiy = [abs(phiym[i] - phiym) for i in range(nelem)]
            rdts["dnux_dJx"] += numpy.sum([
                -b3f[i] / 16 * sqbx[i]  * bbx[i]
                * (3 * numpy.cos(dphix[i] - numpy.pi * nux)
                   / numpy.sin(numpy.pi * nux)
                   + numpy.cos(3 * dphix[i] - 3 * numpy.pi * nux)
                   / numpy.sin(3 * numpy.pi * nux))
                for i in range(nelem)]
            )
            rdts["dnux_dJy"] += numpy.sum([
                b3f[i]  / 8 * sqbx[i]  * betaym[i]
                * (2 * betaxm * numpy.cos(dphix[i]  - numpy.pi * nux)
                   / numpy.sin(numpy.pi * nux)
                   - betaym * numpy.cos(dphix[i]  + 2 * dphiy[i]
                                        - numpy.pi * (nux + 2 * nuy))
                   / numpy.sin(numpy.pi * (nux + 2 * nuy))
                   + betaym * numpy.cos(dphix[i]  - 2 * dphiy[i]
                                        - numpy.pi * (nux - 2 * nuy))
                   / numpy.sin(numpy.pi * (nux - 2 * nuy)))
                for i in range(nelem)]
            )
            rdts["dnuy_dJy"] += numpy.sum([
                -b3f[i]  / 16 * sqbx[i]  * bby[i]
                * (4 * numpy.cos(dphix[i]  - numpy.pi * nux)
                   / numpy.sin(numpy.pi * nux)
                   + numpy.cos(dphix[i]  + 2 * dphiy[i]
                               - numpy.pi * (nux + 2 * nuy))
                   / numpy.sin(numpy.pi * (nux + 2 * nuy))
                   + numpy.cos(dphix[i]  - 2 * dphiy[i]
                               - numpy.pi * (nux - 2 * nuy))
                   / numpy.sin(numpy.pi * (nux - 2 * nuy)))
                for i in range(nelem)]
            )
            # fmt: on
    return rdts


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
    for ii in refpts:
        start_idx = sum(idx_mag < ii)
        beta_rot = numpy.roll(beta, -start_idx, axis=0)
        etax_rot = numpy.roll(etax, -start_idx)
        phi_rot = numpy.roll(phi, -start_idx, axis=0) - avemu[ii]
        if start_idx > 0:
            phi_rot[-start_idx:] += mu
        s_rot = numpy.roll(smag, -start_idx) - sall[ii]
        s_rot[-start_idx:] += sall[-1]
        pols_rot = numpy.roll(pols, -start_idx, axis=0)
        rdtlist.append(
            _computedrivingterms(
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
            )
        )
    return rdtlist


def get_rdts(
    ring: Lattice,
    refpts: Refpts,
    rdt_type: Sequence[RDTType] | RDTType,
    second_order: bool = False,
    use_mp: bool = False,
    pool_size: int = None,
):
    """
    :py:func:`get_rdts` computes the ring RDTs based on the original implementation
    from ELEGANT.
    J.Bengtsson, SLS Note 9 / 97, March 7, 1997, with corrections per W.Guo (NSLS)
    Revised to follow C.X.Wang AOP - TN - 2009 - 020 for second - order terms

    Usage:
      >>> get_rdts(ring, reftps, [RDTType.COUPLING, RDTType.CHROMATIC])

    Parameters:
        ring: :code:`at.Lattice` object
        refpts: Element refpts at which the RDTs are calculated
        rdt_type: Type of RDTs to be calculated. The type can be
        :code:`Sequence[at.RDTType] | at.RDTType`.

    Keyword Args:
        second_order: Compute second order terms (default: False).
          Computation is significantly longer using this method
        use_mp: Activate parallel calculation
        pool_size: Number of processes used for parallelization

    Returns:
        rdts: rdt data (complex) at refpts

        **rdts** is a dataclass with fields:
        =================   ======
        **refts**           location of the rdt

        **h20000**          Normal quadrupole RDTS
        **h00200**

        **h10010**          Coupling RDTs
        **h10100**

        **h11001**          Chromatic RDTs
        **h00111**
        **h20001**
        **h00201**
        **h10002**

        **h21000**          Sextupole geometric RDTs
        **h30000**
        **h10110**
        **h10020**
        **h10200**

        **h22000**          Octupole geometric RDTs
        **h11110**
        **h00220**
        **h31000**
        **h40000**
        **h20110**
        **h11200**
        **h20020**
        **h20200**
        **h00310**
        **h00400**

        **dnux_dJx**        Detuning terms from octupoles
        **dnux_dJy**
        **dnuy_dJy**
        =================   ======
    """
    rdt_type = numpy.atleast_1d(rdt_type)
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
    tune = mu / 2.0 / numpy.pi
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
            pool_size = multiprocessing.cpu_count()
        ctx = multiprocessing.get_context()
        refs = numpy.array_split(refpts, pool_size)
        with ctx.Pool(pool_size) as pool:
            rdtlist = pool.map(fun, refs)
            rdtlist = numpy.concatenate(rdtlist)
    else:
        rdtlist = fun(refpts)
    rdts = _RDT()
    for k in rdts.__annotations__.keys():
        val = [rdt.get(k, None) for rdt in rdtlist]
        if val[0] is not None:
            setattr(rdts, k, val)
    return rdts


Lattice.get_rdts = get_rdts
