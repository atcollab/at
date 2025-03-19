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
        rdts["h20000"] = sum((b2lm / 8) * betaxm * pxm * pxm * pf(2, 0))
        rdts["h00200"] = -sum((b2lm / 8) * betaym * pym * pym * pf(0, 2))

    if (RDTType.COUPLING in rdttype) or (RDTType.ALL in rdttype):
        a2lm = a2l[mask_a2l]
        rbetaxm = rbetax[mask_a2l]
        rbetaym = rbetay[mask_a2l]
        pxm = px[mask_a2l]
        pym = py[mask_a2l]
        rdts["h10010"] = sum((a2lm / 4) * rbetaxm * rbetaym * pxm / pym * pf(1, -1))
        rdts["h10100"] = sum((a2lm / 4) * rbetaxm * rbetaym * pxm * pym * pf(1, 1))

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
            -b3lm * rbetaxm * betaym / 8 * pxm * np.conj(pym * pym) * pf(1, -2)
        )
        rdts["h10200"] = sum(-b3lm * rbetaxm * betaym / 8 * pxm * pym * pym * pf(1, 2))

    if (RDTType.TUNESHIFT in rdttype) or (RDTType.ALL in rdttype):
        b4lm = b4l[mask_b4l]
        betaxm = betax[mask_b4l]
        betaym = betay[mask_b4l]
        rdts["dnux_dJx"] = sum(3 * b4lm * betaxm * betaxm / (8 * np.pi) * nperiods)
        rdts["dnux_dJy"] = sum(-3 * b4lm * betaxm * betaym / (4 * np.pi) * nperiods)
        rdts["dnuy_dJy"] = sum(3 * b4lm * betaym * betaym / (8 * np.pi) * nperiods)

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
            * np.conj(pym * pym) * pf(2, -2)
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
    ardts = np.rec.fromarrays(rdts.values(), names=list(keylist))
    ardts2 = np.rec.fromarrays(rdts2.values(), names=list(keylist))
    ardttot = np.rec.fromarrays(rdttot.values(), names=list(keylist))
    return ardts, ardts2, ardttot


Lattice.get_rdts = get_rdts
