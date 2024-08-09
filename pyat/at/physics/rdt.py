from __future__ import annotations
import numpy
from ..lattice import Lattice, Refpts, Multipole, All
from enum import Enum
from collections.abc import Sequence
from dataclasses import dataclass


__all__ = ['get_rdts', 'RDTType']


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
    #Normal quadrupoles rdts
    h20000: Sequence[complex] | complex = None
    h00200: Sequence[complex] | complex = None
    #linear coupling rdts
    h10010: Sequence[complex] | complex = None
    h10100: Sequence[complex] | complex = None
    #chromatic rdts
    h11001: Sequence[complex] | complex = None
    h00111: Sequence[complex] | complex = None
    h20001: Sequence[complex] | complex = None
    h00201: Sequence[complex] | complex = None
    h10002: Sequence[complex] | complex = None
    #geometric rdts
    h21000: Sequence[complex] | complex = None
    h30000: Sequence[complex] | complex = None
    h10110: Sequence[complex] | complex = None
    h10020: Sequence[complex] | complex = None
    h10200: Sequence[complex] | complex = None
    #Detuning
    dnux_dJx: Sequence[float] | float = None
    dnux_dJy: Sequence[float] | float = None
    dnuy_dJy: Sequence[float] | float = None


def _get_polynom(elem, attr, index):
    try:
        val = getattr(elem, attr)[index]
        return val
    except (IndexError, AttributeError):
        return 0


def _computeDrivingTerms(s, betax, betay, phix, phiy, etax, a2l, b2l, b3l, b4l, tune,
                        rdttype, nperiods, refpt):
    """
    Original implementation from ELEGANT
    Based on J.Bengtsson, SLS Note 9 / 97, March 7, 1997, with corrections per W.Guo (NSLS)
    Revised to follow C.X.Wang AOP - TN - 2009 - 020 for second - order terms
    """
    periodicfactor = numpy.ones((9,9))
    rdts = _RDT()

    def PF(i, j):
        return periodicfactor[4+i][4+j]

    if nperiods != 1:
        for i in range(9):
            for j in range(9):
                a1 = numpy.pi*2 * (tune[0] * (i-4)+tune[1] * (j-4))
                a2 = a1/nperiods
                periodicfactor[i][j] = (numpy.exp(1j * a1)-1.0) / (numpy.exp(1j * a2)-1.0)

    rbetax = numpy.sqrt(betax)
    rbetay = numpy.sqrt(betay)
    px = numpy.exp(1j * phix)
    py = numpy.exp(1j * phiy)

    rdts.refpts = refpt
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
        rdts.h20000 = sum((b2lm / 8) * betaxm * pxm * pxm * PF(2, 0))
        rdts.h00200 = sum((b2lm / 8) * betaym * pym * pym * PF(0, 2))


    if (RDTType.COUPLING in rdttype) or (RDTType.ALL in rdttype):
        a2lm = a2l[mask_a2l]
        rbetaxm = rbetax[mask_a2l]
        rbetaym = rbetay[mask_a2l]
        pxm = px[mask_a2l]
        pym = py[mask_a2l]
        rdts.h10010 = sum((a2lm / 4) * rbetaxm * rbetaym * pxm / pxm * PF(1, -1))
        rdts.h10100 = sum((a2lm / 4) * rbetaxm * rbetaym * pxm / pym * PF(1, 1))

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
        rdts.h11001 = sum((b3lm * betaxm * etaxm / 2 - b2lm * betaxm / 4) * nperiods)
        rdts.h00111 = sum((b2lm * betaym / 4 - b3lm * betaym * etaxm / 2) * nperiods)
        rdts.h20001 = sum((b3lm * betaxm * etaxm / 2 - b2lm * betaxm / 4) / 2 * pxm * pxm * PF(2, 0))
        rdts.h00201 = sum((b2lm * betaym / 4 - b3lm * betaym * etaxm / 2) / 2 * pym * pym * PF(0, 2))
        rdts.h10002 = sum((b3lm * rbetaxm * etaxm * etaxm - b2lm * rbetaxm * etaxm) / 2 * pxm * PF(1, 0))

    if (RDTType.GEOMETRIC1 in rdttype) or (RDTType.ALL in rdttype):
        b3lm = b3l[mask_b3l]
        betaxm = betax[mask_b3l]
        betaym = betay[mask_b3l]
        rbetaxm = rbetax[mask_b3l]
        pxm = px[mask_b3l]
        pym = py[mask_b3l]
        rdts.h21000 = sum(b3lm * rbetaxm * betaxm / 8 * pxm * PF(1, 0))
        rdts.h30000 = sum(b3lm * rbetaxm * betaxm / 24 * pxm*pxm*pxm * PF(3, 0))
        rdts.h10110 = sum(-b3lm * rbetaxm * betaxm / 4 * pxm * PF(1, 0))
        rdts.h10020 = sum(-b3lm * rbetaxm * betaym / 8 * pxm * numpy.conj(pym * pym) * PF(1, -2))
        rdts.h10200 = sum(-b3lm * rbetaxm * betaym / 8 * pxm * pym * pym * PF(1, 2))

    if (RDTType.TUNESHIFT in rdttype) or (RDTType.ALL in rdttype):
        b4lm = b4l[mask_b4l]
        betaxm = betax[mask_b4l]
        betaym = betay[mask_b4l]
        rdts.dnux_dJx = sum(3 * b4lm * betaxm * betaxm / (8 * numpy.pi) * nperiods)
        rdts.dnux_dJy = sum(-3 * b4lm * betaxm * betaym / (4 * numpy.pi) * nperiods)
        rdts.dnuy_dJy = sum(3 * b4lm * betaym * betaym / (8 * numpy.pi) * nperiods)

    if (RDTType.GEOMETRIC2 in rdttype) or (RDTType.ALL in rdttype):
        b4lm = b4l[mask_b4l]
        betaxm = betax[mask_b4l]
        betaym = betay[mask_b4l]
        pxm = px[mask_b4l]
        pym = py[mask_b4l]
        rdts.h22000 = sum(3 * b4lm * betaxm * betaxm / 32 * nperiods)
        rdts.h11110 = sum(-3 * b4lm * betaxm * betaym / 8 * nperiods)
        rdts.h00220 = sum(3 * b4lm * betaym * betaym / 32 * nperiods)
        rdts.h31000 = sum(b4lm * betaxm * betaxm / 16 * pxm * pxm * PF(2, 0))
        rdts.h40000 = sum(b4lm * betaxm * betaxm / 64 * pxm * pxm * pxm * pxm * PF(4, 0))
        rdts.h20110 = sum(-3 * b4lm * betaxm * betaym / 16 * pxm * pxm * PF(2, 0))
        rdts.h11200 = sum(-3 * b4lm * betaxm * betaym / 16 * pym * pym * PF(0, 2))
        rdts.h20020 = sum(-3 * b4lm * betaxm * betaym / 32 * pxm * pxm * numpy.conj(pym * pym) * PF(2, -2))
        rdts.h20200 = sum(-3 * b4lm * betaxm * betaym / 32 * pxm * pxm * pym * pym * PF(2, 2))
        rdts.h00310 = sum(b4lm * betaym * betaym / 16 * pym * pym * PF(0, 2))
        rdts.h00400 = sum(b4lm * betaym * betaym / 64 * pym * pym * pym * pym * PF(0, 4))

    return rdts


def get_rdts(ring: Lattice, refpts: Refpts, rdt_type: Sequence[RDTType] | RDTType, nperiods=1):

    rdt_type = numpy.atleast_1d(rdt_type)

    refpts = ring.uint32_refpts(refpts)
    lo, avebeta, avemu, avedisp, *_ = ring.avlinopt(refpts=All)
    idx_mag = ring.get_uint32_index(Multipole)

    sall =  ring.get_s_pos(All)
    smag = sall[idx_mag]
    betax = avebeta[idx_mag, 0]
    betay = avebeta[idx_mag, 1]
    etax = avedisp[idx_mag, 0]
    phix = avemu[idx_mag, 0]
    phiy = avemu[idx_mag, 1]

    a2l = [_get_polynom(e, 'PolynomA', 1)*e.Length for e in ring[idx_mag]]
    b2l = [_get_polynom(e, 'PolynomB', 1)*e.Length for e in ring[idx_mag]]
    b3l = [_get_polynom(e, 'PolynomB', 2)*e.Length for e in ring[idx_mag]]
    b4l = [_get_polynom(e, 'PolynomB', 3)*e.Length for e in ring[idx_mag]]

    mux = lo[-1].mu[0]
    muy = lo[-1].mu[1]
    tune = numpy.array([mux, muy]) / 2.0 / numpy.pi

    rdtlist = []
    for ii in refpts:
        start_idx = sum(idx_mag < ii)
        betax_rot = numpy.concatenate((betax[start_idx:], betax[:start_idx]))
        betay_rot = numpy.concatenate((betay[start_idx:], betay[:start_idx]))
        etax_rot = numpy.concatenate((etax[start_idx:], etax[:start_idx]))
        phix_rot = numpy.concatenate((phix[start_idx:], phix[:start_idx] + mux)) - avemu[ii, 0]
        phiy_rot = numpy.concatenate((phiy[start_idx:], phiy[:start_idx] + muy)) - avemu[ii, 1]
        s_rot = numpy.concatenate((smag[start_idx:], smag[: start_idx] + sall[-1]))- sall[ii]
        a2l_rot = numpy.concatenate((a2l[start_idx:], a2l[:start_idx]))
        b2l_rot = numpy.concatenate((b2l[start_idx:], b2l[:start_idx]))
        b3l_rot = numpy.concatenate((b3l[start_idx:], b3l[:start_idx]))
        b4l_rot = numpy.concatenate((b4l[start_idx:], b4l[:start_idx]))
        rdtlist.append(_computeDrivingTerms(s_rot, betax_rot, betay_rot, phix_rot, phiy_rot,
                                            etax_rot, a2l_rot, b2l_rot, b3l_rot,
                                            b4l_rot, tune, rdt_type, nperiods, ii))
    rdts = _RDT()
    for k in rdts.__annotations__.keys():
        val = [getattr(rdt, k) for rdt in rdtlist]
        if val[0] != None:
            setattr(rdts, k, val)
    return rdts

Lattice.get_rdts = get_rdts