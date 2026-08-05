"""
Functions relating to fast_ring.
"""

from __future__ import annotations

__all__ = ["fast_ring_new"]

from collections.abc import Sequence
import numpy as np
from ..lattice import Lattice, Refpts, End
from ..lattice import Drift, RFCavity, Element, Marker, SimpleQuantDiff
from ..physics import gen_m66_elem, gen_detuning_elem, gen_quantdiff_elem
from ..physics import ELossMethod


def _replace_cav(cav_element: Element) -> Element:
    name = ("DR_" + cav_element.FamName,)
    elem = Drift(name, cav_element.Length) if cav_element.Length > 0 else Marker(name)
    return elem


def _merge_cavs(all_cavs: Sequence) -> Sequence[RFCavity]:
    m_cavs = []
    freqs = [e.Frequency for e in all_cavs]
    for fr in np.atleast_1d(np.unique(freqs)):
        cavf = [cav for cav in all_cavs if cav.Frequency == fr]
        vol = np.sum([c.Voltage for c in cavf])
        cavl = RFCavity(
            cavf[0].FamName,
            0,
            vol,
            fr,
            cavf[0].HarmNumber,
            cavf[0].Energy,
            TimeLag=cavf[0].TimeLag,
        )
        m_cavs.append(cavl)
    return m_cavs


def _rearrange(ring) -> tuple:
    cav_idx = ring.get_uint32_index(RFCavity)
    mcavs = _merge_cavs(ring[cav_idx])
    ring[cav_idx] = [_replace_cav(ring[i]) for i in cav_idx]
    ring += mcavs


def _split_ring(ring: Lattice, split_inds: Refpts | None = None) -> Sequence:
    inds = ring.get_bool_index(split_inds, endpoint=True)
    inds[[0, -1]] = True
    inds = ring.get_uint32_index(inds)
    split_ring = []
    for b, e in zip(inds[:-1], inds[1:], strict=True):
        split_ring += [Marker("xsplit")]
        split_ring += ring[int(b) : int(e)] 
    split_ring += [Marker("xsplit")]
    split_ring = Lattice(split_ring, **vars(ring)) 
    _rearrange(split_ring)
    return split_ring, split_ring.get_uint32_index("xsplit")


def fast_ring_new(
    ring,
    split_inds: Refpts | None = None,
    qpx: Sequence[float] | None = None,
    qpy: Sequence[float] | None = None,
    detuning_coeff: Sequence[float] | None = None,
    alphac: Sequence[float] | None = None,
    emitx: float = None,
    emity: float = None,
    espread: float = None,
) -> Lattice:
    """
    A fast ring consisting in the following elements.

    * a RF cavity per distinct frequency,
    * a 6x6 linear transfer map,
    * a detuning and chromaticity element,
    * a quantum diffusion element active when radiations are enabled.

    A new lattices is returned with the same attributes (energy, particle,
    circumference, periodicity,…) as the initial one. The radiation state is
    also concerved. Radiations can be turned on or of using ring.disable_6d()
    or ring.enable_6d(). *Rad attribute will be used when radiations are turned
    on.

    It is possible to split the original ring in multiple "fastrings"
    using the ``split_inds`` argument.

    Parameters:
        ring:       Lattice description
        split_inds: List of indexes where to split the ring
        qpx:            horizontal chromatic detuning coefficients. Default None.
          If specified qpy should also be provided, if not first order term computed
          automatically from ring. at.physics.chromaticity can be used to compute
          higher order chromatic coefficients
        qpy:            vertical chromatic detuning coefficients. Default None.
          If specified qpy should also be provided, if not first order term computed
          automatically from ring. at.physics.chromaticity can be used to compute
          higher order chromatic coefficients
        detuning_coeff: First order amplitude detuning coefficients
          [dQx/dJx, dQx/dJy, dQy/dJy]. Default None: coefficient computer from ring.
          at.physics.detuning can be used to compute detuning coefficients
        alphac: Higher order momentum compaction factor. Default None. Can be
          Can be claculated with physics.get_mcf()
        emitx: Equilibrium horizontal emittance. Default is the one given by ring
        emity: Equilibrium vertical emittance. Default is the one given by ring
        espread: Equilibrium energy spread. Default is the one given by ring

    Returns:
        fastring (Lattice):    Fast ring lattice object
    """
    split_ring, split_inds = _split_ring(ring, split_inds)
    _, o4 = split_ring.disable_6d(copy=True).find_orbit(refpts=split_inds)
    _, o6 = split_ring.enable_6d(copy=True).find_orbit(refpts=split_inds)
    _, o4f = split_ring.disable_6d(copy=True).find_orbit(refpts=End)
    _, o6f = split_ring.enable_6d(copy=True).find_orbit(refpts=End)
    fastring = []
    for o4b, o4e, o6b, o6e, sib, sie in zip(o4[:-1], o4[1:], o6[:-1], o6[1:], split_inds[:-1], split_inds[1:]):
        r = split_ring[sib:sie]
        lin_elem = gen_m66_elem(
            r.disable_6d(copy=True), o4b, o4e, r.enable_6d(copy=True), o6b, o6e
        ) 
        fastring = [*fastring, lin_elem]
    fastring = [*fastring, *list(split_ring[split_inds[-1]+1:])]
    detuning_elem = gen_detuning_elem(
        split_ring,
        qpx=qpx,
        qpy=qpy,
        detuning_coeff=detuning_coeff,
        alphac=alphac,
        orbit=o4f[-1],
        orbit6=o6f[-1],
    )
    if emity is None and emitx is None and espread is None:
        qd_elem = gen_quantdiff_elem(split_ring.enable_6d(copy=True), orbit=o6f[-1])
    else:
        params = ring.enable_6d(copy=True).envelope_parameters()
        taux, tauy, tauz = params.Tau
        if emitx is None:
            emitx = params.emittances[0]
        if emity is None:
            emitx = params.emittances[1]
        if espread is None:
            espread = params.sigma_e
        qd_elem = SimpleQuantDiff("Diffusion",
                                  betax=detuning_elem.BetaRad[0],
                                  betay=detuning_elem.BetaRad[1],
                                  emitx=emitx,
                                  emity=emity,
                                  espread=espread,
                                  taux=taux,
                                  tauy=tauy,
                                  tauz=tauz,
        )     
    fastring.append(detuning_elem)
    fastring.append(qd_elem)
    fastring = Lattice(fastring, **vars(ring))  
    if ring.radiation:
        fastring.enable_6d()
    else:
        fastring.disable_6d()
    return fastring, split_ring
