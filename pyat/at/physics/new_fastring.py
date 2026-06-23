"""
Functions relating to fast_ring.
"""

from __future__ import annotations

__all__ = ["fast_ring_new"]

from collections.abc import Sequence
import numpy as np
from ..lattice import Lattice, Refpts
from ..lattice import Drift, RFCavity, Element, Marker, Radiative
from ..physics import gen_m66_elem, gen_detuning_elem, gen_quantdiff_elem
from ..physics import ELossMethod


def _replace_cav(cav_element: Element) -> Element:
    name = ("DR_" + cav_element.FamName,)
    elem = Drift(name, cav_element.Length) if cav_element.Length > 0 else Marker(name)
    return elem


def _merge_cavs(all_cavs: Sequence) -> Sequence[RFCavity]:
    m_cavs = []
    freqs = [e.Frequency for e in all_cavs]
    for i, fr in enumerate(np.atleast_1d(np.unique(freqs))):
        cavf = [cav for cav in all_cavs if cav.Frequency == fr]
        vol = np.sum([c.Voltage for c in cavf])
        cavl = RFCavity(
            f"CAV_{i}",
            0,
            vol,
            fr,
            cavf[0].HarmNumber,
            cavf[0].Energy,
            TimeLag=cavf[0].TimeLag,
        )
        m_cavs.append(cavl)
    return m_cavs


def _split_ring(ring: Lattice, split_inds: Refpts | None = None) -> Sequence:
    inds = ring.get_bool_index(split_inds, endpoint=True)
    inds[[0, -1]] = True
    inds = ring.get_uint32_index(inds)
    all_rings = [ring[int(b) : int(e)] for b, e in zip(inds[:-1], inds[1:], strict=True)]
    split_ring = []
    for r in all_rings:
        r.insert(0, Marker("xsplit"))
        split_ring += r
    split_ring = Lattice(split_ring, **vars(ring))
    split_ring.append(Marker("xsplit"))
    return all_rings, split_ring


def _rearrange(all_rings) -> tuple:
    all_cavs = []
    for r in all_rings:
        cav_idx = r.get_uint32_index(RFCavity)
        mcavs = _merge_cavs(r[cav_idx])
        r[cav_idx] = [_replace_cav(r[i]) for i in cav_idx]
        all_cavs.append(mcavs)
    return all_cavs


def fast_ring_new(
    ring,
    split_inds: Refpts | None = None,
    qpx: Sequence[float] | None = None,
    qpy: Sequence[float] | None = None,
    detuning_coeff: Sequence[float] | None = None,
) -> Lattice:
    """
    A fast ring consisting in:

    * a RF cavity per distinct frequency,
    * a 6x6 linear transfer map,
    * a detuning and chromaticity element,
    * a quantum diffusion element active when radiations are enabled.

    A new lattices is returned with the same attributes (energy, particle,
    circumference, periodicity,…) as the initial one. The radiation state is
    also concerved. Radiations can be turned on or of using ring.disable_6d()
    or ring.enable_6d().

    It is possible to split the original ring in multiple "fastrings"
    using the ``split_inds`` argument.

    Parameters:
        ring:       Lattice description
        split_inds: List of indexes where to split the ring
        qpx:            horizontal chromatic detuning coefficients. Default None.
          If specified qpy should also be provided, if not first order term computed
          automatically from ring.
        qpy:            vertical chromatic detuning coefficients. Default None.
          If specified qpy should also be provided, if not first order term computed
          automatically from ring.
        detuning_coeff: First order amplitude detuning coefficients
          [dQx/dJx, dQx/dJy, dQy/dJy]. Default None: coefficient computer from ring.

    Returns:
        fastring (Lattice):    Fast ring lattice object
    """
    all_rings, split_ring = _split_ring(ring, split_inds)
    all_cavs = _rearrange(all_rings)
    fastring = []
    _, o4 = split_ring.disable_6d(copy=True).find_orbit(refpts="xsplit")
    _, o6 = split_ring.enable_6d(copy=True).find_orbit(refpts="xsplit")
    for r, cav, o4b, o4e, o6b, o6e in zip(
        all_rings, all_cavs, o4[:-1], o4[1:], o6[:-1], o6[1:], strict=True
    ):
        rcav = r.enable_6d(copy=True) + cav
        do6 = np.zeros(6)
        do6[4] = -rcav.get_energy_loss(method=ELossMethod.TRACKING) / r.energy
        lin_elem = gen_m66_elem(
            r.disable_6d(copy=True), o4b, o4e, r.enable_6d(copy=True), o6b, o6e + do6
        )
        fastring = fastring + [lin_elem] + list(np.atleast_1d(cav))
    detuning_elem = gen_detuning_elem(
        ring, qpx=qpx, qpy=qpy, detuning_coeff=detuning_coeff, orbit=o6[-1]
    )
    qd_elem = gen_quantdiff_elem(ring.enable_6d(copy=True), orbit=o6[-1])
    fastring.append(detuning_elem)
    fastring.append(qd_elem)
    fastring = Lattice(fastring, **vars(ring))
    if ring.radiation:
        fastring.enable_6d()
    else:
        fastring.disable_6d()
    return fastring
