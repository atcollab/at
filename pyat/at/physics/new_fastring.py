"""
Functions relating to fast_ring
"""

from __future__ import annotations

__all__ = ["fast_ring", "simple_ring"]

import copy
from collections.abc import Sequence
from functools import reduce

import numpy as np

from ..constants import clight
from ..lattice import Lattice, Particle, Refpts
from ..lattice import Drift, RFCavity, Element, Marker, get_cells, checkname
from ..lattice import get_elements, M66, SimpleQuantDiff, AtError, SimpleRadiation
from ..physics import gen_m66_elem, gen_detuning_elem, gen_quantdiff_elem


def _replace_cav(cav_element):
    name = "DR_"+cav_element.FamName,
    if cav_element.Length > 0:
        elem = Drift(name, cav_element.Length)
    else:
        elem = Marker(name)
    return elem


def _merge_cavs(all_cavs):
    m_cavs= []
    freqs = [e.Frequency for e in all_cavs]
    for i, fr in enumerate(np.atleast_1d(np.unique(freqs))):
        cavf = [cav for cav in all_cavs if cav.Frequency==fr]
        vol = np.sum([c.Voltage for c in cavf])
        cavl = RFCavity(f"CAV_{i}", 0, vol, fr, cavf[0].HarmNumber, cavf[0].Energy,
                        TimeLag=cavf[0].TimeLag)
        m_cavs.append(cavl)
    return m_cavs


def _split_ring(ring: Lattice, split_inds: Refpts = None):
    inds = ring.get_bool_index(split_inds, endpoint=True)
    inds[0] = True
    inds[-1] = True
    inds = ring.get_uint32_index(inds)
    return [ring[int(b) : int(e)] for b, e in zip(inds[:-1], inds[1:])]


def _rearrange(ring: Lattice, split_inds: Refpts = None):
    new_ring = ring.copy()
    cav_idx = new_ring.get_uint32_index(RFCavity)
    all_cavs = []
    for cavi in cav_idx:
        all_cavs.append(new_ring[cavi])
        new_ring[cavi] = _replace_cav(new_ring[cavi])
    all_cavs = _merge_cavs(all_cavs)
    all_rings = _split_ring(new_ring, split_inds)
    return new_ring, all_rings, all_cavs


def _fring(ring, split_inds: Refpts = None, detuning_elem=None)
    

