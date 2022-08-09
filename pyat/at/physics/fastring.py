"""
Functions relating to fast_ring
"""
from functools import reduce
import numpy
from typing import Tuple
from at.lattice import RFCavity, Marker, Lattice, get_cells, checkname
from at.lattice import get_elements
from at.physics import gen_m66_elem, gen_detuning_elem, gen_quantdiff_elem
import copy


__all__ = ['fast_ring']


def _rearrange(ring: Lattice, split_inds=[]):
    inds = numpy.append(split_inds, [0, len(ring)+1])
    inds = numpy.unique(inds)
    all_rings = [ring[int(b):int(e)] for b, e in zip(inds[:-1], inds[1:])]

    ringm = []
    for ring_slice in all_rings:
        ring_slice.insert(0, Marker('xbeg'))
        ring_slice.append(Marker('xend'))
        cavs = [e for e in ring_slice if isinstance(e, RFCavity)]
        newpass = ['IdentityPass' if c.Length == 0
                   else 'DriftPass' for c in cavs]
        for c, pm in zip(cavs, newpass):
            c.PassMethod = pm
        uni_freq = numpy.unique([e.Frequency for e in cavs])
        for fr in numpy.atleast_1d(uni_freq):
            cavf = [c for c in cavs if c.Frequency == fr]
            vol = reduce(lambda x, y: x+y, (c.Voltage for c in cavf))
            cavl = RFCavity('CAVL', 0, vol, fr,
                            cavf[0].HarmNumber, cavf[0].Energy)
            cavl.TimeLag = cavf[0].TimeLag
            ring_slice.append(cavl)
        ringm = ringm + ring_slice
    return all_rings, Lattice(ringm, energy=ring.energy)


def _fring(ring, split_inds=[], detuning_elem=None):
    all_rings, merged_ring = _rearrange(ring, split_inds=split_inds)
    ibegs = get_cells(merged_ring, checkname('xbeg'))
    iends = get_cells(merged_ring, checkname('xend'))
    _, orbit = merged_ring.find_orbit(refpts=ibegs | iends)
    if detuning_elem is None:
        detuning_elem = gen_detuning_elem(merged_ring, orbit[-1])
    else:
        detuning_elem.T1 = -orbit[-1]
        detuning_elem.T2 = orbit[-1]

    fastring = []
    for counter, r in enumerate(all_rings):
        cavs = [e for e in r if e.PassMethod.endswith('CavityPass')]
        [r.remove(c) for c in cavs]
        lin_elem = gen_m66_elem(r, orbit[2*counter],
                                orbit[2*counter+1])
        lin_elem.FamName = lin_elem.FamName + '_' + str(counter)
        [fastring.append(c) for c in cavs]
        fastring.append(lin_elem)
    fastring.append(detuning_elem)
    try:
        qd_elem = gen_quantdiff_elem(merged_ring)
        fastring.append(qd_elem)
    except ValueError:      # No synchrotron radiation => no diffusion element
        pass
    fastring = Lattice(fastring, **vars(ring))
    return fastring


def fast_ring(ring: Lattice, split_inds=[]) -> Tuple[Lattice, Lattice]:
    """Generates a "fast ring"

    A fast ring consisting in:

    * a RF cavity per distinct frequency,
    * a 6x6 linear transfer map,
    * a detuning and chromaticity element,
    * a quantum diffusion element (for radiation ring).

    2 new lattices are returned, one with radiation and one without
    It is possible to split the original ring in multiple "fastrings"
    using the ``split_inds`` argument

    Parameters:
        ring:       Lattice description
        split_inds: List of indexes where to split the ring

    Returns:
        fring (Lattice):    Fast ring without radiation
        fringrad (Lattice): Fast ring with radiation
    """
    ringi = ring.deepcopy()
    fastringnorad = _fring(ringi.radiation_off(copy=True),
                           split_inds=split_inds)
    detuning_elem = copy.deepcopy(get_elements(fastringnorad,
                                               'NonLinear')[0])
    fastringrad = _fring(ringi.radiation_on(copy=True),
                         split_inds=split_inds,
                         detuning_elem=detuning_elem)
    return fastringnorad, fastringrad
