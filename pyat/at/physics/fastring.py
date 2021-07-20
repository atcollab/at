"""
Functions relating to fast_ring
"""
import numpy
from math import sqrt, atan2, pi
from at.lattice import uint32_refpts, get_refpts, get_elements
from at.lattice import Element, RFCavity, Marker, Drift, Lattice
from at.physics import find_orbit6, find_orbit4
from at.physics import gen_m66_elem, gen_detuning_elem, gen_quantdiff_elem


__all__ = ['fast_ring']


def rearrange(ring, split_inds=[]):


    def merge_rings(all_rings):
        ring = all_rings[0]
        if len(all_rings) > 1:
            for ringi in all_rings[1:]:
                ring += ringi
        return ring


    iends = numpy.concatenate((split_inds, [len(ring)+1]))
    iends = iends[iends != 0]
    ibegs = numpy.concatenate(([0], iends[:-1]))
    ibegs = ibegs[ibegs != len(ring) + 1]
    all_rings = [ring[int(ibeg):int(iend)]
                 for ibeg, iend in zip(ibegs, iends)]

    for ring_slice in all_rings:
        # replace cavity with length > 0 with drift
        # set cavity length to 0 and move to start
        icav = get_refpts(ring_slice, RFCavity)
        allcavs= [ring_slice[i] for i in icav]
        for i,c in zip(icav,allcavs):
            if c.Length != 0:
                cavdrift = Drift('CavDrift', c.Length)
                ring_slice[i] = cavdrift
                c.Length = 0.0

        # merge all cavities with the same frequency
        all_freq = numpy.array([e.Frequency for e in allcavs])
        uni_freq = numpy.unique(all_freq)
        for fr in numpy.atleast_1d(uni_freq):
            cavf = [c for c in allcavs if c.Frequency==fr]
            vol = numpy.sum([c.Voltage for c in cavf])
            cavf[0].Voltage=vol
            ring_slice.insert(0,cavf[0])

        ring_slice.insert(len(uni_freq), Marker('xbeg'))
        ring_slice.append(Marker('xend'))
    return all_rings, merge_rings(all_rings)


def fast_ring(ring, split_inds=[]):
    ring = ring.deepcopy()
    ringr = ring.deepcopy()
    ring.radiation_off()
    ringr.radiation_on(quadrupole_pass='auto')
    all_rings, merged_ring = rearrange(ring, split_inds=split_inds)
    all_ringsr, merged_ringr = rearrange(ringr, split_inds=split_inds)

    ibegs = get_refpts(merged_ring, 'xbeg')
    iends = get_refpts(merged_ring, 'xend')
    markers = numpy.sort(numpy.concatenate((ibegs, iends)))
    _, orbit4 = merged_ring.find_sync_orbit(dct=0.0, refpts=markers)
    _, orbit6 = merged_ringr.find_orbit6(refpts=markers)

    detuning_elem = gen_detuning_elem(merged_ring, orbit4[-1])   
    detuning_elem_rad = detuning_elem.deepcopy()
    detuning_elem_rad.T1 = -orbit6[-1]
    detuning_elem_rad.T2 = orbit6[-1]

    qd_elem = gen_quantdiff_elem(merged_ringr)

    fastringnorad, fastringrad= [], []

    cnts = range(len(all_rings))
    for counter, ring_slice, ring_slicer in zip(cnts, all_rings, all_ringsr):
        ring_slice.radiation_off(cavity_pass='CavityPass')
        cavs = get_elements(ring_slice, RFCavity)
        cavs_rad = get_elements(ring_slicer, RFCavity)
        [ring_slice.remove(c) for c in cavs]
        [ring_slicer.remove(c) for c in cavs_rad]
        lin_elem = gen_m66_elem(ring_slice, orbit4[2*counter],
                                orbit4[2*counter+1])
        lin_elem.FamName = lin_elem.FamName + '_' + str(counter)
        lin_elem_rad = gen_m66_elem(ring_slicer, orbit6[2*counter],
                                    orbit6[2*counter+1])
        lin_elem_rad.FamName = lin_elem_rad.FamName + '_' + str(counter)
        [fastringnorad.append(cav) for cav in cavs]
        fastringnorad.append(lin_elem)
        [fastringrad.append(cav) for cav in cavs_rad]
        fastringrad.append(lin_elem_rad)

    fastringnorad.append(detuning_elem)
    fastringrad.append(qd_elem)
    fastringrad.append(detuning_elem_rad)

    fastringnorad = Lattice(fastringnorad, energy=ring.energy)
    fastringrad = Lattice(fastringrad, energy=ring.energy)
    fastringrad._radiation = True
    return fastringnorad, fastringrad
