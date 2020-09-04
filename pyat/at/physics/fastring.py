"""
Functions relating to fast_ring
"""
import numpy
from math import sqrt, atan2, pi
from at.lattice import uint32_refpts, get_refpts
from at.lattice import Element, RFCavity, Marker, Drift
from at.physics import find_orbit6, find_orbit4
from at.physics import gen_m66_elem, gen_detuning_elem, gen_quantdiff_elem

__all__ = ['fast_ring']


def rearrange(ring, split_inds=[]):
    iends = numpy.concatenate((split_inds, [len(ring)+1]))
    ibegs = numpy.concatenate(([0], iends[:-1]))
    iends = iends[iends != 0]
    ibegs = ibegs[ibegs != len(ring) + 1]
    all_rings = [ring[int(ibeg):int(iend)]
                 for ibeg, iend in zip(ibegs, iends)]

    for ring_slice in all_rings:
        # replace cavity with length > 0 with drift
        # set cavity length to 0 and move to start
        icav = get_refpts(ring_slice, RFCavity)
        for i in numpy.arange(len(icav)):
            cav_elem = ring_slice.pop(int(icav[i]))
            if cav_elem.Length != 0:
                cavdrift = Drift('CavDrift', cav_elem.Length)
                ring_slice.insert(icav[i], cavdrift)
                icav = icav + 1
                cav_elem.Length = 0.0
            ring_slice.insert(0, cav_elem)

        # merge all cavities with the same frequency
        icav = get_refpts(ring_slice, RFCavity)
        all_freq = numpy.array([ring_slice[ic].Frequency for ic in icav])
        all_volt = numpy.array([ring_slice[ic].Voltage for ic in icav])
        uni_freq = numpy.unique(all_freq)

        for ii in numpy.arange(len(uni_freq)):
            fr = uni_freq[ii]
            cavmsk = all_freq == fr
            vol = numpy.sum(all_volt[cavmsk])
            ring_slice[ii].Frequency = fr
            ring_slice[ii].Voltage = vol

        for pp in numpy.arange(len(icav)-len(uni_freq)):
            ring_slice.pop(len(uni_freq))

        ring_slice.insert(len(uni_freq), Marker('xbeg'))
        ring_slice.append(Marker('xend'))
    return all_rings


def merge_rings(all_rings):
    ringnorad = all_rings[0]
    if len(all_rings) > 1:
        for ringi in all_rings[1:]:
            ringnorad += ringi
    ringrad = ringnorad.deepcopy()
    if ringnorad.radiation:
        ringnorad.radiation_off(quadrupole_pass='auto')
    else:
        ringrad.radiation_on(quadrupole_pass='auto')
    return ringnorad, ringrad


def fast_ring(ring, split_inds=[]):
    all_rings = rearrange(ring, split_inds=split_inds)
    ringnorad, ringrad = merge_rings(all_rings)

    ibegs = get_refpts(ringnorad, 'xbeg')
    iends = get_refpts(ringnorad, 'xend')
    markers = numpy.sort(numpy.concatenate((ibegs, iends)))

    _, orbit4 = ringnorad.find_sync_orbit(dct=0.0, refpts=markers)
    _, orbit6 = ringrad.find_orbit6(refpts=markers)

    detuning_elem = gen_detuning_elem(ringnorad, orbit4[-1])
    detuning_elem_rad = detuning_elem.deepcopy()
    detuning_elem_rad.T1 = -orbit6[-1]
    detuning_elem_rad.T2 = orbit6[-1]

    for counter, ring_slice in enumerate(all_rings):

        ibeg = get_refpts(ring_slice, 'xbeg')[0]
        iend = get_refpts(ring_slice, 'xend')[0]
        cavs = ring_slice[:ibeg]
        ring_slice = ring_slice[ibeg:iend]
        ring_slice_rad = ring_slice.deepcopy()

        if ring_slice.radiation:
            ring_slice.radiation_off(quadrupole_pass='auto')
        else:
            ring_slice_rad.radiation_on(quadrupole_pass='auto')

        lin_elem = gen_m66_elem(ring_slice, orbit4[2*counter],
                                orbit4[2*counter+1])
        lin_elem.FamName = lin_elem.FamName + '_' + str(counter)

        qd_elem = gen_quantdiff_elem(ring_slice_rad, orbit=orbit6[2*counter])
        qd_elem.FamName = qd_elem.FamName + '_' + str(counter)
        lin_elem_rad = gen_m66_elem(ring_slice_rad, orbit6[2*counter],
                                    orbit6[2*counter+1])
        lin_elem_rad.FamName = lin_elem_rad.FamName + '_' + str(counter)

        [ringnorad.append(cav) for cav in cavs]
        ringnorad.append(lin_elem)
        [ringrad.append(cav) for cav in cavs]
        ringrad.append(lin_elem_rad)
        ringrad.append(qd_elem)

    ringnorad.append(detuning_elem)
    ringrad.append(detuning_elem_rad)
    del ringnorad[:markers[-1]+1]
    del ringrad[:markers[-1]+1]

    return ringnorad, ringrad
