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

    I_cav = get_refpts(ring, RFCavity)

    ring_temp = ring.deepcopy()

    for i in numpy.arange(len(I_cav)):
        if ring_temp[I_cav[i]].Length != 0:
            CavElement = ring_temp[I_cav[i]]
            CavDrift = Drift('CavDrift',CavElement.Length/2)

            ring_temp[I_cav[i]] = CavDrift
            CavElement.Length = 0
            ring_temp.insert(I_cav[i]+1,CavElement)
            ring_temp.insert(I_cav[i]+2,CavDrift)
            I_cav[i+1:] = I_cav[i+1:]+2


    all_rings = []

    iends = numpy.concatenate((split_inds, [len(ring_temp)+1]))
    ibegs = numpy.concatenate(([0], iends[:-1]))

    for ibeg, iend in zip(ibegs, iends):
        ring_slice = ring_temp.deepcopy()
        ring_slice = ring_slice[int(ibeg):int(iend)]
        I_cav1 = get_refpts(ring_slice, RFCavity)
        I_cav2 = get_refpts(ring_slice, 'RingParam')
        I_cav = numpy.array(numpy.concatenate((I_cav1, I_cav2)), dtype=int)

        for it in numpy.arange(len(I_cav)):
            cav_elem = ring_slice.pop(int(I_cav[it]))
            ring_slice.insert(0,cav_elem)
        ring_slice.insert(len(I_cav), Marker('xbeg'))
        ring_slice.append(Marker('xend'))

        all_rings.append(ring_slice)


    return all_rings

def merge_rings(all_rings):

    if len(all_rings)>=2:
        main_ring = all_rings[0].deepcopy()
        for it in numpy.arange(1, len(all_rings)):
            rc = all_rings[it].deepcopy()
            for itt in numpy.arange(len(rc)):
                el = rc.pop(0)
                main_ring.append(el)

    else:
        main_ring = all_rings[0]

    return main_ring


def fast_ring(ring, split_inds):
    all_rings = rearrange(ring, split_inds)
    ring_merge = merge_rings(all_rings.copy())

    ibegs = get_refpts(ring_merge, 'xbeg')
    iends = get_refpts(ring_merge, 'xend')
    markers = numpy.sort(numpy.concatenate((ibegs, iends))) 

    _, orbit4 = ring_merge.find_orbit4(refpts=markers)
    _, orbit6 = ring_merge.find_orbit6(refpts=markers)

    detuning_elem = gen_detuning_elem(ring_merge, orbit4[0])
    detuning_elem_rad = gen_detuning_elem(ring_merge, orbit6[0])

    ring_rg = ring.deepcopy()
    ring_rg.clear()
    ring_rg_rad = ring_rg.deepcopy()
    ring_rg_rad.radiation_on(quadrupole_pass='auto')

    for counter in numpy.arange(len(all_rings)):
        ring_slice = all_rings[counter]
        ibeg = get_refpts(ring_slice, 'xbeg')[0]
        iend = get_refpts(ring_slice, 'xend')[0]

        cavs = ring_slice[:ibeg]
        lin_elem = gen_m66_elem(ring_slice[ibeg:iend], orbit4[2*counter], orbit4[2*counter+1], radiation=False)
        lin_elem.FamName = lin_elem.FamName + '_' + str(counter)
        [ring_rg.append(cav) for cav in cavs]
        ring_rg.append(lin_elem)


        ring_slice_rad = ring_slice.deepcopy()
        ring_slice_rad.radiation_on(quadrupole_pass='auto')
        cavs = ring_slice_rad[:ibeg]
        

        qd_elem = gen_quantdiff_elem(ring_slice_rad, orbit=orbit6[2*counter])
        qd_elem.FamName = qd_elem.FamName + '_' + str(counter)
        lin_elem_rad = gen_m66_elem(ring_slice_rad[ibeg:iend], orbit6[2*counter], orbit6[2*counter+1], radiation=True)
        lin_elem_rad.FamName = lin_elem_rad.FamName + '_' + str(counter)

        [ring_rg_rad.append(cav) for cav in cavs]
        ring_rg_rad.append(lin_elem_rad)
        ring_rg_rad.append(qd_elem)

    ring_rg.append(detuning_elem)
    ring_rg_rad.append(detuning_elem_rad)

    return ring_rg, ring_rg_rad


