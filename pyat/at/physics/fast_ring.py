"""
Functions relating to fast_ring
"""
import numpy
from math import sqrt, atan2, pi
from at.lattice import Lattice, check_radiation, uint32_refpts, get_s_pos, \
    bool_refpts, get_refpts
from at.tracking import lattice_pass
from at.lattice import Element, RFCavity, Marker, Drift, Bend, M66
from at.physics import HarmonicAnalysis, detuning
from at.physics import find_elem_m66, find_orbit6, find_orbit4
from at.physics import find_mpole_raddiff_matrix, linopt
from at.physics import jmat

__all__ = ['fast_ring']

def rearrange(ring):

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

    I_cav = get_refpts(ring_temp, RFCavity)
    for it in numpy.arange(len(I_cav)):
        cav_elem = ring_temp.pop(I_cav[it])
        ring_temp.insert(it,cav_elem)
    ring_temp.insert(len(I_cav), Marker('xbeg'))
    ring_temp.insert(len(ring_temp)+1, Marker('xend'))
    return ring_temp


def symplectify(M):
    '''
    symplectify makes a matrix more symplectic
    follow Healy algorithm as described by McKay
    BNL-75461-2006-CP
    '''
    J = jmat(3)
    V=J*(numpy.identity(6)-M)*numpy.linalg.inv(numpy.identity(6)+M)
    #V should be almost symmetric.  Replace with symmetrized version.
    W=(V+V.T)/2
    #Now reconstruct M from W
    MS=(numpy.identity(6)+J*W)*numpy.linalg.inv(numpy.identity(6)-J*W)
    return MS


def gen_quantdiff_elem(ring_slice):
    
    ring_slice.radiation_on(quadrupole_pass='auto')
    rad_inds = get_rad_indexes(ring_slice)
    quantumDiffMatrix = quantumDiff(ring_slice,rad_inds)

    '''
    #lmat does Cholesky decomp of dmat unless diffusion is 0 in
    #vertical.  Then do chol on 4x4 hor-long matrix and put 0's
    #in vertical
    dmat = quantumDiffMatrix
    try:
        lmat66 = numpy.linalg.cholesky(dmat)
    except:
        lmc = numpy.array([numpy.concatenate((dmat[0][0:2], dmat[0][4:])), 
                        numpy.concatenate((dmat[1][0:2], dmat[1][4:])),
                        numpy.concatenate((dmat[4][0:2], dmat[4][4:])),
                        numpy.concatenate((dmat[5][0:2], dmat[5][4:]))])
        lm=[chol(dmat([1 2 5 6],[1 2 5 6])) zeros(4,2);zeros(2,6)];
        lmat66=lm([1 2 5 6 3 4],[1 2 5 6 3 4]);
    
    lmatp=lmat66';
    '''
    diff_elem = Element('Diffusion',
                        Lmatp=numpy.asfortranarray(quantumDiffMatrix),
                        PassMethod='QuantDiffPass')

    return diff_elem

def quantumDiff(ring, rad_inds):
    NumElements=len(ring);

    #[mring, ms, orbit] = findm66(ring,1:NumElements+1);

    #orb=num2cell(linepass(elems, orb0, 1:NumElements),1)';
    orb0,orb=find_orbit6(ring, refpts=numpy.arange(len(ring)))


    #zr={zeros(6,6)};
    #B=zr(ones(NumElements,1));   # B{i} is the diffusion matrix of the i-th element
    B = []
    for i in numpy.arange(NumElements):
        B.append(numpy.zeros((6,6)))
    B = numpy.array(B)

    # calculate Radiation-Diffusion matrix B for elements with radiation
    for ind in rad_inds:
        B[ind] = find_mpole_raddiff_matrix(ring[ind], orb[ind], 6e9)

    #B(radindex)=cellfun(@findmpoleraddiffmatrix,...
    #    elems(radindex),orb(radindex),'UniformOutput',false);

    # Calculate cumulative Radiation-Diffusion matrix for the ring
    BCUM = numpy.zeros((6,6))
    # Batbeg{i} is the cumulative diffusion matrix from
    # 0 to the beginning of the i-th element
    #Batbeg=[numpy.zeros((6,6))]cellfun(@cumulb,elems,orb,B,'UniformOutput',false)]
    Batbeg = [numpy.zeros((6,6))]
    for i, elem in enumerate(ring):
        m=find_elem_m66(elem,orbit=orb[i]);

        BCUM = m*BCUM*m.T + B[i]

        Batbeg.append(BCUM)


    DiffCum = BCUM

    DiffMat=(DiffCum + DiffCum.T)/2   

    #Lmat=chol((DiffCum+DiffCum')/2);
    return DiffMat



def get_rad_indexes(ring):
    arr = numpy.ones(len(ring),dtype=bool)
    for i in numpy.arange(len(ring)):
        if 'Rad' not in ring[i].PassMethod:
            arr[i] = False
    inds = numpy.arange(len(ring))[arr]
    return inds


def fast_ring(ring):

    [orb04, orb4] = find_orbit4(ring, refpts=[0,len(ring)])   
    [orb06, orb6] = find_orbit6(ring, refpts=[0,len(ring)])   


    ring.radiation_off(quadrupole_pass='auto')

    #[lindata0, tunes, xsi, lindata] = linopt(ring, dp=0, refpts = len(ring), get_chrom=True)


    ring_rg = rearrange(ring)

    ibeg = get_refpts(ring_rg, 'xbeg')[0]
    iend = get_refpts(ring_rg, 'xend')[0]

    ring_slice = ring_rg[ibeg:iend]
    rg = ring_rg[:ibeg]
    rgrad = rg.deepcopy()

    lin_elem = gen_lin_elem(ring_slice, orb4, radiation=False)
    lin_elem_rad = gen_lin_elem(ring_slice, orb6, radiation=True)

    rg.append(lin_elem)
    rgrad.append(lin_elem_rad)

    xm = 1e-4
    zm = 1e-4
    r1 = detuning(ring, xm, zm, orb4)
    nonlin_elem = gen_nonlin_elem(ring, r1, orb4)
    nonlin_elem_rad = gen_nonlin_elem(ring, r1, orb6)
    rg.append(nonlin_elem)
    rgrad.append(nonlin_elem_rad)

    diff_elem = gen_quantdiff_elem(ring) #complaining about no cavities
    #so using the full ring
    rgrad.append(diff_elem)

    return rg, rgrad

Lattice.fast_ring = fast_ring

