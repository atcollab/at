"""
Coupled or non-coupled 4x4 non-linear motion
"""
import numpy
from math import sqrt, atan2, pi
from at.lattice import Lattice, check_radiation, uint32_refpts, get_s_pos, \
    bool_refpts
from at.tracking import lattice_pass
from at.physics import HarmonicAnalysis, get_tune, linopt, find_orbit4, get_tunes_harmonic
from at.lattice import Element

__all__ = ['tune_amplitude', 'detuning', 'gen_nonlin_elem']

def tune_amplitude(ring, ampl, pl, orbit, tunes, betas):
    """
    Compute the change in tune as a function of amplitude

    tunes_h, tunes_v = tune_amplitude(ring, ampl, pl, orbit)

    """
    tuneh, tunev = tunes
    betah, betav = betas

    siza = numpy.size(ampl)
    nampl = numpy.prod(siza)
    p0 = numpy.zeros((6,nampl)) 
    p0[pl,:] = ampl
    p0[0,0] = 1e-9
    p0[2,0] = 1e-9
    if pl==0:
        p0[2,:] = 1e-9 
    else:
        p0[0,:] = 1e-9

    for i in numpy.arange(len(ampl)):
        p0[:,i] += orbit

    p1 = lattice_pass(ring, p0, refpts=len(ring), nturns=512)

    for i in numpy.arange(len(ampl)):
        for ii in numpy.arange(6): #plane
            p1[ii,i,0,:] -= orbit[ii]


    tunetrack_h = numpy.zeros((len(ampl)))
    tunetrack_v = numpy.zeros((len(ampl)))
    for i in numpy.arange(len(ampl)):
        tunetrack_h[i] = detuning_find_tune(p1[0,i,0,:] -1j*p1[1,i,0,:]*betah, tuneh)
        tunetrack_v[i] = detuning_find_tune(p1[2,i,0,:] -1j*p1[3,i,0,:]*betav, tunev)
        tuneh = tunetrack_h[i]
        tunev = tunetrack_v[i]

    return [tunetrack_h, tunetrack_v]

def detuning_find_tune(tbt, tune, window=0.05):

    ha = HarmonicAnalysis(tbt)
    ha_tune, ha_amp = ha.laskar_method(num_harmonics=20)
    ha_amp = numpy.abs(numpy.array(ha_amp))
    ha_tune = numpy.array(ha_tune)

    msk = ha_tune < 0.5
    ha_tune = ha_tune[msk]
    ha_amp = ha_amp[msk]

    msk2 = numpy.logical_and(ha_tune < tune + window, ha_tune > tune - window)
    ha_tune = ha_tune[msk2]
    ha_amp = ha_amp[msk2] 
    return ha_tune[0]

def detuning(ring, xm, zm, orbit):
    #python
    #r1/1e5 gives 
    #array([[ 0.96051517, -1.07237005, -1.07769544,  0.49889636]])

    #matlab
    #    #  1.0e+05 *
    #
    #    0.9887
    #   -1.0518
    #   -1.0537
    #    0.5095

    [lindata0, tunes, xsi, lindata] = ring.linopt(dp=0, refpts=len(ring), get_chrom=True)

    gamma=(1+lindata.alpha*lindata.alpha)/lindata.beta

    betas = lindata.beta[0]

    x2=numpy.linspace(0,xm*xm,11)
    z2=numpy.linspace(0,zm*zm,11)
    [nuxx,nuzx]=tune_amplitude(ring,numpy.sqrt(x2),0, orbit, tunes, betas)
    [nuxz,nuzz]=tune_amplitude(ring,numpy.sqrt(z2),2, orbit, tunes, betas)

    rx1 = numpy.ones((2,len(nuxx)))
    rx1[0,:] = nuxx-nuxx[0]
    rx1[1,:] = nuzx-nuzx[0]

    rz1 = numpy.ones((2,len(nuxz)))
    rz1[0,:] = nuxz-nuxz[0]
    rz1[1,:] = nuzz-nuzz[0]

    z21 = numpy.ones((2,len(z2)))
    z21[0,:] = z2
    z21[1,:] = z2

    x21 = numpy.ones((2,len(x2)))
    x21[0,:] = x2
    x21[1,:] = x2
    
    rx=numpy.array([numpy.dot(rx1,x2)/numpy.dot(x21,x2)/gamma[0,0]])
    rz=numpy.array([numpy.dot(rz1,z2)/numpy.dot(z21,z2)/gamma[0,1]])

    r=2*numpy.concatenate((rx,rz),axis=1)
    return r


def tunes_vs_amp(ring, xamp=None, yamp=None, zamp=None, dp=0, window=0.1, nturns=512):

    def _gen_part(ring,amp,dim,orbit,nturns):
        part = numpy.array([orbit,]*len(amp)).T+1.0e-6
        part[dim,:] += amp
        part = lattice_pass(ring, part, nturns=nturns)
        sh = part.shape
        partx = numpy.reshape(part[0,:],(sh[1],sh[3]))
        party = numpy.reshape(part[2,:],(sh[1],sh[3]))
        return partx,party
        
    orbit,_ = find_orbit4(ring)
    q0 = get_tune(ring) 

    xtunes = [] 
    ytunes = []
    ztunes = []   

    if zamp is not None: 
        for z in zamp:
            ztunes.append(get_tune(ring,dp=dp+z))
    if xamp is not None:
        partx,party = _gen_part(ring,xamp,0,orbit,nturns)
        qx = get_tunes_harmonic(partx, 'laskar', fmin=q0[0]-window, fmax=q0[0]+window)
        qy = get_tunes_harmonic(party, 'laskar', fmin=q0[1]-window, fmax=q0[1]+window)
        xtunes = numpy.vstack((qx,qy)).T
    if yamp is not None:
        partx,party = _gen_part(ring,yamp,2,orbit,nturns)
        qx = get_tunes_harmonic(partx, 'laskar', fmin=q0[0]-window, fmax=q0[0]+window)
        qy = get_tunes_harmonic(party, 'laskar', fmin=q0[1]-window, fmax=q0[1]+window)
        ytunes = numpy.vstack((qx,qy)).T

    return numpy.array(xtunes), numpy.array(ytunes), numpy.array(ztunes)
            


def detuning2(ring, xm=1.0e-4, ym=1.0e-4, npoints=3, dp=0):

    lindata0, _, _, _ = linopt(ring,dp=dp)
    gamma=(1+lindata0.alpha*lindata0.alpha)/lindata0.beta
   
    x=numpy.linspace(-xm, xm, npoints)
    y=numpy.linspace(-ym, ym, npoints)
    x2 = x*x
    y2 = y*y
   
    q_dx,q_dy,_ = tunes_vs_amp(ring, xamp=x, yamp=y)
    fx = numpy.polyfit(x2,q_dx,1)
    fy = numpy.polyfit(y2,q_dy,1)

    q0 = [fx[1,0],fx[1,1],fy[1,0],fy[1,1]]
    q1 = [2*fx[0,0]/gamma[0],2*fx[0,1]/gamma[0],2*fy[0,0]/gamma[1],2*fy[0,1]/gamma[1]]

    return numpy.array(q0),numpy.array(q1)




def gen_nonlin_elem(ring, r1):
    [lindata0, tunes, xsi, lindata] = ring.linopt(dp=0, refpts=len(ring), get_chrom=True)
    [orb04, orb4] = ring.find_orbit4(refpts=[0,len(ring)])   


    nonlin_elem = Element(  'NonLinear',
                            PassMethod='DeltaQPass',
                            Betax=lindata.beta[0][0],
                            Betay=lindata.beta[0][1], 
                            Alphax=lindata.alpha[0][0],
                            Alphay=lindata.alpha[0][1], 
                            Qpx=xsi[0],
                            Qpy=xsi[1], 
                            A1=r1[0][0], 
                            A2=r1[0][1], 
                            A3=r1[0][3], 
                            T1=-orb4[-1],
                            T2= orb4[-1])

    return nonlin_elem

Lattice.tune_amplitude = tune_amplitude
Lattice.detuning = detuning
Lattice.gen_nonlin_elem = gen_nonlin_elem

