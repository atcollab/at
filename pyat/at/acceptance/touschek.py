import at
import numpy
from scipy.constants import physical_constants
from scipy.special import iv
from scipy import integrate
from tlt_da.compute_aperture import get_aperture_array_multi
_qe = physical_constants['elementary charge'][0]
_re = physical_constants['classical electron radius'][0]
_clight = physical_constants['speed of light in vacuum'][0]
_erm = physical_constants['electron mass energy equivalent in MeV'][0]*1.0e6


def blength_espread(ring, zn, bcurr):
    ringtmp = ring.radiation_off(copy=True)
    rad_param = ringtmp.radiation_parameters()
    vrf = ring.get_rf_voltage()
    u0 = rad_param.U0
    e0 = rad_param.E0
    h = ring.get_rf_harmonic_number()
    alphac = rad_param.alphac
    espread = rad_param.sigma_e
    circ = ring.get_s_pos(len(ring))[0]
    phis = rad_param.phi_s
    nus = rad_param.f_s/_clight*circ
    delta = -(2*numpy.pi*bcurr*zn)/(vrf*h*numpy.cos(phis)*(alphac*espread/nus)**3)
    q=delta/(4*numpy.sqrt(numpy.pi))
    blg = (2/3)**(1/3)/(9*q + numpy.sqrt(3)*numpy.sqrt(-4+27*q**2))**(1/3) \
           + (9*q + numpy.sqrt(3)*numpy.sqrt(-4+27*q**2))**(1/3)/(2**(1/3)*3**(2/3))
    bl= espread*(circ * alphac)/(2 * numpy.pi * nus )
    bl *= blg
    return bl, espread


def get_tlt(ring, bcurr, epsx, epsy, sigs, sigp, refpts=None, iends=None, maxdp=0.06, tol=1.0e-3, nturns=1024):
    momap = get_aperture_array_multi(ring, 4, maxdp, tol=tol, refpts=refpts, nturns=nturns)
    return tlt_piwinski(ring, bcurr, epsx, epsy, sigs, sigp, momap, refpts=refpts, iends=iends)/3600


def tlt_intpiw(k, km, B1, B2): 
    t=numpy.tan(k)**2
    tm=numpy.tan(km)**2
    if B2*t<500:
        I=((2*t+1)**2*(t/tm/(1+t)-1)/t + t -numpy.sqrt(t*tm*(1+t)) -(2+1/(2*t))*numpy.log(t/tm/(1+t)))* \
          numpy.exp(-B1*t)*iv(0,B2*t)*numpy.sqrt(1+t)
    else:
        I=((2*t+1)**2*(t/tm/(1+t)-1)/t + t -numpy.sqrt(t*tm*(1+t)) -(2+1/(2*t))*numpy.log(t/tm/(1+t)))* \
          numpy.exp(B2*t-B1*t)/numpy.sqrt(2*numpy.pi*B2*t)*numpy.sqrt(1+t)
    return I


def tlt_piwinski(ring, bcurr, emitx, emity, sigs, sigp, momap, refpts=None, iends=None):
    if refpts is None:
        refpts=range(len(ring))
    if iends is None:
        iends = (0,len(ring))
    ii,ie=iends

    si= ring.get_s_pos(ii)[0]
    se= ring.get_s_pos(ie)[0]
    spos = ring.get_s_pos(refpts)
    l=numpy.zeros(len(refpts))
    l[0]=spos[1]-si
    for ii in range(1,len(refpts)-1):
        l[ii]=spos[ii+1]-spos[ii]
    l[len(refpts)-1]=se-spos[len(refpts)-1] 
  
    tlcol= numpy.zeros(2)

    circ = ring.get_s_pos(len(ring))[0]
    energy = ring.energy
    nc = bcurr/(_clight/circ)/_qe
    gamma = energy/_erm
    beta=numpy.sqrt(1-1/gamma**2)
    beta2 = beta*beta
    gamma2 = gamma*gamma

    ld0,bd,ld = ring.get_optics(refpts=refpts,get_chrom=True)
    bx = ld.beta[:,0]
    by = ld.beta[:,1]
    ax = ld.alpha[:,0]
    ay = ld.alpha[:,1]
    dx = ld.dispersion[:,0]
    dy = ld.dispersion[:,2]
    dpx = ld.dispersion[:,1]
    dpy = ld.dispersion[:,3]

    sigxb=numpy.sqrt(emitx*bx)
    sigyb=numpy.sqrt(emity*by)   
    sigx=numpy.sqrt(emitx*bx+sigp*dx*sigp*dx)
    sigy=numpy.sqrt(emity*by+sigp*dy*sigp*dy) 

    dtx=dx*ax+dpx*bx
    dty=dy*ay+dpy*by

    bx2 = bx*bx
    by2 = by*by
    sigx2 = sigx*sigx
    sigy2 = sigy*sigy
    sigp2=sigp*sigp
    dx2=dx*dx
    dy2=dy*dy
    dtx2=dtx*dtx
    dty2=dty*dty
    sigxb2=sigxb*sigxb
    sigyb2=sigyb*sigyb

    sigtx2=sigx2+sigp2*dtx2
    sigty2=sigy2+sigp2*dty2
    sigtx=numpy.sqrt(sigtx2)
    sigty=numpy.sqrt(sigty2)

    sighinv2=1/(sigp2) +(dx2+dtx2)/(sigxb2) + (dy2+dty2)/(sigyb2)
    sigh2 = 1/sighinv2    
    sigh=numpy.sqrt(sigh2)

    B1=1/(2*beta2*gamma2)*(bx2/sigxb2*(1-sigh2*dtx2/sigxb2) + by2/sigyb2*(1-sigh2*dty2/sigyb2))  
    B2sq=1/(2*beta2*gamma2)**2*(bx2/sigxb2*(1-sigh2*dtx2/sigxb2) - by2/sigyb2*(1-sigh2*dty2/sigyb2))**2+ \
         (sigh2*sigh2*bx2*by2*dtx2*dty2)/(beta2*beta2*gamma2*gamma2*sigxb2*sigxb2*sigyb2*sigyb2)   
    B2=numpy.sqrt(B2sq)
  
    val=numpy.zeros(len(refpts))
 
    for i in range(2): 
        dpp = momap[:,i]          
        um=beta2*dpp*dpp  
        em=bx2*sigx2/(beta2*gamma2*sigxb2*sigtx2)*um  
        km=numpy.arctan(numpy.sqrt(um)) 
        FpiWfact=numpy.sqrt(numpy.pi*(B1**2-B2**2))*um
                
        for ii in range(len(refpts)):
            val[ii],_ = integrate.quad(tlt_intpiw, km[ii], numpy.pi/2,args=(km[ii], B1[ii], B2[ii]),limit=1000,points=range(1000))  
             
        frontfact=(_re**2*_clight*nc)/(8*numpy.pi*(gamma2)*sigs*numpy.sqrt(sigx2*sigy2-sigp2*sigp2*dx2*dy2)*um)*2*FpiWfact
        contributionsTL=frontfact*val          
        tlcol[i]=1/(1/sum(l)*sum(contributionsTL*l.T))
    
    return len(tlcol)/sum(1/tlcol)



