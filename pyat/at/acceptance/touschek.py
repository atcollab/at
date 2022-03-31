import at
import numpy
from ..lattice import Lattice, AtError
from scipy.special import iv
from scipy import integrate
from scipy.optimize import fsolve
from ..lattice.constants import clight, e_mass, qe, clight, _e_radius
from .acceptance import get_momentum_acceptance


__all__ = ['get_bunch_length_espread', 'get_lifetime']


def get_bunch_length_espread(ring, zn=None, bunch_curr=None, espread=None):
    """
    Solves the Haissinski formula and returns the bunch length and energy
    spread for given bunch current and Z/n. If both zn and bunch_curr are None,
    zero current case, otherwise both are needed for the calculation

    PARAMETERS
        ring              ring use for tracking

    KEYWORDS
        zn=None           Z/n for the full ring
        bunch_curr=None   Bunch current
        espread=None      Energy spread, if None use lattice parameter

    OUTPUT
        Bunch length, energy spread
    """
    def haissinski(x, cst):
        return x**3-x-cst

    ep = at.envelope_parameters(ring.radiation_on(copy=True))
    bl0 = ep.sigma_l
    if espread is None:
        espread = ep.sigma_e
    if zn is None and bunch_curr is None:
        bl = bl0
    elif zn is None or bunch_curr is None:
        raise AtError('Please provide both current and Z/n for bunch '
                      'length calculation')
    else:
        vrf = ring.get_rf_voltage()
        h = ring.harmonic_number
        etac = at.get_slip_factor(ring.radiation_off(copy=True))
        cs = numpy.cos(ep.phi_s)
        nus = ep.f_s/ring.revolution_frequency
        cst = (-0.5 * numpy.sqrt(numpy.pi) * bunch_curr * zn /
               (vrf * h * cs * (abs(etac) * espread / nus)**3))
        bl = numpy.float(bl0*fsolve(haissinski, (1,), args=(cst)))
    return bl, espread


def get_beam_sizes(ring, bunch_curr, zn=None, emitx=None,
                   sigs=None, sigp=None):
    if zn is None:
        bunch_curr = None
    if sigs is None or sigp is None:
        sigsi, sigpi = get_bunch_length_espread(ring, zn=zn,
                                                bunch_curr=bunch_curr,
                                                espread=sigp)
        if sigs is None:
            sigs = sigsi
        if sigp is None:
            sigp = sigpi
    if emitx is None:
        _, bd, _ = at.ohmi_envelope(ring.radiation_on(copy=True))
        emitx = bd.mode_emittances[0]
    return emitx, sigs, sigp


def int_piwinski(k, km, B1, B2):
    """
    Integrand of the piwinski formula
    """
    t = numpy.tan(k)**2
    tm = numpy.tan(km)**2
    fact = ((2*t+1)**2*(t/tm/(1+t)-1)/t + t - numpy.sqrt(t*tm*(1+t)) -
            (2+1/(2*t))*numpy.log(t/tm/(1+t)))
    if B2*t < 500:
        intp = fact * numpy.exp(-B1*t)*iv(0, B2*t)*numpy.sqrt(1+t)
    else:
        intp = fact * numpy.exp(B2*t-B1*t) / \
             numpy.sqrt(2*numpy.pi*B2*t)*numpy.sqrt(1+t)
    return intp


def get_lifetime(ring, emity, bunch_curr, emitx=None, sigs=None, sigp=None,
                 zn=None, momap=None, refpts=None, **kwargs):
    """
    Computes the touschek lifetime using the piwinski formula

    PARAMETERS
        ring            ring use for tracking
        emity           verticla emittance
        bunch_curr      bunch current

    KEYWORDS
        emitx=None      horizontal emittance
        sigs=None       rms bunch length
        sigp=None       energy spread
        zn=None         full ring Z/n
        momap=None      momentum aperture, has to be consistent with refpts
                        if provided the momentum aperture is not calculated
        refpts=None     refpts where the momentum aperture is calculated, the
                        default is to compute it for all elements in the ring
                        len(refpts)>2 is required
        resolution      minimum distance between 2 grid points, default=1.0e-3
        amplitude       max. amplitude of the grid or initial step in RECURSIVE
                        default = 0.1
        nturns=1024     Number of turns for the tracking
        dp=None         static momentum offset
        offset=None     initial orbit, default closed orbit
        bounds=None     Allows to define boundaries for the grid default
                        values are:
                        GridMode.GRID: ((-1,1),(0,1))
                        GridMode.RADIAL/RECURSIVE: ((0,1),(pi,0))
        grid_mode       mode for the gird default GridMode.RADIAL
        use_mp=False    Use python multiprocessing (patpass, default use
                        lattice_pass). In case multi-processing is not
                        enabled GridMode is forced to
                        RECURSIVE (most efficient in single core)
        verbose=True    Print out some inform
        epsabs, epsrel  integral absolute and relative tolerances

    OUTPUT
        Returns the touschek lifetime in seconds, the momentum aperture
        and the refpts at which the aperture was computed
    """

    epsabs = kwargs.pop('epsabs', 1.0e-16)
    epsrel = kwargs.pop('epsrel', 1.0e-12)

    emitx, sigs, sigp = get_beam_sizes(ring, bunch_curr, zn=zn,
                                       emitx=emitx, sigs=sigs, sigp=sigp)
    if refpts is None:
        refpts = range(len(ring))
    else:
        refpts = ring.uint32_refpts(refpts)
        assert ring.refcount(refpts) < 3, \
            'len(refpts) > 2 required for lifetime calculation'
            

    if momap is None:
        resolution = kwargs.pop('resolution', 1.0e-3)
        amplitude = kwargs.pop('amplitude', 0.1)
        kwargs.update({'refpts': refpts})
        momap, _, _ = ring.get_momentum_acceptance(resolution,
                                                   amplitude, **kwargs)
    else:
        print('Using user input momentum aperture: refpts and momap '
              'have to be coherent')
        assert len(momap) == len(refpts), \
            'momap and refpts have different lengths'

    spos = numpy.diff(ring.get_s_pos(refpts))
    ma, rp = momap[:-1], refpts[:-1]

    nc = bunch_curr/ring.revolution_frequency/qe
    beta2 = ring.beta*ring.beta
    gamma2 = ring.gamma*ring.gamma

    emit = numpy.array([emitx, emity])
    _, _, ld = ring.get_optics(refpts=rp)
    bxy = ld.beta
    bxy2 = bxy*bxy
    axy = ld.alpha
    dxy = ld.dispersion[:, [0, 2]]
    dxy2 = dxy*dxy
    dpxy = ld.dispersion[:, [1, 3]]
    sigb = numpy.sqrt(bxy*emit)
    sigb2 = sigb*sigb
    sigp2 = sigp*sigp
    sig = numpy.sqrt(sigb2+sigp2*dxy2)
    sig2 = sig*sig
    dt = dxy*axy+dpxy*bxy
    dt2 = dt*dt
    sigt = numpy.sqrt(sig2 + sigp2*dt2)
    sigt2 = sigt*sigt
    sigh2 = 1/(1/sigp2 + numpy.sum((dxy2+dt2)/(sigb2), axis=1))
    sigh = numpy.sqrt(sigh2)

    bs = bxy2/sigb2*(1-(sigh2*(dt2/sigb2).T).T)
    bg2i = 1/(2*beta2*gamma2)
    B1 = bg2i*numpy.sum(bs, axis=1)
    B2sq = bg2i*bg2i*(numpy.diff(bs, axis=1).T**2 +
                      sigh2*sigh2*numpy.prod(bxy2*dt2, axis=1) /
                      numpy.prod(sigb2*sigb2, axis=1))
    B2 = numpy.squeeze(numpy.sqrt(B2sq))

    invtl = numpy.zeros(2)
    for i in range(2):
        dpp = ma[:, i]
        um = beta2*dpp*dpp
        km = numpy.arctan(numpy.sqrt(um))

        val = numpy.zeros(len(rp))
        for ii in range(len(rp)):
            args = (km[ii], B1[ii], B2[ii])
            val[ii], _ = integrate.quad(int_piwinski, args[0], numpy.pi/2,
                                        args=args, epsabs=epsabs,
                                        epsrel=epsrel)

        val *= (_e_radius**2*clight*nc /
                (8*numpy.pi*(gamma2)*sigs *
                 numpy.sqrt(numpy.prod(sig2, axis=1) -
                 sigp2*sigp2*numpy.prod(dxy2, axis=1))) *
                2*numpy.sqrt(numpy.pi*(B1*B1-B2*B2)))

        invtl[i] = sum(val*spos.T)/sum(spos)
    tl = 1/numpy.mean(invtl)

    return tl, momap, refpts


Lattice.get_bunch_length_espread = get_bunch_length_espread
Lattice.get_lifetime = get_lifetime
