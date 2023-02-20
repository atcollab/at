"""
Analytical wake functions
"""
import numpy
from at.constants import clight
from scipy.interpolate import interp1d


def convolve_wakefun_fft(srange, w, sigs, xr=2*0.09, npoints=500000):
    """Convolution of a wake function with a pulse of rms
    length sigs, this is use to generate a wake potential
    that can be added to the output of EM code like GDFIDL"""

    sigt = sigs / clight  # Convert sigmaz to sigmat
    nt = npoints + npoints - 1
    
    # First we have to uniformly sample what is provided

    s_out = numpy.linspace(-xr, xr, npoints)
    min_step = numpy.diff(s_out)[0]
    
    # The actual number of points must be doubled to overlap
    # in the correct position. We fill with zeros everywhere
    # else
    #npoints = len(s_out)

    func = interp1d(srange, w, bounds_error=False, fill_value=0)
    wout = func(s_out)
    wout = numpy.append(wout, numpy.zeros(nt-len(wout)))
    # Perform the fft and multiply with gaussian in freq.
    fftr = numpy.fft.fft(wout)
    
    #f = numpy.fft.fftshift(numpy.linspace(-(npoints-1)/(4*xr/clight),(npoints-1)/(4*xr/clight),nt))
   
    f = numpy.fft.fftfreq(nt, d=min_step/clight)
    fftl = numpy.exp(-(f*2*numpy.pi*sigt)**2/2)
    wout = numpy.fft.ifft(fftr*fftl)

    # some manipulations to resample at the provided
    # srange and recenter the data
    wout = numpy.roll(wout, int(npoints/2))
    s_out = numpy.linspace(-2*xr, 2*xr, nt)
    func = interp1d(s_out, wout, bounds_error=False, fill_value=0)
    wout = func(srange)
    return wout
    
def convolve_wakefun_numpy(srange, w, sigs, gauss_sigma=10):
    """Convolution of a wake function with a pulse of rms
    length sigs, this is use to generate a wake potential
    that can be added to the output of EM code like GDFIDL"""

    def _gauss(s):
        ampl = 1 / (numpy.sqrt(2 * numpy.pi) * sigs)
        expon = numpy.exp(-s**2 / (2 * sigs**2))
        return ampl * expon
        
    ds = numpy.diff(srange)[-1]
    s_gauss = numpy.arange(-gauss_sigma*sigs, gauss_sigma*sigs+1e-15, ds)
    gauss = _gauss(s_gauss)

    conv = numpy.convolve(gauss, w, mode='full')*ds
    s_offset = gauss_sigma * sigs  - numpy.amin(srange)
    s_conv = numpy.arange(len(conv)) * ds - s_offset 

    ifun = interp1d(s_conv, conv, bounds_error=False, fill_value=0)
    conv_wf = ifun(srange)
    return conv_wf
    
    
def long_resonator_wf(srange, frequency, qfactor, rshunt, beta):
    """Define the wake function (longitudinal) of a resonator
    with the given parameters according to Alex Chao's resonator
    model (Eq. 2.82) and definitions of the resonator in HEADTAIL.
    """

    omega = 2 * numpy.pi * frequency
    alpha = omega / (2 * qfactor)
    omegabar = numpy.sqrt(numpy.abs(omega**2 - alpha**2))
    dt = -srange/(beta * clight)
    if qfactor > 0.5:
        wake = (-(numpy.sign(dt) - 1) * rshunt * alpha *
                numpy.exp(alpha * dt) * (numpy.cos(omegabar * dt) +
                alpha / omegabar * numpy.sin(omegabar*dt)))
    elif qfactor == 0.5:
        wake = (-(numpy.sign(dt) - 1) * rshunt * alpha *
                numpy.exp(alpha * dt) * (1. + alpha * dt))
    elif qfactor < 0.5:
        wake = (-(numpy.sign(dt) - 1) * rshunt * alpha *
                numpy.exp(alpha * dt) * (numpy.cosh(omegabar * dt) +
                alpha / omegabar * numpy.sinh(omegabar * dt)))
    return wake


def transverse_resonator_wf(srange, frequency, qfactor, rshunt,
                            yokoya_factor, beta):
    """Define the wake function (transverse) of a resonator
    with the given parameters according to Alex Chao's
    resonator model (Eq. 2.82) and definitions of the resonator
    in HEADTAIL.
    """

    omega = 2 * numpy.pi * frequency
    alpha = omega / (2 * qfactor)
    omegabar = numpy.sqrt(numpy.abs(omega**2 - alpha**2))
    dt = -srange/(beta * clight)
    if qfactor > 0.5:
        wake = (-(numpy.sign(dt) - 1) * yokoya_factor * rshunt * omega**2 /
                (qfactor * omegabar) *
                numpy.exp(alpha*dt) * numpy.sin(omegabar*dt))
    elif qfactor == 0.5:
        wake = (-(numpy.sign(dt) - 1) * yokoya_factor * rshunt * omega**2 /
                qfactor *
                numpy.exp(alpha * dt) * dt)
    else:
        wake = (-(numpy.sign(dt) - 1) * yokoya_factor * rshunt * omega**2 /
                (qfactor * omegabar) *
                numpy.exp(alpha*dt) * numpy.sinh(omegabar*dt))
    return wake


def transverse_reswall_wf(srange, yokoya_factor, length, rvac, conduct, beta):
    """Define the wake function (transverse) of a resistive wall with the given
    parameters according to Alex Chao's RW model (2.53) and definitions used in
    HEADTAIL
    """

    z0 = 119.9169832 * numpy.pi
    dt = -srange/(beta * clight)
    wake = (yokoya_factor * (numpy.sign(dt) - 1) / 2. *
            beta * length / numpy.pi / rvac**3 *
            numpy.sqrt(-z0 * clight / conduct / numpy.pi / dt))

    wake[srange <= 0] = 0.0

    return wake
