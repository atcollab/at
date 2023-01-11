"""
Analytical wake functions
"""
import numpy
from at.constants import clight
from scipy.interpolate import interp1d


def convolve_wakefun(srange, w, sigs):
    """Convolution of a wake function with a pulse of rms
    length sigs, this is use to generate a wake potential
    that can be added to the output of EM code like GDFIDL"""
    sigt = sigs / clight #Convert sigmaz to sigmat
    
    # First we have to uniformly sample what is provided
    min_step = numpy.diff(srange)[0]
    s_out = numpy.arange(srange[0], srange[-1], min_step)
    sdiff = (s_out[-1]-s_out[0])
    
    # The actual number of points must be doubled to overlap
    # in the correct position. We fill with zeros everywhere 
    # else
    npoints = len(s_out)
    nt = npoints + npoints-1
    func = interp1d(srange, w, bounds_error=False, fill_value=0)
    wout = func(s_out)
    wout = numpy.append(wout, numpy.zeros(nt-len(wout)))
    
    # Perform the fft and multiply with gaussian in freq. 
    fftr = numpy.fft.fft(wout)
    f = numpy.fft.fftfreq(nt, d=min_step/clight)
    fftl = numpy.exp(-(f*2*numpy.pi*sigt)**2/2)
    wout = numpy.fft.ifft(fftr*fftl)
    
    # some manipulations to resample at the provided
    # srange and recenter the data
    wout = numpy.roll(wout, int(npoints/2))
    s_out = numpy.linspace(s_out[0], s_out[-1], nt)
    func = interp1d(s_out, wout, bounds_error=False, fill_value=0)
    wout = func(srange)
    return wout



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
