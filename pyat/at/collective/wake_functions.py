"""
Analytical wake functions
"""
import numpy
from at.lattice.constants import clight
from scipy.interpolate import interp1d


def convolve_wakefun(srange, w, sigs):
    """Convolution of a wake function with a pulse of rms
    length sigs, this is use to generate a wake potential
    that can be added to the output of EM code like GDFIDL"""
    sigt = sigs/clight
    min_step = numpy.diff(srange)
    t_out = numpy.arange(srange[0], srange[-1], min_step)
    sdiff = t_out[-1]-t_out[0]
    npoints = len(t_out)
    nt = npoints+npoints-1
    func = interp1d(srange, w, bounds_error=False, fill_value=0)
    wout = func(t_out)
    wout = numpy.append(wout, numpy.zeros(nt-len(wout)))
    fftr = numpy.fft.fft(wout)
    f = numpy.fft.fftshift(numpy.linspace(-(npoints-1)/sdiff,
                           (npoints-1)/sdiff, nt))
    fftl = numpy.exp(-(f*2*numpy.pi*sigt)**2/2)
    wout = numpy.fft.ifft(fftr*fftl)
    wout = numpy.roll(wout, int(npoints/2))
    t_out = numpy.linspace(t_out[0], t_out[-1], nt)
    func = interp1d(t_out, wout, bounds_error=False, fill_value=0)
    wout = func(srange)
    return wout


def long_resonator(srange, frequency, qfactor, rshunt, beta):
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


def transverse_resonator(srange, frequency, qfactor, rshunt,
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
        wake = (-(numpy.sign(dt) - 1) * yokoya_factor * rshunt * omega**2 / (qfactor *
                omegabar) * numpy.exp(alpha*dt) * numpy.sin(omegabar*dt))
    elif qfactor == 0.5:
        wake = (-(numpy.sign(dt) - 1) * yokoya_factor * rshunt * omega**2 / qfactor *
                numpy.exp(alpha * dt) * dt)
    else:
        wake = (-(numpy.sign(dt) - 1) * yokoya_factor * rshunt * omega**2 / (qfactor *
                omegabar) * numpy.exp(alpha*dt) * numpy.sinh(omegabar*dt))
    return wake


def transverse_reswall(srange, yokoya_factor, length, rvac, conduct, beta):
    """Define the wake function (transverse) of a resistive wall with the given
    parameters according to Alex Chao's RW model (2.53) and definitions used in
    HEADTAIL
    """

    if numpy.amin(srange) <= 0:
        raise ValueError("""
                        Provided srange has either negative values or 0s
                        This is not allowed for the transverse resistive wall
                        wake function. Please correct.
                        """)

    z0 = 119.9169832 * numpy.pi
    dt = -srange/(beta * clight)
    wake = (yokoya_factor * (numpy.sign(dt) - 1) / 2. *
            beta * length / numpy.pi / rvac**3 *
            numpy.sqrt(-z0 * clight / conduct / numpy.pi / dt))

    return wake


