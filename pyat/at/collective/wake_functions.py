"""
Analytical wake functions
"""
import numpy
from at.constants import clight
from scipy.interpolate import interp1d


def convolve_wakefun(srange, w, sigs, gauss_sigma=10):
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

    conv = numpy.convolve(gauss, w, mode='full') * ds
    s_offset = gauss_sigma * sigs - numpy.amin(srange)
    s_conv = numpy.arange(len(conv)) * ds - s_offset

    ifun = interp1d(s_conv, conv, bounds_error=False, fill_value=0)
    conv_wf = ifun(srange)
    return conv_wf


def long_resonator_wf(srange, frequency, qfactor, rshunt, beta):
    """
    
    Define the wake function (longitudinal) of a resonator
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
    """
    Define the wake function (transverse) of a resonator
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
    """
    Define the wake function (transverse) of a resistive wall with the given
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
