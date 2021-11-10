"""
Analytical wake functions
"""
import numpy
from at.lattice.constants import clight


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
        wake = (yokoya_factor * rshunt * omega**2 / (qfactor *
                omegabar) * numpy.exp(alpha*dt) * numpy.sin(omegabar*dt))
    elif qfactor == 0.5:
        wake = (yokoya_factor * rshunt * omega**2 / qfactor *
                numpy.exp(alpha * dt) * dt)
    else:
        wake = (yokoya_factor * rshunt * omega**2 / (qfactor *
                omegabar) * numpy.exp(alpha*dt) * numpy.sinh(omegabar*dt))
    return wake


def transverse_reswall(srange, yokoya_factor, length, rvac, conduct, beta):
    """Define the wake function (transverse) of a resistive wall with the given
    parameters according to Alex Chao's RW model (2.53) and definitions used in
    HEADTAIL
    """
    z0 = 119.9169832 * np.pi
    dt = -srange/(beta * clight)
    wake = (yokoya_factor * (numpy.sign(dt) - 1) / 2. *
            beta * length / numpy.pi / rvac**3 *
            numpy.sqrt(-clight * z0 / conduct / numpy.pi / dt))
    return wake
