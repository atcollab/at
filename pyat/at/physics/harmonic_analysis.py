"""
Simplified version of harpy from
Jaime Coello Maria de Portugal - Martinez Vazquez
Github: https://github.com/pylhc/harpy
"""

from __future__ import annotations
import numpy
from warnings import warn
from scipy.fft import fft, fftfreq
from at.lattice import AtWarning
from at.lattice import AtError


__all__ = ['get_spectrum_harmonic', 'get_main_harmonic',
           'get_tunes_harmonic']


def _pad_signal(samples, pad_length):
    """ Pad signal and round to the next power of 2"""
    if pad_length is not None:
        length = len(samples)
        pad_length = (1 << (length + int(pad_length)
                            - 1).bit_length()) - length
        samples = numpy.pad(samples, (0, pad_length))
    return samples


def _interpolate_peak(complex_values):
    """This method interpolates the real frequency of the
       signal using the three highest peaks in the FFT.
    """
    k = numpy.argmax(numpy.abs(complex_values))
    rv = complex_values
    n = len(complex_values)
    delta = numpy.tan(numpy.pi / n) / (numpy.pi / n)
    kp = (k + 1) % n
    km = (k - 1) % n
    dk = rv[km] - rv[kp]
    sk = 2 * rv[k] - rv[km] - rv[kp]
    delta *= numpy.real(dk/sk)
    return (k + delta) / n


def _compute_coef(samples, freq):
    """
    Computes the coefficient of the Discrete Time Fourier
    Transform corresponding to the given frequency
    """
    n = len(samples)
    exponents = numpy.exp(-2j*numpy.pi * freq * numpy.arange(n))
    coef = numpy.sum(exponents * samples)
    return coef


def _interpolated_fft(samples, num_harmonics, fmin, fmax,
                      maxiter):
    """Compute the interpolated FFT of a signal"""
    rn = numpy.arange(len(samples))
    coefficients = numpy.zeros(num_harmonics, dtype=complex)
    frequencies = numpy.zeros(num_harmonics)

    nfound = 0
    niter = 0
    nmax = num_harmonics*maxiter

    while (nfound < num_harmonics) and (niter < nmax):
        fft_data = fft(samples)
        frequency = _interpolate_peak(fft_data)
        coefficient = _compute_coef(samples, frequency)
        if frequency >= fmin and frequency <= fmax:
            frequencies[nfound] = frequency
            coefficients[nfound] = coefficient
            nfound += 1
        samples = samples - coefficient * numpy.exp(2j*numpy.pi*frequency*rn)
        niter += 1

    if nfound < num_harmonics:
        msg = ('{0}/{1} harmonics found in '
               'requested range'.format(nfound, num_harmonics))
        warn(AtWarning(msg))
    if nfound == 0:
        msg = ('No harmonic found within range, '
               'consider extending it or increase maxiter')
        raise AtError(msg)

    coefficients, frequencies = zip(*sorted(zip(coefficients, frequencies),
                                    key=lambda tuple: numpy.abs(tuple[0]),
                                    reverse=True))
    return frequencies, coefficients


def get_spectrum_harmonic(cent: numpy.ndarray, method: str = 'interp_fft',
                          num_harmonics: int = 20,
                          hann: bool = False,
                          fmin: float = 0, fmax: float = 1,
                          maxiter: float = 10,
                          pad_length=None) -> tuple[numpy.ndarray,
                                                    numpy.ndarray]:
    """Frequency analysis of beam motion

    Parameters:
        cent:           Centroid motions of the particle
        method:         ``'interp_fft'`` or ``'interp_fft'``.
                        Default: ``'interp_fft'``
        num_harmonics:  Number of harmonics to search for with interp_fft
        fmin:           Lower bound for spectrum search with interp_fft
        fmax:           Upper bound for spectrum search with interp_fft
        maxiter:        Multiplies ``num_harmonics`` to define the max.
                        number of iteration for the search
        hann:           Turn on Hanning window. Default: :py:obj:`False`.
                        Ignored for interpolated FFT
        pad_length      Zero pad the input signal.
                        Rounded to the higher power of 2
                        Ignored for interpolated FFT

    Returns:
        frequency (ndarray): (num_harmonics,) array of frequencies
        amplitude (ndarray): (num_harmonics,) array of amplitudes
        phase (ndarray): (num_harmonics,) array of phases
    """
    lc = len(cent)
    # laskar kept for backward compatibility
    if method == 'interp_fft' or method == 'laskar':
        if hann:
            warn(AtWarning('Windowing not efficient for'
                           'interpolated FFT: ignored'))
        if pad_length is not None:
            warn(AtWarning('Padding not efficient for'
                           'interpolated FFT: ignored'))
        ha_tune, ha_amp = _interpolated_fft(cent, num_harmonics,
                                            fmin, fmax, maxiter)
    elif method == 'fft':
        if hann:
            cent *= numpy.hanning(lc)
        cent = _pad_signal(cent, pad_length)
        fft_data = fft(cent)
        ha_tune, ha_amp = fftfreq(len(fft_data)), fft_data
    else:
        raise ValueError('The method ' + method + ' is undefined')

    ha_phase = numpy.angle(numpy.array(ha_amp))
    ha_amp = numpy.abs(numpy.array(ha_amp)) / lc
    ha_tune = numpy.array(ha_tune)
    return ha_tune, ha_amp, ha_phase


def get_main_harmonic(cents: numpy.ndarray, method: str = 'interp_fft',
                      hann: bool = False,
                      fmin: float = 0, fmax: float = 1,
                      maxiter: float = 10,
                      pad_length=None) -> numpy.ndarray:
    """Computes tunes, amplitudes and pahses from harmonic analysis

    Parameters:
        cents:          Centroid motions of the particle
        method:         ``'interp_fft'`` or ``'fft'``.
                        Default: ``'interp_fft'``
        fmin:           Lower bound for tune search
        fmax:           Upper bound for tune search
        maxiter:        Maximum number of iterations for the search
        hann:           Turn on Hanning window. Default: :py:obj:`False`.
                        Ignored for interpolated FFT
        pad_length:     Zero pad the input signal.
                        Rounded to the higher power of 2
                        Ignored for interpolated FFT

    Returns:
        tunes (ndarray):    numpy array of length len(cents), max of the
        spectrum within [fmin fmax]
        amplitude (ndarray): (len(cents), ) array of amplitudes
                             corresponding to the tune
        phase (ndarray): (len(cents), ) array of phases
                         corresponding to the tune
    """
    def get_max_spectrum(freq, amp, phase, fmin, fmax, method):
        if method == 'interp_fft':
            return freq[0], amp[0], phase[0]
        msk = numpy.logical_and(freq >= fmin, freq <= fmax)
        amp = amp[msk]
        freq = freq[msk]
        phase = phase[msk]
        freq = freq[numpy.argmax(amp)]
        phase = phase[numpy.argmax(amp)]
        amp = numpy.amax(amp)
        return freq, amp, phase

    cents = numpy.array(cents)
    if cents.ndim > 1:
        npart = cents.shape[0]
    else:
        cents = [cents]
        npart = 1

    tunes = numpy.zeros(npart)
    amps = numpy.zeros(npart)
    phases = numpy.zeros(npart)

    for i in range(npart):
        out = get_spectrum_harmonic(cents[i], num_harmonics=1, method=method,
                                    hann=hann, pad_length=pad_length,
                                    fmin=fmin, fmax=fmax, maxiter=maxiter)
        freq, amp, phase = out
        try:
            tunes[i], amps[i], phases[i] = get_max_spectrum(freq, amp, phase,
                                                            fmin, fmax,
                                                            method)
        except ValueError:
            tunes[i] = numpy.nan
            amps[i] = numpy.nan
            phases[i] = numpy.nan
    return tunes, amps, phases


def get_tunes_harmonic(cents: numpy.ndarray, method: str = 'interp_fft',
                       hann: bool = False,
                       fmin: float = 0, fmax: float = 1,
                       maxiter: float = 10,
                       pad_length=None, **kwargs) -> numpy.ndarray:
    """Computes tunes from harmonic analysis

    Parameters:
        cents:          Centroid motions of the particle
        method:         ``'interp_fft'`` or ``'fft'``.
                        Default: ``'interp_fft'``
        fmin:           Lower bound for tune search
        fmax:           Upper bound for tune search
        maxiter:        Maximum number of iterations for the search
        hann:           Turn on Hanning window. Default: :py:obj:`False`.
                        Ignored for interpolated FFT
        pad_length:     Zero pad the input signal.
                        Rounded to the higher power of 2
                        Ignored for interpolated FFT

    Returns:
        tunes (ndarray):    numpy array of length len(cents), max of the
        spectrum within [fmin fmax]
    """
    num_harmonics = kwargs.pop('num_harmonics', 1)  # Backward compatibility
    if num_harmonics != 1:
        msg = "num_harmonics is deprecated and ignored for tune calculation"
        warn(AtWarning(msg))
    tunes, _, _ = get_main_harmonic(cents, method=method, hann=hann, fmin=fmin,
                                    fmax=fmax, pad_length=pad_length)
    return tunes
