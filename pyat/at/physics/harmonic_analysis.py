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
import multiprocessing
from functools import partial


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
    coef = numpy.sum(exponents * samples)/n
    return coef


def _interpolated_fft(samples, num_harmonics, fmin, fmax,
                      maxiter):
    """Compute the interpolated FFT of a signal"""
    rn = numpy.arange(len(samples))
    coefficients = numpy.zeros(num_harmonics, dtype=complex)
    frequencies = numpy.zeros(num_harmonics)

    nfound = 0
    niter = 0
    nmax = num_harmonics * maxiter

    while (nfound < num_harmonics) and (niter < nmax):
        fft_data = fft(samples)
        frequency = _interpolate_peak(fft_data)
        coefficient = _compute_coef(samples, frequency)
        if (frequency >= fmin) and (frequency <= fmax):
            frequencies[nfound] = frequency
            coefficients[nfound] = coefficient
            nfound += 1
        samples = samples - coefficient * numpy.exp(2j*numpy.pi*frequency*rn)
        niter += 1

    if nfound == 0:
        msg = ('No harmonic found within range, '
               'consider extending it or increase maxiter')
        raise AtError(msg)
    elif nfound < num_harmonics:
        msg = ('{0}/{1} harmonics found in '
               'requested range'.format(nfound, num_harmonics))
        warn(AtWarning(msg))

    coefficients, frequencies = zip(*sorted(zip(coefficients, frequencies),
                                    key=lambda tp: numpy.abs(tp[0]),
                                    reverse=True))
    return frequencies, coefficients


def get_spectrum_harmonic(cent: numpy.ndarray, method: str = 'interp_fft',
                          num_harmonics: int = 20,
                          hann: bool = False,
                          fmin: float = 0, fmax: float = 1,
                          maxiter: float = 10,
                          pad_length: float = None,
                          remove_mean: bool = True) -> tuple[numpy.ndarray,
                                                             numpy.ndarray,
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
        remove_mean:    Remove the mean of the input signal.
                        Default: :py:obj:`True`.

    Returns:
        frequency (ndarray): (num_harmonics,) array of frequencies
        amplitude (ndarray): (num_harmonics,) array of amplitudes
        phase (ndarray): (num_harmonics,) array of phases
    """
    lc = len(cent)
    # laskar kept for backward compatibility
    if remove_mean:
        cent -= numpy.mean(cent)
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


def _get_max_spectrum(freq, amp, phase, fmin, fmax):
    msk = (freq >= fmin) & (freq <= fmax)
    amp = amp[msk]
    freq = freq[msk]
    phase = phase[msk]
    tune = freq[numpy.argmax(amp)]
    phase = phase[numpy.argmax(amp)]
    amp = numpy.amax(amp)
    return tune, amp, phase


def _get_main_single(cents, **kwargs):

    def get_hmain(cents):
        fmin = kwargs.get('fmin', 0)
        fmax = kwargs.get('fmax', 1)
        try:
            out = get_spectrum_harmonic(cents, **kwargs)
            freq, amp, phase = out
            tunes, amps, phases = _get_max_spectrum(freq, amp,
                                                    phase, fmin,
                                                    fmax)
        except AtError:
            msg = ('No harmonic found within range, '
                   'consider extending it or increase maxiter')
            warn(AtWarning(msg))
            tunes = numpy.nan
            amps = numpy.nan
            phases = numpy.nan
        except ValueError:
            msg = ('Invalid input vector provided')
            warn(AtWarning(msg))
            tunes = numpy.nan
            amps = numpy.nan
            phases = numpy.nan
        return tunes, amps, phases

    cents = numpy.atleast_2d(cents)
    results = [get_hmain(c) for c in cents]
    return numpy.transpose(results)


def _get_main_multi(cents, **kwargs):
    cents = numpy.array(cents)
    start_method = kwargs.pop('start_method', None)
    pool_size = kwargs.pop('pool_size', None)
    if cents.ndim > 1:
        npart = cents.shape[0]
    else:
        return _get_main_single(cents, **kwargs)
    if pool_size is None:
        pool_size = min(npart, multiprocessing.cpu_count())
    ctx = multiprocessing.get_context(start_method)
    fun = partial(_get_main_single, **kwargs)
    with ctx.Pool(pool_size) as pool:
        results = pool.map(fun, cents)
    results = numpy.concatenate(results, axis=1)
    return results


def get_main_harmonic(cents: numpy.ndarray, method: str = 'interp_fft',
                      hann: bool = False,
                      fmin: float = 0, fmax: float = 1,
                      num_harmonics: int = 1,
                      maxiter: float = 10,
                      pad_length=None,
                      use_mp: bool = False,
                      pool_size: int = None,
                      start_method: str = None,
                      remove_mean: bool = True) -> tuple[numpy.ndarray,
                                                         numpy.ndarray,
                                                         numpy.ndarray]:
    """Computes tunes, amplitudes and phases from harmonic analysis
    The tune is defined as the harmonic with the maximum amplitude
    within the search range.

    Parameters:
        cents:          Centroid motions of the particle
        method:         ``'interp_fft'`` or ``'fft'``.
                        Default: ``'interp_fft'``
        fmin:           Lower bound for tune search
        fmax:           Upper bound for tune search
        num_harmonics:  Number of harmonics to search for.
                        Default=1.
        maxiter:        Maximum number of iterations for the search
        hann:           Turn on Hanning window. Default: :py:obj:`False`.
                        Ignored for interpolated FFT
        pad_length:     Zero pad the input signal.
                        Rounded to the higher power of 2
                        Ignored for interpolated FFT
        use_mp (bool): Flag to activate multiprocessing (default: False)
        pool_size:     number of processes used when
          *use_mp* is :py:obj:`True`. If None, ``min(npart,nproc)``
          is used. It can be globally set using the variable
          *at.lattice.DConstant.patpass_poolsize*
        start_method:  Python multiprocessing start method. Default None
                       uses the OS default method.
        remove_mean:    Remove the mean of the input signal.
                        Default: :py:obj:`True`.

    Returns:
        tunes (ndarray):    numpy array of length len(cents), max of the
        spectrum within [fmin fmax]
        amplitude (ndarray): (len(cents), ) array of amplitudes
                             corresponding to the tune
        phase (ndarray): (len(cents), ) array of phases
                         corresponding to the tune


    .. note::

       * The tune is defined as the harmonic with the maximum amplitude within
         the search range ``(fmin, fmax)``.
       * In case a ``Nan`` is present in the input vector or the tune cannot
         be found within the range, the function returns ``NaN``.
       * For the method ``'interp_fft'``, harmonics are calculated iteratively
         starting from the maximum peak of the raw FFT. ``num_harmonics=1``
         is the default, only the first harmonic is calculated.
         However, it is possible that the maximum of the interpolated
         FFT does not correspond to the maximum of the raw FFT, in which case
         ``num_harmonics`` has to be increased to get the correct peak.
    """
    if use_mp:
        tunes, amps, phases = _get_main_multi(cents,
                                              num_harmonics=num_harmonics,
                                              method=method, hann=hann,
                                              pad_length=pad_length,
                                              fmin=fmin, fmax=fmax,
                                              maxiter=maxiter,
                                              pool_size=pool_size,
                                              start_method=start_method,
                                              remove_mean=remove_mean)
    else:
        tunes, amps, phases = _get_main_single(cents,
                                               num_harmonics=num_harmonics,
                                               method=method, hann=hann,
                                               pad_length=pad_length,
                                               fmin=fmin, fmax=fmax,
                                               maxiter=maxiter,
                                               remove_mean=remove_mean)
    return tunes, amps, phases


def get_tunes_harmonic(cents: numpy.ndarray, method: str = 'interp_fft',
                       hann: bool = False,
                       fmin: float = 0, fmax: float = 1,
                       num_harmonics: int = 1,
                       maxiter: float = 10,
                       pad_length=None,
                       use_mp: bool = False,
                       pool_size: int = None,
                       start_method: str = None,
                       remove_mean: bool = True,
                       **kwargs) -> numpy.ndarray:
    """Computes tunes from harmonic analysis.

    Parameters:
        cents:          Centroid motions of the particle
        method:         ``'interp_fft'`` or ``'fft'``.
                        Default: ``'interp_fft'``
        fmin:           Lower bound for tune search
        fmax:           Upper bound for tune search
        num_harmonics:  Number of harmonics to search for.
                        Default=1.
        maxiter:        Maximum number of iterations for the search
        hann:           Turn on Hanning window. Default: :py:obj:`False`.
                        Ignored for interpolated FFT.
        pad_length:     Zero pad the input signal.
                        Rounded to the higher power of 2
                        Ignored for interpolated FFT
        use_mp (bool): Flag to activate multiprocessing (default: False)
        pool_size:     number of processes used when
          *use_mp* is :py:obj:`True`. If None, ``min(npart,nproc)``
          is used. It can be globally set using the variable
          *at.lattice.DConstant.patpass_poolsize*
        start_method:  Python multiprocessing start method. Default None
                       uses the OS default method.
        remove_mean:    Remove the mean of the input signal.
                        Default: :py:obj:`True`.

    Returns:
        tunes (ndarray):    numpy array of length len(cents), max of the
        spectrum within [fmin fmax]

    .. note::

       * The tune is defined as the harmonic with the maximum amplitude within
         the search range ``(fmin, fmax)``.
       * In case a ``Nan`` is present in the input vector or the tune cannot
         be found within the range, the function returns ``NaN``.
       * For the method ``'interp_fft'``, harmonics are calculated iteratively
         starting from the maximum peak of the raw FFT. ``num_harmonics=1``
         is the default, only the first harmonic is calculated.
         However, it is possible that the maximum of the interpolated
         FFT does not correspond to the maximum of the raw FFT, in which case
         ``num_harmonics`` has to be increased to get the correct peak.
    """
    tunes, _, _ = get_main_harmonic(cents, method=method, hann=hann, fmin=fmin,
                                    fmax=fmax, pad_length=pad_length,
                                    use_mp=use_mp, pool_size=pool_size,
                                    start_method=start_method, maxiter=maxiter,
                                    remove_mean=remove_mean,
                                    num_harmonics=num_harmonics)
    return tunes
