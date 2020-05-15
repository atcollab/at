'''
Original author 
Jaime Maria Coello De Portugal - Martinez Vazquez


Written while at CERN circa 2016
This is a reduced version keeping only the key components.
Full code can be found at: https://github.com/pylhc/harpy

'''

import numpy as np

PI2I = 2 * np.pi * complex(0, 1)

ZERO_PAD_DEF = False
HANN_DEF = False


class HarmonicAnalysis(object):

    def __init__(self, samples, zero_pad=ZERO_PAD_DEF, hann=HANN_DEF):
        self._samples = samples
        self._compute_orbit()
        if zero_pad:
            self._pad_signal()
        self._length = len(self._samples)
        self._int_range = np.arange(self._length)
        self._hann_window = None
        if hann:
            self._hann_window = np.hanning(self._length)

    def laskar_method(self, num_harmonics):
        samples = self._samples[:]  # Copy the samples array.
        n = self._length
        coefficients = []
        frequencies = []
        for _ in range(num_harmonics):
            # Compute this harmonic frequency and coefficient.
            dft_data = HarmonicAnalysis._fft(samples)
            frequency = self._jacobsen(dft_data)
            coefficient = HarmonicAnalysis._compute_coef(
                samples,
                frequency * n
            ) / n

            # Store frequency and amplitude
            coefficients.append(coefficient)
            frequencies.append(frequency)

            # Subtract the found pure tune from the signal
            new_signal = coefficient * np.exp(PI2I * frequency * self._int_range)
            samples = samples - new_signal

        coefficients, frequencies = zip(*sorted(zip(coefficients, frequencies),
                                                key=lambda tuple: np.abs(tuple[0]),
                                                reverse=True))
        return frequencies, coefficients

    def get_signal(self):
        if self._hann_window is not None:
            return self._samples * self._hann_window
        else:
            return self._samples

    def get_coefficient_for_freq(self, freq):
        return self._compute_coef(self._samples, freq * self._length) / self._length

    def _pad_signal(self):
        """
        Pads the signal with zeros to a "good" FFT size.
        """
        length = len(self._samples)
        # TODO Think proper pad size
        pad_length = (1 << (length - 1).bit_length()) - length
        # pad_length = 6600 - length
        self._samples = np.pad(
            self._samples,
            (0, pad_length),
            'constant'
        )
        # self._samples = self._samples[:6000]

    def _jacobsen(self, dft_values):
        """
        This method interpolates the real frequency of the
        signal using the three highest peaks in the FFT.
        """
        k = np.argmax(np.abs(dft_values))
        n = self._length
        r = dft_values
        delta = np.tan(np.pi / n) / (np.pi / n)
        kp = (k + 1) % n
        km = (k - 1) % n
        delta = delta * np.real((r[km] - r[kp]) / (2 * r[k] - r[km] - r[kp]))
        return (k + delta) / n

    @staticmethod
    def _compute_coef_simple(samples, kprime):
        """
        Computes the coefficient of the Discrete Time Fourier
        Transform corresponding to the given frequency (kprime).
        """
        n = len(samples)
        freq = kprime / n
        exponents = np.exp(-PI2I * freq * np.arange(n))
        coef = np.sum(exponents * samples)
        return coef

    def _compute_coef_goertzel(samples, kprime):
        """
        Computes the coefficient of the Discrete Time Fourier
        Transform corresponding to the given frequency (kprime).
        This function is faster than the previous one if compiled
        with Numba.
        """
        n = len(samples)
        a = 2 * np.pi * (kprime / n)
        b = 2 * np.cos(a)
        c = np.exp(-complex(0, 1) * a)
        d = np.exp(-complex(0, 1) *
                   ((2 * np.pi * kprime) / n) *
                   (n - 1))
        s0 = 0.
        s1 = 0.
        s2 = 0.
        for i in range(n - 1):
            s0 = samples[i] + b * s1 - s2
            s2 = s1
            s1 = s0
        s0 = samples[n - 1] + b * s1 - s2
        y = s0 - s1 * c
        return y * d

    def _compute_orbit(self):
        self.closed_orbit = np.mean(self._samples)
        self.closed_orbit_rms = np.std(self._samples)
        self.peak_to_peak = np.max(self._samples) - np.min(self._samples)

    @staticmethod
    def _conditional_import_compute_coef():
        """
        Checks if Numba is installed.
        If it is, it sets the compiled goertzel algorithm as the
        coefficient function to use. If it isn't, it uses the
        normal Numpy one.
        """
        try:
            from numba import jit
            print("Using compiled Numba functions.")
            return jit(HarmonicAnalysis._compute_coef_goertzel,
                       nopython=True, nogil=True)
        except ImportError:
            print("Numba not found, using numpy functions.")
            return HarmonicAnalysis._compute_coef_simple

    @staticmethod
    def _conditional_import_fft():
        """
        If SciPy is installed, it will set its fft as the one
        to use as it is slightly faster. Otherwise it will use
        the Numpy one.
        """
        try:
            from scipy.fftpack import fft as scipy_fft
            fft = staticmethod(scipy_fft)
            print("Scipy found, using scipy FFT.")
        except ImportError:
            from numpy.fft import fft as numpy_fft
            fft = staticmethod(numpy_fft)
            print("Scipy not found, using numpy FFT.")
        return fft

# Set up conditional functions on load ##############################################
HarmonicAnalysis._compute_coef = staticmethod(HarmonicAnalysis._conditional_import_compute_coef())
HarmonicAnalysis._fft = HarmonicAnalysis._conditional_import_fft()
#####################################################################################


class HarmonicPlotter(object):
    def __init__(self, harmonic_analisys):
        import matplotlib.pyplot as plt
        plt.rcParams['backend'] = "Qt4Agg"
        self._plt = plt
        self._harmonic_analisys = harmonic_analisys

    def plot_laskar(self, num_harmonics):
        (frequencies,
         coefficients) = self._harmonic_analisys.laskar_method(num_harmonics)
        self._plt.scatter(frequencies, coefficients)
        self._show_and_reset()

    def plot_windowed_signal(self):
        signal = self._harmonic_analisys.get_signal()
        self._plt.plot(range(len(signal)), signal)
        self._show_and_reset()

    def plot_fft(self):
        fft_func = HarmonicAnalysis._fft
        signal = self._harmonic_analisys.get_signal()
        fft = fft_func(signal)
        self._plt.plot(range(len(fft)), np.abs(fft))
        self._show_and_reset()

    def _show_and_reset(self):
        self._plt.show()
        self._plt.clf()
