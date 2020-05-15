import numpy as np
import matplotlib.pyplot as plt

from harmonic_analysis import HarmonicAnalysis




def compute_tunes(centroid, num_harmonics=5):
    ''' 
    INPUT
    centroid: turn by turn amplitude data
    num_harmonics: how many lines you want to return

    OUTPUT
    ha_tune, ha_amp: tune and amplitude lines. They are sorted in descending order (wrt. amplitude).

    EXAMPLE

    amplitudes = [3., 0.5]
    tunes = [0.23, 0.16]
    turns = np.arange(100)

    centroid = 0
    for amp, tune in zip(amplitudes, tunes):
	    centroid += amp * np.sin(2. * np.pi * turns * tune)

    '''
    ha = HarmonicAnalysis(centroid)
    ha_tune, ha_amp = ha.laskar_method(num_harmonics=num_harmonics)
    ha_amp = np.abs(np.array(ha_amp))
    ha_tune = np.array(ha_tune)
    return (ha_tune, ha_amp)



