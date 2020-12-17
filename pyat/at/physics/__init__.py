"""
Accelerator physics functions
"""

from math import sqrt, pi
from scipy.constants import physical_constants as cst

_e_radius = cst['classical electron radius'][0]
_hbar = cst['Planck constant over 2 pi times c in MeV fm'][0]

e_mass = 1.0e+06 * cst['electron mass energy equivalent in MeV'][0]  # eV
Cgamma = 4.0 * pi * _e_radius / 3.0 / pow(e_mass, 3)                 # m/eV^3
Cq = 55 / 32 / sqrt(3) * _hbar / e_mass * 1.0e-9                     # m

from .diff_constants import *
from .orbit import *
from .amat import *
from .matrix import *
from .linear import *
from .ring_parameters import *
# noinspection PyUnresolvedReferences
from .diffmatrix import find_mpole_raddiff_matrix
from .radiation import *
from .harmonic_analysis import *
from .nonlinear import *
from .fastring import *
