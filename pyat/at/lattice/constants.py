from math import sqrt as _sqrt, pi
# noinspection PyUnresolvedReferences
from scipy.constants import c as clight
from scipy.constants import physical_constants as _cst
from scipy.constants import e as qe

e_mass = 1.0e+06 * _cst['electron mass energy equivalent in MeV'][0]  # eV
p_mass = 1.0e+06 * _cst['proton mass energy equivalent in MeV'][0]    # eV

_e_radius = _cst['classical electron radius'][0]
_hbar_c = _cst['Planck constant over 2 pi times c in MeV fm'][0]

Cgamma = 4.0 * pi * _e_radius / 3.0 / pow(e_mass, 3)                 # m/eV^3
Cq = 55 / 32 / _sqrt(3) * _hbar_c / e_mass * 1.0e-9                  # m
