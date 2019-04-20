"""
Import of AT lattice from different formats:
- .mat files
"""
from at.load.utils import *
from at.load.utils import _isparam
from at.load.matfile import load_mat


__all__ = ['load_mat']
