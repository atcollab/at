"""Python port of the Accelerator Toolbox"""

import sys
if sys.version_info.minor < 8:
    from importlib_metadata import version, PackageNotFoundError
else:
    from importlib.metadata import version, PackageNotFoundError
try:
    __version__ = version('accelerator-toolbox')
except PackageNotFoundError:
    __version__ = "0.0.0"
# from ._version import version as __version__
# Make all functions visible in the at namespace:
from .lattice import *
from .tracking import *
from .physics import *
from .load import *
from .matching import *
sys.modules['at.constants'] = sys.modules['at.lattice.constants']
