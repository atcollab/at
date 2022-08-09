"""Python port of the Accelerator Toolbox"""

try:
    # setuptools_scm >= 7 for python >= 3.7
    from ._version import __version__, __version_tuple__
except ImportError:
    # setuptools_scm < 7 for python 3.6
    from ._version import version as __version__
    from ._version import version_tuple as __version_tuple__
# Make all functions visible in the at namespace:
from .lattice import *
from .tracking import *
from .physics import *
from .load import *
from .matching import *
from .acceptance import *
from .collective import *
