from setuptools import build_meta as _orig
import os

# Mandatory hooks
build_wheel = _orig.build_wheel
build_sdist = _orig.build_sdist

_mpi = os.environ.get('MPI', None)
if _mpi is not None:
    _addlist = ["mpi4py"]
else:
    _addlist = []


# Optional hook
def get_requires_for_build_wheel(config_settings=None):
    print("config settings:", config_settings)
    print("Additional modules:", _addlist)
    return _orig.get_requires_for_build_wheel(config_settings) + _addlist
