"""
conftest.py is a special pytest file that allows you to share fixtures
between other modules.
"""

import os
import sys

import pytest

if sys.version_info.minor < 9:
    from importlib_resources import files, as_file
else:
    from importlib.resources import files, as_file
import platform
import numpy
import machine_data
import at

try:
    from matlab.engine import connect_matlab, start_matlab, EngineError
except ImportError:
    print("Matlab comparison tests require Matlab Python Engine installed.")
    print("Python will exit.")
    sys.exit()

ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), "../.."))


# noinspection PyUnusedLocal
def pytest_report_header(config):
    try:
        a = platform.uname()
        sysinfo = "system: {} {}".format(a[0], a[2])
    except AttributeError:
        sysinfo = ""
    numpyinfo = "numpy version: {}".format(numpy.__version__)
    return [sysinfo, numpyinfo]


@pytest.fixture(scope="session")
def engine():
    try:
        eng = connect_matlab("pytest")
        # Keep the existing Matlab path
    except EngineError:
        eng = start_matlab()
        # Add the local AT path
        eng.addpath(eng.genpath(os.path.join(ROOT_DIR, "atintegrators/")))
        eng.addpath(eng.genpath(os.path.join(ROOT_DIR, "atmat/")))
    # Add the local test_matlab directory
    eng.addpath(os.path.dirname(__file__))
    yield eng
    eng.quit()


def _load_lattice(engine, name, key):
    with as_file(files(machine_data) / name) as path:
        myl = engine.load(str(path))
        pyl = at.load_lattice(path, key=key, keep_all=True)
    return pyl, myl[key], None


# ----------------------- dba lattice -----------------------
@pytest.fixture(scope="session")
def dba(engine):
    return _load_lattice(engine, "dba.mat", "RING")


# ----------------------- hmba lattice -----------------------
@pytest.fixture(scope="session")
def hmba(engine):
    return _load_lattice(engine, "hmba.mat", "RING")


#                hmba lattice with cavities on
@pytest.fixture(scope="session")
def hmba_cav(engine):
    pyl, myl, _ = _load_lattice(engine, "hmba.mat", "RING")
    myl, radindex = engine.pyproxy("atradon", myl, "BendPass", "", nargout=2)
    return pyl.radiation_on(copy=True, dipole_pass=None), myl, radindex


#               hmba lattice with radiation on
@pytest.fixture(scope="session")
def hmba_rad(engine):
    pyl, myl, _ = _load_lattice(engine, "hmba.mat", "RING")
    myl, radindex = engine.pyproxy("atradon", myl, nargout=2)
    return pyl.radiation_on(copy=True), myl, radindex


# -------------------- lattice with errors --------------------
# Takes too much time for usual tests
@pytest.fixture(scope="session")
def err(engine):
    return _load_lattice(engine, "err.mat", "RING")
