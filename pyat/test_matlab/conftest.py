"""
conftest.py is a special pytest file that allows you to share fixtures
between other modules.
"""
import pytest
import os
import sys
from at.load import load_mat
try:
    from matlab.engine import connect_matlab, start_matlab, EngineError
except ImportError:
    print('Matlab comparison tests require Matlab Python Engine installed.')
    print('Python will exit.')
    sys.exit()

ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '../..'))
dba_ring = os.path.join(ROOT_DIR, 'pyat/test_matlab/dba.mat')
hmba_ring = os.path.join(ROOT_DIR, 'pyat/test_matlab/hmba.mat')
err_ring = os.path.join(ROOT_DIR, 'pyat/test_matlab/err.mat')


@pytest.fixture(scope='session')
def engine():
    try:
        eng = connect_matlab('pytest')
        # Keep the existing Matlab path
    except EngineError:
        eng = start_matlab()
        # Add the local AT path
        eng.addpath(eng.genpath(os.path.join(ROOT_DIR, 'atintegrators/')))
        eng.addpath(eng.genpath(os.path.join(ROOT_DIR, 'atmat/')))
    # Add the local test_matlab directory
    eng.addpath(os.path.dirname(__file__))
    yield eng
    eng.quit()


def _load_lattice(engine, name, key):
    myl = engine.load(name)
    pyl = load_mat(name, key=key, keep_all=True)
    return pyl, myl[key], None


# ----------------------- dba lattice -----------------------
@pytest.fixture(scope='session')
def dba(engine):
    return _load_lattice(engine, dba_ring, 'RING')


# ----------------------- hmba lattice -----------------------
@pytest.fixture(scope='session')
def hmba(engine):
    return _load_lattice(engine, hmba_ring, 'RING')


#                hmba lattice with cavities on
@pytest.fixture(scope='session')
def hmba_cav(engine):
    pyl, myl, _ = _load_lattice(engine, hmba_ring, 'RING')
    myl, radindex = engine.pyproxy('atradon', myl, 'BendPass', '', nargout=2)
    return pyl.radiation_on(copy=True,dipole_pass=None), myl, radindex


#               hmba lattice with radiation on
@pytest.fixture(scope='session')
def hmba_rad(engine):
    pyl, myl, _ = _load_lattice(engine, hmba_ring, 'RING')
    myl, radindex = engine.pyproxy('atradon', myl, nargout=2)
    return pyl.radiation_on(copy=True), myl, radindex


# -------------------- lattice with errors --------------------
# Takes too much time for usual tests
@pytest.fixture(scope='session')
def err(engine):
    return _load_lattice(engine, err_ring, 'RING')
