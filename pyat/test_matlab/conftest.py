"""
conftest.py is a special pytest file that allows you to share fixtures
between other modules.
"""
import pytest
import utils
from at import Lattice
from at.load import load_mat


@pytest.fixture(scope='session')
def engine():
    eng = utils.initialise_matlab()
    yield eng
    eng.quit()


@pytest.fixture
def ml_dba(engine):
    lattice = engine.load(utils.dba_ring)
    return lattice['RING']

@pytest.fixture
def ml_hmba(engine):
    lattice = engine.load(utils.hmba_ring)
    return lattice['RING']


@pytest.fixture
def py_dba():
    return Lattice(load_mat(utils.dba_ring, key='RING'), keep_all=True)


@pytest.fixture
def py_hmba():
    return Lattice(load_mat(utils.hmba_ring, key='RING'), keep_all=True)
