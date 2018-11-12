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
def ml_lattice(engine):
    lattice = engine.load(utils.LATTICE)
    return lattice['RING']


@pytest.fixture
def py_lattice():
    return Lattice(load_mat(utils.LATTICE, key='RING'), keep_all=True)
