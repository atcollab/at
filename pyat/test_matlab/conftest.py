"""
conftest.py is a special pytest file that allows you to share fixtures
between other modules.
"""
import pytest
import utils
from at import load_mat, Lattice


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
    return Lattice(load_mat.load(utils.LATTICE), keep_all=True)
