"""
A special file that contains test fixtures for the other test files to use.
"""
import os
import numpy
import pytest
import at


@pytest.fixture
def rin():
    rin = numpy.array(numpy.zeros((6, 1)), order='F')
    return rin


@pytest.fixture
def hmba_ring():
    path = os.path.realpath(os.path.join(os.path.dirname(__file__),
                                         '../test_matlab/hmba.mat'))
    return at.load.load_mat(path)


@pytest.fixture
def hmba_lattice(hmba_ring):
    return at.lattice.Lattice(hmba_ring)
