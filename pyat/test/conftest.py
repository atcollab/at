"""
A special file that contains test fixtures for the other test files to use.
"""
import os
import numpy
import pytest
from at import elements, load, lattice


@pytest.fixture
def rin():
    rin = numpy.array(numpy.zeros((6, 1)), order='F')
    return rin


@pytest.fixture(scope='session')
def simple_ring():
    ring = [elements.Drift('D1', 1, R1=numpy.eye(6), R2=numpy.eye(6)),
            elements.Marker('M', attr='a_value'), elements.M66('M66'),
            elements.Drift('D2', 1, T1=numpy.zeros(6), T2=numpy.zeros(6)),
            elements.Drift('D3', 1, R1=numpy.eye(6), R2=numpy.eye(6)),
            elements.Drift('D4', 1, T1=numpy.zeros(6), T2=numpy.zeros(6))]
    return ring


@pytest.fixture(scope='session')
def dba_ring():
    path = os.path.realpath(os.path.join(os.path.dirname(__file__),
                                         '../test_matlab/dba.mat'))
    return load.load_mat(path)


@pytest.fixture(scope='session')
def hmba_ring():
    path = os.path.realpath(os.path.join(os.path.dirname(__file__),
                                         '../test_matlab/hmba.mat'))
    return load.load_mat(path)


@pytest.fixture(scope='session')
def hmba_lattice(hmba_ring):
    return lattice.Lattice(hmba_ring)
