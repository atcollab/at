"""
A special file that contains test fixtures for the other test files to use.
"""
import os
import platform
import numpy
import pytest
from at import elements, load, lattice

def pytest_report_header(config):
    try:
        a = platform.uname()
        sysinfo = "system: {} {}".format(a[0], a[2])
    except AttributeError:
        sysinfo = ''
    numpyinfo = "numpy version: {}".format(numpy.__version__)
    return [sysinfo, numpyinfo]

@pytest.fixture
def rin():
    rin = numpy.array(numpy.zeros((6, 1)), order='F')
    return rin


@pytest.fixture(scope='session')
def simple_ring():
    ring = [elements.Drift('D1', 1, R1=numpy.eye(6), R2=numpy.eye(6)),
            elements.Marker('M1', attr='a_value'), elements.M66('M66'),
            elements.Drift('D2', 1, T1=numpy.zeros(6), T2=numpy.zeros(6)),
            elements.Drift('D3', 1, R1=numpy.eye(6), R2=numpy.eye(6)),
            elements.Drift('D4', 1, T1=numpy.zeros(6), T2=numpy.zeros(6))]
    return ring


@pytest.fixture(scope='session')
def simple_lattice(simple_ring):
    return lattice.Lattice(simple_ring, name='lat', energy=5, periodicity=1)


@pytest.fixture(scope='session')
def dba_lattice():
    path = os.path.realpath(os.path.join(os.path.dirname(__file__),
                                         '../test_matlab/dba.mat'))
    return load.load_mat(path)


@pytest.fixture(scope='session')
def hmba_lattice():
    path = os.path.realpath(os.path.join(os.path.dirname(__file__),
                                         '../test_matlab/hmba.mat'))
    return load.load_mat(path)
