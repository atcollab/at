"""
A special file that contains test fixtures for the other test files to use.
"""
import sys
if sys.version_info.minor < 9:
    from importlib_resources import files, as_file
else:
    from importlib.resources import files, as_file
import platform
import numpy
import pytest
import at
import machine_data
from at import elements


# noinspection PyUnusedLocal
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
    return at.Lattice(simple_ring, name='lat', energy=5, periodicity=1)


@pytest.fixture(scope='session')
def dba_lattice():
    with as_file(files(machine_data) / 'dba.mat') as path:
        ring = at.load_lattice(path)
    return ring


@pytest.fixture(scope='session')
def hmba_lattice():
    with as_file(files(machine_data) / 'hmba.mat') as path:
        ring = at.load_lattice(path)
    return ring


@pytest.fixture(scope='session')
def noringparam_lattice():
    with as_file(files(machine_data) / 'noringparam.mat') as path:
        ring = at.load_lattice(path)
    return ring
