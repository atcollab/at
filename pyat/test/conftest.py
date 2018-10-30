"""
A special file that contains test fixtures for the other test files to use.
"""
import numpy
import pytest
from at import load_mat


@pytest.fixture
def rin():
    rin = numpy.array(numpy.zeros((6, 1)), order='F')
    return rin

@pytest.fixture(scope='session')
def class_list():
    return list(load_mat.CLASSES)

@pytest.fixture(scope='session')
def class_mapped():
    return load_mat.CLASS_MAPPING

@pytest.fixture(scope='session')
def famname_mapped():
    return load_mat.FAMILY_MAPPING

@pytest.fixture(scope='session')
def passmethod_mapped():
    return load_mat.PASSMETHOD_MAPPING
