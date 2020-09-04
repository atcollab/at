import at
import numpy
import pytest


@pytest.fixture
def r_in(engine):
    r_in = engine.zeros(6, 1)
    return r_in


@pytest.mark.parametrize('ml_lattice, py_lattice',
                         [(pytest.lazy_fixture('ml_dba'),
                           pytest.lazy_fixture('py_dba')),
                          (pytest.lazy_fixture('ml_hmba'),
                           pytest.lazy_fixture('py_hmba'))])
def test_one_turn_for_demo_lattice(r_in, engine, ml_lattice, py_lattice):
    for i in range(6):
        # Change each item in r_in before calling.
        r_in[i][0] = 1e-5
        # Matlab call
        r_out = engine.atpass(ml_lattice, r_in, 1, 1)

        # Python setup
        py_r_in = numpy.asfortranarray(r_in).reshape(6, 1)
        py_r_out = numpy.asfortranarray(r_out).reshape(6, 1)

        # Python call; py_r_in modified in place
        at.atpass(py_lattice, py_r_in, 1)

        numpy.testing.assert_almost_equal(py_r_in, py_r_out)
