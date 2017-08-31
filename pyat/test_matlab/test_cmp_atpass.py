import at
import numpy
import pytest


@pytest.fixture
def r_in(engine):
    r_in = engine.zeros(6, 1)
    return r_in


def test_one_turn_for_demo_lattice(r_in, engine, ml_lattice, py_lattice):
    for i in range(6):
        # Change each item in r_in before calling.
        r_in[i][0] = 1e-5
        # Matlab call
        r_out = engine.atpass(ml_lattice, r_in, 1, 1)

        # Python setup
        py_r_in = numpy.asarray(r_in).reshape(1, 6)
        py_r_out = numpy.asarray(r_out).reshape(1, 6)

        # Python call; py_r_in modified in place
        at.atpass(py_lattice, py_r_in, 1)

        numpy.testing.assert_almost_equal(py_r_in, py_r_out)
