from at import elements, track_function
import numpy
import pytest


@pytest.mark.parametrize('input_dim', [(0,), (5,), (7,), (1, 1), (6, 1, 1)])
def test_track_function_raises_AssertionError_if_rin_incorrect_shape(input_dim):
    rin = numpy.zeros(input_dim)
    lattice = []
    with pytest.raises(AssertionError):
        track_function(lattice, rin)


def test_multiple_particles_track_function():
    lattice = [elements.Drift('Drift', 1.0)]
    rin = numpy.zeros((6, 2))
    rin[0, 0] = 1e-6  # particle one offset in x
    rin[2, 1] = 1e-6  # particle two offset in y
    r_original = numpy.copy(rin)
    r_out, *_ = track_function(lattice, rin, nturns=2)
    # particle position is not changed passing through the drift
    numpy.testing.assert_equal(r_original[:, 0], r_out[:, 0, 0, 0])
    numpy.testing.assert_equal(r_original[:, 0], r_out[:, 0, 0, 1])
    numpy.testing.assert_equal(r_original[:, 1], r_out[:, 1, 0, 0])
    numpy.testing.assert_equal(r_original[:, 1], r_out[:, 1, 0, 1])


def test_lattice_convert_to_list_if_incorrect_type():
    lattice = numpy.array([elements.Drift('Drift', 1.0)])
    rin = numpy.zeros((6, 2))
    rin[0, 0] = 1e-6
    r_original = numpy.copy(rin)
    r_out, *_ = track_function(lattice, rin, 1)
    numpy.testing.assert_equal(r_original, r_out.reshape(6, 2))
