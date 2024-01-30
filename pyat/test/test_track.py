import numpy
import pytest

from at import elements
from at import lattice_pass, internal_lpass
from at import lattice_track


@pytest.mark.parametrize("func", (lattice_track, lattice_pass, internal_lpass))
@pytest.mark.parametrize("input_dim", [(0,), (5,), (7,), (1, 1), (6, 1, 1)])
def test_lattice_track_raises_AssertionError_if_rin_incorrect_shape(input_dim, func):
    rin = numpy.zeros(input_dim)
    lattice = []
    with pytest.raises(AssertionError):
        func(lattice, rin)


@pytest.mark.parametrize("func", (lattice_track, lattice_pass, internal_lpass))
def test_multiple_particles_lattice_track(func):
    lattice = [elements.Drift("Drift", 1.0)]
    rin = numpy.zeros((6, 2))
    rin[0, 0] = 1e-6  # particle one offset in x
    rin[2, 1] = 1e-6  # particle two offset in y
    r_original = numpy.copy(rin)
    r_out = func(lattice, rin, nturns=2)
    if isinstance(r_out, tuple):
        r_out, *_ = r_out
    # particle position is not changed passing through the drift
    numpy.testing.assert_equal(r_original[:, 0], r_out[:, 0, 0, 0])
    numpy.testing.assert_equal(r_original[:, 0], r_out[:, 0, 0, 1])
    numpy.testing.assert_equal(r_original[:, 1], r_out[:, 1, 0, 0])
    numpy.testing.assert_equal(r_original[:, 1], r_out[:, 1, 0, 1])


@pytest.mark.parametrize("func", (lattice_track, lattice_pass))
def test_lattice_convert_to_list_if_incorrect_type(func):
    lattice = numpy.array([elements.Drift("Drift", 1.0)])
    rin = numpy.zeros((6, 2))
    rin[0, 0] = 1e-6
    r_original = numpy.copy(rin)
    r_out = func(lattice, rin, 1)
    if isinstance(r_out, tuple):
        r_out, *_ = r_out
    numpy.testing.assert_equal(r_original, r_out.reshape(6, 2))


def test_in_place_flag_tracking():
    lattice = numpy.array([elements.Drift("Drift", 1.0)])
    rin = numpy.zeros((6, 1))
    rin[1] = 1e-6
    r_original = numpy.copy(rin)
    lattice_track(lattice, rin)
    numpy.testing.assert_equal(r_original, rin.reshape(6, 1))
    rout, *_ = lattice_track(lattice, rin, in_place=True)
    numpy.testing.assert_equal(rin, rout.reshape(6, 1))
