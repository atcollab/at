from at import track
import numpy
import pytest


@pytest.mark.parametrize('input_dim', [(0,), (5,), (7,), (1,1), (6,1,1)])
def test_lattice_pass_raises_AssertionError_if_rin_incorrect_shape(input_dim):
    r_in = numpy.zeros(input_dim)
    lattice = []
    with pytest.raises(AssertionError):
        track.lattice_pass(lattice, r_in)
