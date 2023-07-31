from at import elements
from at import lattice_track
from at import patpass, internal_plpass
import numpy
import pytest
# import sys


# @pytest.mark.skipif(not sys.platform.startswith("win"),
#                    reason="May hang on linux and MacOS")
@pytest.mark.parametrize('func', (lattice_track, patpass, internal_plpass))
def test_patpass_multiple_particles_and_turns(func):
    nturns = 10
    nparticles = 10
    rin = numpy.zeros((6, nparticles))
    d = elements.Drift('drift', 1.0)
    lattice = [d]
    rin[1, 0] = 1e-6
    rin[3, 0] = -2e-6
    if func == lattice_track:
        rout, *_ = func(lattice, rin, nturns, use_mp=True)
    else:
        rout = func(lattice, rin, nturns)
    # results from Matlab
    assert rout.shape == (6, nparticles, 1, nturns)
    rout_expected = numpy.array([1e-6, 1e-6, -2e-6, -2e-6, 0,
                                 2.5e-12]).reshape(6, 1)
    zeros = numpy.zeros((6,))
    # The fourth index is for nturns; the first element is after one turn.
    numpy.testing.assert_equal(rout[:, 0, :, 0], rout_expected)
    for i in range(nturns):
        # The first particle is not all zeros.
        with pytest.raises(AssertionError):
            numpy.testing.assert_equal(rout[:, 0, 0, i], zeros)
        # All other particles are all zeros.
        for j in range(1, nparticles):
            numpy.testing.assert_equal(rout[:, j, 0, i], zeros)
