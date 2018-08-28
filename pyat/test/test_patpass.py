from at import elements, atpass, patpass
import numpy
import pytest


def test_patpass_raises_ValueError_if_reuse_False(rin):
    with pytest.raises(ValueError):
        patpass([], rin, 1, reuse=False)


def test_patpass_multiple_particles_and_turns():
    nturns = 10
    nparticles = 10
    rin = numpy.zeros((6, nparticles))
    d = elements.Drift('drift', 1.0)
    assert d.Length == 1
    lattice = [d]
    rin[1, 0] = 1e-6
    rin[3, 0] = -2e-6
    rout = patpass(lattice, rin, nturns)
    # results from Matlab
    assert rout.shape == (6, nparticles*nturns)
    rout_expected = numpy.array([1e-6, 1e-6, -2e-6, -2e-6, 0, 2.5e-12]).reshape(6,)
    zeros = numpy.zeros((6,))
    numpy.testing.assert_equal(rout[:, 0], rout_expected)
    # only the first particle is not all zeros
    for i in range(nturns):
        with pytest.raises(AssertionError):
            numpy.testing.assert_equal(rout[:, i * nturns], zeros)
        for j in range(1, nparticles):
            numpy.testing.assert_equal(rout[:, i * nturns + j], zeros)
