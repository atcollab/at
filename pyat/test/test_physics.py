import numpy
import pytest
from at import physics
from at import load_mat
from at import atpass


LATTICE_FILE = '../atmat/atdemos/atmatchExamples/ExampleATMATCH/dba.mat'


DP = 1e-5
M44_MATLAB = numpy.array([[-0.66380202, 2.23414498, 0, 0],
                          [-0.25037182, -0.66380182, 0, 0],
                          [0, 0, -0.99921978, 0.262170798],
                          [0, 0, -0.0059496965, -0.99921979]])


@pytest.fixture
def ring():
    ring = load_mat.load(LATTICE_FILE)
    return ring


def test_find_orbit4(ring):
    orbit4 = physics.find_orbit4(ring, DP)
    expected = numpy.array([1.091636e-7, 1.276747e-15, 0, 0, DP, 0])
    numpy.testing.assert_allclose(orbit4, expected, atol=1e-12)


def test_find_orbit4_finds_zeros_if_dp_zero(ring):
    orbit4 = physics.find_orbit4(ring, 0)
    expected = numpy.zeros((6,))
    numpy.testing.assert_allclose(orbit4, expected)


def test_find_orbit4_result_unchanged_by_atpass(ring):
    orbit = physics.find_orbit4(ring, DP)
    orbit_copy = numpy.copy(orbit)
    orbit[4] = DP
    atpass(ring, orbit, 1)
    numpy.testing.assert_allclose(orbit[:4], orbit_copy[:4], atol=1e-12)


def test_find_orbit4_with_two_refpts(ring):
    _, all_points = physics.find_orbit4(ring, DP, [49, 99])
    expected = numpy.array(
        [[8.148212e-6, 1.0993354e-5, 0, 0, DP, 2.963929e-6],
         [3.0422808e-8, 9.1635269e-8, 0, 0, DP, 5.9280346e-6]]
    ).T
    numpy.testing.assert_allclose(all_points, expected, atol=1e-12)


@pytest.mark.parametrize('refpts', ([145], [20], [1, 2, 3]))
def test_find_m44_returns_same_answer_as_matlab(ring, refpts):
    m44, mstack = physics.find_m44(ring, dp=DP, refpts=refpts)

    numpy.testing.assert_allclose(m44[:4], M44_MATLAB[:4], rtol=1e-5, atol=1e-7)
    stack_size = 0 if refpts is None else len(refpts)
    assert mstack.shape == (4, 4, stack_size)


@pytest.mark.parametrize('refpts', ([145], [1, 2, 3, 145]))
def test_get_twiss(ring, refpts):
    twiss, tune, chrom = physics.get_twiss(ring, DP, refpts, get_chrom=True)
    numpy.testing.assert_allclose(twiss['s_pos'][-1], 56.209377216)
    numpy.testing.assert_allclose(twiss['closed_orbit'][0][:5],
                                  [1.0916359e-7, 0, 0, 0, DP], atol=1e-12)
    numpy.testing.assert_allclose(twiss['m44'][-1, :, :],
                                  M44_MATLAB, rtol=1e-5, atol=1e-7)
    numpy.testing.assert_almost_equal(twiss['beta'][-1, :],
                                      (2.9872, 6.6381), decimal=4)

    numpy.testing.assert_allclose(tune, (0.3655291, 0.4937126),
                                  rtol=1e-5, atol=1e-12)
    numpy.testing.assert_allclose(chrom, (-0.3090409, -0.44186077),
                                  rtol=1e-5)
