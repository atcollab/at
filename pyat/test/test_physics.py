import numpy
import pytest
from at import physics, load_mat, atpass, elements


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
    orbit4, _ = physics.find_orbit4(ring, DP)
    expected = numpy.array([1.091636e-7, 1.276747e-15, 0, 0, DP, 0])
    numpy.testing.assert_allclose(orbit4, expected, atol=1e-12)


def test_find_orbit4_finds_zeros_if_dp_zero(ring):
    orbit4, _ = physics.find_orbit4(ring, 0)
    expected = numpy.zeros((6,))
    numpy.testing.assert_allclose(orbit4, expected)


def test_find_orbit4_result_unchanged_by_atpass(ring):
    orbit, _ = physics.find_orbit4(ring, DP)
    orbit_copy = numpy.copy(orbit)
    orbit[4] = DP
    atpass(ring, orbit, 1)
    numpy.testing.assert_allclose(orbit[:4], orbit_copy[:4], atol=1e-12)


def test_find_orbit4_with_two_refpts_with_and_without_guess(ring):
    expected = numpy.array(
        [[8.148212e-6, 1.0993354e-5, 0, 0, DP, 2.963929e-6],
         [3.0422808e-8, 9.1635269e-8, 0, 0, DP, 5.9280346e-6]]
    )
    _, all_points = physics.find_orbit4(ring, DP, [49, 99])
    numpy.testing.assert_allclose(all_points, expected, atol=1e-12)
    _, all_points = physics.find_orbit4(ring, DP, [49, 99],
                                        numpy.array([0., 0., 0., 0., DP, 0.]))
    numpy.testing.assert_allclose(all_points, expected, atol=1e-12)


@pytest.mark.parametrize('refpts', ([145], [20], [1, 2, 3]))
def test_find_m44_returns_same_answer_as_matlab(ring, refpts):
    m44, mstack = physics.find_m44(ring, dp=DP, refpts=refpts)

    numpy.testing.assert_allclose(m44[:4], M44_MATLAB[:4], rtol=1e-5, atol=1e-7)
    stack_size = 0 if refpts is None else len(refpts)
    assert mstack.shape == (stack_size, 4, 4)


def test_find_sync_orbit_finds_zeros(ring):
    sync_orbit = physics.find_sync_orbit(ring)[0]
    numpy.testing.assert_equal(sync_orbit, numpy.zeros(6))


def test_find_m44_no_refpts(ring):
    m44 = physics.find_m44(ring, dp=DP)[0]
    expected = numpy.array([[-0.66380, 2.23415, 0., 0.],
                            [-0.25037, -0.66380, 0.,0.],
                            [-1.45698e-31, -1.15008e-30, -0.99922, 0.26217],
                            [6.57748e-33, 8.75482e-32, -5.94970e-3, -0.99922]])
    numpy.testing.assert_allclose(m44, expected, rtol=1e-5, atol=1e-7)


@pytest.mark.parametrize('refpts', ([145], [1, 2, 3, 145]))
def test_get_twiss(ring, refpts):
    twiss0, tune, chrom, twiss = physics.get_twiss(ring, DP, refpts,
                                                   get_chrom=True)
    numpy.testing.assert_allclose(twiss['s_pos'][-1], 56.209377216)
    numpy.testing.assert_allclose(twiss['closed_orbit'][0][:5],
                                  [1.0916359e-7, 0, 0, 0, DP], atol=1e-12)
    numpy.testing.assert_allclose(twiss['m44'][-1, :, :], M44_MATLAB, rtol=1e-5,
                                  atol=1e-7)
    numpy.testing.assert_almost_equal(twiss['beta'][-1, :], [2.9872, 6.6381],
                                      decimal=4)
    numpy.testing.assert_allclose(tune, [0.3655291, 0.4937126], rtol=1e-5,
                                  atol=1e-12)
    numpy.testing.assert_allclose(chrom, [-0.30903657, -0.4418593], rtol=1e-5)


def test_get_twiss_no_refpts(ring):
    twiss0, tune, chrom, twiss = physics.get_twiss(ring, DP, get_chrom=True)
    assert list(twiss) == []
    assert len(physics.get_twiss(ring, DP, get_chrom=True)) is 4


@pytest.mark.parametrize('refpts', ([145], [1, 2, 3, 145]))
def test_linopt(ring, refpts):
    lindata0, tune, chrom, lindata = physics.linopt(ring, DP, refpts,
                                                    get_chrom=True)
    numpy.testing.assert_allclose(tune, [0.365529, 0.493713], rtol=1e-5)
    numpy.testing.assert_allclose(chrom, [-0.309037, -0.441859], rtol=1e-5)
    numpy.testing.assert_allclose(lindata['s_pos'][-1], 56.209377216)
    numpy.testing.assert_allclose(lindata['closed_orbit'][-1][:5],
                                  [1.091636e-7, 1.276757e-15, 4.238871e-33,
                                   1.117703e-33, DP], atol=1e-12)
    numpy.testing.assert_allclose(lindata['dispersion'][-1],
                                  [1.107402e-2, 1.262031e-10, -2.139355e-25,
                                   3.757804e-25], rtol=1e-5)
    expected = [[-0.663802,  2.234145,  0,  0], [-0.250372, -0.663802,  0,  0],
                [-1.456977e-31, -1.150075e-30, -0.99922,  0.262171],
                [6.577482e-33,  8.75482e-32, -5.949696e-03, -0.99922]]
    numpy.testing.assert_allclose(lindata['m44'][-1], expected, rtol=1e-5)
    numpy.testing.assert_allclose(lindata['alpha'][-1], [-1.32787e-7,
                                                         1.85909e-7], rtol=1e-5)
    numpy.testing.assert_almost_equal(lindata['beta'][-1], [2.98719, 6.638115],
                                      decimal=4)
    numpy.testing.assert_almost_equal(lindata['mu'][-1], [2.296687, 3.102088],
                                      decimal=4)
    numpy.testing.assert_almost_equal(lindata['gamma'][-1], 1)
    numpy.testing.assert_allclose(lindata['A'][-1],
                                  [[-0.6638, 2.23415],
                                   [-0.25037, -0.6638]], rtol=1e-5)
    numpy.testing.assert_allclose(lindata['B'][-1],
                                  [[-0.99922, 0.262171],
                                   [-0.00595, -0.99922]], rtol=1e-4)
    numpy.testing.assert_allclose(lindata['C'][-1],
                                  [[-9.87933e-32, -1.65044e-30],
                                   [-2.44501e-32, -2.91703e-31]], rtol=1e-5)


def test_linopt_no_refpts(ring):
    lindata0, tune, chrom, lindata = physics.linopt(ring, DP, get_chrom=True)
    assert list(lindata) == []
    assert len(physics.linopt(ring, DP, get_chrom=True)) is 4
