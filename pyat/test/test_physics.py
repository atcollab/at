import at
import numpy
import pytest
from at import physics, load, atpass, elements
from at.lattice import AtWarning, AtError


DP = 1e-5
M44_MATLAB = numpy.array([[-0.66380202, 2.23414498, 0, 0],
                          [-0.25037182, -0.66380182, 0, 0],
                          [0, 0, -0.99921978, 0.262170798],
                          [0, 0, -0.0059496965, -0.99921979]])


def test_find_orbit4(dba_ring):
    orbit4, _ = physics.find_orbit4(dba_ring, DP)
    expected = numpy.array([1.091636e-7, 1.276747e-15, 0, 0, DP, 0])
    numpy.testing.assert_allclose(orbit4, expected, atol=1e-12)


def test_find_orbit4_finds_zeros_if_dp_zero(dba_ring):
    orbit4, _ = physics.find_orbit4(dba_ring, 0)
    expected = numpy.zeros((6,))
    numpy.testing.assert_allclose(orbit4, expected, atol=1e-7)


def test_find_orbit4_result_unchanged_by_atpass(dba_ring):
    orbit, _ = physics.find_orbit4(dba_ring, DP)
    orbit_copy = numpy.copy(orbit)
    orbit[4] = DP
    atpass(dba_ring, orbit, 1)
    numpy.testing.assert_allclose(orbit[:4], orbit_copy[:4], atol=1e-12)


def test_find_orbit4_with_two_refpts_with_and_without_guess(dba_ring):
    expected = numpy.array(
        [[8.148212e-6, 1.0993354e-5, 0, 0, DP, 2.963929e-6],
         [3.0422808e-8, 9.1635269e-8, 0, 0, DP, 5.9280346e-6]]
    )
    _, all_points = physics.find_orbit4(dba_ring, DP, [49, 99])
    numpy.testing.assert_allclose(all_points, expected, atol=1e-12)
    _, all_points = physics.find_orbit4(dba_ring, DP, [49, 99],
                                        numpy.array([0., 0., 0., 0., DP, 0.]))
    numpy.testing.assert_allclose(all_points, expected, atol=1e-12)


def test_orbit_maxiter_warnings(hmba_ring):
    with pytest.warns(AtWarning):
        physics.find_orbit4(hmba_ring, max_iterations=1)
    with pytest.warns(AtWarning):
        physics.find_sync_orbit(hmba_ring, max_iterations=1)
    with pytest.warns(AtWarning):
        physics.find_orbit6(hmba_ring, max_iterations=1)


@pytest.mark.parametrize('refpts', ([145], [20], [1, 2, 3]))
def test_find_m44_returns_same_answer_as_matlab(dba_ring, refpts):
    m44, mstack = physics.find_m44(dba_ring, dp=DP, refpts=refpts)
    numpy.testing.assert_allclose(m44, M44_MATLAB, rtol=1e-5, atol=1e-7)
    assert mstack.shape == (len(refpts), 4, 4)
    m44, mstack = physics.find_m44(dba_ring, dp=DP, refpts=refpts, full=True)
    numpy.testing.assert_allclose(m44, M44_MATLAB, rtol=1e-5, atol=1e-7)
    assert mstack.shape == (len(refpts), 4, 4)


@pytest.mark.parametrize('refpts', ([145], [20], [1, 2, 3]))
def test_find_m66(hmba_ring, refpts):
    m66, mstack = physics.find_m66(hmba_ring, refpts=refpts)
    expected = numpy.array([[-0.735654, 4.673766, 0., 0., 2.997161e-3, 0.],
                            [-9.816788e-2, -0.735654, 0., 0., 1.695263e-4, 0.],
                            [0., 0., 0.609804, -2.096051, 0., 0.],
                            [0., 0., 0.299679, 0.609799, 0., 0.],
                            [0., 0., 0., 0., 1., 0.],
                            [1.695128e-4, 2.997255e-3, 0., 0., 2.243281e-3, 1.]])
    numpy.testing.assert_allclose(m66, expected, rtol=1e-5, atol=1e-7)
    stack_size = 0 if refpts is None else len(refpts)
    assert mstack.shape == (stack_size, 6, 6)


@pytest.mark.parametrize('index', (20, 1, 2))
def test_find_elem_m66(hmba_ring, index):
    m66 = physics.find_elem_m66(hmba_ring[index])
    if index is 20:
        expected = numpy.array([[1.0386, 0.180911, 0., 0., 0., 0.],
                                [0.434959, 1.0386, 0., 0., 0., 0.],
                                [0., 0., 0.961891, 0.176344, 0., 0.],
                                [0., 0., -0.423978, 0.961891, 0., 0.],
                                [0., 0., 0., 0., 1., 0.],
                                [0., 0., 0., 0., 0., 1.]])
    else:
        expected = numpy.eye(6)
    numpy.testing.assert_allclose(m66, expected, rtol=1e-5, atol=1e-7)


def test_find_sync_orbit(dba_ring):
    expected = numpy.array([[1.030844e-5, 1.390795e-5, -2.439041e-30,
                             4.701621e-30, 1.265181e-5, 3.749859e-6],
                            [3.86388e-8, 1.163782e-7, -9.671192e-30,
                             3.567819e-30, 1.265181e-5, 7.5e-6]])
    _, all_points = physics.find_sync_orbit(dba_ring, DP, [49, 99])
    numpy.testing.assert_allclose(all_points, expected, rtol=1e-5, atol=1e-7)


def test_find_sync_orbit_finds_zeros(dba_ring):
    sync_orbit = physics.find_sync_orbit(dba_ring)[0]
    numpy.testing.assert_equal(sync_orbit, numpy.zeros(6))


def test_find_orbit6(hmba_ring):
    expected = numpy.zeros((len(hmba_ring), 6))
    refpts = numpy.ones(len(hmba_ring), dtype=bool)
    _, all_points = physics.find_orbit6(hmba_ring, refpts)
    numpy.testing.assert_allclose(all_points, expected, atol=1e-12)

def test_find_orbit6_raises_AtError_if_there_is_no_cavity(dba_ring):
    with pytest.raises(at.lattice.utils.AtError):
        physics.find_orbit6(dba_ring)

def test_find_m44_no_refpts(dba_ring):
    m44 = physics.find_m44(dba_ring, dp=DP)[0]
    expected = numpy.array([[-0.66380, 2.23415, 0., 0.],
                            [-0.25037, -0.66380, 0.,0.],
                            [-1.45698e-31, -1.15008e-30, -0.99922, 0.26217],
                            [6.57748e-33, 8.75482e-32, -5.9497e-3, -0.99922]])
    numpy.testing.assert_allclose(m44, expected, rtol=1e-5, atol=1e-7)


@pytest.mark.parametrize('refpts', ([145], [1, 2, 3, 145]))
def test_get_twiss(dba_ring, refpts):
    twiss0, tune, chrom, twiss = physics.get_twiss(dba_ring, DP, refpts,
                                                   get_chrom=True)
    numpy.testing.assert_allclose(twiss['s_pos'][-1], 56.209377216, atol=1e-9)
    numpy.testing.assert_allclose(twiss['closed_orbit'][0][:5],
                                  [1.0916359e-7, 0, 0, 0, DP], atol=1e-12)
    numpy.testing.assert_allclose(twiss['m44'][-1, :, :], M44_MATLAB,
                                  rtol=1e-5, atol=1e-7)
    numpy.testing.assert_almost_equal(twiss['beta'][-1, :], [2.9872, 6.6381],
                                      decimal=4)
    numpy.testing.assert_allclose(tune, [0.3655291, 0.4937126], rtol=1e-5,
                                  atol=1e-7)
    numpy.testing.assert_allclose(chrom, [-0.30903657, -0.4418593], rtol=1e-5,
                                  atol=1e-7)


def test_get_twiss_no_refpts(dba_ring):
    twiss0, tune, chrom, twiss = physics.get_twiss(dba_ring, DP, get_chrom=True)
    assert list(twiss) == []
    assert len(physics.get_twiss(dba_ring, DP, get_chrom=True)) is 4


@pytest.mark.parametrize('refpts', ([145], [1, 2, 3, 145]))
def test_linopt(dba_ring, refpts):
    lindata0, tune, chrom, lindata = physics.linopt(dba_ring, DP, refpts,
                                                    get_chrom=True)
    numpy.testing.assert_allclose(tune, [0.365529, 0.493713], rtol=1e-5)
    numpy.testing.assert_allclose(chrom, [-0.309037, -0.441859], rtol=1e-5)
    numpy.testing.assert_allclose(lindata['s_pos'][-1], 56.209377216, atol=1e-9)
    numpy.testing.assert_allclose(lindata['closed_orbit'][-1][:5],
                                  [1.091636e-7, 1.276757e-15, 4.238871e-33,
                                   1.117703e-33, DP], atol=1e-12)
    numpy.testing.assert_allclose(lindata['dispersion'][-1],
                                  [1.107402e-2, 1.262031e-10, -2.139355e-25,
                                   3.757804e-25], rtol=1e-5, atol=1e-7)
    expected = [[-0.663802, 2.234145, 0, 0], [-0.250372, -0.663802, 0, 0],
                [-1.456977e-31, -1.150075e-30, -0.99922, 0.262171],
                [6.577482e-33, 8.75482e-32, -5.949696e-3, -0.99922]]
    numpy.testing.assert_allclose(lindata['m44'][-1], expected, rtol=1e-5,
                                  atol=1e-7)
    numpy.testing.assert_allclose(lindata['alpha'][-1], [-1.32787e-7,
                                                         1.85909e-7],
                                  rtol=1e-5, atol=1e-7)
    numpy.testing.assert_almost_equal(lindata['beta'][-1], [2.98719, 6.638115],
                                      decimal=4)
    numpy.testing.assert_almost_equal(lindata['mu'][-1], [2.296687, 3.102088],
                                      decimal=4)
    numpy.testing.assert_almost_equal(lindata['gamma'][-1], 1)
    numpy.testing.assert_allclose(lindata['A'][-1],
                                  [[-0.6638, 2.23415],
                                   [-0.25037, -0.6638]], rtol=1e-5, atol=1e-7)
    numpy.testing.assert_allclose(lindata['B'][-1],
                                  [[-0.99922, 0.262171],
                                   [-0.00595, -0.99922]], rtol=1e-4, atol=1e-7)
    numpy.testing.assert_allclose(lindata['C'][-1], [[-9.87933e-32, -1.65044e-30],
                                                     [-2.44501e-32, -2.91703e-31]],
                                  rtol=1e-5, atol=1e-7)


@pytest.mark.parametrize('refpts', ([145], [1, 2, 3, 145]))
def test_linopt_uncoupled(dba_ring, refpts):
    lindata0, tune, chrom, lindata = physics.linopt(dba_ring, DP, refpts,
                                                    coupled=False)
    numpy.testing.assert_allclose(tune, [0.365529, 0.493713], rtol=1e-5)
    numpy.testing.assert_allclose(lindata['s_pos'][-1], 56.209377216, atol=1e-9)
    numpy.testing.assert_allclose(lindata['closed_orbit'][-1][:5],
                                  [1.091636e-7, 1.276757e-15, 4.238871e-33,
                                   1.117703e-33, DP], atol=1e-12)
    expected_m44 = [[-0.663802, 2.234145, 0, 0], [-0.250372, -0.663802, 0, 0],
                    [-1.456977e-31, -1.150075e-30, -0.99922, 0.262171],
                    [6.577482e-33, 8.75482e-32, -5.949696e-3, -0.99922]]
    numpy.testing.assert_allclose(lindata['m44'][-1], expected_m44, rtol=1e-5,
                                  atol=1e-7)
    numpy.testing.assert_allclose(lindata['alpha'][-1], [-1.32787e-7,
                                                         1.85909e-7],
                                  rtol=1e-5, atol=1e-7)
    numpy.testing.assert_almost_equal(lindata['beta'][-1], [2.98719, 6.638115],
                                      decimal=4)
    numpy.testing.assert_almost_equal(lindata['mu'][-1], [2.296687, 3.102088],
                                      decimal=4)
    numpy.testing.assert_almost_equal(lindata['gamma'][-1], 1)
    numpy.testing.assert_allclose(lindata['A'][-1],
                                  [[-0.6638, 2.23415],
                                   [-0.25037, -0.6638]], rtol=1e-5, atol=1e-7)
    numpy.testing.assert_allclose(lindata['B'][-1],
                                  [[-0.99922, 0.262171],
                                   [-0.00595, -0.99922]], rtol=1e-4, atol=1e-7)
    numpy.testing.assert_allclose(lindata['C'][-1], [[0., 0.], [0., 0.,]],
                                  rtol=1e-5, atol=1e-7)


def test_linopt_no_refpts(dba_ring):
    lindata0, tune, chrom, lindata = physics.linopt(dba_ring, DP, get_chrom=True)
    assert list(lindata) == []
    assert len(physics.linopt(dba_ring, DP, get_chrom=True)) is 4


@pytest.mark.parametrize('refpts', ([145], [1, 2, 3, 145]))
@pytest.mark.parametrize('ring_test', (False, True))
def test_ohmi_envelope(hmba_lattice, refpts, ring_test):
    hmba_lattice.radiation_on()
    if ring_test:
        hmba_lattice = hmba_lattice[:]
    emit0, beamdata, emit = physics.ohmi_envelope(hmba_lattice, refpts)
    expected_beamdata = [([0.38156302, 0.85437641, 1.0906073e-4]),
                         ([1.0044543e-5, 6.6238162e-6, 9.6533473e-6]),
                         ([[[6.9000153, -2.6064253e-5, 1.643376e-25,
                             -5.6606657e-26, 6.7852489e-8, -2.1641716e-6],
                            [-2.6064253e-5, 0.14492721, -2.5861685e-26,
                             8.9100196e-27, -1.8308682e-6, -2.5027415e-4],
                            [1.643376e-25, -2.5861685e-26, 8.5287203e-51,
                             -2.9380814e-51, 3.2831931e-31, 4.4607812e-29],
                            [-5.6606657e-26, 8.9100196e-27, -2.9380814e-51,
                             1.0121475e-51, -1.1311438e-31, -1.5368549e-29],
                            [6.7852489e-8, -1.8308682e-6, 3.2831931e-31,
                             -1.1311438e-31, 2.3130054e-11, 3.1616964e-9],
                            [-2.1641716e-6, -2.5027415e-4, 4.4607812e-29,
                             -1.5368549e-29, 3.1616964e-9, 4.32198e-7]],
                           [[1.3961211e-38, -1.855871e-40, -5.7793903e-20,
                             6.9292292e-20, 5.2371579e-36, -7.8318756e-36],
                            [-1.855871e-40, 2.0480679e-40, 2.2829887e-20,
                             1.7097162e-21, 1.4101058e-36, 2.2039876e-37],
                            [-5.7793903e-20, 2.2829887e-20, 2.6446806,
                             2.6295361e-6, 1.3965834e-16, 4.510019e-17],
                            [6.9292292e-20, 1.7097162e-21, 2.6295361e-6,
                             0.3781175, 4.5232436e-17, -3.7359176e-17],
                            [5.2371579e-36, 1.4101058e-36, 1.3965833e-16,
                             4.5232435e-17, 1.2785886e-32, -2.0874784e-33],
                            [-7.8318756e-36, 2.2039875e-37, 4.5100189e-17,
                             -3.7359176e-17, -2.0874784e-33, 4.460312e-33]],
                           [[9.1097876e-7, -1.5294479e-10, -1.8583128e-19,
                             6.4012846e-20, 5.2753238e-4, -2.4327187e-5],
                            [-1.5294479e-10, 2.6020095e-14, 3.1563872e-23,
                             -1.0872729e-23, -8.8571243e-8, -2.9379751e-8],
                            [-1.8583128e-19, 3.1563871e-23, 3.8296329e-32,
                             -1.3191842e-32, -1.0761548e-16, -3.0696962e-17],
                            [6.4012845e-20, -1.0872729e-23, -1.3191842e-32,
                             4.5441615e-33, 3.7070041e-17, 1.0574103e-17],
                            [5.2753238e-4, -8.8571243e-8, -1.0761548e-16,
                             3.7070041e-17, 0.30548511, -1.3744688e-2],
                            [-2.4327187e-5, -2.9379751e-8, -3.0696963e-17,
                             1.0574104e-17, -1.3744688e-2, 3.2741004]]]),
                         ([1.3252781e-10, -1.0946339e-38, 2.8586822e-6])]
    numpy.testing.assert_almost_equal(beamdata[0], expected_beamdata[0])
    numpy.testing.assert_almost_equal(beamdata[1], expected_beamdata[1])
    numpy.testing.assert_almost_equal(beamdata[2], expected_beamdata[2])
    numpy.testing.assert_almost_equal(beamdata[3], expected_beamdata[3])
    numpy.testing.assert_almost_equal(emit['r66'][-1], [[6.9056838e-9, 4.0964602e-9, 3.2992917e-24,
                                                         -3.3789475e-24, 7.0173936e-8, -5.9491604e-11],
                                                        [4.0964602e-9, 2.4440257e-9, 1.9637594e-24,
                                                         -2.0111714e-24, 4.1767993e-8, -3.6520695e-11],
                                                        [3.2992917e-24, 1.9637594e-24, -3.529181e-37,
                                                         2.1693618e-37, 4.1050357e-23, 1.362656e-23],
                                                        [-3.3789475e-24, -2.0111714e-24, 2.1693618e-37,
                                                         -1.4562392e-37, -4.2041462e-23, -1.374786e-23],
                                                        [7.0173936e-8, 4.1767993e-8, 4.1050357e-23,
                                                         -4.2041462e-23, 8.7311835e-7, -7.274907e-10],
                                                        [-5.9491604e-11, -3.6520695e-11, 1.362656e-23,
                                                         -1.374786e-23, -7.274907e-10, 9.3577536e-6]])
    numpy.testing.assert_almost_equal(emit['r44'][-1], numpy.zeros((4,4)))
    numpy.testing.assert_almost_equal(emit['m66'][-1], [[-1.0819411, 3.1880957, 0.,
                                                         0., 8.2240779e-2, -1.7215898e-5],
                                                        [-0.68052274, 1.0809957, 0.,
                                                         0., 4.9013119e-2, -1.0260176e-5],
                                                        [-3.7152887e-30, 1.4294634e-29, 0.75592965,
                                                         3.8705927, -1.614862e-31, 3.3803709e-35],
                                                        [-5.926115e-30, 2.2741196e-29, -0.67927929,
                                                         -2.1552475, 5.3473212e-31, -1.1193798e-34],
                                                        [-1.2076084e-8, -1.0553822e-7, 0.,
                                                         0.,0.99999591, -2.0933444e-4],
                                                        [2.9374173e-3, 6.7356792e-2, 0.,
                                                         0., 2.835826e-4, 0.99999994]])
    numpy.testing.assert_almost_equal(emit['orbit6'][-1], [4.5069513e-7, 2.6511223e-7, -1.7884602e-31,
                                                           2.2368375e-31, 4.4139269e-6, -5.8843664e-2])
    numpy.testing.assert_almost_equal(emit['emitXY'][-1], [1.3252791e-10, 0.])
    numpy.testing.assert_almost_equal(emit['emitXYZ'][-1], [3.10938742e-10, 6.58179984e-38, 2.85839567e-6])
