import at
import numpy
from numpy.testing import assert_allclose as assert_close
from numpy.testing import assert_equal
import pytest
from at import AtWarning, physics
from at import lattice_track
from at import lattice_pass, internal_lpass


DP = 1e-5
DP2 = 0.005

M44_MATLAB = numpy.array([[-0.66380202, 2.23414498, 0, 0],
                          [-0.25037182, -0.66380182, 0, 0],
                          [0, 0, -0.99921978, 0.262170798],
                          [0, 0, -0.0059496965, -0.99921979]])

# New values after update of the C constants
orbit6_MATLAB = numpy.array(
    [-2.635267925399e-09, -1.458487743354e-10, 0.000000000000e+00,
     0.000000000000e+00, -6.696944241670e-06, -5.884521162896e-02]
)

# New values after update of the C constants
M66_MATLAB = numpy.array([
    [-0.735631089984354,   4.673714021050702,   0.0,
      0.0,                 0.002998514873422,  -0.000000627691677],
    [-0.098166216409437,  -0.735666318495682,   0.0,
      0.0,                 0.000169015216087,  -0.000000035380663],
    [ 0.0,                 0.0,                 0.609804485536081,
     -2.096029241592979,   0.0,                 0.0],
    [ 0.0,                 0.0,                 0.299675179703767,
      0.609800209788848,   0.0,                 0.0],
    [ 0.000001246879424,   0.000021544877348,   0.0,
      0.0,                 0.999980690841359,  -0.000209330146676],
    [ 0.000170098657382,   0.002995796897591,   0.0,
      0.0,                 0.002243258562703,   0.999999530409042]
])


def test_find_orbit4(dba_lattice):
    orbit4, _ = physics.find_orbit4(dba_lattice, DP)
    expected = numpy.array([1.091636e-7, 1.276747e-15, 0, 0, DP, 0])
    assert_close(orbit4, expected, atol=1e-12)


def test_find_orbit4_finds_zeros_if_dp_zero(dba_lattice):
    orbit4, _ = physics.find_orbit4(dba_lattice, 0)
    expected = numpy.zeros((6,))
    assert_close(orbit4, expected, atol=1e-7)


@pytest.mark.parametrize('func', (lattice_track, lattice_pass, internal_lpass))
def test_find_orbit4_result_unchanged_by_atpass(dba_lattice, func):
    orbit, _ = physics.find_orbit4(dba_lattice, DP)
    orbit_copy = numpy.copy(orbit)
    orbit[4] = DP
    func(dba_lattice, orbit, 1)
    assert_close(orbit[:4], orbit_copy[:4], atol=1e-12)


def test_find_orbit4_produces_same_result_with_keep_lattice_True(dba_lattice):
    orbit0, _ = physics.find_orbit4(dba_lattice)
    orbit1, _ = physics.find_orbit4(dba_lattice, keep_lattice=True)
    assert_close(orbit0, orbit1, rtol=0, atol=1e-12)


def test_find_orbit4_with_two_refpts_with_and_without_guess(dba_lattice):
    expected = numpy.array(
        [[8.148212e-6, 1.0993354e-5, 0, 0, DP, 2.963929e-6],
         [3.0422808e-8, 9.1635269e-8, 0, 0, DP, 5.9280346e-6]]
    )
    _, all_points = physics.find_orbit4(dba_lattice, DP, [49, 99])
    assert_close(all_points, expected, atol=1e-12)
    _, all_points = physics.find_orbit4(dba_lattice, DP, [49, 99],
                                        guess=numpy.array([0., 0.,
                                                           0., 0.,
                                                           DP, 0.]))
    assert_close(all_points, expected, atol=1e-12)


def test_orbit_maxiter_warnings(hmba_lattice):
    hmba_lattice_rad = hmba_lattice.radiation_on(copy=True)
    with pytest.warns(AtWarning):
        physics.find_orbit4(hmba_lattice, max_iterations=1)
    with pytest.warns(AtWarning):
        physics.find_sync_orbit(hmba_lattice, max_iterations=1)
    with pytest.warns(AtWarning):
        physics.find_orbit6(hmba_lattice_rad, max_iterations=1)


@pytest.mark.parametrize('refpts', ([145], [20], [1, 2, 3]))
def test_find_m44_returns_same_answer_as_matlab(dba_lattice, refpts):
    m44, mstack = physics.find_m44(dba_lattice, dp=DP, refpts=refpts)
    assert_close(m44, M44_MATLAB, rtol=1e-5, atol=1e-7)
    assert mstack.shape == (len(refpts), 4, 4)
    m44, mstack = physics.find_m44(dba_lattice, dp=DP, refpts=refpts,
                                   full=True)
    assert_close(m44, M44_MATLAB, rtol=1e-5, atol=1e-7)
    assert mstack.shape == (len(refpts), 4, 4)


@pytest.mark.parametrize('refpts', ([145], [20], [1, 2, 3]))
def test_find_m66(hmba_lattice, refpts):
    hmba_lattice = hmba_lattice.radiation_on(copy=True)
    m66, mstack = physics.find_m66(hmba_lattice, refpts=refpts)
    assert_close(m66, M66_MATLAB, rtol=0, atol=1e-8)
    stack_size = 0 if refpts is None else len(refpts)
    assert mstack.shape == (stack_size, 6, 6)


@pytest.mark.parametrize('index', (19, 0, 1))
def test_find_elem_m66(hmba_lattice, index):
    m66 = physics.find_elem_m66(hmba_lattice[index])
    if index == 19:
        expected = numpy.array([[1.0386, 0.180911, 0., 0., 0., 0.],
                                [0.434959, 1.0386, 0., 0., 0., 0.],
                                [0., 0., 0.961891, 0.176344, 0., 0.],
                                [0., 0., -0.423978, 0.961891, 0., 0.],
                                [0., 0., 0., 0., 1., 0.],
                                [0., 0., 0., 0., 0., 1.]])
    else:
        expected = numpy.eye(6)
    assert_close(m66, expected, rtol=1e-5, atol=1e-7)


def test_find_sync_orbit(dba_lattice):
    expected = numpy.array([[1.030844e-5, 1.390795e-5, -2.439041e-30,
                             4.701621e-30, 1.265181e-5, 3.749859e-6],
                            [3.86388e-8, 1.163782e-7, -9.671192e-30,
                             3.567819e-30, 1.265181e-5, 7.5e-6]])
    _, all_points = physics.find_sync_orbit(dba_lattice, DP, [49, 99])
    assert_close(all_points, expected, rtol=1e-5, atol=1e-7)


def test_find_sync_orbit_finds_zeros(dba_lattice):
    sync_orbit = physics.find_sync_orbit(dba_lattice)[0]
    numpy.testing.assert_equal(sync_orbit, numpy.zeros(6))


def test_find_orbit6(hmba_lattice):
    hmba_lattice = hmba_lattice.radiation_on(copy=True)
    refpts = numpy.ones(len(hmba_lattice), dtype=bool)
    orbit6, all_points = physics.find_orbit6(hmba_lattice, refpts)
    assert_close(orbit6, orbit6_MATLAB, rtol=0, atol=1e-12)


def test_find_orbit6_produces_same_result_with_keep_lattice_True(hmba_lattice):
    hmba_lattice = hmba_lattice.radiation_on(quadrupole_pass=None, copy=True)
    orbit0, _ = physics.find_orbit6(hmba_lattice)
    # Technicality - the default arguments to find_orbit6 mean that
    # keep_lattice argument is always false.
    orbit1, _ = physics.find_orbit6(hmba_lattice, keep_lattice=True)
    # With INTEGRAL keep_lattice does take effect.
    orbit2, _ = physics.find_orbit6(
        hmba_lattice, keep_lattice=True, method=physics.ELossMethod.INTEGRAL
    )
    assert_close(orbit0, orbit1, rtol=0, atol=1e-12)
    assert_close(orbit0, orbit2, rtol=0, atol=1e-12)


def test_find_orbit6_raises_AtError_if_there_is_no_cavity(dba_lattice):
    with pytest.raises(at.lattice.utils.AtError):
        physics.find_orbit6(dba_lattice)


def test_find_m44_no_refpts(dba_lattice):
    m44 = physics.find_m44(dba_lattice, dp=DP)[0]
    expected = numpy.array(
        [[-0.66380, 2.23415, 0., 0.],
         [-0.25037, -0.66380, 0., 0.],
         [-1.45698e-31, -1.15008e-30, -0.99922, 0.26217],
         [6.57748e-33, 8.75482e-32, -5.9497e-3, -0.99922]])
    assert_close(m44, expected, rtol=1e-5, atol=1e-7)


@pytest.mark.parametrize('refpts', ([145], [1, 2, 3, 145]))
def test_linopt(dba_lattice, refpts):
    """Compare with Matlab results"""
    lindata0, tune, chrom, lindata = physics.linopt(dba_lattice, DP2, refpts,
                                                    get_chrom=True)
    obs = lindata[-1]
    assert_close(tune, [0.355804633927360, 0.488487169156598],
                 rtol=1e-8)
    assert_close(chrom, [-3.428312247995742, -1.597924047969101],
                 rtol=2e-4)
    assert_close(obs['s_pos'], 56.209377216,  atol=1e-9)
    assert_close(obs['closed_orbit'][:5],
                 [0.000426438389644, -0.000000000287482, 0, 0, DP2],
                 atol=1e-12)
    assert_close(obs['dispersion'],
                 [0.156822576442091, -0.000000162902610, 0, 0],
                 rtol=1e-7, atol=2e-10)
    expected = [[-0.616893565445970, 2.191800192047084, 0, 0],
                [-0.282617257709865, -0.616893094967279, 0, 0],
                [0, 0, -0.997384485665288, 0.466909288129080],
                [0, 0, -0.011187515624062, -0.997384713772755]]
    assert_close(obs['m44'], expected, rtol=1e-7, atol=1e-12)
    assert_close(obs['alpha'], [-2.988886505944512e-07, 1.578070569086581e-06],
                 rtol=1e-8, atol=1e-8)
    assert_close(obs['beta'],
                 [2.784841119739221, 6.460251554763623], rtol=5e-8, atol=1e-12)
    assert_close(obs['mu'],
                 [2.235586452367286, 3.069255403991328], rtol=1e-8, atol=1e-12)
    assert_close(obs['gamma'], 1.0, rtol=1e-6, atol=1e-20)
    assert_close(obs['A'],
                 [[-0.616892545020345, 2.191796566512299],
                  [-0.282616790222590, -0.616892074542433]],
                 rtol=1e-7, atol=1e-12)
    assert_close(obs['B'],
                 [[-0.997384485665288, 0.466909288129080],
                  [-0.011187515624062, -0.997384713772754]],
                 rtol=1e-7, atol=1e-12)
    assert_close(obs['C'], [[0, 0], [0, 0]], rtol=1e-5, atol=1e-7)


@pytest.mark.parametrize('refpts', ([145], [1, 2, 3, 145]))
def test_linopt_uncoupled(dba_lattice, refpts):
    """Compare with Matlab results"""
    lindata0, tune, chrom, lindata = physics.linopt(dba_lattice, DP2, refpts,
                                                    coupled=False)
    obs = lindata[-1]
    assert_close(tune, [0.355804634603528, 0.488487169156732], rtol=1e-8)
    assert_close(obs['s_pos'], 56.209377216,  atol=1e-9)
    assert_close(obs['closed_orbit'][:5],
                 [0.000426438389644, -0.000000000287482, 0, 0, DP2],
                 atol=1e-12)
    expected = [[-0.616893565445970, 2.191800192047084, 0, 0],
                [-0.282617257709865, -0.616893094967279, 0, 0],
                [0, 0, -0.997384485665288, 0.466909288129080],
                [0, 0, -0.011187515624062, -0.997384713772755]]
    assert_close(obs['m44'], expected, rtol=1e-7, atol=1e-12)
    assert_close(obs['alpha'],
                 [-2.988885294797512e-07, 1.578069929495847e-06],
                 rtol=1e-8, atol=1e-8)
    assert_close(obs['beta'],
                 [2.784839991078270, 6.460248936505217], rtol=5e-8, atol=1e-12)
    assert_close(obs['mu'],
                 [2.235586452367286, 3.069255403991328], rtol=1e-8, atol=1e-12)


def test_linopt_no_refpts(dba_lattice):
    lindata0, tune, chrom, lindata = physics.linopt(dba_lattice, DP,
                                                    get_chrom=True)
    assert list(lindata) == []
    assert len(physics.linopt(dba_lattice, DP, get_chrom=True)) == 4


@pytest.mark.parametrize('refpts', ([121], [1, 2, 3, 121]))
def test_linopt_line(hmba_lattice, refpts):
#    refpts.append(len(hmba_lattice))
    l0, q, qp, ld = at.linopt(hmba_lattice, refpts=refpts)
    lt0, qt, qpt, ltd = at.linopt(hmba_lattice, refpts=refpts, twiss_in=l0)
    assert_close(ld['beta'], ltd['beta'], rtol=1e-12)
    assert_close(ld['s_pos'], ltd['s_pos'], rtol=1e-12)
    assert_close(ld['closed_orbit'], ltd['closed_orbit'], rtol=1e-12)
    assert_close(ld['alpha'], ltd['alpha'], rtol=1e-12)
    assert_close(ld['dispersion'], ltd['dispersion'], rtol=1e-7, atol=1e-12)


def test_get_tune_chrom(hmba_lattice):
    qlin = hmba_lattice.get_tune()
    qplin = hmba_lattice.get_chrom()
    qharm = hmba_lattice.get_tune(method='interp_fft')
    qpharm = hmba_lattice.get_chrom(method='interp_fft')
    assert_close(qlin, [0.38156245, 0.85437541], rtol=1e-8)
    assert_close(qharm, [0.38156245, 0.85437541], rtol=1e-8)
    assert_close(qplin, [0.1791909, 0.12242558], rtol=1e-5)
    assert_close(qpharm, [0.17919145, 0.12242622], rtol=1e-5)


def test_nl_detuning_chromaticity(hmba_lattice):
    nlqplin, _, _ = at.nonlinear.chromaticity(hmba_lattice, npoints=11)
    nlqpharm, _, _ = at.nonlinear.chromaticity(hmba_lattice,
                                               method='interp_fft', npoints=11)
    q0, q1, _, _, _, _ = at.nonlinear.detuning(hmba_lattice,
                                               npoints=11)
    assert_close(nlqplin, [[0.38156741, 0.17908231, 1.18656034, -16.47368694],
                           [0.85437409, 0.1224062, 2.01744075, -3.06407764]],
                 atol=1e-12, rtol=1e-5)
    assert_close(nlqpharm, [[0.38156741, 0.17908228, 1.18656178, -16.47370342],
                            [0.85437409, 0.12240619, 2.01744051, -3.06407046]],
                 atol=1e-12, rtol=1e-5)
    assert_close(q0, [[0.38156263, 0.85437553], [0.38156263, 0.85437553]],
                 atol=1e-12, rtol=1e-5)
    assert_close(q1, [[3005.74776344, -3256.81838517],
                      [-3258.24669916,  1615.13729938]],
                 atol=1e-12, rtol=1e-5)


def test_quantdiff(hmba_lattice):
    hmba_lattice = hmba_lattice.radiation_on(copy=True)
    dmat = physics.radiation.quantdiffmat(hmba_lattice)
    lmat = physics.radiation._lmat(dmat)
    assert_close(lmat,
                 [[1.45502934e-07, 0.00000000e+00, 0.00000000e+00,
                   0.00000000e+00, 0.00000000e+00, 0.00000000e+00],
                  [2.40963396e-09, 1.79735260e-08, 0.00000000e+00,
                   0.00000000e+00, 0.00000000e+00, 0.00000000e+00],
                  [0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
                   0.00000000e+00, 0.00000000e+00, 0.00000000e+00],
                  [0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
                   0.00000000e+00, 0.00000000e+00, 0.00000000e+00],
                  [3.72874832e-07, 2.37718999e-07, 0.00000000e+00,
                   0.00000000e+00, 5.78954180e-06, 0.00000000e+00],
                  [-1.72955964e-09, -5.42857509e-11, 0.00000000e+00,
                   0.00000000e+00, 6.52385922e-09, 3.25943528e-09]],
                 rtol=1e-4, atol=1e-20)
    assert_close(dmat,
                 [[2.11711037e-14, 3.50608810e-16, 0.00000000e+00,
                   0.00000000e+00, 5.42543819e-14, -2.51656002e-16],
                  [3.50608810e-16, 3.28853971e-16, -0.00000000e+00,
                   0.00000000e+00, 5.17114045e-15, -5.14331200e-18],
                  [0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
                   0.00000000e+00, 0.00000000e+00, 0.00000000e+00],
                  [0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
                   0.00000000e+00, 0.00000000e+00, 0.00000000e+00],
                  [5.42543819e-14, 5.17114045e-15, 0.00000000e+00,
                   0.00000000e+00, 3.37143402e-11, 3.71123417e-14],
                  [-2.51656002e-16, -5.14331200e-18, 0.00000000e+00,
                   0.00000000e+00, 3.71123417e-14, 5.61789810e-17]],
                 rtol=1e-5, atol=1e-20)


def test_simple_ring():
    ring = physics.simple_ring(6e9, 844, 992, 0.1, 0.2, 6e6, 8.5e-5)
    assert_equal(len(ring), 5)
    assert_equal(ring[-1].PassMethod, 'SimpleQuantDiffPass')
    ring.disable_6d()
    assert_equal(ring[-1].PassMethod, 'IdentityPass')
    assert_close(ring.get_tune(), [0.1, 0.2], atol=1e-10)
    

@pytest.mark.parametrize('refpts', ([121], [0, 40, 121]))
def test_ohmi_envelope(hmba_lattice, refpts):
    hmba_lattice = hmba_lattice.radiation_on(copy=True)
    emit0, beamdata, emit = hmba_lattice.ohmi_envelope(refpts)
    obs = emit[-1]

    # All expected values are Matlab results

    assert_close(beamdata['tunes'],
                 [0.381563018594173, 0.854376396669209, 1.090604388784677e-04],
                 rtol=2e-6)

    assert_close(beamdata['damping_rates'],
                 [1.00820161236330e-05, 6.5784570179901e-06, 9.6533371920752e-06],
                 rtol=2e-4)

    assert_close(
        beamdata['mode_matrices'], [
            [[ 6.900015076697874, -0.000026006243611,  0.000000000000000,  0.000000000000000,  0.000000070332460, -0.000001958895231],
             [-0.000026003127402,  0.144927219661856,  0.000000000000000,  0.000000000000000, -0.000000026206635, -0.000250268964394],
             [ 0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000],
             [ 0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000],
             [-0.000000027098267, -0.000001830741714,  0.000000000000000,  0.000000000000000,  0.000000000000331,  0.000000003161442],
             [-0.000002243290569, -0.000250274565508,  0.000000000000000,  0.000000000000000,  0.000000000045233,  0.000000432189707]],
            [[ 0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000],
             [ 0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000],
             [ 0.000000000000000,  0.000000000000000,  2.644680916022736,  0.000002697478376,  0.000000000000000,  0.000000000000000],
             [ 0.000000000000000,  0.000000000000000,  0.000002697478376,  0.378117448478870,  0.000000000000000,  0.000000000000000],
             [ 0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000],
             [ 0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000]],
            [[ 0.000000911299729, -0.000000000175037,  0.000000000000000,  0.000000000000000,  0.000527531006343, -0.000024327856198],
             [-0.000000000152588,  0.000000000000029,  0.000000000000000,  0.000000000000000, -0.000000088574316, -0.000000029379401],
             [ 0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000],
             [ 0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000,  0.000000000000000],
             [ 0.000527713959257, -0.000000101359885,  0.000000000000000,  0.000000000000000,  0.305484300478136, -0.013745076014696],
             [-0.000065095373481,  0.000000017572362,  0.000000000000000,  0.000000000000000, -0.013745076037145,  3.274109096814317]]
         ],
        rtol=0, atol=1e-8)

    assert_close(beamdata['mode_emittances'],
                 [1.320573957833556e-10, 0.0, 2.858758561755633e-06],
                 atol=1e-9)

    assert_close(obs['r66'], [
        [9.136524587096e-10, -3.482504553178e-15, -1.204862185189e-25, 9.079259501530e-26, 1.507759938693e-09, -1.184593870253e-15],
        [-3.482505991775e-15, 1.913557299753e-11, -4.047269952244e-29, -3.082870375266e-29, -2.511080835175e-13, -1.287047621221e-13],
        [-4.558591483169e-25, 7.673242260804e-29, 5.698353128109e-38, -4.743039437447e-38, -2.639781522189e-22, 3.360592031329e-23],
        [-1.793031624988e-25, 2.991500424099e-29, 5.352246900428e-38, -2.885933430723e-39, -1.038305919889e-22, 1.321918444089e-23],
        [1.507759938693e-09, -2.511080835186e-13, -6.905920845369e-23, 5.268875066089e-23, 8.731196550253e-07, 9.794607684193e-10],
        [-1.184856850504e-15, -1.287047620838e-13, 6.688601516891e-21, 1.180736964716e-21, 9.794606160964e-10, 9.357798394974e-06]
    ], atol=1E-10)

    assert_close(obs['r44'], [
        [9.110487602509e-10, -3.048897973283e-15, 0.000000000000e+00, 0.000000000000e+00],
        [-3.048899411882e-15, 1.913557292355e-11, 0.000000000000e+00, 0.000000000000e+00],
        [0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00],
        [0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00]
    ], atol=1E-20)

    assert_close(obs['m66'], [
        [-7.356310899844e-01, 4.673714021051e+00, 0.000000000000e+00, 0.000000000000e+00, 2.998514873422e-03, -6.276916774085e-07],
        [-9.816621640944e-02, -7.356663184957e-01, 0.000000000000e+00, 0.000000000000e+00, 1.690152160869e-04, -3.538066334557e-08],
        [0.000000000000e+00, 0.000000000000e+00, 6.098044855361e-01, -2.096029241593e+00, 0.000000000000e+00, 0.000000000000e+00],
        [0.000000000000e+00, 0.000000000000e+00, 2.996751797038e-01, 6.098002097888e-01, 0.000000000000e+00, 0.000000000000e+00],
        [1.246879423984e-06, 2.154487734820e-05, 0.000000000000e+00, 0.000000000000e+00, 9.999806908414e-01, -2.093301466764e-04],
        [1.700986573816e-04, 2.995796897591e-03, 0.000000000000e+00, 0.000000000000e+00, 2.243258562703e-03, 9.999995304090e-01]],
                 atol=1E-8)

    assert_close(obs['orbit6'],
                 [-2.635267925409e-09, -1.458487743359e-10, 0.000000000000e+00,
                  0.000000000000e+00, -6.696944241670e-06, -5.884521162896e-02],
                 atol=1E-20)

    assert_close(obs['emitXY'],
                 [1.320357526559e-10, 0.000000000000e+00],
                 atol=3e-12)
    assert_close(obs['emitXYZ'],
                 [1.322242916634e-10, 4.872515915668e-38, 2.858404580719e-06],
                 atol=3e-12)
