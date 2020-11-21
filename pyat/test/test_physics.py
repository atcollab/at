import at
import numpy
from numpy.testing import assert_allclose as assert_close
import pytest
from at import AtWarning, physics, atpass


DP = 1e-5
DP2 = 0.005

M44_MATLAB = numpy.array([[-0.66380202, 2.23414498, 0, 0],
                          [-0.25037182, -0.66380182, 0, 0],
                          [0, 0, -0.99921978, 0.262170798],
                          [0, 0, -0.0059496965, -0.99921979]])


def test_find_orbit4(dba_lattice):
    orbit4, _ = physics.find_orbit4(dba_lattice, DP)
    expected = numpy.array([1.091636e-7, 1.276747e-15, 0, 0, DP, 0])
    assert_close(orbit4, expected, atol=1e-12)


def test_find_orbit4_finds_zeros_if_dp_zero(dba_lattice):
    orbit4, _ = physics.find_orbit4(dba_lattice, 0)
    expected = numpy.zeros((6,))
    assert_close(orbit4, expected, atol=1e-7)


def test_find_orbit4_result_unchanged_by_atpass(dba_lattice):
    orbit, _ = physics.find_orbit4(dba_lattice, DP)
    orbit_copy = numpy.copy(orbit)
    orbit[4] = DP
    atpass(dba_lattice, orbit, 1)
    assert_close(orbit[:4], orbit_copy[:4], atol=1e-12)


def test_find_orbit4_with_two_refpts_with_and_without_guess(dba_lattice):
    expected = numpy.array(
        [[8.148212e-6, 1.0993354e-5, 0, 0, DP, 2.963929e-6],
         [3.0422808e-8, 9.1635269e-8, 0, 0, DP, 5.9280346e-6]]
    )
    _, all_points = physics.find_orbit4(dba_lattice, DP, [49, 99])
    assert_close(all_points, expected, atol=1e-12)
    _, all_points = physics.find_orbit4(dba_lattice, DP, [49, 99],
                                        numpy.array([0., 0., 0., 0., DP, 0.]))
    assert_close(all_points, expected, atol=1e-12)


def test_orbit_maxiter_warnings(hmba_lattice):
    with pytest.warns(AtWarning):
        physics.find_orbit4(hmba_lattice, max_iterations=1)
    with pytest.warns(AtWarning):
        physics.find_sync_orbit(hmba_lattice, max_iterations=1)
    with pytest.warns(AtWarning):
        physics.find_orbit6(hmba_lattice, max_iterations=1)


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
    m66, mstack = physics.find_m66(hmba_lattice, refpts=refpts)
    expected = numpy.array(
        [[-0.735654, 4.673766, 0., 0., 2.997161e-3, 0.],
         [-9.816788e-2, -0.735654, 0., 0., 1.695263e-4, 0.],
         [0., 0., 0.609804, -2.096051, 0., 0.],
         [0., 0., 0.299679, 0.609799, 0., 0.],
         [0., 0., 0., 0., 1., 0.],
         [1.695128e-4, 2.997255e-3, 0., 0., 2.243281e-3, 1.]])
    assert_close(m66, expected, rtol=1e-5, atol=1e-7)
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
    expected = numpy.zeros((len(hmba_lattice), 6))
    refpts = numpy.ones(len(hmba_lattice), dtype=bool)
    _, all_points = physics.find_orbit6(hmba_lattice, refpts)
    assert_close(all_points, expected, atol=1e-12)


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
    assert_close(tune, [0.355804634603528, 0.488487169156732], rtol=1e-8)
    assert_close(chrom, [-3.429156145429157, -1.597924048635235], rtol=2e-4)
    assert_close(obs['s_pos'], 56.209377216,  atol=1e-9)
    assert_close(obs['closed_orbit'][:5],
                 [0.000426438389644, -0.000000000287482, 0, 0, DP2], atol=1e-12)
    assert_close(obs['dispersion'],
                 [0.156822576442091, -0.000000162902610, 0, 0], rtol=1e-7, atol=2e-10)
    expected = [[-0.616893565445970, 2.191800192047084, 0, 0],
                [-0.282617257709865, -0.616893094967279, 0, 0],
                [0, 0, -0.997384485665288, 0.466909288129080],
                [0, 0, -0.011187515624062, -0.997384713772755]]
    assert_close(obs['m44'], expected, rtol=1e-7, atol=1e-12)
    assert_close(obs['alpha'], [-2.988886505944512e-07, 1.578070569086581e-06],
                 rtol=1e-8, atol=1e-8)
    assert_close(obs['beta'],
                 [2.784841119739221, 6.460251554763623], rtol=1e-8, atol=1e-12)
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
                 [0.000426438389644, -0.000000000287482, 0, 0, DP2], atol=1e-12)
    expected = [[-0.616893565445970, 2.191800192047084, 0, 0],
                [-0.282617257709865, -0.616893094967279, 0, 0],
                [0, 0, -0.997384485665288, 0.466909288129080],
                [0, 0, -0.011187515624062, -0.997384713772755]]
    assert_close(obs['m44'], expected, rtol=1e-7, atol=1e-12)
    assert_close(obs['alpha'],
                 [-2.988885294797512e-07, 1.578069929495847e-06],
                 rtol=1e-8, atol=1e-8)
    assert_close(obs['beta'],
                 [2.784839991078270, 6.460248936505217], rtol=1e-8, atol=1e-12)
    assert_close(obs['mu'],
                 [2.235586452367286, 3.069255403991328], rtol=1e-8, atol=1e-12)


def test_linopt_no_refpts(dba_lattice):
    lindata0, tune, chrom, lindata = physics.linopt(dba_lattice, DP,
                                                    get_chrom=True)
    assert list(lindata) == []
    assert len(physics.linopt(dba_lattice, DP, get_chrom=True)) == 4


@pytest.mark.parametrize('refpts', ([145], [1, 2, 3, 145]))
def test_linopt_line(hmba_lattice, refpts):
    refpts.append(len(hmba_lattice))
    l0, q, qp, ld = at.linopt(hmba_lattice, refpts=refpts, get_chrom=True)
    lt0, qt, qpt, ltd = at.linopt(hmba_lattice, refpts=refpts, twiss_in=l0, get_chrom=True)
    assert_close(ld['beta'], ltd['beta'], rtol=1e-12)
    assert_close(ld['s_pos'], ltd['s_pos'], rtol=1e-12)
    assert_close(ld['closed_orbit'], ltd['closed_orbit'], rtol=1e-12)
    assert_close(ld['alpha'], ltd['alpha'], rtol=1e-12)
    assert_close(ld['dispersion'], ltd['dispersion'], rtol=1e-7, atol=1e-12)
    assert_close(q, qt, rtol=1e-12)


def test_get_tune_chrom(hmba_lattice):
    qlin = hmba_lattice.get_tune()
    qplin = hmba_lattice.get_chrom()
    qharm = hmba_lattice.get_tune(method='laskar')
    qpharm = hmba_lattice.get_chrom(method='laskar')
    assert_close(qlin, [0.38156245, 0.85437541], rtol=1e-8)
    assert_close(qharm, [0.38156245, 0.85437541], rtol=1e-8)
    assert_close(qplin, [0.1791909, 0.12242558], rtol=1e-5)
    assert_close(qpharm, [0.17919145, 0.12242622], rtol=1e-5)


def test_nl_detuning_chromaticity(hmba_lattice):
    nlqplin, _, _ = at.nonlinear.chromaticity(hmba_lattice, npoints=11)
    nlqpharm, _, _ = at.nonlinear.chromaticity(hmba_lattice,
                                               method='laskar', npoints=11)
    q0, q1, _, _, _, _ = at.nonlinear.detuning(hmba_lattice,
                                               npoints=11, window=1)
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
    hmba_lattice = hmba_lattice.radiation_on(quadrupole_pass='auto', copy=True)
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


@pytest.mark.parametrize('refpts', ([121], [0, 40, 121]))
def test_ohmi_envelope(hmba_lattice, refpts):
    hmba_lattice = hmba_lattice.radiation_on(copy=True)
    emit0, beamdata, emit = hmba_lattice.ohmi_envelope(refpts)
    obs = emit[-1]
    assert_close(beamdata['tunes'], [3.81563019e-01, 8.54376397e-01, 1.09060761e-04], rtol=2e-6)
    assert_close(beamdata['damping_rates'], [1.00820168e-05, 6.57845603e-06, 9.64745588e-06])
    assert_close(
        beamdata['mode_matrices'], [
            [[6.90001508e+00, -2.60056091e-05, 0.0, 0.0, 6.78601607e-08, -2.24127548e-06],
             [-2.60056091e-05, 1.44927220e-01, 0.0, 0.0, -1.83087793e-06, -2.50275358e-04],
             [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
             [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
             [6.78601607e-08, -1.83087793e-06, 0.0, 0.0, 2.31302981e-11, 3.16172758e-09],
             [-2.24127548e-06, -2.50275358e-04, 0.0, 0.0, 3.16172758e-09, 4.32202209e-07]],
            [[0.0, 0.0, -1.76566316e-20, 7.71623094e-20, 0.0, 0.0],
             [0.0, 0.0, 2.12004507e-20, -6.21315614e-21, 0.0, 0.0],
             [-1.76566316e-20, 2.12004507e-20, 2.64468092e+00, 2.69748826e-06, 1.32045066e-16, 4.46538091e-17],
             [7.71623094e-20, -6.21315614e-21, 2.69748826e-06, 3.78117448e-01, 4.41614321e-17, 8.30871957e-17],
             [0.0, 0.0, 1.32045066e-16, 4.41614321e-17, 0.0, 0.0],
             [0.0, 0.0, 4.46538091e-17, 8.30871957e-17, 0.0, 0.0]],
            [[9.10978501e-07, -1.52946455e-10, -2.16234823e-19, -2.45466411e-20, 5.27532228e-04, -2.43420242e-05],
             [-1.52946455e-10, 2.60206584e-14, 3.67142758e-23, 4.16774785e-24, -8.85722067e-08, -2.93771029e-08],
             [-2.16234823e-19, 3.67142758e-23, 5.18182300e-32, 5.88232532e-33, -1.25222125e-16, -3.43353554e-17],
             [-2.45466411e-20, 4.16774785e-24, 5.88232532e-33, 6.67752471e-34, -1.42150211e-17, -3.89772629e-18],
             [5.27532228e-04, -8.85722067e-08, -1.25222125e-16, -1.42150211e-17, 3.05485018e-01, -1.37532801e-02],
             [-2.43420242e-05, -2.93771029e-08, -3.43353554e-17, -3.89772629e-18, -1.37532801e-02, 3.27410214e+00]],
        ],
        rtol=2e-6, atol=1e-20)
    assert_close(beamdata['mode_emittances'], [1.32038990e-10, 0.0, 2.86050277e-06], atol=1e-20)
    assert_close(obs['r66'],
                 [[9.13675728e-10, -3.48276991e-15, 0.0, 0.0, 1.50872017e-09, -8.62044057e-14],
                  [-3.48276991e-15, 1.91360269e-11, 0.0, 0.0, -2.51262711e-13, -1.28752116e-13],
                  [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                  [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                  [1.50872017e-09, -2.51262711e-13, 0.0, 0.0, 8.73675730e-07, 9.30851497e-10],
                  [-8.62044057e-14, -1.28752116e-13, 0.0, 0.0, 9.30851497e-10, 9.36372058e-06]],
                 atol=1E-20)
    assert_close(obs['r44'],
                 [[9.11070371e-10, -3.04889633e-15, 0.00000000e+00, 0.00000000e+00],
                  [-3.04889633e-15, 1.91360268e-11, 0.00000000e+00, 0.00000000e+00],
                  [0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00],
                  [0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00]],
                 atol=1E-20)
    assert_close(obs['m66'],
                 [[-7.35631091e-01, 4.67371402e+00, 0.00000000e+00, 0.00000000e+00, 2.99851473e-03, -6.27694955e-07],
                  [-9.81662164e-02, -7.35666318e-01, 0.00000000e+00, 0.00000000e+00, 1.69015223e-04, -3.53807757e-08],
                  [0.00000000e+00, 0.00000000e+00, 6.09804486e-01, -2.09602924e+00, 0.00000000e+00, 0.00000000e+00],
                  [0.00000000e+00, 0.00000000e+00, 2.99675180e-01, 6.09800210e-01, 0.00000000e+00, 0.00000000e+00],
                  [1.22308593e-06, 2.16108213e-05, 0.00000000e+00, 0.00000000e+00, 9.99980691e-01, -2.09331258e-04],
                  [1.70098657e-04, 2.99580777e-03, 0.00000000e+00, 0.00000000e+00, 2.24326044e-03, 9.99999542e-01]],
                 atol=1E-20)
    assert_close(obs['orbit6'], [-2.63520321e-09, -1.45845186e-10, 0.00000000e+00, 0.00000000e+00, -6.69677956e-06, -5.88436653e-02],
                 atol=1E-20)
    assert_close(obs['emitXY'], [1.32038885e-10, 0.00000000e+00], atol=1e-20)
    assert_close(obs['emitXYZ'], [1.32227544e-10, 0.00000000e+00, 2.86021932e-06], atol=1e-20)
