import numpy as np
import pytest
from numpy.testing import assert_allclose as assert_close
from numpy.testing import assert_equal

import at
from at import AtWarning, physics
from at import lattice_pass, internal_lpass
from at import lattice_track

DP = 1e-5
DP2 = 0.005

# fmt: off
M44_MATLAB = np.array(
    [[-0.66380202, 2.23414498, 0, 0],
     [-0.25037182, -0.66380182, 0, 0],
     [0, 0, -0.99921978, 0.262170798],
     [0, 0, -0.0059496965, -0.99921979]]
)

# New values with CODATA 2022 constants
orbit6_MATLAB = np.array(
    [-2.6352679435692906e-09, -1.4584877344311391e-10,  0.0000000000000000e+00,
      0.0000000000000000e+00, -6.6969442207863660e-06, -5.8845211247215173e-02]
)
# New values after update of the C constants
M66_MATLAB = np.array(
    [[-0.735631089984354,   4.673714021050702,   0.0,
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
       0.0,                 0.002243258562703,   0.999999530409042]]
)
# fmt: on


def test_find_orbit4(dba_lattice):
    orbit4, _ = physics.find_orbit4(dba_lattice, DP)
    expected = np.array([1.091636e-7, 1.276747e-15, 0, 0, DP, 0])
    assert_close(orbit4, expected, atol=1e-12)


def test_find_orbit4_finds_zeros_if_dp_zero(dba_lattice):
    orbit4, _ = physics.find_orbit4(dba_lattice, 0)
    expected = np.zeros((6,))
    assert_close(orbit4, expected, atol=1e-7)


@pytest.mark.parametrize("func", (lattice_track, lattice_pass, internal_lpass))
def test_find_orbit4_result_unchanged_by_atpass(dba_lattice, func):
    orbit, _ = physics.find_orbit4(dba_lattice, DP)
    orbit_copy = np.copy(orbit)
    orbit[4] = DP
    func(dba_lattice, orbit, 1)
    assert_close(orbit[:4], orbit_copy[:4], atol=1e-12)


def test_find_orbit4_produces_same_result_with_keep_lattice_True(dba_lattice):
    orbit0, _ = physics.find_orbit4(dba_lattice)
    orbit1, _ = physics.find_orbit4(dba_lattice, keep_lattice=True)
    assert_close(orbit0, orbit1, rtol=0, atol=1e-12)


def test_find_orbit4_with_two_refpts_with_and_without_guess(dba_lattice):
    expected = np.array(
        [
            [8.148212e-6, 1.0993354e-5, 0, 0, DP, 2.963929e-6],
            [3.0422808e-8, 9.1635269e-8, 0, 0, DP, 5.9280346e-6],
        ]
    )
    _, all_points = physics.find_orbit4(dba_lattice, DP, [49, 99])
    assert_close(all_points, expected, atol=1e-12)
    _, all_points = physics.find_orbit4(
        dba_lattice, DP, [49, 99], guess=np.array([0.0, 0.0, 0.0, 0.0, DP, 0.0])
    )
    assert_close(all_points, expected, atol=1e-12)


def test_orbit_maxiter_warnings(hmba_lattice):
    hmba_lattice_rad = hmba_lattice.radiation_on(copy=True)
    with pytest.warns(AtWarning):
        physics.find_orbit4(hmba_lattice, max_iterations=1)
    with pytest.warns(AtWarning):
        physics.find_sync_orbit(hmba_lattice, max_iterations=1)
    with pytest.warns(AtWarning):
        physics.find_orbit6(hmba_lattice_rad, max_iterations=1)


@pytest.mark.parametrize("refpts", ([145], [20], [1, 2, 3]))
def test_find_m44_returns_same_answer_as_matlab(dba_lattice, refpts):
    m44, mstack = physics.find_m44(dba_lattice, dp=DP, refpts=refpts)
    assert_close(m44, M44_MATLAB, rtol=1e-5, atol=1e-7)
    assert mstack.shape == (len(refpts), 4, 4)
    m44, mstack = physics.find_m44(dba_lattice, dp=DP, refpts=refpts, full=True)
    assert_close(m44, M44_MATLAB, rtol=1e-5, atol=1e-7)
    assert mstack.shape == (len(refpts), 4, 4)


@pytest.mark.parametrize(
    "lattice",
    ["hmba_lattice", "noringparam_lattice"],
)
@pytest.mark.parametrize("refpts", ([145], [20], [1, 2, 3]))
def test_find_m66(request, lattice, refpts):
    lattice = request.getfixturevalue(lattice).enable_6d(copy=True)
    m66, mstack = lattice.find_m66(refpts=refpts)
    assert_close(m66, M66_MATLAB, rtol=0, atol=1e-8)
    stack_size = 0 if refpts is None else len(refpts)
    assert mstack.shape == (stack_size, 6, 6)


@pytest.mark.parametrize("index", (19, 0, 1))
def test_find_elem_m66(hmba_lattice, index):
    m66 = physics.find_elem_m66(hmba_lattice[index])
    if index == 19:
        expected = np.array(
            [
                [1.0386, 0.180911, 0.0, 0.0, 0.0, 0.0],
                [0.434959, 1.0386, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.961891, 0.176344, 0.0, 0.0],
                [0.0, 0.0, -0.423978, 0.961891, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
            ]
        )
    else:
        expected = np.eye(6)
    assert_close(m66, expected, rtol=1e-5, atol=1e-7)


def test_find_sync_orbit(dba_lattice):
    # fmt: off
    expected = np.array(
        [[1.030844e-5, 1.390795e-5, -2.439041e-30,
          4.701621e-30, 1.265181e-5, 3.749859e-6],
         [3.86388e-8, 1.163782e-7, -9.671192e-30,
          3.567819e-30, 1.265181e-5, 7.5e-6]]
    )
    # fmt: on
    _, all_points = physics.find_sync_orbit(dba_lattice, DP, [49, 99])
    assert_close(all_points, expected, rtol=1e-5, atol=1e-7)


def test_find_sync_orbit_finds_zeros(dba_lattice):
    sync_orbit = physics.find_sync_orbit(dba_lattice)[0]
    np.testing.assert_equal(sync_orbit, np.zeros(6))


@pytest.mark.parametrize(
    "lattice",
    ["hmba_lattice", "noringparam_lattice"],
)
def test_find_orbit6(request, lattice):
    lattice = request.getfixturevalue(lattice).enable_6d(copy=True)
    refpts = np.ones(len(lattice), dtype=bool)
    orbit6, all_points = lattice.find_orbit6(refpts)
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
    expected = np.array(
        [
            [-0.66380, 2.23415, 0.0, 0.0],
            [-0.25037, -0.66380, 0.0, 0.0],
            [-1.45698e-31, -1.15008e-30, -0.99922, 0.26217],
            [6.57748e-33, 8.75482e-32, -5.9497e-3, -0.99922],
        ]
    )
    assert_close(m44, expected, rtol=1e-5, atol=1e-7)


@pytest.mark.parametrize("refpts", ([145], [1, 2, 3, 145]))
def test_linopt(dba_lattice, refpts):
    """Compare with Matlab results"""
    lindata0, tune, chrom, lindata = physics.linopt(
        dba_lattice, DP2, refpts, get_chrom=True
    )
    obs = lindata[-1]
    assert_close(tune, [0.355804633927360, 0.488487169156598], rtol=1e-8)
    assert_close(chrom, [-3.428312247995742, -1.597924047969101], rtol=2e-4)
    assert_close(obs["s_pos"], 56.209377216, atol=1e-9)
    assert_close(
        obs["closed_orbit"][:5],
        [0.000426438389644, -0.000000000287482, 0, 0, DP2],
        atol=1e-12,
    )
    assert_close(
        obs["dispersion"],
        [0.156822576442091, -0.000000162902610, 0, 0],
        rtol=1e-7,
        atol=2e-10,
    )
    expected = np.array(
        [
            [-0.616893565445970, 2.191800192047084, 0, 0],
            [-0.282617257709865, -0.616893094967279, 0, 0],
            [0, 0, -0.997384485665288, 0.466909288129080],
            [0, 0, -0.011187515624062, -0.997384713772755],
        ]
    )
    assert_close(obs["m44"], expected, rtol=1e-7, atol=1e-12)
    assert_close(
        obs["alpha"],
        [-2.988886505944512e-07, 1.578070569086581e-06],
        rtol=1e-8,
        atol=1e-8,
    )
    assert_close(
        obs["beta"], [2.784841119739221, 6.460251554763623], rtol=5e-8, atol=1e-12
    )
    assert_close(
        obs["mu"], [2.235586452367286, 3.069255403991328], rtol=1e-8, atol=1e-12
    )
    assert_close(obs["gamma"], 1.0, rtol=1e-6, atol=1e-20)
    assert_close(
        obs["A"],
        np.array(
            [
                [-0.616892545020345, 2.191796566512299],
                [-0.282616790222590, -0.616892074542433],
            ]
        ),
        rtol=1e-7,
        atol=1e-12,
    )
    assert_close(
        obs["B"],
        np.array(
            [
                [-0.997384485665288, 0.466909288129080],
                [-0.011187515624062, -0.997384713772754],
            ]
        ),
        rtol=1e-7,
        atol=1e-12,
    )
    assert_close(obs["C"], np.array([[0, 0], [0, 0]]), rtol=1e-5, atol=1e-7)


@pytest.mark.parametrize("refpts", ([145], [1, 2, 3, 145]))
def test_linopt_uncoupled(dba_lattice, refpts):
    """Compare with Matlab results"""
    lindata0, tune, chrom, lindata = physics.linopt(
        dba_lattice, DP2, refpts, coupled=False
    )
    obs = lindata[-1]
    assert_close(tune, [0.355804634603528, 0.488487169156732], rtol=1e-8)
    assert_close(obs["s_pos"], 56.209377216, atol=1e-9)
    assert_close(
        obs["closed_orbit"][:5],
        [0.000426438389644, -0.000000000287482, 0, 0, DP2],
        atol=1e-12,
    )
    expected = np.array(
        [
            [-0.616893565445970, 2.191800192047084, 0, 0],
            [-0.282617257709865, -0.616893094967279, 0, 0],
            [0, 0, -0.997384485665288, 0.466909288129080],
            [0, 0, -0.011187515624062, -0.997384713772755],
        ]
    )
    assert_close(obs["m44"], expected, rtol=1e-7, atol=1e-12)
    assert_close(
        obs["alpha"],
        [-2.988885294797512e-07, 1.578069929495847e-06],
        rtol=1e-8,
        atol=1e-8,
    )
    assert_close(
        obs["beta"], [2.784839991078270, 6.460248936505217], rtol=5e-8, atol=1e-12
    )
    assert_close(
        obs["mu"], [2.235586452367286, 3.069255403991328], rtol=1e-8, atol=1e-12
    )


def test_linopt_no_refpts(dba_lattice):
    lindata0, tune, chrom, lindata = physics.linopt(dba_lattice, DP, get_chrom=True)
    assert list(lindata) == []
    assert len(physics.linopt(dba_lattice, DP, get_chrom=True)) == 4


@pytest.mark.parametrize("refpts", ([121], [1, 2, 3, 121]))
def test_linopt_line(hmba_lattice, refpts):
    #    refpts.append(len(hmba_lattice))
    l0, q, qp, ld = at.linopt(hmba_lattice, refpts=refpts)
    lt0, qt, qpt, ltd = at.linopt(hmba_lattice, refpts=refpts, twiss_in=l0)
    assert_close(ld["beta"], ltd["beta"], rtol=1e-12)
    assert_close(ld["s_pos"], ltd["s_pos"], rtol=1e-12)
    assert_close(ld["closed_orbit"], ltd["closed_orbit"], rtol=1e-12)
    assert_close(ld["alpha"], ltd["alpha"], rtol=1e-12)
    assert_close(ld["dispersion"], ltd["dispersion"], rtol=1e-7, atol=1e-12)


def test_get_tune_chrom(hmba_lattice):
    qlin = hmba_lattice.get_tune()
    qplin = hmba_lattice.get_chrom()
    qharm = hmba_lattice.get_tune(method="interp_fft")
    qpharm = hmba_lattice.get_chrom(method="interp_fft")
    print(qlin, qharm)
    assert_close(qlin, [0.2099983, 0.34001317], atol=1e-8)
    assert_close(qharm, [0.20999833, 0.34001324], atol=1e-8)
    assert_close(qplin, [5.734099, 3.917612], atol=1e-8)
    assert_close(qpharm, [5.734123, 3.917639], atol=1e-8)


def test_nl_detuning_chromaticity(hmba_lattice):
    nlqplin, _, _ = at.nonlinear.chromaticity(hmba_lattice, npoints=11)
    nlqpharm, _, _ = at.nonlinear.chromaticity(
        hmba_lattice, method="interp_fft", npoints=11
    )
    q0, q1, _, _, _, _ = at.nonlinear.detuning(hmba_lattice, npoints=11)
    assert_close(
        nlqplin,
        np.array(
            [
                [0.2101570, 5.730634, 151.87972, -18977.6808],
                [0.3399707, 3.916998, 258.2324, -3529.81728],
            ]
        ),
        atol=1e-12,
        rtol=1e-5,
    )
    assert_close(
        nlqpharm,
        np.array(
            [
                [0.2101570, 5.730630, 151.87968, -18977.7132],
                [0.3399708, 3.916997, 258.23236, -3529.8072],
            ]
        ),
        atol=1e-12,
        rtol=1e-5,
    )
    assert_close(
        q0,
        np.array([[0.210004, 0.340017], [0.210004, 0.340017]]),
        atol=1e-12,
        rtol=1e-5,
    )
    assert_close(
        q1,
        np.array([[96183.925683, -104218.18371], [-104263.908197, 51684.400417]]),
        atol=1e-12,
        rtol=1e-5,
    )


def test_periodicity(hmba_lattice):
    q32 = at.linear.get_tune(hmba_lattice)
    qp32 = at.linear.get_chrom(hmba_lattice)
    dq32, *_ = at.nonlinear.detuning(hmba_lattice)
    dqp32, *_ = at.nonlinear.chromaticity(hmba_lattice, npoints=11)
    hmba_lattice.periodicity = 1
    hmba_lattice = hmba_lattice.repeat(32)
    q1 = at.linear.get_tune(hmba_lattice)
    qp1 = at.linear.get_chrom(hmba_lattice)
    dq1, *_ = at.nonlinear.detuning(hmba_lattice)
    dqp1, *_ = at.nonlinear.chromaticity(hmba_lattice, npoints=11)
    assert_close(q32, q1, atol=1e-12)
    assert_close(qp32, qp1, atol=1e-12)
    assert_close(dq32, dq1, atol=1e-7)
    assert_close(dqp32, dqp1, atol=1e-7)


def test_quantdiff(hmba_lattice):
    hmba_lattice = hmba_lattice.radiation_on(copy=True)
    dmat = physics.radiation.quantdiffmat(hmba_lattice)
    lmat = physics.radiation._lmat(dmat)
    # fmt: off
    assert_close(lmat, np.array(
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
                   0.00000000e+00, 6.52385922e-09, 3.25943528e-09]]),
                 rtol=1e-4, atol=1e-20)
    assert_close(dmat, np.array(
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
                   0.00000000e+00, 3.71123417e-14, 5.61789810e-17]]),
                 rtol=1e-5, atol=1e-20)
    # fmt: on


def test_simple_ring():
    ring = physics.simple_ring(6e9, 844, 992, 0.1, 0.2, 6e6, 8.5e-5)
    assert_equal(len(ring), 5)
    assert_equal(ring[-1].PassMethod, "SimpleQuantDiffPass")
    ring.disable_6d()
    assert_equal(ring[-1].PassMethod, "IdentityPass")
    assert_close(ring.get_tune(), [0.1, 0.2], atol=1e-10)


@pytest.mark.parametrize(
    "lattice",
    ["hmba_lattice", "noringparam_lattice"],
)
@pytest.mark.parametrize("refpts", ([121], [0, 40, 121]))
def test_ohmi_envelope(request, lattice, refpts):
    lattice = request.getfixturevalue(lattice).enable_6d(copy=True)
    emit0, beamdata, emit = lattice.ohmi_envelope(refpts)
    obs = emit[-1]

    # All expected values are Matlab results

    assert_close(
        beamdata["tunes"],
        [3.8156301859417241e-01, 8.5437639666919862e-01, 1.0906043886603917e-04],
        rtol=2e-6,
    )

    assert_close(
        beamdata["damping_rates"],
        [1.0082255849735316e-05, 6.5786188107455773e-06, 9.6535944754221195e-06],
        rtol=2e-6,
    )
    # fmt: off
    assert_close(
        beamdata["mode_matrices"], np.array([
        [[ 6.9000150766975468e+00, -2.6006243587486107e-05,  3.2653578639235223e-19,
          -3.7967745269216019e-19,  7.0332461709808689e-08, -1.9588952305233929e-06],
         [-2.6003127378572992e-05,  1.4492721966186267e-01, -1.9213610720373289e-20,
           2.0416359702093117e-20, -2.6206635396777689e-08, -2.5026896439390100e-04],
         [ 5.6444262224190266e-26, -8.8841529063567627e-27,  3.8488698134843282e-45,
          -4.3573005366593365e-45,  2.1817740491973322e-33,  1.5325293690413872e-29],
         [-1.9696012530780710e-26,  3.1000916669610156e-27, -1.3430486127116434e-45,
           1.5204635969825715e-45, -7.6132183004964948e-34, -5.3477034630266425e-30],
         [-2.7096903009017986e-08, -1.8307416373635866e-06,  2.4141125168786456e-25,
          -2.5639333062849341e-25,  3.3076648482142071e-13,  3.1614417843383483e-09],
         [-2.2432905706817794e-06, -2.5027456550832711e-04,  3.3071665955236401e-23,
          -3.5131068233223940e-23,  4.5232871220285155e-11,  4.3218970718496246e-07]],

        [[-2.2928294576283527e-40,  1.5052548450645591e-57, -2.9188149330766947e-19,
          -1.5140778203379946e-19, -3.5051801544373919e-37, -2.9766163106410513e-35],
         [-6.2160514262232196e-41,  3.7882279090924725e-58, -1.4553702864044294e-19,
          -1.7875951010571767e-20, -9.5028350345329565e-38, -8.0698544767557654e-36],
         [ 8.3548833738190147e-22, -4.7882541912711381e-39,  2.6446809160227356e+00,
           2.6974783759959333e-06,  1.2772590345573058e-18,  1.0846546846854328e-16],
         [ 3.4232225141711210e-22, -2.4394082787927380e-39,  2.6974783759959333e-06,
           3.7811744847886963e-01,  5.2332769805336290e-19,  4.4441246760564963e-17],
         [-3.8353443543983292e-38,  2.9389076758041900e-55,  4.6703386529315050e-17,
          -5.8660977848835226e-17, -5.8633113211900941e-35, -4.9791529519320920e-33],
         [-9.9107960119587158e-38,  6.8871348009966869e-55, -3.9791944217300570e-17,
          -9.5586285962631032e-17, -1.5151203409488728e-34, -1.2866476816442893e-32]],

        [[ 9.1129972881601156e-07, -1.7503747040171271e-10,  1.3713152284793508e-19,
          -7.5075423107737713e-21,  5.2753100634258230e-04, -2.4327856196424929e-05],
         [-1.5258771338022758e-10,  2.9256399010740371e-14, -1.7833250548974794e-23,
           1.8515034658507326e-24, -8.8574315741309644e-08, -2.9379401411635980e-08],
         [-6.2113700707272211e-20,  1.1937412750029960e-23, -1.0034900096364489e-32,
           4.3194611901458771e-34, -3.5923395227888331e-17,  6.1469644287517179e-18],
         [ 2.1674341717930857e-20, -4.1655151798392467e-24,  3.5016405610563869e-33,
           -1.5072596996559408e-34,  1.2535333347903046e-17, -2.1449600658093055e-18],
         [ 5.2771395924912238e-04, -1.0135988186593673e-07,  7.9357371970076673e-17,
          -4.3535442566604639e-18,  3.0548430047810055e-01, -1.3745076013796528e-02],
         [-6.5095371736554014e-05,  1.7571715133671204e-08, -5.1141897662730421e-16,
          -5.7612397975267779e-17, -1.3745076036246284e-02,  3.2741090968146191e+00]]
        ]),
        rtol=0, atol=1e-8)

    assert_close(
        beamdata["mode_emittances"],
        [1.3203595267442241e-10, 1.1573828900796540e-37, 2.8586919200065274e-06],
        atol=1e-9,
    )

    assert_close(
        obs["r66"],
        np.array([
        [ 9.1365311679111472e-10, -3.4825080173542211e-15,  6.9592739907066549e-25,
          3.2842455917050560e-25,  1.5077619241445821e-09, -1.2683061573113966e-15],
        [-3.4825080173566153e-15,  1.9135586746581878e-11, -1.1673711355982052e-28,
         -5.5030817148786995e-29, -2.5110841375848474e-13, -1.2870496457572302e-13],
        [ 6.9592739907066539e-25, -1.1673711355981774e-28,  1.0481996276492947e-37,
          1.2657948880467825e-37,  4.0299728048696993e-22, -4.9337080646554969e-23],
        [ 3.2842455917050565e-25, -5.5030817148784842e-29,  1.2657948880467825e-37,
          5.5172109613132283e-38,  1.9018302161174277e-22, -2.3925455405223881e-23],
        [ 1.5077619241445821e-09, -2.5110841375848878e-13,  4.0299728048696993e-22,
          1.9018302161174277e-22,  8.7312080470949738e-07,  9.7941378315054599e-10],
        [-1.2683061573110081e-15, -1.2870496457572292e-13, -4.9337080646554964e-23,
         -2.3925455405223884e-23,  9.7941378315054641e-10,  9.3578103418469929e-06]
        ]),
        rtol=1.0E-12, atol=2.0E-12
    )

    assert_close(
        obs["r44"],
        np.array(
        [[ 9.1104941490365636e-10, -3.0489008671947536e-15, 0.0000000000000000e+00, 0.0000000000000000e+00],
         [-3.0489008671971892e-15,  1.9135586672600989e-11, 0.0000000000000000e+00, 0.0000000000000000e+00],
         [ 0.0000000000000000e+00,  0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00],
         [ 0.0000000000000000e+00,  0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00]],
        ),
        rtol=1E-12, atol=1.0E-15
    )

    assert_close(
        obs["m66"],
        np.array(
        [[-7.3563108998454019e-01,  4.6737140210495083e+00,  0.0000000000000000e+00,
           0.0000000000000000e+00,  2.9985148734129178e-03, -6.2769167863546892e-07],
         [-9.8166216409474760e-02, -7.3566631849557673e-01,  0.0000000000000000e+00,
           0.0000000000000000e+00,  1.6901521608948762e-04, -3.5380663193032396e-08],
         [ 0.0000000000000000e+00,  0.0000000000000000e+00,  6.0980448553604349e-01,
          -2.0960292415936137e+00,  0.0000000000000000e+00,  0.0000000000000000e+00],
         [ 0.0000000000000000e+00,  0.0000000000000000e+00,  2.9967517970374508e-01,
           6.0980020978885463e-01,  0.0000000000000000e+00,  0.0000000000000000e+00],
         [ 1.2468793675147450e-06,  2.1544876952915747e-05,  0.0000000000000000e+00,
           0.0000000000000000e+00,  9.9998069084148322e-01, -2.0933014695165365e-04],
         [ 1.7009958256745489e-04,  2.9957945846259553e-03,  0.0000000000000000e+00,
           0.0000000000000000e+00,  2.2432585580767217e-03,  9.9999953041829404e-01]],
        ),
        atol=5.0E10
    )
    # fmt: on

    assert_close(
        obs["orbit6"],
        [
            -2.6352679435773159e-09,
            -1.4584877344363697e-10,
            0.0000000000000000e00,
            0.0000000000000000e00,
            -6.6969442207863711e-06,
            -5.8845211247215208e-02,
        ],
        atol=1e-20,
    )

    assert_close(
        obs["emitXY"],
        [1.320358475286751e-10, 0.000000000000000e00],
        rtol=1.0e-9,
        atol=1e-18,
    )
    assert_close(
        obs["emitXYZ"],
        [1.3222438678441061e-10, 0.0000000000000000e00, 2.8584082872712468e-06],
        rtol=1.0e-9,
        atol=1e-12,
    )


def test_rdt(hmba_lattice):
    rdtsd = {
        "refpts": [5],
        "h20000": [-2.95967705 - 1.77129145j],
        "h00200": [1.91618925 + 1.03505564j],
        "h10010": [0.0 + 0.0j],
        "h10100": [0.0 + 0.0j],
        "h11001": [18.80466492 + 0.0j],
        "h00111": [13.52540184 + 0.0j],
        "h20001": [-4.13510682 - 2.41357511j],
        "h00201": [1.71499168 + 0.93163935j],
        "h10002": [-0.00846011 - 0.00231724j],
        "h21000": [-0.82318652 - 0.22468352j],
        "h30000": [-1.93582223 - 1.99073744j],
        "h10110": [5.86104821 + 1.60003465j],
        "h10020": [7.9571528 - 1.86220774j],
        "h10200": [-3.24150177 - 3.09881147j],
        "dnux_dJx": [-86180.58479521 + 0.0j],
        "dnux_dJy": [116574.53372438 + 0.0j],
        "dnuy_dJy": [-39421.93576833 + 0.0j],
        "h22000": [-67686.07301868 + 0.0j],
        "h11110": [183114.84937209 + 0.0j],
        "h00220": [-30961.91595002 + 0.0j],
        "h31000": [1726.39070848 + 1018.3135644j],
        "h40000": [80.47712691 + 145.59516188j],
        "h20110": [-3502.874286 - 2066.17535059j],
        "h11200": [2564.73587249 + 1389.09487007j],
        "h20020": [6026.4984792 + 220.32116222j],
        "h20200": [-210.39145057 - 349.80230878j],
        "h00310": [-578.20990151 - 313.16651706j],
        "h00400": [-68.59147195 - 105.14365305j],
    }

    ring = at.Lattice(hmba_lattice, periodicity=32)
    ring.disable_6d()
    mult0 = ring.get_uint32_index(at.Multipole)[0]
    _, _, rdt = ring.get_rdts(mult0, at.RDTType.ALL)

    ring = at.Lattice(hmba_lattice, periodicity=1)
    ring.disable_6d()
    ring = ring.repeat(32)
    mult0 = ring.get_uint32_index(at.Multipole)[0]
    _, _, rdt32 = ring.get_rdts(mult0, at.RDTType.ALL)

    for k in rdt.dtype.names:
        assert_close(np.absolute(rdt[k]), np.absolute(rdt32[k]), rtol=1.0e-8)
        assert_close(np.absolute(rdt[k]), np.absolute(rdtsd[k]), atol=1.0e-8)
