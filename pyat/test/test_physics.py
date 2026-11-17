import at
import numpy as np
import pytest
from at import AtWarning, internal_lpass, lattice_pass, lattice_track, physics
from numpy.testing import assert_allclose as assert_close
from numpy.testing import assert_equal

DP = 1e-5
DP2 = 0.005

# fmt: off
M44_MATLAB = np.array(
    [[-0.66380202, 2.23414498, 0, 0],
     [-0.25037182, -0.66380182, 0, 0],
     [0, 0, -0.99921978, 0.262170798],
     [0, 0, -0.0059496965, -0.99921979]]
)

# New values with the symmetric radiation split and the default number of
# integration steps (same C integrators as Matlab). The horizontal closed orbit
# at the symmetry point of the HMBA cell is now zero to numerical precision
# (was -2.6e-9 m from the asymmetric radiation kick).
orbit6_MATLAB = np.array(
    [-2.4900915770091941e-13, -1.8787867374938079e-11,  0.0000000000000000e+00,
      0.0000000000000000e+00, -6.5785970288250395e-06, -5.8845229684990161e-02]
)
M66_MATLAB = np.array(
    [[-7.3562865717966230e-01,  4.6737183891599381e+00,  0.0000000000000000e+00,
       0.0000000000000000e+00,  2.9984670980781362e-03, -6.2768163690922679e-07],
     [-9.8166890600206669e-02, -7.3566388512222480e-01,  0.0000000000000000e+00,
       0.0000000000000000e+00,  1.6901135876663662e-04, -3.5379853102774782e-08],
     [ 0.0000000000000000e+00,  0.0000000000000000e+00,  6.0979972432272289e-01,
      -2.0960380737551336e+00,  0.0000000000000000e+00,  0.0000000000000000e+00],
     [ 0.0000000000000000e+00,  0.0000000000000000e+00,  2.9967677140040944e-01,
       6.0979515952440000e-01,  0.0000000000000000e+00,  0.0000000000000000e+00],
     [ 1.2471352561681105e-06,  2.1547275298471464e-05,  0.0000000000000000e+00,
       0.0000000000000000e+00,  9.9998069084755870e-01, -2.0933013372890466e-04],
     [ 1.7009541923111254e-04,  2.9957173316071581e-03,  0.0000000000000000e+00,
       0.0000000000000000e+00,  2.2432589304640280e-03,  9.9999953034890510e-01]]
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
    orbit6, _ = lattice.find_orbit6(refpts)
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
    """Compare with Matlab results."""
    _, tune, chrom, lindata = physics.linopt(dba_lattice, DP2, refpts, get_chrom=True)
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
    """Compare with Matlab results."""
    _, tune, _, lindata = physics.linopt(dba_lattice, DP2, refpts, coupled=False)
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
    _, _, _, lindata = physics.linopt(dba_lattice, DP, get_chrom=True)
    assert list(lindata) == []
    assert len(physics.linopt(dba_lattice, DP, get_chrom=True)) == 4


@pytest.mark.parametrize("refpts", ([121], [1, 2, 3, 121]))
def test_linopt_line(hmba_lattice, refpts):
    #    refpts.append(len(hmba_lattice))
    l0, _, _, ld = at.linopt(hmba_lattice, refpts=refpts)
    _, _, _, ltd = at.linopt(hmba_lattice, refpts=refpts, twiss_in=l0)
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
    assert_close(qlin, [0.209998303584, 0.340013166682], atol=1e-8)
    assert_close(qharm, [0.209998327527, 0.340013235807], atol=1e-8)
    assert_close(qplin, [5.729114185134, 3.931703139652], atol=1e-8)
    assert_close(qpharm, [5.729138523591, 3.931730165145], atol=1e-8)


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
                [2.10157015e-01, 5.72564645e00, 1.51815354e02, -1.89809140e04],
                [3.39970733e-01, 3.93108893e00, 2.58206212e02, -3.53120138e03],
            ]
        ),
        atol=1e-12,
        rtol=1e-5,
    )
    assert_close(
        nlqpharm,
        np.array(
            [
                [2.10156985e-01, 5.72564318e00, 1.51815337e02, -1.89809436e04],
                [3.39970750e-01, 3.93108797e00, 2.58206152e02, -3.53119072e03],
            ]
        ),
        atol=1e-12,
        rtol=1e-5,
    )
    assert_close(
        q0,
        np.array([[0.21000425, 0.34001688], [0.21000425, 0.34001688]]),
        atol=1e-12,
        rtol=1e-5,
    )
    assert_close(
        q1,
        np.array(
            [[96172.00066329, -103658.3658934], [-103704.24788915, 51570.58307421]]),
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
                [[ 1.4541985719833904e-07,  0.0000000000000000e+00,  0.0000000000000000e+00,
                   0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00],
                 [ 2.6087354817708038e-09,  1.7858729798055611e-08,  0.0000000000000000e+00,
                   0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00],
                 [-1.1817600965048644e-27, -2.1165670749147745e-24,  8.6522374271460576e-14,
                   0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00],
                 [-2.1913364288070925e-28, -5.1802969563265213e-27,  3.1110495235333243e-30,
                   8.6522374271460589e-14,  0.0000000000000000e+00,  0.0000000000000000e+00],
                 [ 3.2473092926023373e-07,  1.8127786848211402e-07,  4.4398079685742852e-18,
                   1.1707260764900171e-20,  5.7944863395245962e-06,  0.0000000000000000e+00],
                 [-1.6900856598354796e-09, -4.7225349451423462e-10, -7.0595588115280499e-19,
                  -2.7354874313493366e-20,  6.6352838696557402e-09,  3.2670509467417113e-09]]),
                rtol=1e-7, atol=1e-15)
    assert_close(dmat, np.array(
                [[ 2.1146934867581783e-14,  3.7936194122737923e-16,  0.0000000000000000e+00,
                   0.0000000000000000e+00,  4.7222325360910154e-14, -2.4577201530625062e-16],
                 [ 3.7936194122737923e-16,  3.2573973082191703e-16,  0.0000000000000000e+00,
                   0.0000000000000000e+00,  4.0845295687868189e-15, -1.2842833982781298e-17],
                 [ 0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
                   0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00],
                 [ 0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
                   0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00],
                 [ 4.7222325360910154e-14,  4.0845295687868189e-15,  0.0000000000000000e+00,
                   0.0000000000000000e+00,  3.3714383780956780e-11,  3.7813629547871473e-14],
                 [-2.4577201530625062e-16, -1.2842833982781298e-17,  0.0000000000000000e+00,
                   0.0000000000000000e+00,  3.7813629547871473e-14,  5.7780026820182457e-17]]),
                 rtol=1e-7, atol=1e-18)
    # fmt: on


def test_simple_ring():
    ring = physics.simple_ring(6e9, 844, 992, 0.1, 0.2, 6e6, 8.5e-5)
    assert_equal(len(ring), 5)
    num_cavities = len(ring[at.RFCavity])
    assert_equal(ring[num_cavities + 1].PassMethod, "SimpleQuantDiffPass")
    ring.disable_6d()
    assert_equal(ring[num_cavities + 1].PassMethod, "IdentityPass")
    assert_close(ring.get_tune(), [0.1, 0.2], atol=1e-10)


@pytest.mark.parametrize(
    "lattice",
    ["hmba_lattice", "noringparam_lattice"],
)
@pytest.mark.parametrize("refpts", ([121], [0, 40, 121]))
def test_ohmi_envelope(request, lattice, refpts):
    lattice = request.getfixturevalue(lattice).enable_6d(copy=True)
    _, beamdata, emit = lattice.ohmi_envelope(refpts)
    obs = emit[-1]

    # All expected values are Matlab results

    assert_close(
        beamdata["tunes"],
        [3.8156301859417241e-01, 8.5437639666919862e-01, 1.0906043886603917e-04],
        rtol=2e-6,
    )

    assert_close(
        beamdata["damping_rates"],
        [1.0082262545446862e-05, 6.5786199873676790e-06, 9.6536259190180764e-06],
        rtol=2e-6,
    )
    # fmt: off
    assert_close(
        beamdata["mode_matrices"], np.array([
        [[ 6.8999946070918163e+00, -2.6005722348292342e-05, -2.9781938608550528e-19,
          -7.3260921126587428e-20,  7.0330738927780742e-08, -1.9677147044391813e-06],
         [-2.6002605794226621e-05,  1.4492764960465360e-01, -2.4873095525343931e-20,
           2.2553584119619260e-21, -2.6206330535317807e-08, -2.5026556907180037e-04],
         [ 2.4211916511318089e-25, -3.8108648720349019e-26, -3.9099054296362707e-45,
          -3.1636667550371526e-45,  9.3585948985348524e-33,  6.5736578426430154e-29],
         [-1.1839508615720857e-25,  1.8634942617910471e-26,  1.9119246094869474e-45,
           1.5470175517094127e-45, -4.5763070792206616e-33, -3.2144864958703742e-29],
         [-2.7479256026558483e-08, -1.8309707338625978e-06,  3.1543919701949437e-25,
          -2.8198246679555716e-26,  3.3079918639866208e-13,  3.1617852731079379e-09],
         [-2.2428057790265535e-06, -2.5026898780137964e-04,  4.3050965197100805e-23,
          -3.8703865082922794e-24,  4.5231207383147443e-11,  4.3217293271967485e-07]],

        [[-4.7367812726960025e-40,  1.8761008670049967e-56, -1.0165704162815024e-20,
          -1.1254496537509481e-19, -7.0783579135145311e-37, -5.9255895791438120e-35],
         [ 1.0460698462254075e-41, -2.2583645417320225e-58, -2.8085171729165267e-20,
           1.6328731271644136e-20,  1.5631831718302292e-38,  1.3086060392065374e-36],
         [ 5.2127278835978783e-21, -2.2332397701745144e-37,  2.6446794645679828e+00,
           2.8798208715449149e-06,  7.7895835888706291e-18,  6.5209863507965734e-16],
         [ 1.5241057741222565e-21, -6.0147630004642340e-38,  2.8798208715449149e-06,
           3.7811765599792518e-01,  2.2775309954625032e-18,  1.9066164918169257e-16],
         [-7.5848689378629678e-38,  2.8738230111440134e-54,  1.7946451564768057e-17,
          -2.7593210563594824e-17, -1.1334366942118844e-34, -9.4884728151679796e-33],
         [ 2.5644421782631650e-38, -1.1718635439348159e-54,  2.4005833304205087e-17,
          -5.3765487081233685e-18,  3.8321464600639427e-35,  3.2080501448156755e-33]],

        [[ 9.1126159329021369e-07, -1.7499920516562724e-10,  9.2762623917747402e-21,
           4.9534175662593902e-21,  5.2752222637296464e-04, -2.4327369982396695e-05],
         [-1.5325863172415792e-10,  2.9378202781601049e-14,  8.7762086951757147e-25,
          -6.0829425161822012e-25, -8.8964965524663717e-08, -2.9361123715391454e-08],
         [-2.6647832080472849e-19,  5.1205999174948942e-23, -4.1363776326084372e-33,
          -1.5798022582209344e-33, -1.5411928235937421e-16,  2.6651751181524822e-17],
         [ 1.3030659401107661e-19, -2.5039482856510091e-23,  2.0226684077755680e-33,
           7.7251556846392663e-34,  7.5363574398976395e-17, -1.3032575822331510e-17],
         [ 5.2770059968619174e-04, -1.0133937959074853e-07,  5.3468016411928593e-18,
           2.8661618612050439e-18,  3.0548426675696755e-01, -1.3745026670191198e-02],
         [-6.5099744556956936e-05,  1.7752228401209417e-08, -2.3912439480121196e-16,
          -2.2342705476507698e-17, -1.3745026691931771e-02,  3.2741094537891997e+00]]
        ]),
        rtol=0, atol=1e-8)

    assert_close(
        beamdata["mode_emittances"],
        [1.3150325093260280e-10, 3.6844529119012627e-38, 2.8584014570799761e-06],
        atol=1e-9,
    )

    assert_close(
        obs["r66"],
        np.array([
        [ 9.0997531253512183e-10, -3.4499451305465016e-15, -3.1135748172010471e-25,
          1.3748944167543841e-25,  1.5077279804344444e-09, -7.7400205830988580e-16],
        [-3.4499451305553384e-15,  1.9058457172336212e-11,  5.2408777633896871e-29,
         -2.3107022771456250e-29, -2.5256779653254439e-13, -1.2857330566507051e-13],
        [-3.1135748172010476e-25,  5.2408777633899976e-29, -3.2419983706288599e-37,
         -6.5443500553074712e-38, -1.8030254035480886e-22,  2.2173818462329910e-23],
        [ 1.3748944167543841e-25, -2.3107022771456699e-29, -6.5443500553074702e-38,
         -4.9079439697095379e-38,  7.9618582969445665e-23, -1.0408662107264176e-23],
        [ 1.5077279804344442e-09, -2.5256779653257100e-13, -1.8030254035480889e-22,
          7.9618582969445677e-23,  8.7311846685946889e-07,  9.7970088027940138e-10],
        [-7.7400205830952845e-16, -1.2857330566507046e-13,  2.2173818462329907e-23,
         -1.0408662107264176e-23,  9.7970088027940158e-10,  9.3577906773634131e-06]
        ]),
        rtol=1.0E-12, atol=3.0E-12
    )

    assert_close(
        obs["r44"],
        np.array(
        [[ 9.0737172090745627e-10, -3.0138264496580735e-15, 0.0000000000000000e+00, 0.0000000000000000e+00],
         [-3.0138264496669135e-15,  1.9058457097516912e-11, 0.0000000000000000e+00, 0.0000000000000000e+00],
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
            -2.4900916101289919e-13,
            -1.8787867374667303e-11,
            0.0000000000000000e00,
            0.0000000000000000e00,
            -6.5785970288250514e-06,
            -5.8845229684990140e-02,
        ],
        # x, x' are zero to numerical precision at the symmetry point
        atol=1e-12,
    )

    assert_close(
        obs["emitXY"],
        [1.3150325093065836e-10, 0.000000000000000e00],
        rtol=1.0e-9,
        atol=2e-18,
    )
    assert_close(
        obs["emitXYZ"],
        [1.3169178224145862e-10, 0.0000000000000000e00, 2.8584014570381047e-06],
        rtol=1.0e-9,
        atol=1e-12,
    )


def test_rdt(hmba_lattice):
    rdtsd = {
        "refpts": [5],
        "h20000": [-2.95967705 - 1.77129145j],
        "h00200": [1.91618925 + 1.03505564j],
        "h11000": [342.6894347 + 0.0j],
        "h00110": [257.16817596 + 0.0j],
        "h10001": [0.00644487 + 0.00176723j],
        "h00002": [0.24292177 + 0.0j],
        "h10010": [0.0 + 0.0j],
        "h10100": [0.0 + 0.0j],
        "h11001": [18.80466492 + 0.0j],
        "h00111": [13.52540184 + 0.0j],
        "h20001": [-4.13510682 - 2.41357511j],
        "h00201": [1.71499168 + 0.93163935j],
        "h10002": [-0.00846011 - 0.00231724j],
        "h00003": [-0.07697747 + 0.0j],
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
        "h21001": [8.49759816 + 2.31940215j],
        "h30001": [19.48566545 + 20.03848294j],
        "h10021": [32.25603524 - 7.54967556j],
        "h10111": [-16.24214009 - 4.43318979j],
        "h10201": [47.9919823 + 45.87430765j],
        "h11002": [-388.82676151 + 0.0j],
        "h20002": [11.21340231 + 6.58871615j],
        "h00112": [236.73487128 + 0.0j],
        "h00202": [1.79019363 + 0.96681565j],
        "h10003": [0.01545246 + 0.00422576j],
        "h00004": [-0.00730645 + 0.0j],
        "h10400": [0.0 + 0.0j],
        "h10040": [0.0 + 0.0j],
        "h10310": [0.0 + 0.0j],
        "h10130": [0.0 + 0.0j],
        "h10220": [0.0 + 0.0j],
        "h30200": [0.0 + 0.0j],
        "h30020": [0.0 + 0.0j],
        "h21200": [0.0 + 0.0j],
        "h21020": [0.0 + 0.0j],
        "h30110": [0.0 + 0.0j],
        "h21110": [0.0 + 0.0j],
        "h50000": [0.0 + 0.0j],
        "h41000": [0.0 + 0.0j],
        "h32000": [0.0 + 0.0j],
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
