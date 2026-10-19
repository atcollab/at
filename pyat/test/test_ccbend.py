"""CCBEND native interface, import/export, caching and geometry tests."""
import json
from pathlib import Path

import numpy as np
import pytest

import at
from at.load.elegant import CCBEND, ElegantParser, at2elegant


def test_ccbend_model0_reference():
    bend = at.CCBend("CB", 1.2, .12, -.8, FringeModel=0)
    np.testing.assert_allclose(at.element_pass(bend, np.zeros((6, 1))), 0, atol=3e-14)


def test_ccbend_import_export():
    parser = ElegantParser()
    parser.parse_lines("""CB: CCBEND, L=1.2, ANGLE=-0.12, K1=-0.8, K2=3, &
        FRINGEMODEL=1, FRINGE1K0=0.00001, FRINGE1I1=0.0002, &
        FRINGE2K7=0.0003, N_SLICES=24, INTEGRATION_ORDER=6, &
        OPTIMIZE_FSE=0, OPTIMIZE_DX=0
        RING: LINE=(CB,CB)""".splitlines())
    ring = parser.lattice(use="RING", energy=3e9)
    bend = ring[0]
    assert isinstance(bend, at.CCBend)
    assert bend.PassMethod == "CCBendPass"
    assert bend.NumIntSteps == 24
    assert bend.IntegrationOrder == 6
    assert bend.PolynomB[2] == 1.5
    np.testing.assert_array_equal(bend.fringeInt1, [1e-5, 0, 0, 0, 0, 0, 0, 2e-4])
    assert bend.fringeInt2[5] == 3e-4
    descr = at2elegant(at.CCBend)(bend.to_dict())
    reconstructed = descr.to_at(**descr)[0]
    np.testing.assert_array_equal(reconstructed.PolynomB, bend.PolynomB)
    np.testing.assert_array_equal(reconstructed.fringeInt1, bend.fringeInt1)
    assert reconstructed.BendingAngle == bend.BendingAngle


@pytest.mark.parametrize("fields", ({"NumIntSteps": 0}, {"IntegrationOrder": 5},
                                     {"FringeModel": 2}, {"Length": 0},
                                     {"BendingAngle": 0}, {"PolynomB": np.array([0.])},
                                     {"fringeInt1": np.zeros(7)},
                                     {"fringeInt2": np.full(8, np.nan)}))
def test_invalid_native_parameters(fields):
    bend = at.CCBend("CB", 1.2, .12, -.8)
    # Bypass class validation to exercise native array bounds checks as well.
    vars(bend).update(fields)
    with pytest.raises((RuntimeError, ValueError)):
        at.element_pass(bend, np.zeros((6, 1), order="F"))


def test_ccbend_geometry_protected():
    bend = at.CCBend("CB", 1.2, .12)
    with pytest.raises(at.AtError, match="split"):
        bend.divide([.5, .5])
    bend.set_longt_motion(True, new_pass="auto", energy=3e9)
    assert bend.PassMethod == "CCBendRadPass"
    bend.set_longt_motion(False, new_pass="auto")
    assert bend.PassMethod == "CCBendPass"


def test_lost_particles_and_delta():
    bend = at.CCBend("CB", 1.2, .12)
    pin = np.zeros((6, 3), order="F")
    pin[0, 0] = np.nan
    pin[1, 1] = 2.0  # no real longitudinal momentum
    pin[4, 2] = .02
    out = at.element_pass(bend, pin)
    assert np.isnan(out[0, 0]) and np.isnan(out[0, 1])
    assert np.isfinite(out[:, 2]).all()
    assert out[4, 2] == .02


@pytest.mark.parametrize("method", ("CCBendRadPass", "CCBendQuantPass"))
def test_radiation_energy_required(method):
    bend = at.CCBend("CB", 1.2, .12, PassMethod=method)
    with pytest.raises((RuntimeError, ValueError), match="Energy"):
        at.element_pass(bend, np.zeros((6, 1)))
    out = at.element_pass(bend, np.zeros((6, 1)), energy=3e9, seed=42)
    assert np.isfinite(out).all()
    if method == "CCBendRadPass":
        assert out[4, 0] < 0


def test_radiation_energy_scaling():
    bend = at.CCBend("CB", 1.2, .12, PassMethod="CCBendRadPass", NumIntSteps=24)
    low = at.element_pass(bend, np.zeros((6, 1)), energy=3e9)
    high = at.element_pass(bend, np.zeros((6, 1)), energy=6e9)
    # Classical fractional momentum loss scales as E^3 at fixed curvature.
    np.testing.assert_allclose(high[4, 0]/low[4, 0], 8, rtol=5e-4)


def test_quantum_statistics_and_diffusion():
    bend = at.CCBend("CB", 1.2, .12, PassMethod="CCBendRadPass", NumIntSteps=12)
    deterministic = at.element_pass(bend, np.zeros((6, 1)), energy=3e9)
    diffusion = at.diffusion_matrix(bend, np.zeros((6, 1), order="F"), energy=3e9)
    assert diffusion[4, 4] > 0
    np.testing.assert_allclose(diffusion, diffusion.T, atol=1e-25)
    assert np.linalg.eigvalsh(diffusion).min() > -1e-24
    # Analytic uniform-dipole ISR variance, before the tiny damping correction.
    gamma = 3e9/510998.95069
    p0 = np.sqrt(gamma**2-1)
    radius = 2.8179403205e-15
    variance = radius**2*55/(24*np.sqrt(3))*p0**5*137.0359895*(.12/1.2)**3*1.2
    np.testing.assert_allclose(diffusion[4, 4], variance, rtol=.003)
    bend.PassMethod = "CCBendQuantPass"
    pin = np.zeros((6, 20000), order="F")
    out = at.element_pass(bend, pin, energy=3e9, seed=42)
    # elegant's default +/-3-sigma rejection reduces sampled variance slightly.
    cutoff = 3.0
    import math
    accepted = math.erf(cutoff/np.sqrt(2))
    truncated_variance = 1-2*cutoff*np.exp(-cutoff**2/2)/np.sqrt(2*np.pi)/accepted
    np.testing.assert_allclose(np.var(out[4]), diffusion[4, 4]*truncated_variance, rtol=.035)
    assert abs(np.mean(out[4])-deterministic[4, 0]) < 5*np.sqrt(variance/len(out[4]))
    assert np.std(out[0]) > 0 and np.std(out[1]) > 0
    assert np.isfinite(out).all()


def test_quantum_rng_and_single_particle_option():
    bend = at.CCBend("CB", 1.2, .12, PassMethod="CCBendQuantPass")
    a = at.element_pass(bend, np.zeros((6, 5)), energy=3e9, seed=123)
    b = at.element_pass(bend, np.zeros((6, 5)), energy=3e9, seed=123)
    c = at.element_pass(bend, np.zeros((6, 5)), energy=3e9, seed=124)
    np.testing.assert_array_equal(a, b)
    assert not np.array_equal(a, c)
    bend.QuantumSingleParticle = 0
    quiet = at.element_pass(bend, np.zeros((6, 1)), energy=3e9, seed=123)
    bend.PassMethod = "CCBendRadPass"
    loss = at.element_pass(bend, np.zeros((6, 1)), energy=3e9)
    np.testing.assert_array_equal(quiet, loss)


@pytest.mark.parametrize("quantum", (False, True))
def test_radiative_import_export(quantum):
    descr = CCBEND(name="CB")
    bend = descr.to_at(l=1.2, angle=.12, synch_rad=1, isr=int(quantum), isr1part=0)[0]
    assert bend.PassMethod == ("CCBendQuantPass" if quantum else "CCBendRadPass")
    assert bend.QuantumSingleParticle == 0
    exported = at2elegant(at.CCBend)(bend.to_dict())
    assert exported["synch_rad"] == 1 and exported["isr"] == int(quantum)


def test_ccbend1_elegant_matrices():
    """Q5: all 22 parameters and full-precision R/T from elegant ccbend1."""
    ref = json.loads(Path(__file__).with_name("ccbend1_reference.json").read_text())
    p = ref["parameters"]
    params = {key.lower(): value for key, value in p.items()}
    bend = CCBEND(name=ref["family_name"]).to_at(**params)[0]
    bend.Energy = ref["energy_eV"]
    np.testing.assert_array_equal(bend.fringeInt1, ref["fringeInt1"])
    np.testing.assert_array_equal(bend.fringeInt2, ref["fringeInt2"])
    exported = at2elegant(at.CCBend)(bend.to_dict())
    for key, value in params.items():
        assert exported[key] == value
    # Independent four-degree fits on the example's 9-point sampling grid.
    n = np.arange(-4, 5)
    h = ref["matrix_step"]
    vectors = [np.zeros((6, 1))]
    for j in range(6):
        q = np.zeros((6, 9)); q[j] = n*h; vectors.append(q)
        for k in range(j):
            q = np.zeros((6, 81))
            q[j] = np.repeat(n, 9)*h; q[k] = np.tile(n, 9)*h
            vectors.append(q)
    q = np.concatenate(vectors, axis=1)
    rin = q[[0, 1, 2, 3, 5, 4]].copy(order="F")
    norm = (1+q[5])/np.sqrt(1+q[1]**2+q[3]**2)
    rin[1] = norm*q[1]; rin[3] = norm*q[3]
    out = at.element_pass(bend, rin.copy(order="F"))
    cached = at.lattice_pass([bend], rin.copy(order="F"), refpts=at.End)
    np.testing.assert_allclose(cached[:, :, 0, 0], out, atol=1e-14, rtol=0)
    pz = np.sqrt((1+out[4])**2-out[1]**2-out[3]**2)
    f = np.array([out[0], out[1]/pz, out[2], out[3]/pz, out[5]+bend.Length, out[4]])
    w = np.linalg.pinv(np.vander(n, 5, increasing=True))
    C = f[:, 0]; R = np.zeros((6, 6)); T = np.zeros((6, 6, 6)); offset = 1
    for j in range(6):
        axis = f[:, offset:offset+9]-C[:, None]; offset += 9
        R[:, j] = axis@w[1]/h; T[:, j, j] = axis@w[2]/h**2
        for k in range(j):
            pair = (f[:, offset:offset+81]-C[:, None]).reshape(6, 9, 9); offset += 81
            T[:, j, k] = np.einsum('ijk,j,k->i', pair, w[1], w[1])/h**2
    np.testing.assert_allclose(C, ref["C"], atol=1e-12, rtol=0)
    np.testing.assert_allclose(R, ref["R"], atol=2e-11, rtol=0)
    np.testing.assert_allclose(T, ref["T"], atol=2e-8, rtol=0)
