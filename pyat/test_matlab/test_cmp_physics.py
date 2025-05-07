from at import physics, get_s_pos, uint32_refpts
from scipy.constants import speed_of_light
import matlab
import numpy as np
from numpy.testing import assert_allclose as assert_close
import pytest


def _py_data(ml_data):
    """Convert a Matlab vector to a numpy vector"""
    return np.squeeze(np.asarray(ml_data))


def _ml_refs(refpts, nelems):
    """Convert refpoints to Matlab"""
    uintrefs = uint32_refpts(refpts, nelems)
    return matlab.double([float(ref + 1) for ref in uintrefs])


def _compare_physdata(py_data, ml_data, fields, atol=12, rtol=1.0e-7):
    for ml_key, py_key in fields:
        ml_val = _py_data(ml_data[ml_key])
        if py_key == "closed_orbit":
            py_val = np.squeeze(py_data[py_key][:, :4])
        else:
            py_val = np.squeeze(py_data[py_key])
        assert_close(py_val, ml_val, atol=atol, rtol=rtol)


@pytest.mark.parametrize("dp", (-0.01, 0.0, 0.01))
@pytest.mark.parametrize("lattices", ["dba", "hmba"])
def test_linear_analysis(engine, request, lattices, dp):
    py_lattice, ml_lattice, _ = request.getfixturevalue(lattices)
    nelems = len(py_lattice)
    fields = [
        ("SPos", "s_pos"),
        ("ClosedOrbit", "closed_orbit"),
        ("Dispersion", "dispersion"),
        ("alpha", "alpha"),
        ("beta", "beta"),
        ("mu", "mu"),
        ("M44", "m44"),
        ("A", "A"),
        ("B", "B"),
        ("C", "C"),
        ("gamma", "gamma"),
    ]
    refpts = range(nelems + 1)
    py_data0, py_tune, py_chrom, py_data = physics.linopt(py_lattice, dp, refpts, True)
    # Matlab call
    ml_data, ml_tune, ml_chrom = engine.pyproxy(
        "atlinopt", ml_lattice, dp, _ml_refs(refpts, nelems), nargout=3
    )
    ml_data0 = engine.pyproxy(
        "atlinopt", ml_lattice, dp, _ml_refs(nelems, nelems), nargout=3
    )[0]
    # Comparison
    assert_close(py_tune, _py_data(ml_tune), atol=1.0e-8, rtol=0)
    assert_close(py_chrom, _py_data(ml_chrom), atol=3.0e-3, rtol=0)
    _compare_physdata(np.expand_dims(py_data0, 0), ml_data0, fields, rtol=1e-8)
    _compare_physdata(py_data, ml_data, fields, rtol=1e-8)


@pytest.mark.parametrize("lattices", ["hmba"])
def test_ohmi_envelope(engine, request, lattices):
    py_lattice, ml_lattice, _ = request.getfixturevalue(lattices)
    fields = [
        ("beam66", "r66"),
        ("beam44", "r44"),
        ("emit66", "emitXYZ"),
        ("emit44", "emitXY"),
    ]
    nelems = len(py_lattice)
    refpts = range(nelems + 1)

    # Python call
    py_lattice = py_lattice.radiation_on(copy=True)
    py_emit0, py_beamdata, py_emit = py_lattice.ohmi_envelope(refpts)
    # Matlab call
    ml_emit = engine.pyproxy("atx", ml_lattice, 0.0, _ml_refs(refpts, nelems))
    ml_emit0, ml_params = engine.pyproxy(
        "atx", ml_lattice, 0.0, _ml_refs(0, nelems), nargout=2
    )
    revolution_period = get_s_pos(py_lattice, nelems)[0] / speed_of_light
    damping_times = revolution_period / py_beamdata.damping_rates
    # Comparison
    assert_close(
        damping_times, _py_data(ml_params["dampingtime"]), atol=1.0e-6, rtol=1.0e-3
    )
    assert_close(
        py_beamdata.mode_emittances,
        _py_data(ml_emit0["modemit"]),
        atol=1.0e-20,
        rtol=2.0e-3,
    )
    _compare_physdata(np.expand_dims(py_emit0, 0), ml_emit0, fields)
    _compare_physdata(py_emit, ml_emit, fields)


@pytest.mark.parametrize("lattices", ["hmba_rad"])
def test_quantdiff(engine, request, lattices):
    py_lattice, ml_lattice, radindex = request.getfixturevalue(lattices)
    # Python call
    dmat = physics.radiation.quantdiffmat(py_lattice)
    # Matlab call
    dmat_ml = engine.pyproxy("quantumDiff", ml_lattice, radindex)
    # Comparison
    assert_close(dmat, dmat_ml, rtol=1.0e-8, atol=1.0e-18)


@pytest.mark.parametrize("lattices", ["hmba"])
def test_fastring(engine, request, lattices):
    py_lattice, ml_lattice, _ = request.getfixturevalue(lattices)
    # Python call
    ring, ringrad = physics.fastring.fast_ring(py_lattice)
    # Matlab call
    ring_ml, ringrad_ml = engine.pyproxy("atfastring", ml_lattice, nargout=2)
    # Comparison
    for r, rml in zip([ring, ringrad], [ring_ml, ringrad_ml]):
        assert_close(r[0].Frequency, rml[1]["Frequency"], rtol=1.0e-20)
        assert_close(r[0].Voltage, rml[1]["Voltage"], rtol=1.0e-20)
        assert_close(r[1].I2, rml[2]["I2"], rtol=1.0e-20)
        assert_close(r[1].Length, rml[2]["Length"], rtol=1.0e-20)
        assert_close(r[1].M66, rml[2]["M66"], atol=1.0e-7)
        assert_close(r[1].T1, np.squeeze(rml[2]["T1"]), rtol=1.0e-8, atol=1.0e-11)
        assert_close(r[1].T2, np.squeeze(rml[2]["T2"]), rtol=1.0e-8, atol=1.0e-11)
        assert_close(r[2].A1, rml[-1]["A1"], rtol=0.02)
        assert_close(r[2].A2, rml[-1]["A2"], rtol=0.02)
        assert_close(r[2].A3, rml[-1]["A3"], rtol=0.02)
        assert_close(r[2].Alphax, rml[-1]["Alphax"], rtol=1.0e-5)
        assert_close(r[2].Alphay, rml[-1]["Alphay"], rtol=1.0e-6)
        assert_close(r[2].Betax, rml[-1]["Betax"], rtol=1.0e-10)
        assert_close(r[2].Betay, rml[-1]["Betay"], rtol=1.0e-10)
        assert_close(r[2].chromx_arr, rml[-1]["chromx_arr"], rtol=1.0e-5)
        assert_close(r[2].chromy_arr, rml[-1]["chromy_arr"], rtol=1.0e-5)
        assert_close(r[2].T1, np.squeeze(rml[-1]["T1"]), rtol=1.0e-8, atol=1.0e-11)
        assert_close(r[2].T2, np.squeeze(rml[-1]["T2"]), rtol=1.0e-8, atol=1.0e-11)
        if len(r) >= 4:
            assert_close(r[3].Lmatp, rml[-2]["Lmatp"], rtol=0.02, atol=1.0e-10)


@pytest.mark.parametrize("dp", (0.00, 0.01, -0.01))
@pytest.mark.parametrize("lattices", ["hmba"])
def test_parameters(engine, request, lattices, dp):
    py_lattice, ml_lattice, _ = request.getfixturevalue(lattices)

    # Test perimeter
    py_length = py_lattice.get_s_pos(len(py_lattice))
    ml_length = engine.findspos(ml_lattice, len(ml_lattice) + 1)
    assert_close(py_length, ml_length, rtol=1.0e-8)

    # test energy loss
    ml_energy, ml_periods, ml_voltage, ml_harms, ml_eloss = engine.pyproxy(
        "atenergy", ml_lattice, nargout=5
    )
    assert_close(py_lattice.rf_voltage, ml_voltage, rtol=1.0e-8)
    assert_close(py_lattice.energy_loss, ml_eloss, rtol=1.0e-6)
    assert py_lattice.energy == ml_energy
    assert py_lattice.periodicity == int(ml_periods)
    assert py_lattice.harmonic_number == int(ml_harms)

    # test momentum compaction factor
    py_mcf = py_lattice.get_mcf(dp)
    ml_mcf = engine.mcf(ml_lattice, dp)
    assert_close(py_mcf, ml_mcf, rtol=1.0e-8)
