import matlab
import numpy as np
import pytest
from numpy.testing import assert_allclose as assert_close

from at import linopt6


def _py_data(ml_data):
    """Convert a Matlab vector to a numpy vector"""
    return np.squeeze(np.asarray(ml_data))


def _ml_refs(ring, py_refs):
    """Convert refpoints to Matlab"""
    uintrefs = ring.uint32_refpts(py_refs)
    return matlab.double([float(ref + 1) for ref in uintrefs])


def _compare_4d(py_data, ml_data, fields, atol=1.0e-12, rtol=1.0e-7):
    for ml_key, py_key in fields:
        ml_val = _py_data(ml_data[ml_key])
        if py_key == "closed_orbit":
            py_val = np.squeeze(py_data[py_key][:, :4])
        else:
            py_val = np.squeeze(py_data[py_key])
        assert_close(py_val, ml_val, atol=atol, rtol=rtol)


def _compare_6d(py_data, ml_data, fields, atol=1.0e-12, rtol=1.0e-7):
    for ml_key, py_key in fields:
        ml_val = _py_data(ml_data[ml_key])
        py_val = np.squeeze(py_data[py_key])
        assert_close(py_val, ml_val, atol=atol, rtol=rtol)


@pytest.mark.parametrize("dp", (-0.01, 0.0, 0.01))
@pytest.mark.parametrize("lattices", ["dba", "hmba"])
def test_4d_analysis(engine, request, lattices, dp):
    """Compare linopt6 in 4D"""
    py_lattice, ml_lattice, _ = request.getfixturevalue(lattices)
    fields = [
        ("SPos", "s_pos"),
        ("ClosedOrbit", "closed_orbit"),
        ("M", "M"),
        ("A", "A"),
        ("Dispersion", "dispersion"),
        ("alpha", "alpha"),
        ("beta", "beta"),
        ("mu", "mu"),
    ]
    pypts = range(10)
    mlpts = _ml_refs(py_lattice, pypts)

    # python call
    py_data0, py_ringdata, py_data = linopt6(
        py_lattice, refpts=pypts, dp=dp, get_w=True
    )
    # Matlab call
    ml_ringdata, ml_data = engine.pyproxy(
        "atlinopt6", ml_lattice, mlpts, "dp", dp, "get_w", nargout=2
    )
    # Comparison
    assert_close(py_ringdata.tune, _py_data(ml_ringdata["tune"]), atol=1.0e-8, rtol=0)
    assert_close(
        py_ringdata.chromaticity,
        _py_data(ml_ringdata["chromaticity"]),
        atol=5.0e-4,
        rtol=0,
    )

    _compare_6d(py_data, ml_data, fields, atol=1.0e-6, rtol=1e-8)
    _compare_6d(py_data, ml_data, [("W", "W")], atol=1.0e-6, rtol=4.0e-4)


@pytest.mark.parametrize("lattices", ["hmba"])
def test_6d_analysis(engine, request, lattices):
    """Compare linopt6 in 6D"""
    py_lattice, ml_lattice, _ = request.getfixturevalue(lattices)
    fields = [
        ("SPos", "s_pos"),
        ("ClosedOrbit", "closed_orbit"),
        ("M", "M"),
        ("A", "A"),
        ("Dispersion", "dispersion"),
        ("alpha", "alpha"),
        ("beta", "beta"),
        ("mu", "mu"),
        ("W", "W"),
    ]
    ml_lattice = engine.pyproxy("atradon", ml_lattice)
    py_lattice = py_lattice.radiation_on(copy=True)
    pypts = range(10)
    mlpts = _ml_refs(py_lattice, pypts)

    # python call
    py_data0, py_ringdata, py_data = linopt6(py_lattice, refpts=pypts, get_w=True)
    # Matlab call
    ml_ringdata, ml_data = engine.pyproxy(
        "atlinopt6", ml_lattice, mlpts, "get_w", nargout=2
    )
    # Comparison
    assert_close(py_ringdata.tune, _py_data(ml_ringdata["tune"]), atol=1.0e-12, rtol=0)
    assert_close(
        py_ringdata.chromaticity,
        _py_data(ml_ringdata["chromaticity"]),
        atol=1.0e-5,
        rtol=0,
    )

    _compare_6d(py_data, ml_data, fields, atol=1.0e-6, rtol=1e-8)
