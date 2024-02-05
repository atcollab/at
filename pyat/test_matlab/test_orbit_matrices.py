import matlab
import numpy as np
import pytest
from numpy.testing import assert_allclose as assert_close

from at import find_orbit4, find_sync_orbit, find_orbit6
from at import uint32_refpts, find_m44, find_m66


def _py_data(ml_data):
    """Convert a Matlab vector to a numpy vector"""
    return np.squeeze(np.asarray(ml_data))


def _ml_refs(refpts, nelems):
    """Convert refpoints to Matlab"""
    uintrefs = uint32_refpts(refpts, nelems)
    return matlab.double([float(ref + 1) for ref in uintrefs])


@pytest.mark.parametrize("ct", (-0.0002, 0.0002))
@pytest.mark.parametrize("lattices", ["dba", "hmba"])
def test_find_syncorbit(engine, request, lattices, ct):
    py_lattice, ml_lattice, _ = request.getfixturevalue(lattices)
    nelems = len(py_lattice)
    nrefs = nelems + 1
    refpts = range(nrefs)

    # Python call
    py_orb5, py_orbit5 = find_sync_orbit(py_lattice, ct, refpts)
    # Matlab call
    ml_orbit5, ml_orb5 = engine.findsyncorbit(
        ml_lattice, ct, _ml_refs(refpts, nelems), nargout=2
    )
    ml_orbit5 = np.rollaxis(np.asarray(ml_orbit5), -1)
    # Comparison
    assert_close(py_orb5, _py_data(ml_orb5), atol=1.0e-12, rtol=0)
    assert_close(py_orbit5[:, :5], ml_orbit5, atol=1.0e-12, rtol=0)


@pytest.mark.parametrize("dp", (-0.01, 0.01))
@pytest.mark.parametrize("lattices", ["dba", "hmba"])
def test_find_orbit4(engine, request, lattices, dp):
    py_lattice, ml_lattice, _ = request.getfixturevalue(lattices)
    nelems = len(py_lattice)
    nrefs = nelems + 1
    refpts = range(nrefs)

    # Python call
    py_orb4, py_orbit4 = find_orbit4(py_lattice, dp, refpts)
    # Matlab call
    ml_orbit4, ml_orb4 = engine.findorbit4(
        ml_lattice, dp, _ml_refs(refpts, nelems), nargout=2
    )
    ml_orbit4 = np.rollaxis(np.asarray(ml_orbit4), -1)
    # Comparison
    assert_close(py_orb4, _py_data(ml_orb4), atol=1e-15, rtol=0)
    assert_close(py_orbit4[:, :4], ml_orbit4, atol=1.5e-15, rtol=0)


@pytest.mark.parametrize("lattices", ["hmba_rad"])
def test_find_orbit6(engine, request, lattices):
    py_lattice, ml_lattice, _ = request.getfixturevalue(lattices)
    nelems = len(py_lattice)
    nrefs = nelems + 1
    refpts = range(nrefs)

    # Python call
    py_orb6, py_orbit6 = find_orbit6(py_lattice, refpts)
    # Matlab call
    ml_orbit6, ml_orb6 = engine.findorbit6(
        ml_lattice, _ml_refs(refpts, nelems), nargout=2
    )
    ml_orbit6 = np.rollaxis(np.asarray(ml_orbit6), -1)
    # Comparison
    assert_close(py_orb6, _py_data(ml_orb6), atol=5e-15, rtol=0)
    assert_close(py_orbit6, ml_orbit6, atol=5e-15, rtol=0)


@pytest.mark.parametrize("dp", (-0.01, 0.0, 0.01))
@pytest.mark.parametrize("lattices", ["dba", "hmba"])
def test_find_m44(engine, request, lattices, dp):
    py_lattice, ml_lattice, _ = request.getfixturevalue(lattices)
    nelems = len(py_lattice)
    nrefs = nelems + 1
    refpts = range(nrefs)
    morbit = engine.zeros(6, 1)
    porbit = np.zeros((6, 1))

    # Python call
    py_m44, py_mstack = find_m44(py_lattice, dp, refpts, orbit=porbit)
    # Matlab call
    ml_m44, ml_mstack = engine.findm44(
        ml_lattice, dp, _ml_refs(refpts, nelems), "orbit", morbit, nargout=2
    )
    ml_mstack = np.rollaxis(np.asarray(ml_mstack).reshape((4, 4, nrefs)), -1)
    # Comparison
    assert_close(py_m44, np.asarray(ml_m44), atol=1.0e-11, rtol=0)
    assert_close(py_mstack, ml_mstack, atol=1.0e-11, rtol=0)


@pytest.mark.parametrize("lattices", ["hmba_rad"])
def test_find_m66_orb(engine, request, lattices):
    py_lattice, ml_lattice, _ = request.getfixturevalue(lattices)
    nelems = len(py_lattice)
    nrefs = nelems + 1
    refpts = range(nrefs)
    morbit = engine.zeros(6, 1)
    porbit = np.zeros((6, 1))

    # Python call
    py_m66, py_mstack = find_m66(py_lattice, refpts, orbit=porbit)
    # Matlab call
    ml_m66, ml_mstack = engine.findm66(
        ml_lattice, _ml_refs(refpts, nelems), "orbit", morbit, nargout=2
    )
    assert_close(py_m66, np.asarray(ml_m66), atol=1.0e-11, rtol=0)

    ml_mstack = np.rollaxis(np.asarray(ml_mstack).reshape((6, 6, nrefs)), -1)
    assert_close(py_mstack, ml_mstack, atol=1.0e-11, rtol=0)


@pytest.mark.parametrize("lattices", ["hmba_rad"])
def test_find_m66(engine, request, lattices):
    py_lattice, ml_lattice, _ = request.getfixturevalue(lattices)
    nelems = len(py_lattice)
    nrefs = nelems + 1
    refpts = range(nrefs)

    # Python call
    py_m66, py_mstack = find_m66(py_lattice, refpts)
    # Matlab call
    ml_m66, ml_mstack = engine.findm66(ml_lattice, _ml_refs(refpts, nelems), nargout=2)
    assert_close(py_m66, np.asarray(ml_m66), atol=1e-7, rtol=0)

    ml_mstack = np.rollaxis(np.asarray(ml_mstack).reshape((6, 6, nrefs)), -1)
    assert_close(py_mstack, ml_mstack, atol=1e-7, rtol=0)
