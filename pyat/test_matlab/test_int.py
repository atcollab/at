import numpy as np
import pytest

from at.physics import get_radiation_integrals


@pytest.mark.parametrize("lattices", ["hmba"])
def test_radiation_integrals(engine, request, lattices):
    py_lattice, ml_lattice, _ = request.getfixturevalue(lattices)

    # Python call
    py_integrals = get_radiation_integrals(py_lattice)
    # Matlab call
    results = engine.pyproxy("ringpara", ml_lattice)
    ml_integrals = np.squeeze(results["integrals"])
    # Comparison
    np.testing.assert_almost_equal(py_integrals, ml_integrals[:5])
