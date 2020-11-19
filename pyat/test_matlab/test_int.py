import numpy
import pytest
from at.physics import get_radiation_integrals


@pytest.mark.parametrize('lattices',
                         [pytest.lazy_fixture('hmba')])
def test_radiation_integrals(engine, lattices):
    py_lattice, ml_lattice, _ = lattices

    # Python call
    py_integrals = get_radiation_integrals(py_lattice)
    # Matlab call
    results = engine.pyproxy('ringpara', ml_lattice)
    ml_integrals = numpy.squeeze(results['integrals'])
    # Comparison
    numpy.testing.assert_almost_equal(py_integrals, ml_integrals[:5])
