import numpy
import pytest
from at.physics import get_radiation_integrals


@pytest.mark.parametrize('ml_lattice, py_lattice',
                         [(pytest.lazy_fixture('ml_hmba'),
                           pytest.lazy_fixture('py_hmba'))])
def test_radiation_integrals(engine, ml_lattice, py_lattice):

    # Python call
    py_integrals = get_radiation_integrals(py_lattice)
    # Matlab call
    results = engine.pyproxy('ringpara',ml_lattice)
    ml_integrals = numpy.squeeze(results['integrals'])
    # Comparison
    numpy.testing.assert_almost_equal(py_integrals, ml_integrals[:5])
