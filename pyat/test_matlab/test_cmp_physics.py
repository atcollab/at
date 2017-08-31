from at import physics
import matlab
import numpy
import pytest


@pytest.mark.parametrize('dp', (0, 1e-8, 1e-7, 1e-6))
@pytest.mark.parametrize('refpts', (None, [1], [1, 2, 3]))
def test_find_orbit4(engine, ml_lattice, py_lattice, dp, refpts):
    # Matlab call
    ml_refpts = (matlab.double([]) if refpts is None else
                 matlab.double(list(r + 1 for r in refpts)))
    ml_orbit4 = engine.findorbit4(ml_lattice, dp, ml_refpts)
    py_ml_orbit4 = numpy.asarray(ml_orbit4)

    # Python call
    py_orbit4 = physics.find_orbit4(py_lattice, dp, refpts)

    numpy.testing.assert_almost_equal(py_ml_orbit4, py_orbit4.T)


@pytest.mark.parametrize('dp', (0.0, 1e-8, 1e-7, 1e-6))
@pytest.mark.parametrize('refpts', (None, [1], [1, 2, 3], [145]))
def test_find_m44(engine, ml_lattice, py_lattice, dp, refpts):
    # Matlab call
    ml_refpts = (matlab.double([]) if refpts is None else
                 matlab.double(list(r + 1 for r in refpts)))
    ml_m44, ml_mstack = engine.findm44(ml_lattice, dp, ml_refpts, nargout=2)
    py_ml_m44 = numpy.asarray(ml_m44)

    # Python call
    py_m44, py_mstack = physics.find_m44(py_lattice, dp, refpts)

    py_mstack = numpy.squeeze(py_mstack)
    # Matches to 5 d.p.
    numpy.testing.assert_almost_equal(py_ml_m44, py_m44.T, decimal=5)
    assert py_mstack.T.shape == tuple(numpy.asarray(ml_mstack).shape)
    # Matches to 5 d.p.
    numpy.testing.assert_almost_equal(py_mstack.T, numpy.asarray(ml_mstack), decimal=5)
