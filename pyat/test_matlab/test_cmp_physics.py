from at import physics
import matlab
import numpy
import pytest


@pytest.mark.parametrize('dp', (0, 1e-8, 1e-7, 1e-6))
@pytest.mark.parametrize('refpts', ([1], [1, 2, 3]))
def test_find_orbit4(engine, ml_lattice, py_lattice, dp, refpts):
    # Matlab call
    ml_refpts = (matlab.double([]) if refpts is None else
                 matlab.double(list(r + 1 for r in refpts)))
    ml_orbit4 = engine.findorbit4(ml_lattice, dp, ml_refpts)
    py_ml_orbit4 = numpy.asarray(ml_orbit4)

    # Python call
    orbs, py_orbit4 = physics.find_orbit4(py_lattice, dp, refpts)

    numpy.testing.assert_almost_equal(py_ml_orbit4, py_orbit4[:4,:])


@pytest.mark.parametrize('dp', (0.0, 1e-8, 1e-7, 1e-6))
@pytest.mark.parametrize('refpts', ([1], [1, 2, 3], [145]))
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
    numpy.testing.assert_almost_equal(py_ml_m44, py_m44, decimal=5)
    assert py_mstack.shape == tuple(numpy.asarray(ml_mstack).shape)
    # Matches to 5 d.p.
    numpy.testing.assert_almost_equal(py_mstack, numpy.asarray(ml_mstack), decimal=5)

@pytest.mark.parametrize('dp', (0.0, 1e-8, 1e-7, 1e-6))
@pytest.mark.parametrize('refpts', ([1], [1, 2, 3], [145]))
def test_get_twiss(engine, ml_lattice, py_lattice, dp, refpts):
    # Matlab call
    ml_refpts = (matlab.double([]) if refpts is None else
                 matlab.double(list(r + 1 for r in refpts)))
    ml_twiss = engine.twiss_transfer(ml_lattice, dp, ml_refpts)
    # Python call
    py_twiss = physics.get_twiss(py_lattice, dp, refpts)[3]
    # Local assignment
    ml_orbit = numpy.array(ml_twiss['O']).T
    ml_m44 = numpy.moveaxis(numpy.array(ml_twiss['M']), [0, 1], [-2, -1])
    ml_alpha = numpy.array(ml_twiss['A']).T
    ml_beta = numpy.array(ml_twiss['B']).T
    ml_mu = numpy.array(ml_twiss['U']).T  # Broken
    py_orbit = numpy.array(py_twiss['closed_orbit'])[:, :4]
    py_alpha = numpy.array(py_twiss['alpha'])
    py_beta = numpy.array(py_twiss['beta'])
    py_mu = numpy.array(py_twiss['mu'])  # Broken
    py_m44 = numpy.squeeze(numpy.array(py_twiss['m44']))  # Squeeze to resolve difference in array shape
    #  between the output of find_m44 (Python) and findm44 (Matlab) when len(refpts) is 1.
    # Matches to 10 d.p.
    numpy.testing.assert_almost_equal(py_orbit, ml_orbit, decimal=10)
    numpy.testing.assert_almost_equal(py_m44, ml_m44, decimal=10)
    numpy.testing.assert_almost_equal(py_alpha, ml_alpha, decimal=10)
    numpy.testing.assert_almost_equal(py_beta, ml_beta, decimal=10)
