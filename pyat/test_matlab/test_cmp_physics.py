from at import physics
import matlab
import numpy
import pytest


@pytest.mark.parametrize('dp', (0, 1e-8, 1e-7, 1e-6))
@pytest.mark.parametrize('refpts', ([1], [1, 2, 3], None))
def test_find_orbit4(engine, ml_lattice, py_lattice, dp, refpts):

    def ml_convert(orbit):
        return numpy.rollaxis(numpy.asarray(orbit), -1)

    py_orb4, py_orbit4 = physics.find_orbit4(py_lattice, dp, refpts)
    if refpts is None:
        ml_orbit4, ml_orb4 = engine.findorbit4(ml_lattice, dp, matlab.double([]), nargout=2)
    else:
        ml_orbit4, ml_orb4 = engine.findorbit4(ml_lattice, dp, matlab.double(list(r + 1 for r in refpts)), nargout=2)
        # test refpoints
        numpy.testing.assert_almost_equal(ml_convert(ml_orbit4), py_orbit4[:, :4], decimal=8)
    # test initial orbit
    numpy.testing.assert_almost_equal(numpy.squeeze(numpy.asarray(ml_orb4)), py_orb4, decimal=8)


@pytest.mark.parametrize('dp', (0.0, 1e-8, 1e-7, 1e-6))
@pytest.mark.parametrize('refpts', ([1], [1, 2, 3], None))
def test_find_m44(engine, ml_lattice, py_lattice, dp, refpts):

    def ml_convert(matrix, nrefs):
        return numpy.rollaxis(numpy.asarray(matrix).reshape((4, 4, nrefs)), -1)

    py_m44, py_mstack = physics.find_m44(py_lattice, dp, refpts)
    if refpts is None:
        ml_m44 = engine.findm44(ml_lattice, dp)
    else:
        ml_m44, ml_mstack = engine.findm44(ml_lattice, dp, matlab.double(list(r + 1 for r in refpts)), nargout=2)
        # test refpoints: Matlab wrongly squeezes ml_stack if len(refpts)==1
        numpy.testing.assert_almost_equal(py_mstack, ml_convert(ml_mstack, len(refpts)), decimal=8)
    # test 1-turn matrix
    numpy.testing.assert_almost_equal(py_m44, numpy.asarray(ml_m44), decimal=8)


@pytest.mark.parametrize('dp', (0.0, 0.01))
@pytest.mark.parametrize('refpts', ([0], [0,1,2], None))
def test_linopt(engine, ml_lattice, py_lattice, dp, refpts):

    def compare_lindata(py_data, ml_data, decimal=8):
        for (ml_key, py_key) in [('SPos', 's_pos'), ('ClosedOrbit', 'closed_orbit'), ('Dispersion', 'dispersion'),
                                 ('beta', 'beta'), ('alpha', 'alpha'), ('mu', 'mu'), ('M44', 'm44'),
                                 ('A', 'A'), ('B', 'B'), ('C', 'C'), ('gamma', 'gamma')]:
            ml_val = numpy.squeeze(numpy.asarray(ml_data[ml_key]))
            if py_key == 'closed_orbit':
                py_val = numpy.squeeze(py_data[py_key][:, :4])
                numpy.testing.assert_almost_equal(py_val, ml_val, decimal=decimal)
            else:
                py_val = numpy.squeeze(py_data[py_key])
                numpy.testing.assert_almost_equal(py_val, ml_val, decimal=decimal)

    # Python call
    py_lindata0, py_nu, py_xsi, py_lindata = physics.linopt(py_lattice, 0.0, refpts, get_chrom=True)

    if refpts is None:
        # Matlab call
        ml_lindata, ml_nu, ml_xsi = engine.pyproxy('atlinopt', ml_lattice, 0.0, matlab.double([len(ml_lattice)+1]), nargout=3)
        # test initial values
        numpy.testing.assert_almost_equal(py_nu, numpy.squeeze(numpy.asarray(ml_nu)), decimal=6)
        numpy.testing.assert_almost_equal(py_xsi, numpy.squeeze(numpy.asarray(ml_xsi)), decimal=5)
        compare_lindata(numpy.expand_dims(py_lindata0, 0), ml_lindata, decimal=6)
    else:
        # Matlab call
        ml_lindata, ml_nu, ml_xsi = engine.pyproxy('atlinopt', ml_lattice, 0.0, matlab.double([r + 1 for r in refpts]), nargout=3)
        # test refpoints
        compare_lindata(py_lindata, ml_lindata)
