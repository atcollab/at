from at import physics, lattice
import matlab
import numpy
import pytest

num_elems = 145


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

    numpy.testing.assert_almost_equal(py_ml_orbit4, py_orbit4[:4, :])


@pytest.mark.parametrize('dp', (0.0, 1e-8, 1e-7, 1e-6))
@pytest.mark.parametrize('refpts', ([1], [1, 2, 3], [num_elems]))
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
    numpy.testing.assert_almost_equal(py_ml_m44, py_m44, decimal=8)
    assert py_mstack.shape == tuple(numpy.asarray(ml_mstack).shape)
    # Matches to 5 d.p.
    numpy.testing.assert_almost_equal(py_mstack, numpy.asarray(ml_mstack), decimal=8)


# Because of a limitation of the Matlab interface, on ly a songle reference point is allowed in functions returning
# structures

@pytest.mark.parametrize('dp', (0.0, 0.01))
@pytest.mark.parametrize('refpts', ([0], [10], [num_elems]))
def test_linopt(engine, ml_lattice, py_lattice, dp, refpts):

    def compare_lindata(py_data, ml_data, decimal=8):
        for (ml_key, py_key) in [('SPos', 's_pos'), ('ClosedOrbit', 'closed_orbit'), ('Dispersion', 'dispersion'),
                                 ('beta', 'beta'), ('alpha', 'alpha'), ('mu', 'mu'), ('M44', 'm44'),
                                 ('A', 'A'), ('B', 'B'), ('C', 'C'), ('gamma', 'gamma')]:
            ml_val = numpy.squeeze(numpy.asarray(ml_data[ml_key]))
            if py_key == 'closed_orbit':
                py_val = py_data[py_key][:4]
            else:
                py_val = py_data[py_key]
            numpy.testing.assert_almost_equal(py_val, ml_val, decimal=decimal)

        ml_refpts = matlab.double([r + 1 for r in refpts])
        # Matlab call
        ml_lindata, ml_nu, ml_xsi = engine.atlinopt(ml_lattice, 0.0, ml_refpts, nargout=3)
        # Python call
        py_lindata0, py_nu, py_xsi, py_lindata = physics.linopt(py_lattice, 0.0, refpts, get_chrom=True)

        if refpts[0] == len(py_lattice):
            numpy.testing.assert_almost_equal(py_nu, numpy.squeeze(numpy.asarray(ml_nu)), decimal=6)
            numpy.testing.assert_almost_equal(py_xsi, numpy.squeeze(numpy.asarray(ml_xsi)), decimal=5)
            compare_lindata(py_lindata0, ml_lindata)

        compare_lindata(py_lindata[0], ml_lindata)
