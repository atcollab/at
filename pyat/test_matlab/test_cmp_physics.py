import math
from at import physics, get_s_pos, uint32_refpts
from scipy.constants import speed_of_light
import matlab
import numpy
import pytest


def _py_data(ml_data):
    """Convert a Matlab vector to a numpy vector"""
    return numpy.squeeze(numpy.asarray(ml_data))


def _ml_refs(refpts, nelems):
    """Convert refpoints to Matlab"""
    uintrefs = uint32_refpts(refpts, nelems)
    return matlab.double([ref+1 for ref in uintrefs])


@pytest.mark.parametrize('dp', (0, 1e-8, 1e-7, 1e-6))
@pytest.mark.parametrize('refpts', (0, [0, 1, 2, -1], None, '-'))
def test_find_orbit4(engine, ml_lattice, py_lattice, dp, refpts):

    def ml_convert(orbit):
        return numpy.rollaxis(numpy.asarray(orbit), -1)

    nelems = len(py_lattice)
    if isinstance(refpts, str):
        refpts = range(nelems+1)

    # Python call
    py_orb4, py_orbit4 = physics.find_orbit4(py_lattice, dp, refpts)
    ml_orbit4, ml_orb4 = engine.findorbit4(ml_lattice, dp, _ml_refs(refpts, nelems), nargout=2)

    if refpts is None:
        # test initial orbit
        numpy.testing.assert_almost_equal(py_orb4, _py_data(ml_orb4), decimal=8)
    else:
        # test refpoints
        numpy.testing.assert_almost_equal(py_orbit4[:, :4], ml_convert(ml_orbit4), decimal=8)


@pytest.mark.parametrize('dp', (0.0, 1e-8, 1e-7, 1e-6))
@pytest.mark.parametrize('refpts', (0, [0, 1, 2, -1], None, '-'))
def test_find_m44(engine, ml_lattice, py_lattice, dp, refpts):

    def ml_convert(matrix, nrefs):
        return numpy.rollaxis(numpy.asarray(matrix).reshape((4, 4, nrefs)), -1)

    nelems = len(py_lattice)
    if isinstance(refpts, str):
        refpts = range(nelems+1)

    # Python call
    py_m44, py_mstack = physics.find_m44(py_lattice, dp, refpts)

    if refpts is None:
        ml_m44 = engine.findm44(ml_lattice, dp)
        # test 1-turn matrix
        numpy.testing.assert_almost_equal(py_m44, numpy.asarray(ml_m44), decimal=8)
    else:
        ml_m44, ml_mstack = engine.findm44(ml_lattice, dp, _ml_refs(refpts, nelems), nargout=2)
        # test refpoints: Matlab wrongly squeezes ml_stack if len(refpts)==1
        nrefs = len(uint32_refpts(refpts, nelems))
        numpy.testing.assert_almost_equal(py_mstack, ml_convert(ml_mstack, nrefs), decimal=8)


@pytest.mark.parametrize('dp', (0.0, 0.01))
@pytest.mark.parametrize('refpts', (0, [0, 1, 2, -1], None, '-'))
def test_linopt(engine, ml_lattice, py_lattice, dp, refpts):

    def compare_lindata(py_data, ml_data, decimal=8):
        for (ml_key, py_key) in [('SPos', 's_pos'), ('ClosedOrbit', 'closed_orbit'), ('Dispersion', 'dispersion'),
                                 ('beta', 'beta'), ('alpha', 'alpha'), ('mu', 'mu'), ('M44', 'm44'),
                                 ('A', 'A'), ('B', 'B'), ('C', 'C'), ('gamma', 'gamma')]:
            ml_val = _py_data(ml_data[ml_key])
            if py_key == 'closed_orbit':
                py_val = numpy.squeeze(py_data[py_key][:, :4])
            else:
                py_val = numpy.squeeze(py_data[py_key])
            numpy.testing.assert_almost_equal(py_val, ml_val, decimal=decimal)

    nelems = len(py_lattice)
    if isinstance(refpts, str):
        refpts = range(nelems+1)

    # Python call
    py_lindata0, py_nu, py_xsi, py_lindata = physics.linopt(py_lattice, 0.0, refpts, get_chrom=True)

    if refpts is None:
        # Matlab call for last element
        ml_lindata, ml_nu, ml_xsi = engine.pyproxy('atlinopt', ml_lattice, 0.0, _ml_refs(nelems, nelems), nargout=3)
        # test initial values
        numpy.testing.assert_almost_equal(py_nu, _py_data(ml_nu), decimal=6)
        numpy.testing.assert_almost_equal(py_xsi, _py_data(ml_xsi), decimal=5)
        compare_lindata(numpy.expand_dims(py_lindata0, 0), ml_lindata, decimal=6)
    else:
        # Matlab call
        ml_lindata, ml_nu, ml_xsi = engine.pyproxy('atlinopt', ml_lattice, 0.0, _ml_refs(refpts, nelems), nargout=3)
        # test refpoints
        compare_lindata(py_lindata, ml_lindata)


@pytest.mark.parametrize('refpts', (0, [0, 1, 2, -1], None, '-'))
def test_ohmi_envelope(engine, ml_lattice, py_lattice, refpts):

    def compare_lindata(py_data, ml_data, decimal=8):
        for (ml_key, py_key) in [('beam66', 'r66'), ('beam44', 'r44'), ('emit66', 'emitXYZ'), ('emit44', 'emitXY')]:
            ml_val = _py_data(ml_data[ml_key])
            py_val = numpy.squeeze(py_data[py_key])
            numpy.testing.assert_almost_equal(py_val, ml_val, decimal=decimal)

    nelems = len(py_lattice)
    if isinstance(refpts, str):
        refpts = range(nelems+1)

    # Python call
    py_lattice.radiation_on()
    emit0, beamdata, emit = physics.ohmi_envelope(py_lattice, refpts)

    if refpts is None:
        # Matlab call at initial point
        ml_data, ml_params = engine.pyproxy('atx', ml_lattice, 0.0, _ml_refs(0, nelems), nargout=2)
        # test initial values
        revolution_period = get_s_pos(py_lattice, nelems) / speed_of_light
        damping_times = revolution_period / beamdata.damping_rates
        numpy.testing.assert_almost_equal(damping_times, _py_data(ml_params['dampingtime']), decimal=8)
        numpy.testing.assert_almost_equal(beamdata.mode_emittances, _py_data(ml_data['modemit']), decimal=8)
        compare_lindata(numpy.expand_dims(emit0, 0), ml_data)
    else:
        # Matlab call
        ml_data = engine.pyproxy('atx', ml_lattice, 0.0, _ml_refs(refpts, nelems))
        # test refpoints
        compare_lindata(emit, ml_data)
