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


def _compare_physdata(py_data, ml_data, fields, decimal=8):
    for (ml_key, py_key) in fields:
        ml_val = _py_data(ml_data[ml_key])
        if py_key == 'closed_orbit':
            py_val = numpy.squeeze(py_data[py_key][:, :4])
        else:
            py_val = numpy.squeeze(py_data[py_key])
        numpy.testing.assert_almost_equal(py_val, ml_val, decimal=decimal)


@pytest.mark.parametrize('dp', (-0.01, 0.0, 0.01))
@pytest.mark.parametrize('refpts', (0, [0, 1, 2, -1], None))
@pytest.mark.parametrize('ml_lattice, py_lattice',
                         [(pytest.lazy_fixture('ml_dba'), pytest.lazy_fixture('py_dba')),
                          (pytest.lazy_fixture('ml_hmba'), pytest.lazy_fixture('py_hmba'))])
def test_find_orbit4(engine, ml_lattice, py_lattice, dp, refpts):
    nelems = len(py_lattice)
    refpts = range(nelems + 1) if refpts is None else refpts

    # Python call
    py_orb4, py_orbit4 = physics.find_orbit4(py_lattice, dp, refpts)
    # Matlab call
    ml_orbit4, ml_orb4 = engine.findorbit4(ml_lattice, dp, _ml_refs(refpts, nelems), nargout=2)
    ml_orbit4 = numpy.rollaxis(numpy.asarray(ml_orbit4), -1)
    # Comparison
    numpy.testing.assert_almost_equal(py_orb4, _py_data(ml_orb4), decimal=8)
    numpy.testing.assert_almost_equal(py_orbit4[:, :4], ml_orbit4, decimal=8)


@pytest.mark.parametrize('dp', (-0.01, 0.0, 0.01))
@pytest.mark.parametrize('refpts', (0, [0, 1, 2, -1], None))
@pytest.mark.parametrize('ml_lattice, py_lattice',
                         [(pytest.lazy_fixture('ml_dba'), pytest.lazy_fixture('py_dba')),
                          (pytest.lazy_fixture('ml_hmba'), pytest.lazy_fixture('py_hmba'))])
def test_find_m44(engine, ml_lattice, py_lattice, dp, refpts):
    nelems = len(py_lattice)
    refpts = range(nelems + 1) if refpts is None else refpts
    nrefs = len(uint32_refpts(refpts, nelems))

    # Python call
    py_m44, py_mstack = physics.find_m44(py_lattice, dp, refpts)
    # Matlab call
    ml_m44, ml_mstack = engine.findm44(ml_lattice, dp, _ml_refs(refpts, nelems), nargout=2)
    ml_mstack = numpy.rollaxis(numpy.asarray(ml_mstack).reshape((4, 4, nrefs)), -1)
    # Comparison
    numpy.testing.assert_almost_equal(py_m44, numpy.asarray(ml_m44), decimal=8)
    numpy.testing.assert_almost_equal(py_mstack, ml_mstack, decimal=8)


@pytest.mark.parametrize('dp', (-0.01, 0.0, 0.01))
@pytest.mark.parametrize('refpts', (0, [0, 1, 2, -1], None))
@pytest.mark.parametrize('func_data', (('twissring', [('SPos', 's_pos'),
    ('ClosedOrbit', 'closed_orbit'), ('Dispersion', 'dispersion'),
    ('alpha', 'alpha'), ('beta', 'beta'), ('M44', 'm44')]),
    ('atlinopt', [('SPos', 's_pos'), ('ClosedOrbit', 'closed_orbit'),
    ('Dispersion', 'dispersion'), ('alpha', 'alpha'), ('beta', 'beta'),
    ('mu', 'mu'), ('M44', 'm44'), ('A', 'A'), ('B', 'B'), ('C', 'C'),
    ('gamma', 'gamma')])))
@pytest.mark.parametrize('ml_lattice, py_lattice', [(pytest.lazy_fixture('ml_hmba'),
                                                     pytest.lazy_fixture('py_hmba'))])
def test_linear_analysis(engine, ml_lattice, py_lattice, dp, refpts, func_data):
    """N.B. a 'mu' comparison is left out for twiss data as the values for 'mu'
        returned by 'twissring' in Matlab are inconsistent with those from
        'get_twiss' and 'linopt' in Python as well as those returned from
        'atlinopt' in Matlab.
    """
    nelems = len(py_lattice)
    refpts = range(nelems + 1) if refpts is None else refpts

    # Python call
    if func_data[0] == 'twissring':
        py_data0, py_tune, py_chrom, py_data = physics.get_twiss(py_lattice, dp, refpts, True)
    else:
        py_data0, py_tune, py_chrom, py_data = physics.linopt(py_lattice, dp, refpts, True)
    # Matlab call
    ml_data, ml_tune, ml_chrom = engine.pyproxy(func_data[0], ml_lattice, dp, _ml_refs(refpts, nelems), nargout=3)
    ml_data0 = engine.pyproxy(func_data[0], ml_lattice, dp, _ml_refs(nelems, nelems), nargout=3)[0]
    # Comparison
    numpy.testing.assert_almost_equal(py_tune, _py_data(ml_tune), decimal=6)
    numpy.testing.assert_almost_equal(py_chrom, _py_data(ml_chrom), decimal=4)
    _compare_physdata(numpy.expand_dims(py_data0, 0), ml_data0, func_data[1], decimal=5)
    _compare_physdata(py_data, ml_data, func_data[1], decimal=5)


@pytest.mark.parametrize('refpts', (0, [0, 1, 2, -1], None))
@pytest.mark.parametrize('ml_lattice, py_lattice', [(pytest.lazy_fixture('ml_hmba'),
                                                     pytest.lazy_fixture('py_hmba'))])
def test_ohmi_envelope(engine, ml_lattice, py_lattice, refpts):
    fields = [('beam66', 'r66'), ('beam44', 'r44'), ('emit66', 'emitXYZ'), ('emit44', 'emitXY')]
    nelems = len(py_lattice)
    refpts = range(nelems + 1) if refpts is None else refpts

    # Python call
    py_lattice.radiation_on()
    py_emit0, py_beamdata, py_emit = physics.ohmi_envelope(py_lattice, refpts)
    # Matlab call
    ml_emit = engine.pyproxy('atx', ml_lattice, 0.0, _ml_refs(refpts, nelems))
    ml_emit0, ml_params = engine.pyproxy('atx', ml_lattice, 0.0, _ml_refs(0, nelems), nargout=2)
    revolution_period = get_s_pos(py_lattice, nelems) / speed_of_light
    damping_times = revolution_period / py_beamdata.damping_rates
    # Comparison
    numpy.testing.assert_almost_equal(damping_times, _py_data(ml_params['dampingtime']), decimal=8)
    numpy.testing.assert_almost_equal(py_beamdata.mode_emittances, _py_data(ml_emit0['modemit']), decimal=8)
    _compare_physdata(numpy.expand_dims(py_emit0, 0), ml_emit0, fields)
    _compare_physdata(py_emit, ml_emit, fields)
