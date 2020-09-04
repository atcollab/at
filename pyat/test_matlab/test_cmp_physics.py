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
    print(fields)
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
                         [(pytest.lazy_fixture('ml_dba'),
                           pytest.lazy_fixture('py_dba')),
                          (pytest.lazy_fixture('ml_hmba'),
                           pytest.lazy_fixture('py_hmba'))])
def test_find_orbit4(engine, ml_lattice, py_lattice, dp, refpts):
    nelems = len(py_lattice)
    refpts = range(nelems + 1) if refpts is None else refpts

    # Python call
    py_orb4, py_orbit4 = physics.find_orbit4(py_lattice, dp, refpts)
    # Matlab call
    ml_orbit4, ml_orb4 = engine.findorbit4(ml_lattice, dp,
                                           _ml_refs(refpts, nelems), nargout=2)
    ml_orbit4 = numpy.rollaxis(numpy.asarray(ml_orbit4), -1)
    # Comparison
    numpy.testing.assert_almost_equal(py_orb4, _py_data(ml_orb4), decimal=8)
    numpy.testing.assert_almost_equal(py_orbit4[:, :4], ml_orbit4, decimal=8)


@pytest.mark.parametrize('dp', (-0.01, 0.0, 0.01))
@pytest.mark.parametrize('refpts', (0, [0, 1, 2, -1], None))
@pytest.mark.parametrize('ml_lattice, py_lattice',
                         [(pytest.lazy_fixture('ml_dba'),
                           pytest.lazy_fixture('py_dba')),
                          (pytest.lazy_fixture('ml_hmba'),
                           pytest.lazy_fixture('py_hmba'))])
def test_find_m44(engine, ml_lattice, py_lattice, dp, refpts):
    nelems = len(py_lattice)
    refpts = range(nelems + 1) if refpts is None else refpts
    nrefs = len(uint32_refpts(refpts, nelems))

    # Python call
    py_m44, py_mstack = physics.find_m44(py_lattice, dp, refpts)
    # Matlab call
    ml_m44, ml_mstack = engine.findm44(ml_lattice, dp,
                                       _ml_refs(refpts, nelems), nargout=2)
    ml_mstack = numpy.rollaxis(numpy.asarray(ml_mstack).reshape((4, 4,
                                                                 nrefs)), -1)
    # Comparison
    numpy.testing.assert_almost_equal(py_m44, numpy.asarray(ml_m44), decimal=8)
    numpy.testing.assert_almost_equal(py_mstack, ml_mstack, decimal=8)


@pytest.mark.parametrize('dp', (-0.01, 0.0, 0.01))
@pytest.mark.parametrize('refpts', (0, [0, 1, 2, -1], None))
@pytest.mark.parametrize('ml_lattice, py_lattice',
                         [(pytest.lazy_fixture('ml_hmba'),
                           pytest.lazy_fixture('py_hmba'))])
def test_linear_analysis(engine, ml_lattice, py_lattice, dp, refpts):
    nelems = len(py_lattice)
    fields = [('SPos', 's_pos'),
              ('ClosedOrbit', 'closed_orbit'),
              ('Dispersion', 'dispersion'),
              ('alpha', 'alpha'), ('beta', 'beta'),
              ('mu', 'mu'), ('M44', 'm44'),
              ('A', 'A'), ('B', 'B'), ('C', 'C'),
              ('gamma', 'gamma')]
    refpts = range(nelems + 1) if refpts is None else refpts
    py_data0, py_tune, py_chrom, py_data = physics.linopt(py_lattice, dp,
                                                          refpts, True,
                                                          ddp=1.E-6)
    # Matlab call
    ml_data, ml_tune, ml_chrom = engine.pyproxy('atlinopt', ml_lattice, dp,
                                                _ml_refs(refpts, nelems),
                                                nargout=3)
    ml_data0 = engine.pyproxy('atlinopt', ml_lattice, dp,
                              _ml_refs(nelems, nelems), nargout=3)[0]
    # Comparison
    numpy.testing.assert_almost_equal(py_tune, _py_data(ml_tune), decimal=6)
    numpy.testing.assert_almost_equal(py_chrom, _py_data(ml_chrom), decimal=4)
    _compare_physdata(numpy.expand_dims(py_data0, 0), ml_data0, fields,
                      decimal=5)
    _compare_physdata(py_data, ml_data, fields, decimal=6)


@pytest.mark.parametrize('refpts', (0, [0, 1, 2, -1], None))
@pytest.mark.parametrize('ml_lattice, py_lattice',
                         [(pytest.lazy_fixture('ml_hmba'),
                           pytest.lazy_fixture('py_hmba'))])
def test_ohmi_envelope(engine, ml_lattice, py_lattice, refpts):
    fields = [('beam66', 'r66'), ('beam44', 'r44'), ('emit66', 'emitXYZ'),
              ('emit44', 'emitXY')]
    nelems = len(py_lattice)
    refpts = range(nelems + 1) if refpts is None else refpts

    # Python call
    py_lattice = py_lattice.radiation_on(copy=True)
    py_emit0, py_beamdata, py_emit = py_lattice.ohmi_envelope(refpts)
    # Matlab call
    ml_emit = engine.pyproxy('atx', ml_lattice, 0.0, _ml_refs(refpts, nelems))
    ml_emit0, ml_params = engine.pyproxy('atx', ml_lattice, 0.0,
                                         _ml_refs(0, nelems), nargout=2)
    revolution_period = get_s_pos(py_lattice, nelems) / speed_of_light
    damping_times = revolution_period / py_beamdata.damping_rates
    # Comparison
    numpy.testing.assert_almost_equal(damping_times,
                                      _py_data(ml_params['dampingtime']),
                                      decimal=8)
    numpy.testing.assert_almost_equal(py_beamdata.mode_emittances,
                                      _py_data(ml_emit0['modemit']), decimal=8)
    _compare_physdata(numpy.expand_dims(py_emit0, 0), ml_emit0, fields)
    _compare_physdata(py_emit, ml_emit, fields)


@pytest.mark.parametrize('ml_lattice, py_lattice',
                         [(pytest.lazy_fixture('ml_hmba'),
                           pytest.lazy_fixture('py_hmba'))])
def test_quantdiff(engine, ml_lattice, py_lattice):
    py_lattice = py_lattice.radiation_on(copy=True)
    # Python call
    dmat = physics.radiation.quantdiffmat(py_lattice)
    lmat = physics.radiation._lmat(dmat)
    # Matlab call
    radring, ind = engine.pyproxy('atradon', ml_lattice, nargout=2)
    dmat_ml = engine.pyproxy('quantumDiff', radring, ind)
    # Comparison
    numpy.testing.assert_allclose(dmat, dmat_ml,
                                  rtol=1.0e-8, atol=1.0e-20)


@pytest.mark.parametrize('ml_lattice, py_lattice',
                         [(pytest.lazy_fixture('ml_hmba'),
                           pytest.lazy_fixture('py_hmba'))])
def test_fastring(engine, ml_lattice, py_lattice):
    ml_lattice[0]['Periodicity'] = 1.0
    # Python call
    ring, ringrad = physics.fastring.fast_ring(py_lattice)
    # Matlab call
    ring_ml, ringrad_ml = engine.pyproxy('atfastring',
                                         ml_lattice, nargout=2)
    # Comparison
    for r, rml, idq in zip([ring, ringrad], [ring_ml, ringrad_ml], [3, 4]):
        numpy.testing.assert_allclose(r[0].Frequency,
                                      rml[1]['Frequency'], rtol=1.0e-20)
        numpy.testing.assert_allclose(r[0].Voltage,
                                      rml[1]['Voltage'], rtol=1.0e-20)
        numpy.testing.assert_allclose(r[1].I2,
                                      rml[2]['I2'], rtol=1.0e-20)
        numpy.testing.assert_allclose(r[1].Length,
                                      rml[2]['Length'], rtol=1.0e-20)
        numpy.testing.assert_allclose(r[1].M66,
                                      rml[2]['M66'], atol=1.0e-7)
        numpy.testing.assert_allclose(r[1].T1,
                                      numpy.squeeze(rml[2]['T1']),
                                      rtol=1.0e-8, atol=1.0e-11)
        numpy.testing.assert_allclose(r[1].T2,
                                      numpy.squeeze(rml[2]['T2']),
                                      rtol=1.0e-8, atol=1.0e-11)
        numpy.testing.assert_allclose(r[idq-1].A1,
                                      rml[idq]['A1'], rtol=0.02)
        numpy.testing.assert_allclose(r[idq-1].A2,
                                      rml[idq]['A2'], rtol=0.02)
        numpy.testing.assert_allclose(r[idq-1].A3,
                                      rml[idq]['A3'], rtol=0.02)
        numpy.testing.assert_allclose(r[idq-1].Alphax,
                                      rml[idq]['Alphax'], rtol=1.0e-15)
        numpy.testing.assert_allclose(r[idq-1].Alphay,
                                      rml[idq]['Alphay'], rtol=1.0e-15)
        numpy.testing.assert_allclose(r[idq-1].Betax,
                                      rml[idq]['Betax'], rtol=1.0e-15)
        numpy.testing.assert_allclose(r[idq-1].Betay,
                                      rml[idq]['Betay'], rtol=1.0e-15)
        numpy.testing.assert_allclose(r[idq-1].Qpx,
                                      rml[idq]['Qpx'], rtol=1.0e-5)
        numpy.testing.assert_allclose(r[idq-1].Qpy,
                                      rml[idq]['Qpy'], rtol=1.0e-5)
        numpy.testing.assert_allclose(r[idq-1].T1,
                                      numpy.squeeze(rml[idq]['T1']),
                                      rtol=1.0e-8, atol=1.0e-11)
        numpy.testing.assert_allclose(r[idq-1].T2,
                                      numpy.squeeze(rml[idq]['T2']),
                                      rtol=1.0e-8, atol=1.0e-11)
        if idq == 4:
            numpy.testing.assert_allclose(r[2].Lmatp,
                                          rml[3]['Lmatp'],
                                          rtol=0.02, atol=1.0e-10)
    ml_lattice[0]['Periodicity'] = 32.0


@pytest.mark.parametrize('dp', (0.00, 0.01, -0.01))
@pytest.mark.parametrize('ml_lattice, py_lattice',
                         [(pytest.lazy_fixture('ml_hmba'),
                           pytest.lazy_fixture('py_hmba'))])
def test_parameters(engine, ml_lattice, py_lattice, dp):

    # Test perimeter
    py_length = py_lattice.get_s_pos(len(py_lattice))
    ml_length = engine.findspos(ml_lattice, len(ml_lattice)+1)
    numpy.testing.assert_allclose(py_length, ml_length, rtol=1.E-8)

    # test energy loss
    ml_energy, ml_periods, ml_voltage, \
        ml_harms, ml_eloss = engine.pyproxy('atenergy', ml_lattice, nargout=5)
    numpy.testing.assert_allclose(py_lattice.voltage, ml_voltage, rtol=1.E-8)
    numpy.testing.assert_allclose(py_lattice.energy_loss, ml_eloss, rtol=1.E-6)
    assert py_lattice.energy == ml_energy
    assert py_lattice.periodicity == int(ml_periods)
    assert py_lattice.harmonic_number == int(ml_harms)

    # test momentum compaction factor
    py_mcf = py_lattice.get_mcf(dp, ddp=1.E-6)  # Matlab uses ddp=1.E-6
    ml_mcf = engine.mcf(ml_lattice, dp)
    numpy.testing.assert_allclose(py_mcf, ml_mcf, rtol=1.E-8)
