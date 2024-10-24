import pytest
import numpy
from at import element_track, lattice_track
from at import lattice_pass, internal_lpass
from at import element_pass, internal_epass
from at import elements
from numpy.testing import assert_equal


def test_data_checks():
    val = numpy.zeros([6, 6])
    assert elements._array(val).shape == (36,)
    assert elements._array66(val).shape == (6, 6)


def test_element_string_ordering():
    d = elements.Drift('D0', 1, attr=numpy.array(0))
    assert d.__str__() == ("Drift:\n       FamName: D0\n        Length: 1.0\n"
                           "    PassMethod: DriftPass\n          attr: 0")
    assert d.__repr__() == "Drift('D0', 1.0, attr=array(0))"


def test_element_creation_raises_exception():
    with pytest.raises(ValueError):
        elements.Element('family_name', R1='not_an_array')


def test_base_element_methods():
    e = elements.Element('family_name')
    assert e.divide([0.2, 0.5, 0.3]) == [e]
    assert id(e.copy()) != id(e)


def test_argument_checks():
    q = elements.Quadrupole('quad', 1.0, 0.5)
    # Test type
    with pytest.raises(ValueError):
        q.Length = 'a'
    # Test shape
    with pytest.raises(ValueError):
        q.T1 = [0.0, 0.0]
    # Test coherence of polynoms
    with pytest.raises(ValueError):
        q.MaxOrder = 2
    with pytest.raises(ValueError):
        q.PolynomA = [0.0]
    with pytest.raises(ValueError):
        q.PolynomB = [0.0]


def test_dipole():
    d = elements.Dipole('dipole', 1.0, 0.01)
    assert d.MaxOrder == 0
    assert len(d.PolynomA) == 2
    assert d.K == 0.0
    d = elements.Dipole('dipole', 1.0, 0.01, -0.5)
    assert d.MaxOrder == 1
    assert len(d.PolynomA) == 2
    assert d.K == -0.5
    d = elements.Dipole('dipole', 1.0, 0.01, PolynomB=[0.0, 0.1, 0.0])
    assert d.MaxOrder == 1
    assert len(d.PolynomA) == 3
    assert d.K == 0.1
    d = elements.Dipole('dipole', 1.0, 0.01, PolynomB=[0.0, 0.0, 0.005])
    assert d.MaxOrder == 2
    assert len(d.PolynomA) == 3
    assert d.K == 0.0
    d = elements.Dipole('dipole', 1.0, 0.01, PolynomB=[0.0, 0.0, 0.005],
                        MaxOrder=0)
    assert d.MaxOrder == 0
    assert len(d.PolynomA) == 3
    assert d.K == 0.0


def test_quadrupole():
    q = elements.Quadrupole('quadrupole', 1.0)
    assert q.MaxOrder == 1
    assert len(q.PolynomA) == 2
    assert q.K == 0.0
    q = elements.Quadrupole('quadrupole', 1.0, -0.5)
    assert q.MaxOrder == 1
    assert len(q.PolynomA) == 2
    assert q.K == -0.5
    q = elements.Quadrupole('quadrupole', 1.0, PolynomB=[0.0, 0.0, 0.005])
    assert q.MaxOrder == 2
    assert len(q.PolynomA) == 3
    assert q.K == 0.0
    q = elements.Quadrupole('quadrupole', 1.0, PolynomB=[0.0, 0.5, 0.005],
                            MaxOrder=1)
    assert q.MaxOrder == 1
    assert len(q.PolynomA) == 3
    assert q.K == 0.5


def test_sextupole():
    s = elements.Sextupole('sextupole', 1.0)
    assert s.MaxOrder == 2
    assert len(s.PolynomA) == 3
    assert s.H == 0.0
    s = elements.Sextupole('sextupole', 1.0, -0.5)
    assert s.MaxOrder == 2
    assert len(s.PolynomA) == 3
    assert s.H == -0.5
    s = elements.Sextupole('sextupole', 1.0, PolynomB=[0.0, 0.0, 0.005, 0.0])
    assert s.MaxOrder == 2
    assert len(s.PolynomA) == 4
    assert s.H == 0.005
    s = elements.Sextupole('sextupole', 1.0, PolynomB=[0.0, 0.0, 0.005, 0.001])
    assert s.MaxOrder == 3
    assert len(s.PolynomA) == 4
    assert s.H == 0.005
    s = elements.Sextupole('sextupole', 1.0, PolynomB=[0.0, 0.5, 0.005, 0.001],
                           MaxOrder=2)
    assert s.MaxOrder == 2
    assert len(s.PolynomA) == 4
    assert s.H == 0.005


def test_octupole():
    o = elements.Octupole('octupole', 1.0, [], [0.0, 0.0, 0.0, 0.0])
    assert o.MaxOrder == 3
    assert len(o.PolynomA) == 4


def test_thinmultipole():
    m = elements.ThinMultipole('thin', [], [0.0, 0.0, 0.0, 0.0])
    assert m.MaxOrder == 0
    assert len(m.PolynomA) == 4
    m = elements.ThinMultipole('thin', [], [0.0, 0.0, 1.0, 0.0])
    assert m.MaxOrder == 2
    assert len(m.PolynomA) == 4


def test_multipole():
    m = elements.Multipole('multi', 1.0, [], [0.0, 0.0, 0.0, 0.0])
    assert m.Length == 1.0
    assert m.MaxOrder == 0
    assert m.NumIntSteps == 10
    assert m.PassMethod == 'StrMPoleSymplectic4Pass'


def test_divide_splits_attributes_correctly():
    pre = elements.Drift('drift', 1)
    post = pre.divide([0.2, 0.5, 0.3])
    assert len(post) == 3
    assert sum([e.Length for e in post]) == pre.Length
    pre = elements.Dipole('dipole', 1, KickAngle=[0.5, -0.5], BendingAngle=0.2)
    post = pre.divide([0.2, 0.5, 0.3])
    assert len(post) == 3
    assert sum([e.Length for e in post]) == pre.Length
    assert sum([e.KickAngle[0] for e in post]) == pre.KickAngle[0]
    assert sum([e.KickAngle[1] for e in post]) == pre.KickAngle[1]
    assert sum([e.BendingAngle for e in post]) == pre.BendingAngle
    pre = elements.RFCavity('rfc', 1, voltage=187500, frequency=3.5237e+8,
                            harmonic_number=31, energy=6.e+9)
    post = pre.divide([0.2, 0.5, 0.3])
    assert len(post) == 3
    assert sum([e.Length for e in post]) == pre.Length
    assert sum([e.Voltage for e in post]) == pre.Voltage


def test_insert_into_drift():
    # Create elements
    drift = elements.Drift('drift', 1)
    monitor = elements.Monitor('bpm')
    quad = elements.Quadrupole('quad', 0.3)
    # Test None splitting behaviour
    el_list = drift.insert([(0., None), (0.3, None), (0.7, None), (1., None)])
    assert len(el_list) == 3
    numpy.testing.assert_almost_equal([e.Length for e in el_list],
                                      [0.3, 0.4, 0.3])
    # Test normal insertion
    el_list = drift.insert([(0.3, monitor), (0.7, quad)])
    assert len(el_list) == 5
    numpy.testing.assert_almost_equal([e.Length for e in el_list],
                                      [0.3, 0.0, 0.25, 0.3, 0.15])
    # Test insertion at either end produces -ve length drifts
    el_list = drift.insert([(0.0, quad), (1.0, quad)])
    assert len(el_list) == 5
    numpy.testing.assert_almost_equal([e.Length for e in el_list],
                                      [-0.15, 0.3, 0.7, 0.3, -0.15])


@pytest.mark.parametrize('func', (lattice_track, lattice_pass, internal_lpass))
def test_correct_dimensions_does_not_raise_error(rin, func):
    func([], rin, 1)
    rin = numpy.zeros((6,))
    func([], rin, 1)
    rin = numpy.array(numpy.zeros((6, 2), order='F'))
    func([], rin, 1)


@pytest.mark.parametrize("dipole_class", (elements.Dipole, elements.Bend))
@pytest.mark.parametrize('func', (element_track, element_pass, internal_epass))
def test_dipole_bend_synonym(rin, dipole_class, func):
    b = dipole_class('dipole', 1.0, 0.1, EntranceAngle=0.05, ExitAngle=0.05)
    rin[0, 0] = 1e-6
    if func == element_track:
        func(b, rin, in_place=True)
    else:
        func(b, rin)
    rin_expected = numpy.array([1e-6, 0, 0, 0, 0, 1e-7]).reshape((6, 1))
    numpy.testing.assert_almost_equal(rin, rin_expected)
    assert b.K == 0.0
    b.PolynomB[1] = 0.2
    assert b.K == 0.2
    b.K = 0.1
    assert b.PolynomB[1] == 0.1


@pytest.mark.parametrize('func', (element_track, element_pass, internal_epass))
def test_marker(rin, func):
    m = elements.Marker('marker')
    assert m.Length == 0
    rin = numpy.array(numpy.random.rand(*rin.shape), order='F')
    rin_orig = numpy.array(rin, copy=True, order='F')
    if func == element_track:
        func(m, rin, in_place=True)
    else:
        func(m, rin)
    numpy.testing.assert_equal(rin, rin_orig)


@pytest.mark.parametrize('func', (element_track, element_pass, internal_epass))
def test_monitor(rin, func):
    mon = elements.Monitor('monitor')
    assert mon.Length == 0
    rin = numpy.array(numpy.random.rand(*rin.shape), order='F')
    rin_orig = rin.copy()
    if func == element_track:
        func(mon, rin, in_place=True)
    else:
        func(mon, rin)
    numpy.testing.assert_equal(rin, rin_orig)


@pytest.mark.parametrize('func', (element_track, element_pass, internal_epass))
def test_aperture_inside_limits(rin, func):
    a = elements.Aperture('aperture', [-1e-3, 1e-3, -1e-4, 1e-4])
    assert a.Length == 0
    rin[0, 0] = 1e-5
    rin[2, 0] = -1e-5
    rin_orig = rin.copy()
    if func == element_track:
        func(a, rin, in_place=True)
    else:
        func(a, rin)
    numpy.testing.assert_equal(rin, rin_orig)


@pytest.mark.parametrize('func', (lattice_track, lattice_pass, internal_lpass))
def test_aperture_outside_limits(rin, func):
    a = elements.Aperture('aperture', [-1e-3, 1e-3, -1e-4, 1e-4])
    assert a.Length == 0
    lattice = [a]
    rin[0, 0] = 1e-2
    rin[2, 0] = -1e-2
    if func == lattice_track:
        func(lattice, rin, in_place=True)
    else:
        func(lattice, rin)
    assert numpy.isnan(rin[0, 0])
    assert rin[2, 0] == 0.0  # Only the 1st coordinate is nan, the rest is zero


@pytest.mark.parametrize('func', (element_track, element_pass, internal_epass))
def test_drift_offset(rin, func):
    d = elements.Drift('drift', 1)
    rin[0, 0] = 1e-6
    rin[2, 0] = 2e-6
    rin_orig = rin.copy()
    if func == element_track:
        func(d, rin, in_place=True)
    else:
        func(d, rin)
    numpy.testing.assert_equal(rin, rin_orig)


@pytest.mark.parametrize('func', (element_track, element_pass, internal_epass))
def test_drift_divergence(rin, func):
    d = elements.Drift('drift', 1.0)
    assert d.Length == 1
    rin[1, 0] = 1e-6
    rin[3, 0] = -2e-6
    if func == element_track:
        func(d, rin, in_place=True)
    else:
        func(d, rin)
    # results from Matlab
    rin_expected = numpy.array([1e-6, 1e-6, -2e-6, -2e-6, 0,
                                2.5e-12]).reshape(6, 1)
    numpy.testing.assert_equal(rin, rin_expected)


@pytest.mark.parametrize('func', (element_track, element_pass, internal_epass))
def test_drift_two_particles(rin, func):
    d = elements.Drift('drift', 1.0)
    assert d.Length == 1
    two_rin = numpy.array(numpy.concatenate((rin, rin), axis=1), order='F')
    # particle one is offset
    two_rin[0, 0] = 1e-6
    two_rin[2, 0] = 2e-6
    # particle two has divergence
    two_rin[1, 1] = 1e-6
    two_rin[3, 1] = -2e-6
    two_rin_orig = two_rin.copy()
    if func == element_track:
        func(d, two_rin, in_place=True)
    else:
        func(d, two_rin)
    # results from Matlab
    p1_expected = numpy.array(two_rin_orig[:, 0]).reshape(6, 1)
    p2_expected = numpy.array([1e-6, 1e-6, -2e-6, -2e-6, 0,
                               2.5e-12]).reshape(6, 1)
    two_rin_expected = numpy.concatenate((p1_expected, p2_expected), axis=1)
    numpy.testing.assert_equal(two_rin, two_rin_expected)


@pytest.mark.parametrize('func', (element_track, element_pass, internal_epass))
def test_quad(rin, func):
    q = elements.Quadrupole('quad', 0.4, k=1)
    rin[0, 0] = 1e-6
    if func == element_track:
        func(q, rin, in_place=True)
    else:
        func(q, rin)
    expected = numpy.array([0.9210610203854122, -0.3894182419439, 0,
                            0, 0, 0.0000000103303797478]).reshape(6, 1) * 1e-6
    numpy.testing.assert_allclose(rin, expected)
    assert q.K == 1
    q.PolynomB[1] = 0.2
    assert q.K == 0.2
    q.K = 0.1
    assert q.PolynomB[1] == 0.1


@pytest.mark.parametrize('func', (lattice_track, lattice_pass, internal_lpass))
def test_rfcavity(rin, func):
    rf = elements.RFCavity('rfcavity', 0.0, 187500, 3.5237e+8, 31, 6.e+9)
    lattice = [rf, rf, rf, rf]
    rin[4, 0] = 1e-6
    rin[5, 0] = 1e-6
    if func == lattice_track:
        func(lattice, rin, in_place=True)
    else:
        func(lattice, rin)
    expected = numpy.array([0., 0., 0., 0., 9.990769e-7, 1.e-6]).reshape(6, 1)
    numpy.testing.assert_allclose(rin, expected, atol=1e-12)


@pytest.mark.parametrize('func', (element_track, element_pass, internal_epass))
@pytest.mark.parametrize("n", (0, 1, 2, 3, 4, 5))
def test_m66(rin, n, func):
    m = numpy.random.rand(6, 6)
    m66 = elements.M66('m66', m)
    assert m66.Length == 0
    rin[n, 0] = 1e-6
    if func == element_track:
        func(m66, rin, in_place=True)
    else:
        func(m66, rin)
    expected = numpy.array([m[0, n], m[1, n], m[2, n], m[3, n], m[4, n],
                            m[5, n]]).reshape(6, 1) * 1e-6
    numpy.testing.assert_equal(rin, expected)


@pytest.mark.parametrize('func', (element_track, element_pass, internal_epass))
def test_corrector(rin, func):
    c = elements.Corrector('corrector', 0.0, numpy.array([0.9, 0.5],
                                                         dtype=numpy.float64))
    assert c.Length == 0
    rin[0, 0] = 1e-6
    rin_orig = rin.copy()
    rin_orig[1] = 0.9
    rin_orig[3] = 0.5
    if func == element_track:
        func(c, rin, in_place=True)
    else:
        func(c, rin)
    numpy.testing.assert_equal(rin, rin_orig)


@pytest.mark.parametrize('func', (element_track, element_pass, internal_epass))
def test_wiggler(rin, func):
    period = 0.05
    periods = 23
    bmax = 1
    by = numpy.array([1, 1, 0, 1, 1, 0], dtype=numpy.float64)
    c = elements.Wiggler('wiggler', period * periods, period, bmax, 3e9, By=by)
    assert abs(c.Length - 1.15) < 1e-10
    # Expected value from Matlab AT.
    expected = numpy.array(rin, copy=True)
    expected[5] = 0.000000181809691064259
    if func == element_track:
        func(c, rin, energy=3e9, in_place=True)
    else:
        func(c, rin, energy=3e9)
    numpy.testing.assert_allclose(rin, expected, atol=1e-12)


def test_exit_entrance():
    q = elements.Quadrupole('quad', 0.4, k=1)
    for kin, kout in zip(q._entrance_fields, q._exit_fields):
        assert_equal(kin.replace('Entrance', ''). replace('1', ''),
                     kout.replace('Exit', '').replace('2', ''))
