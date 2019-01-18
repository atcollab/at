import pytest
import numpy
from at import atpass
from at import elements


def test_data_checks():
    val = numpy.zeros([6,6])
    assert elements._array(val).shape == (36,)
    assert elements._array66(val).shape == (6, 6)


def test_element_string_ordering():
    d = elements.Drift('D0', 1, attr=numpy.array(0))
    assert d.__str__() == ("Drift:\nFamName : D0\nLength : 1.0\nPassMethod :"
                           " DriftPass\nattr : 0")
    assert d.__repr__() == "Drift('D0', 1.0, attr=array(0))"


def test_element_creation_raises_exception():
    with pytest.raises(ValueError):
        elements.Element('family_name', R1='not_an_array')


def test_correct_dimensions_does_not_raise_error(rin):
    l = []
    atpass(l, rin, 1)
    rin = numpy.zeros((6,))
    atpass(l, rin, 1)
    rin = numpy.array(numpy.zeros((6, 2), order='F'))
    atpass(l, rin, 1)


@pytest.mark.parametrize("dipole_class", (elements.Dipole, elements.Bend))
def test_dipole(rin, dipole_class):
    b = dipole_class('dipole', 1.0, 0.1, EntranceAngle=0.05, ExitAngle=0.05)
    l = [b]
    rin[0, 0] = 1e-6
    rin_orig = numpy.copy(rin)
    atpass(l, rin, 1)
    rin_expected = numpy.array([1e-6, 0, 0, 0, 0, 1e-7]).reshape((6, 1))
    numpy.testing.assert_almost_equal(rin_orig, rin_expected)


def test_marker(rin):
    m = elements.Marker('marker')
    assert m.Length == 0
    lattice = [m]
    rin = numpy.array(numpy.random.rand(*rin.shape), order='F')
    rin_orig = numpy.array(rin, copy=True, order='F')
    atpass(lattice, rin, 1)
    numpy.testing.assert_equal(rin, rin_orig)


def test_monitor(rin):
    mon = elements.Monitor('monitor')
    assert mon.Length == 0
    lattice = [mon]
    rin = numpy.array(numpy.random.rand(*rin.shape), order='F')
    rin_orig = numpy.array(rin, copy=True, order='F')
    atpass(lattice, rin, 1)
    numpy.testing.assert_equal(rin, rin_orig)


def test_aperture_inside_limits(rin):
    a = elements.Aperture('aperture', [-1e-3, 1e-3, -1e-4, 1e-4])
    assert a.Length == 0
    lattice = [a]
    rin[0, 0] = 1e-5
    rin[2, 0] = -1e-5
    rin_orig = numpy.array(rin, copy=True)
    atpass(lattice, rin, 1)
    numpy.testing.assert_equal(rin, rin_orig)


def test_aperture_outside_limits(rin):
    a = elements.Aperture('aperture', [-1e-3, 1e-3, -1e-4, 1e-4])
    assert a.Length == 0
    lattice = [a]
    rin[0, 0] = 1e-2
    rin[2, 0] = -1e-2
    atpass(lattice, rin, 1)
    assert numpy.isinf(rin[5, 0])
    assert rin[2, 0] == -1e-2  # Only the 6th coordinate is marked as infinity


def test_drift_offset(rin):
    d = elements.Drift('drift', 1)
    lattice = [d]
    rin[0, 0] = 1e-6
    rin[2, 0] = 2e-6
    rin_orig = numpy.array(rin, copy=True)
    atpass(lattice, rin, 1)
    numpy.testing.assert_equal(rin, rin_orig)


def test_drift_divergence(rin):
    d = elements.Drift('drift', 1.0)
    assert d.Length == 1
    lattice = [d]
    rin[1, 0] = 1e-6
    rin[3, 0] = -2e-6
    atpass(lattice, rin, 1)
    # results from Matlab
    rin_expected = numpy.array([1e-6, 1e-6, -2e-6, -2e-6, 0, 2.5e-12]).reshape(6, 1)
    numpy.testing.assert_equal(rin, rin_expected)


def test_drift_two_particles(rin):
    d = elements.Drift('drift', 1.0)
    assert d.Length == 1
    lattice = [d]
    two_rin = numpy.array(numpy.concatenate((rin, rin), axis=1), order='F')
    # particle one is offset
    two_rin[0, 0] = 1e-6
    two_rin[2, 0] = 2e-6
    # particle two has divergence
    two_rin[1, 1] = 1e-6
    two_rin[3, 1] = -2e-6
    two_rin_orig = numpy.array(two_rin, copy=True)
    atpass(lattice, two_rin, 1)
    # results from Matlab
    p1_expected = numpy.array(two_rin_orig[:, 0]).reshape(6, 1)
    p2_expected = numpy.array([1e-6, 1e-6, -2e-6, -2e-6, 0, 2.5e-12]).reshape(6, 1)
    two_rin_expected = numpy.concatenate((p1_expected, p2_expected), axis=1)
    numpy.testing.assert_equal(two_rin, two_rin_expected)


def test_quad(rin):
    q = elements.Quadrupole('quad', 0.4, k=1)
    lattice = [q]
    rin[0, 0] = 1e-6
    atpass(lattice, rin, 1)
    expected = numpy.array([0.921060994002885, -0.389418342308651, 0,
                            0, 0, 0.000000010330489]).reshape(6, 1) * 1e-6
    numpy.testing.assert_allclose(rin, expected)


def test_quad_incorrect_array(rin):
    q = elements.Quadrupole('quad', 0.4, k=1)
    q.PolynomB = 'a'
    lattice = [q]
    with pytest.raises(RuntimeError):
        atpass(lattice, rin, 1)


def test_rfcavity(rin):
    rf = elements.RFCavity('rfcavity', 0.0, 187500, 3.5237e+8, 31, 6.e+9)
    lattice = [rf, rf, rf, rf]
    rin[4, 0] = 1e-6
    rin[5, 0] = 1e-6
    atpass(lattice, rin, 1)
    expected = numpy.array([0., 0., 0., 0., 9.990769e-7, 1.e-6]).reshape(6, 1)
    numpy.testing.assert_allclose(rin, expected, atol=1e-12)


def test_ringparam(rin):
    rp = elements.RingParam('ringparam', 1.e+09)
    assert rp.Length == 0
    lattice = [rp]
    rin = numpy.array(numpy.random.rand(*rin.shape), order='F')
    rin_orig = numpy.array(rin, copy=True, order='F')
    atpass(lattice, rin, 1)
    numpy.testing.assert_equal(rin, rin_orig)

@pytest.mark.parametrize("n", (0, 1, 2, 3, 4, 5))
def test_m66(rin, n):
    m = numpy.random.rand(6, 6)
    m66 = elements.M66('m66', m)
    assert m66.Length == 0
    lattice = [m66]
    rin[n, 0] = 1e-6
    atpass(lattice, rin, 1)
    expected = numpy.array([m[0, n], m[1, n], m[2, n], m[3, n], m[4, n],
                            m[5, n]]).reshape(6, 1) * 1e-6
    numpy.testing.assert_equal(rin, expected)


def test_corrector(rin):
    c = elements.Corrector('corrector', 0.0, numpy.array([0.9, 0.5], dtype=numpy.float64))
    assert c.Length == 0
    lattice = [c]
    rin[0, 0] = 1e-6
    rin_orig = numpy.array(rin, copy=True)
    rin_orig[1] = 0.9
    rin_orig[3] = 0.5
    atpass(lattice, rin, 1)
    numpy.testing.assert_equal(rin, rin_orig)
