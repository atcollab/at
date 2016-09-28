import pytest
import numpy
import at
import elements


@pytest.fixture
def rin():
    rin = numpy.array(numpy.zeros((1,6)))
    return rin


def test_correct_dimensions_does_not_raise_error(rin):
    l = []
    at.atpass(l, rin, 1)
    rin = numpy.zeros((6,))
    at.atpass(l, rin, 1)
    rin = numpy.zeros((2,6))


def test_incorrect_types_raises_value_error(rin):
    l = []
    with pytest.raises(ValueError):
        at.atpass(1, rin, 1)
    with pytest.raises(ValueError):
        at.atpass(l, 1, 1)
    with pytest.raises(ValueError):
        at.atpass(l, rin, 'a')


def test_incorrect_dimensions_raises_value_error():
    l = []
    rin = numpy.array(numpy.zeros((1,7)))
    with pytest.raises(ValueError):
        at.atpass(l, rin, 1)
    rin = numpy.array(numpy.zeros((6,1)))
    with pytest.raises(ValueError):
        at.atpass(l, rin, 1)


def test_fortran_aligned_array_raises_value_error():
    rin = numpy.asfortranarray(numpy.zeros((2,6)))
    l = []
    with pytest.raises(ValueError):
        at.atpass(l, rin, 1)


def test_missing_pass_method_raises_attribute_error(rin):
    m = elements.Marker('marker')
    l = [m]
    del m.PassMethod
    with pytest.raises(AttributeError):
        at.atpass(l, rin, 1)


def test_missing_length_raises_attribute_error(rin):
    m = elements.Drift('drift', 1.0)
    l = [m]
    del m.Length
    with pytest.raises(AttributeError):
        at.atpass(l, rin, 1)


def test_dipole(rin):
    print(elements.__file__)
    b = elements.Dipole('dipole', 1.0, 0.1, EntranceAngle=0.05, ExitAngle=0.05)
    l = [b]
    rin[0,0] = 1e-6
    rin_orig = numpy.copy(rin)
    at.atpass(l, rin, 1)
    rin_expected = numpy.array([1e-6, 0, 0, 0, 0, 1e-7]).reshape((1,6))
    numpy.testing.assert_almost_equal(rin_orig, rin_expected)


def test_marker(rin):
    m = elements.Marker('marker')
    assert m.Length == 0
    lattice = [m]
    rin = numpy.random.rand(*rin.shape)
    rin_orig = numpy.array(rin, copy=True)
    at.atpass(lattice, rin, 1)
    numpy.testing.assert_equal(rin, rin_orig)


def test_aperture_inside_limits(rin):
    a = elements.Aperture('aperture', [-1e-3, 1e-3, -1e-4, 1e-4])
    assert a.Length == 0
    lattice = [a]
    rin[0][0] = 1e-5
    rin[0][2] = -1e-5
    rin_orig = numpy.array(rin, copy=True)
    at.atpass(lattice, rin, 1)
    numpy.testing.assert_equal(rin, rin_orig)


def test_aperture_outside_limits(rin):
    a = elements.Aperture('aperture', [-1e-3, 1e-3, -1e-4, 1e-4])
    assert a.Length == 0
    lattice = [a]
    rin[0][0] = 1e-2
    rin[0][2] = -1e-2
    at.atpass(lattice, rin, 1)
    assert numpy.isinf(rin[0][0])
    assert rin[0][2] == -1e-2  # Only the first coordinate is marked as infinity


def test_drift_offset(rin):
    d = elements.Drift('drift', 1)
    lattice = [d]
    rin[0][0] = 1e-6
    rin[0][2] = 2e-6
    rin_orig = numpy.array(rin, copy=True)
    at.atpass(lattice, rin, 1)
    numpy.testing.assert_equal(rin, rin_orig)


def test_drift_divergence(rin):
    d = elements.Drift('drift', 1.0)
    assert d.Length == 1
    lattice = [d]
    rin[0][1] = 1e-6
    rin[0][3] = -2e-6
    at.atpass(lattice, rin, 1)
    # results from Matlab
    rin_expected = numpy.array([1e-6, 1e-6, -2e-6, -2e-6, 0, 2.5e-12]).reshape(1,6)
    numpy.testing.assert_equal(rin, rin_expected)


def test_drift_two_particles(rin):
    d = elements.Drift('drift', 1.0)
    assert d.Length == 1
    lattice = [d]
    two_rin = numpy.concatenate((rin, rin), axis=0)
    # particle one is offset
    two_rin[0][0] = 1e-6
    two_rin[0][2] = 2e-6
    # particle two has divergence
    two_rin[1][1] = 1e-6
    two_rin[1][3] = -2e-6
    two_rin_orig = numpy.array(two_rin, copy=True)
    at.atpass(lattice, two_rin, 1)
    # results from Matlab
    p1_expected = numpy.array(two_rin_orig[0,:]).reshape(1,6)
    p2_expected = numpy.array([1e-6, 1e-6, -2e-6, -2e-6, 0, 2.5e-12]).reshape(1,6)
    two_rin_expected = numpy.concatenate((p1_expected, p2_expected), axis=0)
    numpy.testing.assert_equal(two_rin, two_rin_expected)
