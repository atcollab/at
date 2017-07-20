import pytest
import numpy
from at import elements, atpass


def test_incorrect_types_raises_value_error(rin):
    l = []
    with pytest.raises(ValueError):
        atpass(1, rin, 1)
    with pytest.raises(ValueError):
        atpass(l, 1, 1)
    with pytest.raises(ValueError):
        atpass(l, rin, 'a')


def test_incorrect_dimensions_raises_value_error():
    l = []
    rin = numpy.array(numpy.zeros((1,7)))
    with pytest.raises(ValueError):
        atpass(l, rin, 1)
    rin = numpy.array(numpy.zeros((6, 1)))
    with pytest.raises(ValueError):
        atpass(l, rin, 1)


def test_fortran_aligned_array_raises_value_error():
    rin = numpy.asfortranarray(numpy.zeros((2, 6)))
    l = []
    with pytest.raises(ValueError):
        atpass(l, rin, 1)


def test_transposed_array_raises_value_error():
    rin_transposed = numpy.zeros((6, 2))
    with pytest.raises(ValueError):
        atpass([], rin_transposed.T, 1)


def test_transposed_fortran_array_does_not_raise_value_error():
    rin_fortran = numpy.asfortranarray(numpy.zeros((6, 2)))
    atpass([], rin_fortran.T, 1)


def test_transposed_fortran_array_gives_same_result_as_c_array():
    """
    This test explores the behaviour of C- and Fortran- aligned arrays.
    Since we use the data from the numpy array directly, some of the
    behaviour we see does not honour what we do in numpy.
    """
    # Standard numpy array (2, 6) as used in pyAT.
    rin_c = numpy.arange(12).reshape(2, 6) * 1e-5
    # Taking an actual copy of the transpose and making sure it is
    # Fortran-aligned should give us the same layout in memory.
    rin_fortran = numpy.copy(rin_c.T, order='F')
    e = elements.Drift('drift', 1.0)
    l = [e]
    rout_c = atpass(l, rin_c, 1)
    # at.c does not currently accept (6, x) arrays.  This transpose
    # passes the dimension check, but does NOT change the layout in memory
    # since it returns a 'view' on the array.  The AT C code would give the
    # same result without the transpose.
    rin_fortran_transposed = rin_fortran.T
    rout_fortran = atpass(l, rin_fortran_transposed, 1)
    numpy.testing.assert_equal(rout_c, rout_fortran)


def test_missing_pass_method_raises_attribute_error(rin):
    m = elements.Marker('marker')
    l = [m]
    del m.PassMethod
    with pytest.raises(AttributeError):
        atpass(l, rin, 1)


def test_missing_length_raises_attribute_error(rin):
    m = elements.Drift('drift', 1.0)
    l = [m]
    del m.Length
    with pytest.raises(AttributeError):
        atpass(l, rin, 1)


@pytest.mark.parametrize("reuse", (True, False))
def test_reuse_attributes(rin, reuse):
    m = elements.Drift('drift', 1.0)
    l = [m]
    rin[0,0] = 1e-6
    rin[0,1] = 1e-6
    rin_copy = numpy.copy(rin)
    # two turns with original lattice
    atpass(l, rin, 2)
    # one turn with original lattice
    atpass(l, rin_copy, 1)
    # change an attribute
    m.Length = 2
    # one turn with altered lattice
    atpass(l, rin_copy, 1, reuse=reuse)
    if reuse:
        numpy.testing.assert_equal(rin, rin_copy)
    else:
        with pytest.raises(AssertionError):
            numpy.testing.assert_equal(rin, rin_copy)


def test_two_particles_for_two_turns():
    rin = numpy.zeros((2, 6))
    d = elements.Drift('drift', 1.0)
    lattice = [d]
    rin[0][1] = 1e-6
    rin[0][3] = -2e-6
    rout = atpass(lattice, rin, 2)
    # results from Matlab
    rout1_expected = numpy.array([1e-6, 1e-6, -2e-6, -2e-6, 0, 2.5e-12]).reshape(6)
    rout2_expected = numpy.array([2e-6, 1e-6, -4e-6, -2e-6, 0, 5e-12]).reshape(6)
    # the second particle doesn't change
    rout3_expected = numpy.zeros((6,))
    # the first two 1x6 rows are the two particles after the first turn
    numpy.testing.assert_equal(rout[0,:], rout1_expected)
    numpy.testing.assert_equal(rout[1,:], rout3_expected)
    # the second two 1x6 rows are the two particles after the second turn
    numpy.testing.assert_equal(rout[2,:], rout2_expected)
    numpy.testing.assert_equal(rout[3,:], rout3_expected)


def test_one_particle_for_two_turns_with_empty_refpts(rin):
    d = elements.Drift('drift', 1.0)
    lattice = [d]
    rin[0][1] = 1e-6
    rin[0][3] = -2e-6
    # an empty refpts returns only the value at the end of the last turn
    rout = atpass(lattice, rin, 2, refpts=numpy.zeros((0,), dtype=numpy.uint32))
    assert rout.shape == (1, 6)
    rout_expected = numpy.array([2e-6, 1e-6, -4e-6, -2e-6, 0, 5e-12]).reshape(1, 6)
    # rin is changed in place
    numpy.testing.assert_equal(rin, rout_expected)
    numpy.testing.assert_equal(rout, rout_expected)


def test_1d_particle():
    # this is just a demonstration that while you can pass a 1d particle
    # (shape (6,), you get back a 2d array (1, 6)
    d = elements.Drift('drift', 1.0)
    lattice = [d]
    rin = numpy.zeros(6,)
    rin[1] = 1e-6
    # an empty refpts returns only the value at the end of the last turn
    rout = atpass(lattice, rin, 1)
    # the input may be 1d but the output is 2d
    assert rout.shape == (1, 6)
    rout_expected = numpy.array([1e-6, 1e-6, 0, 0, 0, 5e-13])
    # rin is changed in place
    numpy.testing.assert_equal(rin, rout_expected)
    numpy.testing.assert_equal(rout, rout_expected.reshape(1, 6))

