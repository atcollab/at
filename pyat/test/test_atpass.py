import numpy
import pytest

from at import elements, uint32_refpts
from at.tracking.atpass import atpass  # Direct access to the C function


def test_incorrect_types_raises_value_error(rin):
    with pytest.raises(TypeError):
        atpass(1, rin, 1)
    with pytest.raises(TypeError):
        atpass([], 1, 1)
    with pytest.raises(TypeError):
        atpass([], rin, "a")


def test_incorrect_dimensions_raises_value_error():
    rin = numpy.array(numpy.zeros((7, 1)))
    with pytest.raises(ValueError):
        atpass([], rin, 1)
    rin = numpy.array(numpy.zeros((1, 6)))
    with pytest.raises(ValueError):
        atpass([], rin, 1)


def test_c_aligned_array_raises_value_error():
    rin = numpy.zeros((2, 6))
    with pytest.raises(ValueError):
        atpass([], rin, 1)


def test_transposed_array_raises_value_error():
    rin_transposed = numpy.asfortranarray(numpy.zeros((2, 6)))
    with pytest.raises(ValueError):
        atpass([], rin_transposed.T, 1)


def test_transposed_c_array_does_not_raise_value_error():
    rin_c = numpy.zeros((2, 6))
    atpass([], rin_c.T, 1)


def test_transposed_c_array_gives_same_result_as_fortran_array():
    """
    This test explores the behaviour of C- and Fortran- aligned arrays.
    Since we use the data from the numpy array directly, some of the
    behaviour we see is not what we'd expect if we used numpy as intended.
    """
    # Standard numpy array (6, 2) as used in pyAT.
    rin_fortran = numpy.array(numpy.arange(12).reshape(6, 2), order="F") * 1e-5
    # Taking a copy of the transpose and making sure it is C-aligned should
    # give us the same layout of data in memory, but with a (2, 6) array.
    rin_c = numpy.copy(rin_fortran.T, order="C")
    lat = [elements.Drift("drift", 1.0)]
    rout_fortran = atpass(lat, rin_fortran, 1)
    # at.c does not accept (x, 6) arrays.  This transpose allows rin_c
    # to pass the dimension check, but does NOT change the layout in memory
    # since in Python it returns a 'view' on the array.
    # The AT C code would give the same result for rin_c without the
    # transpose if the dimension check were not there.
    rin_c_transposed = rin_c.T
    rout_c = atpass(lat, rin_c_transposed, 1)
    numpy.testing.assert_equal(rout_c, rout_fortran)


def test_missing_pass_method_raises_attribute_error(rin):
    lat = [elements.Marker("marker")]
    del lat[0].PassMethod
    with pytest.raises(AttributeError):
        atpass(lat, rin, 1)


def test_missing_integrator_raises_runtime_error(rin):
    lat = [elements.Drift("drift", 1.0, PassMethod="UnknownPass")]
    with pytest.raises(RuntimeError):
        atpass(lat, rin, 1)


def test_integrator_error(rin):
    lat = [elements.Drift("drift", 1.0)]
    del lat[0].Length
    with pytest.raises(AttributeError):
        atpass(lat, rin, 1)


@pytest.mark.parametrize("reuse", (True, False))
def test_reuse_attributes(rin, reuse):
    lat = [elements.Drift("drift", 1.0)]
    rin[0, 0] = 1e-6
    rin[1, 0] = 1e-6
    rin_copy = numpy.copy(rin)
    # two turns with original lattice
    atpass(lat, rin, 2)
    # one turn with original lattice
    atpass(lat, rin_copy, 1)
    # change an attribute
    lat[0].Length = 2
    # one turn with altered lattice
    atpass(lat, rin_copy, 1, reuse=reuse)
    if reuse:
        numpy.testing.assert_equal(rin, rin_copy)
    else:
        with pytest.raises(AssertionError):
            numpy.testing.assert_equal(rin, rin_copy)


def test_two_particles_for_two_turns():
    rin = numpy.asfortranarray(numpy.zeros((6, 2)))
    lat = [elements.Drift("drift", 1.0)]
    rin[1][0] = 1e-6
    rin[3][0] = -2e-6
    rout = atpass(lat, rin, 2, refpts=uint32_refpts([1], 1))
    # results from Matlab
    rout_particle1_turn1 = numpy.array([1e-6, 1e-6, -2e-6, -2e-6, 0, 2.5e-12]).reshape(
        6, 1
    )
    rout_particle1_turn2 = numpy.array([2e-6, 1e-6, -4e-6, -2e-6, 0, 5e-12]).reshape(
        6, 1
    )
    # the second particle doesn't change
    rout_particle2 = numpy.zeros((6, 1))
    # the second index is particle number
    numpy.testing.assert_equal(rout[:, 0, :, 0], rout_particle1_turn1)
    numpy.testing.assert_equal(rout[:, 1, :, 0], rout_particle2)
    # the fourth index is turn number
    numpy.testing.assert_equal(rout[:, 0, :, 1], rout_particle1_turn2)
    numpy.testing.assert_equal(rout[:, 1, :, 1], rout_particle2)


def test_one_particle_for_two_turns_with_no_refpts(rin):
    lat = [elements.Drift("drift", 1.0)]
    rin[1][0] = 1e-6
    rin[3][0] = -2e-6
    atpass(lat, rin, 2)
    rout_expected = numpy.array([2e-6, 1e-6, -4e-6, -2e-6, 0, 5e-12]).reshape(6, 1)
    # rin is changed in place
    numpy.testing.assert_equal(rin, rout_expected)


def test_1d_particle():
    # This is just a demonstration that if you pass no refpts you get back
    # a (6, *, 0, *) array. You may do this if you want only to operate
    # on rin in-place.
    lat = [elements.Drift("drift", 1.0)]
    rin = numpy.zeros(6)
    rin[1] = 1e-6
    # an empty refpts returns only the value at the end of the last turn
    rout = atpass(lat, rin, 1)
    # output shape: (dimensions, nparticles, refpts, nturns)
    # if no refpts are supplied, the output is (6, 1, 0, 1) and contains
    # no data
    assert rout.shape == (6, 1, 0, 1)
    rout_expected = numpy.array([1e-6, 1e-6, 0, 0, 0, 5e-13])
    # rin is changed in place
    numpy.testing.assert_equal(rin, rout_expected)
