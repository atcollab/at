import numpy
from at import lattice, elements
import pytest


@pytest.mark.parametrize('ref_in, expected', (
    [2, numpy.array([2], dtype=numpy.uint32)],
    [-1, numpy.array([4], dtype=numpy.uint32)],
    [7, numpy.array([2], dtype=numpy.uint32)],
    [[1, 2, 3], numpy.array([1, 2, 3], dtype=numpy.uint32)],
    [[0, 1, 2, 3, 4, 5], numpy.array([0, 1, 2, 3, 4, 5], dtype=numpy.uint32)],
    [[0, 6, 2, -2, 4, 5], numpy.array([0, 1, 2, 3, 4, 5], dtype=numpy.uint32)],
))
def test_uint32_refpts_converts_numerical_inputs_correctly(ref_in, expected):
    numpy.testing.assert_equal(lattice.uint32_refpts(ref_in, 5), expected)
    ref_in2 = numpy.asarray(ref_in).astype(numpy.float64)
    numpy.testing.assert_equal(lattice.uint32_refpts(ref_in2, 5), expected)
    ref_in2 = numpy.asarray(ref_in).astype(numpy.int64)
    numpy.testing.assert_equal(lattice.uint32_refpts(ref_in2, 5), expected)


@pytest.mark.parametrize('ref_in, expected', (
    [None, numpy.array([], dtype=numpy.uint32)],
    [[], numpy.array([], dtype=numpy.uint32)],
    [True, numpy.array([0], dtype=numpy.uint32)],
    [False, numpy.array([], dtype=numpy.uint32)],
    [[True, False, True], numpy.array([0, 2], dtype=numpy.uint32)],
    [[False, False, False, False, False, False], numpy.array([], dtype=numpy.uint32)],
    [[True, True, True, True, True, True], numpy.array([0, 1, 2, 3, 4, 5], dtype=numpy.uint32)]
))
def test_uint32_refpts_converts_other_input_types_correctly(ref_in, expected):
    numpy.testing.assert_equal(lattice.uint32_refpts(ref_in, 5), expected)


# too long, misordered, duplicate, -ve indexing misordered, -ve indexing
# duplicate, wrap-around indexing misordered, wrap-around indexing duplicate,
# -ve indexing and wrap-around indexing duplicate, -ve indexing and wrap-around
# indexing misordered
@pytest.mark.parametrize('ref_in', ([0, 1, 2, 3], [2, 1], [0, 0], [-1, 0],
                                    [0, -2], [3, 0], [1, 3], [-1, 3], [3, -2]))
def test_uint32_refpts_throws_ValueError_correctly(ref_in):
    with pytest.raises(ValueError):
        r = lattice.uint32_refpts(ref_in, 2)


def test_get_s_pos_returns_zero_for_empty_lattice():
    numpy.testing.assert_equal(lattice.get_s_pos([], None), numpy.array((0,)))


def test_get_s_pos_returns_length_for_lattice_with_one_element():
    e = elements.Element('e', 0.1)
    assert lattice.get_s_pos([e], [1]) == numpy.array([0.1])


def test_get_s_pos_returns_all_points_for_lattice_with_two_elements_and_refpts_None():
    e = elements.Element('e', 0.1)
    f = elements.Element('e', 2.1)
    print(lattice.get_s_pos([e, f], None))
    numpy.testing.assert_equal(lattice.get_s_pos([e, f], None),
                               numpy.array([0, 0.1, 2.2]))


def test_get_s_pos_returns_all_points_for_lattice_with_two_elements_using_int_refpts():
    e = elements.Element('e', 0.1)
    f = elements.Element('e', 2.1)
    lat = [e, f]
    numpy.testing.assert_equal(lattice.get_s_pos(lat, range(len(lat) + 1)),
                               numpy.array([0, 0.1, 2.2]))

def test_get_s_pos_returns_all_points_for_lattice_with_two_elements_using_bool_refpts():
    e = elements.Element('e', 0.1)
    f = elements.Element('e', 2.1)
    lat = [e, f]
    numpy.testing.assert_equal(lattice.get_s_pos(lat, numpy.array((True, True, True))),
                               numpy.array([0, 0.1, 2.2]))
