import numpy
from at import lattice, elements
import pytest


def test_uint32_refpts_handles_array_as_input():
    expected = numpy.array([1, 2, 3], dtype=numpy.uint32)
    requested = numpy.array([1, 2, 3], dtype=numpy.uint32)
    numpy.testing.assert_equal(lattice.uint32_refpts(requested, 5),
                               expected)
    requested = numpy.array([1, 2, 3], dtype=numpy.float64)
    numpy.testing.assert_equal(lattice.uint32_refpts(requested, 5),
                               expected)


@pytest.mark.parametrize('input', ([-1, 0], [0, 0], [2, 1], [0, 1, 2, 3]))
def test_uint32_refpts_throws_ValueError_if_input_invalid(input):
    with pytest.raises(ValueError):
        r = lattice.uint32_refpts(input, 2)


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
