import numpy
import pytest

import math
from at.lattice import elements as elt
from at.lattice import checktype, get_cells, refpts_iterator, get_elements
from at.lattice import elements, uint32_refpts, bool_refpts, checkattr
from at.lattice import get_s_pos, tilt_elem, shift_elem, set_tilt, set_shift


@pytest.mark.parametrize(
    "ref_in, expected",
    (
        [2, numpy.array([2], dtype=numpy.uint32)],
        [-1, numpy.array([4], dtype=numpy.uint32)],
        [7, numpy.array([2], dtype=numpy.uint32)],
        [[1, 2, 3], numpy.array([1, 2, 3], dtype=numpy.uint32)],
        [[0, 1, 2, 3, 4, 5], numpy.array([0, 1, 2, 3, 4, 5], dtype=numpy.uint32)],
        [[0, 6, 2, -2, 4, 5], numpy.array([0, 1, 2, 3, 4, 5], dtype=numpy.uint32)],
    ),
)
def test_uint32_refpts_converts_numerical_inputs_correctly(ref_in, expected):
    numpy.testing.assert_equal(uint32_refpts(ref_in, 5), expected)
    # float indexes are deprecated
    # ref_in2 = numpy.asarray(ref_in).astype(numpy.float64)
    # numpy.testing.assert_equal(uint32_refpts(ref_in2, 5), expected)
    ref_in2 = numpy.asarray(ref_in).astype(numpy.int64)
    numpy.testing.assert_equal(uint32_refpts(ref_in2, 5), expected)


@pytest.mark.parametrize(
    "ref_in, expected",
    (
        [None, numpy.array([], dtype=numpy.uint32)],
        [[], numpy.array([], dtype=numpy.uint32)],
        [True, numpy.array([0], dtype=numpy.uint32)],
        [False, numpy.array([], dtype=numpy.uint32)],
        [[True, False, True], numpy.array([0, 2], dtype=numpy.uint32)],
        [
            [False, False, False, False, False, False],
            numpy.array([], dtype=numpy.uint32),
        ],
        [
            [True, True, True, True, True, True],
            numpy.array([0, 1, 2, 3, 4, 5], dtype=numpy.uint32),
        ],
    ),
)
def test_uint32_refpts_converts_other_input_types_correctly(ref_in, expected):
    numpy.testing.assert_equal(uint32_refpts(ref_in, 5), expected)


# too long, misordered, duplicate, -ve indexing misordered, -ve indexing
# duplicate, wrap-around indexing misordered, wrap-around indexing duplicate,
# -ve indexing and wrap-around indexing duplicate, -ve indexing and wrap-around
# indexing misordered
@pytest.mark.parametrize(
    "ref_in",
    ([0, 1, 2, 3], [2, 1], [0, 0], [-1, 0], [0, -2], [3, 0], [1, 3], [-1, 3], [3, -2]),
)
def test_uint32_refpts_throws_IndexError_correctly(ref_in):
    with pytest.raises(IndexError):
        uint32_refpts(ref_in, 2)


@pytest.mark.parametrize(
    "ref_in, expected",
    (
        [[True], numpy.array([True, False, False, False], dtype=bool)],
        [[], numpy.array([False, False, False, False], dtype=bool)],
        [None, numpy.array([False, False, False, False], dtype=bool)],
        [2, numpy.array([False, False, True, False], dtype=bool)],
        [[1, 2], numpy.array([False, True, True, False], dtype=bool)],
    ),
)
def test_bool_refpts(ref_in, expected):
    numpy.testing.assert_equal(bool_refpts(ref_in, 3), expected)


@pytest.mark.parametrize("ref_in", (10, range(100)))
def test_bool_refpts_throws_IndexErrorrs(ref_in):
    with pytest.raises(IndexError):
        bool_refpts(ref_in, 3)


def test_checkattr(simple_ring):
    assert checkattr("Length")(simple_ring[0]) is True
    assert checkattr("not_an_attr")(simple_ring[0]) is False
    assert list(filter(checkattr("Length", 1), simple_ring)) == [
        simple_ring[0],
        simple_ring[3],
        simple_ring[4],
        simple_ring[5],
    ]
    assert list(filter(checkattr("Length", 2), simple_ring)) == []
    assert list(filter(checkattr("not_an_attr"), simple_ring)) == []


def test_checktype(simple_ring):
    assert checktype(elements.Drift)(simple_ring[0]) is True
    assert checktype(elements.Marker)(simple_ring[0]) is False
    assert list(filter(checktype(elements.Drift), simple_ring)) == [
        simple_ring[0],
        simple_ring[3],
        simple_ring[4],
        simple_ring[5],
    ]
    assert list(filter(checktype(elements.Monitor), simple_ring)) == []


def test_get_cells(simple_ring):
    a = numpy.array([True, True, True, True, True, True, False])
    numpy.testing.assert_equal(get_cells(simple_ring, checkattr("Length")), a)
    a = numpy.array([False, True, False, False, False, False, False])
    numpy.testing.assert_equal(get_cells(simple_ring, "attr"), a)
    a = numpy.array([True, False, False, False, False, False, False])
    numpy.testing.assert_equal(get_cells(simple_ring, "FamName", "D1"), a)


def test_refpts_iterator(simple_ring):
    assert list(refpts_iterator(simple_ring, [0, 1, 2, 3, 4, 5])) == simple_ring
    assert list(refpts_iterator(simple_ring, numpy.ones(6, dtype=bool))) == simple_ring
    assert list(refpts_iterator(simple_ring, [1])) == [simple_ring[1]]
    a = numpy.array([False, True, False, False, False, False])
    assert list(refpts_iterator(simple_ring, a)) == [simple_ring[1]]


def test_get_elements(hmba_lattice):
    # test FamName direct match
    assert get_elements(hmba_lattice, "BPM_06") == [hmba_lattice[65]]
    # test FamName wildcard matching
    assert get_elements(hmba_lattice, "QD2?") == list(hmba_lattice[9, 113])
    assert get_elements(hmba_lattice, "QD3*") == list(hmba_lattice[19, 105])
    assert get_elements(hmba_lattice, "S*H2B") == list([hmba_lattice[55]])
    assert get_elements(hmba_lattice, "*C_1") == list(hmba_lattice[59, 60])
    assert get_elements(hmba_lattice, "DR_2[1-3]") == list(hmba_lattice[54, 56, 58])
    assert get_elements(hmba_lattice, "DR_2[!1-7]") == list(hmba_lattice[52, 78, 80])
    # test element instance
    marker = elements.Marker("M1")
    assert get_elements(hmba_lattice, marker) == list(hmba_lattice[1, 12, 61, 67, 73])
    # test element type
    assert get_elements(hmba_lattice, elements.RFCavity) == [hmba_lattice[0]]
    # test invalid key raises TypeError
    with pytest.raises(TypeError):
        get_elements(hmba_lattice, 1.0)


def test_get_s_pos_returns_zero_for_empty_lattice():
    numpy.testing.assert_equal(get_s_pos([]), numpy.array((0,)))


def test_get_s_pos_returns_length_for_lattice_with_one_element():
    e = elements.LongElement("e", 0.1)
    assert get_s_pos([e], [1]) == numpy.array([0.1])


def test_get_s_pos_returns_all_pts_for_lat_with_2_elements():
    e = elements.LongElement("e", 0.1)
    f = elements.LongElement("e", 2.1)
    numpy.testing.assert_equal(get_s_pos([e, f]), numpy.array([0, 0.1, 2.2]))


def test_get_s_pos_returns_all_pts_for_lat_with_2_elements_using_int_refpts():
    e = elements.LongElement("e", 0.1)
    f = elements.LongElement("e", 2.1)
    lat = [e, f]
    numpy.testing.assert_equal(
        get_s_pos(lat, range(len(lat) + 1)), numpy.array([0, 0.1, 2.2])
    )


def test_get_s_pos_returns_all_pts_for_lat_with_2_elements_using_bool_refpts():
    e = elements.LongElement("e", 0.1)
    f = elements.LongElement("e", 2.1)
    lat = [e, f]
    numpy.testing.assert_equal(
        get_s_pos(lat, numpy.ones(3, dtype=bool)), numpy.array([0, 0.1, 2.2])
    )


def test_tilt_elem():
    elem = elt.Drift("Drift", 1.0)
    # Test tilt_elem function
    tilt_elem(elem, math.pi / 4.0)
    v = math.sqrt(2.0) / 2.0
    a = numpy.diag([v, v, v, v, 1.0, 1.0])
    a[0, 2], a[1, 3], a[2, 0], a[3, 1] = v, v, -v, -v
    numpy.testing.assert_allclose(elem.R1, a, atol=1.0e-15)
    numpy.testing.assert_allclose(elem.R2, a.T, atol=1.0e-15)
    numpy.testing.assert_allclose(elem.tilt, numpy.pi / 4.0)
    # Test tilt property
    elem.tilt = math.pi / 2.0
    a = numpy.diag([0.0, 0.0, 0.0, 0.0, 1.0, 1.0])
    a[0, 2], a[1, 3], a[2, 0], a[3, 1] = 1.0, 1.0, -1.0, -1.0
    numpy.testing.assert_allclose(elem.R1, a, atol=1.0e-15)
    numpy.testing.assert_allclose(elem.R2, a.T, atol=1.0e-15)


def test_shift_elem():
    elem = elt.Drift("Drift", 1.0)
    # Test shift_elem function
    shift_elem(elem, 1.0, 0.5)
    a = numpy.array([1.0, 0.0, 0.5, 0.0, 0.0, 0.0])
    numpy.testing.assert_equal(elem.T1, -a)
    numpy.testing.assert_equal(elem.T2, a)
    numpy.testing.assert_equal(elem.dx, a[0])
    numpy.testing.assert_equal(elem.dy, a[2])
    # Test dx, dy properties
    elem.dx = -2.0
    elem.dy = -1.0
    numpy.testing.assert_equal(elem.T1, 2 * a)
    numpy.testing.assert_equal(elem.T2, -2 * a)


def test_set_tilt(simple_ring):
    set_tilt(simple_ring, [(numpy.pi / 4), 0, 0, 0, (numpy.pi / 4), 0])
    v = 1 / 2**0.5
    a = numpy.diag([v, v, v, v, 1.0, 1.0])
    a[0, 2], a[1, 3], a[2, 0], a[3, 1] = v, v, -v, -v
    numpy.testing.assert_allclose(simple_ring[0].R1, a)
    numpy.testing.assert_allclose(simple_ring[0].R2, a.T)
    numpy.testing.assert_allclose(simple_ring[0].R1, a)
    numpy.testing.assert_allclose(simple_ring[0].R2, a.T)
    ring = [simple_ring[0]]
    set_tilt(ring, 0)
    numpy.testing.assert_allclose(ring[0].R1, numpy.eye(6))
    numpy.testing.assert_allclose(ring[0].R2, numpy.eye(6))


def test_set_shift(simple_ring):
    set_shift(
        simple_ring,
        numpy.array([0.0, 0.0, 0.0, 1.0, 0.0, 0.5]),
        numpy.array([0.0, 0.0, 0.0, 2.0, 0.0, 1.0]),
    )
    a = numpy.array([0.5, 0.0, 1.0, 0.0, 0.0, 0.0])
    numpy.testing.assert_equal(simple_ring[3].T1, -a * 2)
    numpy.testing.assert_equal(simple_ring[3].T2, a * 2)
    numpy.testing.assert_equal(simple_ring[5].T1, -a)
    numpy.testing.assert_equal(simple_ring[5].T2, a)
    ring = [simple_ring[3]]
    set_shift(ring, 3, 5)
    a = numpy.array([3.0, 0.0, 5.0, 0.0, 0.0, 0.0])
    numpy.testing.assert_equal(simple_ring[3].T1, -a)
    numpy.testing.assert_equal(simple_ring[3].T2, a)
