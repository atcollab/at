import numpy
import pytest

from at.load.tracy import (
    tokenise_expression,
    parse_float,
    expand_tracy,
    parse_lines,
    tracy_element_from_string,
)
from at.lattice.elements import (
    Dipole,
    Drift,
    Marker,
    Multipole,
    Quadrupole,
    RFCavity,
    Sextupole,
)


@pytest.fixture
def defaults():
    return {"energy": 3e9, "harmonic_number": 936}


@pytest.mark.parametrize(
    "exp,target",
    [
        ["3", ["3"]],
        ["-3", ["-", "3"]],
        ["-4.85d-3", ["-", "4.85d-3"]],
        ["3+4", ["3", "+", "4"]],
        ["3 + 4", ["3", "+", "4"]],
        ["3+4e3", ["3", "+", "4e3"]],
        ["1e3+4e-3", ["1e3", "+", "4e-3"]],
        ["1 - (3 / 2)", ["1", "-", "(", "3", "/", "2", ")"]],
        ["(1 - 1e-1) + 3", ["(", "1", "-", "1e-1", ")", "+", "3"]],
        ["0.485-0.002", ["0.485", "-", "0.002"]],
        [" (3.482-0.02)-1.85", ["(", "3.482", "-", "0.02", ")", "-", "1.85"]],
        ["tbx-tbev-tbes", ["tbx", "-", "tbev", "-", "tbes"]],
    ],
)
def test_tokenise_expression(exp, target):
    assert tokenise_expression(exp) == target


@pytest.mark.parametrize(
    "exp,target",
    [
        ["3.4", 3.4],
        ["-3", -3],
        ["3+4", 7],
        ["3.8*1.1", 4.18],
        ["3+4e3", 4003],
        ["1e3+4e-3", 1000.004],
        ["1 + 12 / (a + 1)", 5],
        ["0.485-0.002", 0.483],
        ["4.85d3", 4850],
        ["-4.85d-3", -0.00485],
    ],
)
def test_parse_float(exp, target):
    assert parse_float(exp, {"a": 2}) == target


def test_parse_lines_removes_comments():
    assert parse_lines("a{}b") == ["ab"]


def test_parse_lines_removes_whitespace():
    assert parse_lines("a\n{a\nb}b") == ["ab"]


def test_parse_line_splits_on_semicolons():
    assert parse_lines("a;b") == ["a", "b"]


def test_tracy_element_from_string_handles_drift(defaults):
    drift = "drift,l=0.0450000"
    assert tracy_element_from_string("d1", drift, defaults).equals(
        Drift("d1", 0.0450000)
    )


def test_tracy_element_from_string_handles_marker(defaults):
    marker = "marker"
    expected = Marker("m1")
    assert tracy_element_from_string("m1", marker, defaults).equals(expected)


def test_tracy_element_from_string_handles_quadrupole(defaults):
    quad = "quadrupole,l=0.15,k=10.24,n=nquad,method=4"
    defaults["nquad"] = 10
    expected = Quadrupole(
        "q1",
        0.15,
        k=10.24,
        NumIntSteps=10,
        method="4",
        PassMethod="StrMPoleSymplectic4Pass",
    )
    assert tracy_element_from_string("q1", quad, defaults).equals(expected)


def test_tracy_element_from_string_handles_sextupole(defaults):
    sext = "sextupole,l=0.14,k=174.4,n=nsext,method=4"
    defaults["nsext"] = 2
    expected = Sextupole("s1", 0.14, h=174.4, NumIntSteps=2, method="4")
    assert tracy_element_from_string("s1", sext, defaults).equals(expected)


def test_tracy_element_from_string_handles_hom(defaults):
    oct = "multipole,l=0.0,hom=(4,1.0,0.3)"
    expected = Multipole("m1", 0.0, poly_a=[0, 0, 0, 0.3], poly_b=[0, 0, 0, 1])
    assert tracy_element_from_string("m1", oct, defaults).equals(expected)


def test_tracy_element_from_string_handles_cavity(defaults):
    cavity = "cavity,l=0.0,frequency=499.654e6,voltage=2.2e6,phi=0.0"
    expected = RFCavity("c1", 0.0, 2.2e6, 4.99654e8, 936, 3e9, Phi=0.0)
    constructed = tracy_element_from_string("c1", cavity, defaults)
    assert constructed.equals(expected)


def test_tracy_element_from_string_handles_bending(defaults):
    bend = "bending,l= 0.20000000,t=0.32969999,t1=0.00000000," \
           "t2=0.32969999,k=-0.124,gap=0.06,n=nbend,method=4"
    defaults["nbend"] = 2
    expected = Dipole(
        "b1",
        0.2,
        BendingAngle=0.00575435036929238,
        EntranceAngle=0,
        ExitAngle=0.00575435036929238,
        k=-0.12411107,
        FullGap=0.03,
        FringeInt1=1,
        FringeInt2=1,
        NumIntSteps=2,
        PolynomA=numpy.array([0, 0, 0, 0], dtype=numpy.float64),
        PolynomB=numpy.array([0, -0.124, 0, 0], dtype=numpy.float64),
        method="4",
        PassMethod="BndMPoleSymplectic4Pass",
    )
    assert tracy_element_from_string("b1", bend, defaults).equals(expected)


def test_tracy_element_from_string_handles_variable(defaults):
    drift = "drift,l=a"
    defaults["a"] = 1
    correct_drift = Drift("d1", 1)
    constructed_drift = tracy_element_from_string("d1", drift, defaults)
    assert correct_drift.equals(constructed_drift)


def test_expand_tracy():
    contents = "define lattice;energy=1;dmult:drift,l=1;cell:dmult;end;"
    elements, energy = expand_tracy(contents, "cell", 936)
    assert len(elements) == 1
    assert energy == 1e9
    assert elements[0].equals(Drift("dmult", 1.0))
