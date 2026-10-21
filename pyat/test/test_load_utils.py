import numpy
import pytest
from at import AtWarning
from at.lattice import Lattice, elements, params_filter, no_filter
from at.load.utils import find_class

# noinspection PyProtectedMember
from at.load.utils import _CLASS_MAP, _PASS_MAP
from at.load.utils import RingParam, split_ignoring_parentheses

# noinspection PyProtectedMember
from at.load.matfile import ringparam_filter


def _matlab_scanner(element_list, **kwargs):
    """This function simulates a mat-file reading but replaces the mat-file by
    a list of elements."""
    latt = Lattice(
        ringparam_filter, no_filter, element_list, iterator=params_filter, **kwargs
    )
    return vars(latt)


# noinspection pep8-naming
def test_lattice_gets_attributes_from_RingParam():
    p = RingParam("lattice_name", 3.0e6, Periodicity=32)
    params = _matlab_scanner([p])
    assert params["name"] == "lattice_name"
    assert params["_energy"] == 3.0e6
    assert params["periodicity"] == 32
    assert params["_radiation"] is False


def test_lattice_gets_attributes_from_elements():
    d = elements.Dipole(
        "d1",
        1,
        BendingAngle=numpy.pi,
        Energy=3.0e6,
        PassMethod="BndMPoleSymplectic4RadPass",
    )
    params = _matlab_scanner([d])
    assert params["name"] == ""
    assert params["_energy"] == 3.0e6
    assert params["periodicity"] == 2
    assert params["_radiation"] is True


def test_no_bending_in_the_cell_warns_correctly():
    with pytest.warns(AtWarning):
        params = _matlab_scanner([], energy=3.0e9)
    assert params["periodicity"] == 1


def test_non_integer_number_of_cells_warns_correctly():
    d = elements.Dipole("d1", 1, BendingAngle=0.195)
    with pytest.warns(AtWarning):
        params = _matlab_scanner([d], energy=3.0e9)
    assert params["periodicity"] == 32


def test_inconsistent_energy_values_warns_correctly():
    d1 = elements.Dipole(
        "d1",
        1,
        BendingAngle=numpy.pi,
        Energy=5.0e9,
        PassMethod="BndMPoleSymplectic4RadPass",
    )
    d2 = elements.Dipole(
        "d2",
        1,
        BendingAngle=numpy.pi,
        Energy=3.0e9,
        PassMethod="BndMPoleSymplectic4RadPass",
    )
    m1 = elements.Marker("m1", Energy=5.0e9)
    m2 = elements.Marker("m2", Energy=3.0e9)
    with pytest.warns(AtWarning):
        params = _matlab_scanner([m1, m2])
    assert params["_energy"] == 5.0e9
    with pytest.warns(AtWarning):
        params = _matlab_scanner([d1, d2])
    assert params["_energy"] == 5.0e9


def test_more_than_one_RingParam_in_ring_raises_warning():
    p1 = RingParam("rp1", 5.0e6, 1)
    p2 = RingParam("rp2", 3.0e6, 1)
    with pytest.warns(AtWarning):
        params = _matlab_scanner([p1, p2])
    assert params["_energy"] == 5.0e6


def test_family_mapping():
    for class_name, cls in _CLASS_MAP.items():
        with pytest.warns(AtWarning):
            cls_name = find_class({"FamName": class_name})
        assert cls_name == cls.__name__


def test_PassMethod_mapping():
    for pass_method, class_name in _PASS_MAP.items():
        with pytest.warns(AtWarning):
            cls_name = find_class({"FamName": "fam", "PassMethod": pass_method})
        assert cls_name == class_name


def test_find_Aperture():
    elem_kwargs = {"Limits": [-0.5, 0.5, -0.5, 0.5], "FamName": "fam"}
    with pytest.warns(AtWarning):
        cls_name = find_class(elem_kwargs)
    assert cls_name == "Aperture"


@pytest.mark.parametrize(
    "elem_kwargs",
    (
        {"Voltage": 2.5e6, "FamName": "fam"},
        {"Frequency": 5.0e8, "FamName": "fam"},
        {"HarmNumber": 1000, "FamName": "fam"},
    ),
)
def test_find_RFCavity(elem_kwargs):
    with pytest.warns(AtWarning):
        cls_name = find_class(elem_kwargs)
    assert cls_name == "RFCavity"


def test_find_Monitor():
    elem_kwargs = {"GCR": [1, 1, 0, 0], "FamName": "fam"}
    with pytest.warns(AtWarning):
        cls_name = find_class(elem_kwargs)
    assert cls_name == "Monitor"


@pytest.mark.parametrize(
    "elem_kwargs",
    (
        {"BendingAngle": 0.05, "FamName": "fam"},
        {"EntranceAngle": 0.05, "FamName": "fam"},
        {"ExitAngle": 0.05, "FamName": "fam"},
    ),
)
def test_find_Dipole(elem_kwargs):
    with pytest.warns(AtWarning):
        cls_name = find_class(elem_kwargs)
    assert cls_name == "Dipole"


def test_find_Corrector():
    elem_kwargs = {"KickAngle": [0, 0], "FamName": "fam"}
    with pytest.warns(AtWarning):
        cls_name = find_class(elem_kwargs)
    assert cls_name == "Corrector"


def test_find_RingParam():
    elem_kwargs = {"Periodicity": 1, "FamName": "fam"}
    with pytest.warns(AtWarning):
        cls_name = find_class(elem_kwargs)
    assert cls_name == "RingParam"


def test_find_M66():
    elem_kwargs = {"M66": numpy.eye(6), "FamName": "fam"}
    with pytest.warns(AtWarning):
        cls_name = find_class(elem_kwargs)
    assert cls_name == "M66"


@pytest.mark.parametrize(
    "elem_kwargs",
    ({"K": -0.5, "FamName": "fam"}, {"PolynomB": [0, 1, 0, 0], "FamName": "fam"}),
)
def test_find_Quadrupole(elem_kwargs):
    with pytest.warns(AtWarning):
        cls_name = find_class(elem_kwargs)
    assert cls_name == "Quadrupole"


@pytest.mark.parametrize(
    "elem_kwargs",
    (
        {
            "PolynomB": [0, 0, 0, 0],
            "PassMethod": "StrMPoleSymplectic4Pass",
            "FamName": "fam",
        },
        {"PolynomB": [0, 0, 0, 0], "Length": 1, "FamName": "fam"},
    ),
)
def test_find_Multipole(elem_kwargs):
    with pytest.warns(AtWarning):
        cls_name = find_class(elem_kwargs)
    assert cls_name == "Multipole"


def test_find_Drift():
    elem_kwargs = {"Length": 1, "FamName": "fam"}
    with pytest.warns(AtWarning):
        cls_name = find_class(elem_kwargs)
    assert cls_name == "Drift"


def test_find_Sextupole():
    elem_kwargs = {"PolynomB": [0, 0, 1, 0], "FamName": "fam"}
    with pytest.warns(AtWarning):
        cls_name = find_class(elem_kwargs)
    assert cls_name == "Sextupole"


def test_find_Octupole():
    elem_kwargs = {"PolynomB": [0, 0, 0, 1], "FamName": "fam"}
    with pytest.warns(AtWarning):
        cls_name = find_class(elem_kwargs)
    assert cls_name == "Octupole"


@pytest.mark.parametrize(
    "elem_kwargs",
    (
        {"PolynomB": [0, 0, 0, 0], "FamName": "fam"},
        {"PolynomB": [0, 0, 0, 0, 1], "FamName": "fam"},
        {"PolynomB": [0, 0, 0, 0, 1], "Length": 0, "FamName": "fam"},
    ),
)
def test_find_ThinMultipole(elem_kwargs):
    with pytest.warns(AtWarning):
        cls_name = find_class(elem_kwargs)
    assert cls_name == "ThinMultipole"


@pytest.mark.parametrize(
    "elem_kwargs",
    (
        {"FamName": "fam", "PassMethod": "IdentityPass"},
        {"Length": 0.0, "FamName": "fam", "PassMethod": "IdentityPass"},
    ),
)
def test_find_Marker(elem_kwargs):
    with pytest.warns(AtWarning):
        cls_name = find_class(elem_kwargs)
    assert cls_name == "Marker"


@pytest.mark.parametrize(
    "string,delimiter,target",
    [
        ["a,b", ",", ["a", "b"]],
        ["a,b(c,d)", ",", ["a", "b(c,d)"]],
        ["l=0,hom(4,0.0,0)", ",", ["l=0", "hom(4,0.0,0)"]],
        ["inv(arca_c1r),3*(ms,arca_c2)", ",", ["inv(arca_c1r)", "3*(ms,arca_c2)"]],
    ],
)
def test_split_ignoring_parentheses(string, delimiter, target):
    assert split_ignoring_parentheses(string, delimiter) == target
