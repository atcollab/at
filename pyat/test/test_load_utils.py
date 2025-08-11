import at
import numpy
import pytest
import warnings
from at import AtWarning
from at.lattice import Lattice, elements, params_filter, no_filter
from at.load.utils import find_class, element_from_dict
# noinspection PyProtectedMember
from at.load.utils import _CLASS_MAP, _PASS_MAP
from at.load.utils import RingParam, split_ignoring_parentheses
# noinspection PyProtectedMember
from at.load.matfile import ringparam_filter


def _matlab_scanner(element_list, **kwargs):
    """This function simulates a mat-file reading but replaces the mat-file by
    a list of elements"""
    latt = Lattice(ringparam_filter, no_filter, element_list,
                   iterator=params_filter, **kwargs)
    return vars(latt)


def test_lattice_gets_attributes_from_RingParam():
    p = RingParam('lattice_name', 3.e+6, Periodicity=32)
    params = _matlab_scanner([p])
    assert params['name'] == 'lattice_name'
    assert params['_energy'] == 3.e+6
    assert params['periodicity'] == 32
    assert params['_radiation'] is False


def test_lattice_gets_attributes_from_elements():
    d = elements.Dipole('d1', 1, BendingAngle=numpy.pi, Energy=3.e+6,
                        PassMethod='BndMPoleSymplectic4RadPass')
    params = _matlab_scanner([d])
    assert params['name'] == ''
    assert params['_energy'] == 3.e+6
    assert params['periodicity'] == 2
    assert params['_radiation'] is True


def test_no_bending_in_the_cell_warns_correctly():
    with pytest.warns(AtWarning):
        params = _matlab_scanner([], energy=3.e+9)
        assert params['periodicity'] == 1


def test_non_integer_number_of_cells_warns_correctly():
    d = elements.Dipole('d1', 1, BendingAngle=0.195)
    with pytest.warns(AtWarning):
        params = _matlab_scanner([d], energy=3.e+9)
        assert params['periodicity'] == 32


def test_inconsistent_energy_values_warns_correctly():
    d1 = elements.Dipole('d1', 1, BendingAngle=numpy.pi, Energy=5.e+9,
                         PassMethod='BndMPoleSymplectic4RadPass')
    d2 = elements.Dipole('d2', 1, BendingAngle=numpy.pi, Energy=3.e+9,
                         PassMethod='BndMPoleSymplectic4RadPass')
    m1 = elements.Marker('m1', Energy=5.e+9)
    m2 = elements.Marker('m2', Energy=3.e+9)
    with pytest.warns(AtWarning):
        params = _matlab_scanner([m1, m2])
        assert params['_energy'] == 5.e+9
    with pytest.warns(AtWarning):
        params = _matlab_scanner([d1, d2])
        assert params['_energy'] == 5.e+9


def test_more_than_one_RingParam_in_ring_raises_warning():
    p1 = RingParam('rp1', 5.e+6, 1)
    p2 = RingParam('rp2', 3.e+6, 1)
    with pytest.warns(AtWarning):
        params = _matlab_scanner([p1, p2])
        assert params['_energy'] == 5.e+6


def test_invalid_class_warns_correctly():
    elem_kwargs = {'Class': 'Invalid'}
    with pytest.warns(at.AtWarning):
        find_class(elem_kwargs, quiet=False)
    with warnings.catch_warnings(record=True) as record:
        find_class(elem_kwargs, quiet=True)
    assert len(record) == 0


def test_no_pass_method_warns_correctly():
    elem_kwargs = {}
    with pytest.warns(at.AtWarning):
        find_class(elem_kwargs, quiet=False)
    with warnings.catch_warnings(record=True) as record:
        find_class(elem_kwargs, quiet=True)
    assert len(record) == 0


def test_invalid_pass_method_warns_correctly():
    elem_kwargs = {'PassMethod': 'invalid'}
    with pytest.warns(at.AtWarning):
        find_class(elem_kwargs, quiet=False)
    with warnings.catch_warnings(record=True) as record:
        find_class(elem_kwargs, quiet=True)
    assert len(record) == 0


def test_class_mapping():
    elem_kwargs = {'FamName': 'fam'}
    for class_name in _CLASS_MAP:
        elem_kwargs['Class'] = class_name
        assert find_class(elem_kwargs) is _CLASS_MAP[class_name]
        elem_kwargs['Class'] = class_name.upper()
        assert find_class(elem_kwargs) is _CLASS_MAP[class_name]


def test_family_mapping():
    elem_kwargs = {}
    for family_name in _CLASS_MAP:
        elem_kwargs['FamName'] = family_name
        assert find_class(elem_kwargs) is _CLASS_MAP[family_name]
        elem_kwargs['FamName'] = family_name.upper()
        assert find_class(elem_kwargs) is _CLASS_MAP[family_name]


def test_PassMethod_mapping():
    elem_kwargs = {'FamName': 'fam'}
    for pass_method in _PASS_MAP.keys():
        elem_kwargs['PassMethod'] = pass_method
        assert find_class(elem_kwargs) == _PASS_MAP[pass_method]


def test_find_Aperture():
    elem_kwargs = {'Limits': [-0.5, 0.5, -0.5, 0.5], 'FamName': 'fam'}
    assert find_class(elem_kwargs, True) is elements.Aperture


@pytest.mark.parametrize('elem_kwargs', (
        {'Voltage': 2.5e+6, 'FamName': 'fam'},
        {'Frequency': 5.e+8, 'FamName': 'fam'},
        {'HarmNumber': 1000, 'FamName': 'fam'},
        {'PhaseLag': 0, 'FamName': 'fam'},
        {'TimeLag': 0.0, 'FamName': 'fam'}))
def test_find_RFCavity(elem_kwargs):
    assert find_class(elem_kwargs, True) is elements.RFCavity


def test_find_Monitor():
    elem_kwargs = {'GCR': [1, 1, 0, 0], 'FamName': 'fam'}
    assert find_class(elem_kwargs, True) is elements.Monitor


@pytest.mark.parametrize('elem_kwargs', (
        {'FullGap': 0.05, 'FamName': 'fam'},
        {'FringeInt1': 0.5, 'FamName': 'fam'},
        {'FringeInt2': 0.5, 'FamName': 'fam'},
        {'gK': 0.03, 'FamName': 'fam'},
        {'EntranceAngle': 0.05, 'FamName': 'fam'},
        {'ExitAngle': 0.05, 'FamName': 'fam'}))
def test_find_Dipole(elem_kwargs):
    assert find_class(elem_kwargs, True) is elements.Dipole


def test_find_Corrector():
    elem_kwargs = {'KickAngle': [0, 0], 'FamName': 'fam'}
    assert find_class(elem_kwargs, True) is elements.Corrector


def test_find_RingParam():
    elem_kwargs = {'Periodicity': 1, 'FamName': 'fam'}
    assert find_class(elem_kwargs, True) is RingParam


def test_find_M66():
    elem_kwargs = {'M66': numpy.eye(6), 'FamName': 'fam'}
    assert find_class(elem_kwargs, True) is elements.M66


@pytest.mark.parametrize('elem_kwargs', (
        {'K': -0.5, 'FamName': 'fam'},
        {'PolynomB': [0, 1, 0, 0], 'FamName': 'fam'}))
def test_find_Quadrupole(elem_kwargs):
    assert find_class(elem_kwargs, True) is elements.Quadrupole


@pytest.mark.parametrize('elem_kwargs', (
        {'PolynomB': [0, 0, 0, 0], 'PassMethod': 'StrMPoleSymplectic4Pass',
         'FamName': 'fam'},
        {'PolynomB': [0, 0, 0, 0], 'Length': 1, 'FamName': 'fam'}))
def test_find_Multipole(elem_kwargs):
    assert find_class(elem_kwargs, True) is elements.Multipole


def test_find_Drift():
    elem_kwargs = {'Length': 1, 'FamName': 'fam'}
    assert find_class(elem_kwargs, True) is elements.Drift


def test_find_Sextupole():
    elem_kwargs = {'PolynomB': [0, 0, 1, 0], 'FamName': 'fam'}
    assert find_class(elem_kwargs, True) is elements.Sextupole


def test_find_Octupole():
    elem_kwargs = {'PolynomB': [0, 0, 0, 1], 'FamName': 'fam'}
    assert find_class(elem_kwargs, True) is elements.Octupole


@pytest.mark.parametrize('elem_kwargs', (
        {'PolynomB': [0, 0, 0, 0], 'FamName': 'fam'},
        {'PolynomB': [0, 0, 0, 0, 1], 'FamName': 'fam'},
        {'PolynomB': [0, 0, 0, 0], 'FamName': 'fam'},
        {'PolynomB': [0, 0, 0, 0, 1], 'Length': 0, 'FamName': 'fam'}))
def test_find_ThinMultipole(elem_kwargs):
    assert find_class(elem_kwargs, True) is elements.ThinMultipole


@pytest.mark.parametrize('elem_kwargs', (
        {'FamName': 'fam', 'PassMethod': 'IdentityPass'},
        {'Length': 0.0, 'FamName': 'fam', 'PassMethod': 'IdentityPass'}))
def test_find_Marker(elem_kwargs):
    assert find_class(elem_kwargs, True) is elements.Marker


@pytest.mark.parametrize('elem_kwargs', (
        {'FamName': 'fam'},
        {'Length': 0.0, 'FamName': 'fam'}))
def test_find_Element(elem_kwargs):
    assert find_class(elem_kwargs, True) is elements.Element


@pytest.mark.parametrize('elem_kwargs', (
        {'Class': 'Marker', 'PassMethod': 'IdentityPass', 'Length': 1,
         'FamName': 'fam'},
        {'Class': 'Quadrupole', 'PassMethod': 'CavityPass', 'Length': 1,
         'FamName': 'fam'},
        {'Class': 'Marker', 'PassMethod': 'StrMPoleSymplectic4Pass',
         'FamName': 'fam'},
        {'Class': 'Drift', 'PassMethod': 'StrMPoleSymplectic4Pass',
         'Length': 1.0, 'FamName': 'fam'},
        {'Class': 'Drift', 'PassMethod': 'InvalidPass',
         'Length': 1.0, 'FamName': 'fam'}))
def test_sanitise_class_error(elem_kwargs):
    with pytest.warns(AtWarning):
        element_from_dict(elem_kwargs)


@pytest.mark.parametrize(
    "string,delimiter,target",
    [
        ["a,b", ",", ["a", "b"]],
        ["a,b(c,d)", ",", ["a", "b(c,d)"]],
        ["l=0,hom(4,0.0,0)", ",", ["l=0", "hom(4,0.0,0)"]],
        ["inv(arca_c1r),3*(ms,arca_c2)", ",",
         ["inv(arca_c1r)", "3*(ms,arca_c2)"]]
    ],
)
def test_split_ignoring_parentheses(string, delimiter, target):
    assert split_ignoring_parentheses(string, delimiter) == target
