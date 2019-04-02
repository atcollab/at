import at
import numpy
import pytest
from at.load import find_class_name, element_from_dict
from at.load import CLASS_MAPPING, PASS_MAPPING


def test_invalid_class_warns_correctly():
    elem_kwargs = {'Class': 'Invalid'}
    with pytest.warns(at.AtWarning):
        find_class_name(elem_kwargs, quiet=False)
    with pytest.warns(None) as record:
        find_class_name(elem_kwargs, quiet=True)
    assert len(record) is 0


def test_no_pass_method_warns_correctly():
    elem_kwargs = {}
    with pytest.warns(at.AtWarning):
        find_class_name(elem_kwargs, quiet=False)
    with pytest.warns(None) as record:
        find_class_name(elem_kwargs, quiet=True)
    assert len(record) is 0


def test_invalid_pass_method_warns_correctly():
    elem_kwargs = {'PassMethod': 'invalid'}
    with pytest.warns(at.AtWarning):
        find_class_name(elem_kwargs, quiet=False)
    with pytest.warns(None) as record:
        find_class_name(elem_kwargs, quiet=True)
    assert len(record) is 0


def test_class_mapping():
    elem_kwargs = {'FamName': 'fam'}
    for class_name in CLASS_MAPPING:
        elem_kwargs['Class'] = class_name
        assert find_class_name(elem_kwargs) == CLASS_MAPPING[class_name]
        elem_kwargs['Class'] = class_name.upper()
        assert find_class_name(elem_kwargs) == CLASS_MAPPING[class_name]


def test_family_mapping():
    elem_kwargs = {}
    for family_name in CLASS_MAPPING:
        elem_kwargs['FamName'] = family_name
        assert find_class_name(elem_kwargs) == CLASS_MAPPING[family_name]
        elem_kwargs['FamName'] = family_name.upper()
        assert find_class_name(elem_kwargs) == CLASS_MAPPING[family_name]


def test_PassMethod_mapping():
    elem_kwargs = {'FamName': 'fam'}
    for pass_method in PASS_MAPPING.keys():
        elem_kwargs['PassMethod'] = pass_method
        assert find_class_name(elem_kwargs) == PASS_MAPPING[pass_method]


def test_find_Aperture():
    elem_kwargs = {'Limits': [-0.5, 0.5, -0.5, 0.5], 'FamName': 'fam'}
    assert find_class_name(elem_kwargs, True) == 'Aperture'


@pytest.mark.parametrize('elem_kwargs', (
        {'Voltage': 2.5e+6, 'FamName': 'fam'},
        {'Frequency': 5.e+8, 'FamName': 'fam'},
        {'HarmNumber': 1000, 'FamName': 'fam'},
        {'PhaseLag': 0, 'FamName': 'fam'},
        {'TimeLag': 0.0, 'FamName': 'fam'}))
def test_find_RFCavity(elem_kwargs):
    assert find_class_name(elem_kwargs, True) == 'RFCavity'


def test_find_Monitor():
    elem_kwargs = {'GCR': [1, 1, 0, 0], 'FamName': 'fam'}
    assert find_class_name(elem_kwargs, True) == 'Monitor'


@pytest.mark.parametrize('elem_kwargs', (
        {'FullGap': 0.05, 'FamName': 'fam'},
        {'FringeInt1': 0.5, 'FamName': 'fam'},
        {'FringeInt2': 0.5, 'FamName': 'fam'},
        {'gK': 0.03, 'FamName': 'fam'},
        {'EntranceAngle': 0.05, 'FamName': 'fam'},
        {'ExitAngle': 0.05, 'FamName': 'fam'}))
def test_find_Dipole(elem_kwargs):
    assert find_class_name(elem_kwargs, True) == 'Dipole'


def test_find_Corrector():
    elem_kwargs = {'KickAngle': [0, 0], 'FamName': 'fam'}
    assert find_class_name(elem_kwargs, True) == 'Corrector'


def test_find_RingParam():
    elem_kwargs = {'Periodicity': 1, 'FamName': 'fam'}
    assert find_class_name(elem_kwargs, True) == 'RingParam'


def test_find_M66():
    elem_kwargs = {'M66': numpy.eye(6), 'FamName': 'fam'}
    assert find_class_name(elem_kwargs, True) == 'M66'


@pytest.mark.parametrize('elem_kwargs', (
        {'K': -0.5, 'FamName': 'fam'},
        {'PolynomB': [0, 1, 0, 0], 'FamName': 'fam'}))
def test_find_Quadrupole(elem_kwargs):
    assert find_class_name(elem_kwargs, True) == 'Quadrupole'


@pytest.mark.parametrize('elem_kwargs', (
        {'PolynomB': [0, 0, 0, 0], 'PassMethod': 'StrMPoleSymplectic4Pass',
         'FamName': 'fam'},
        {'PolynomB': [0, 0, 0, 0], 'Length': 1, 'FamName': 'fam'}))
def test_find_Multipole(elem_kwargs):
    assert find_class_name(elem_kwargs, True) == 'Multipole'


def test_find_Drift():
    elem_kwargs = {'Length': 1, 'FamName': 'fam'}
    assert find_class_name(elem_kwargs, True) == 'Drift'


def test_find_Sextupole():
    elem_kwargs = {'PolynomB': [0, 0, 1, 0], 'FamName': 'fam'}
    assert find_class_name(elem_kwargs, True) == 'Sextupole'


def test_find_Octupole():
    elem_kwargs = {'PolynomB': [0, 0, 0, 1], 'FamName': 'fam'}
    assert find_class_name(elem_kwargs, True) == 'Octupole'


@pytest.mark.parametrize('elem_kwargs', (
        {'PolynomB': [0, 0, 0, 0], 'FamName': 'fam'},
        {'PolynomB': [0, 0, 0, 0, 1], 'FamName': 'fam'},
        {'PolynomB': [0, 0, 0, 0], 'FamName': 'fam'},
        {'PolynomB': [0, 0, 0, 0, 1], 'Length': 0, 'FamName': 'fam'}))
def test_find_ThinMultipole(elem_kwargs):
    assert find_class_name(elem_kwargs, True) == 'ThinMultipole'


@pytest.mark.parametrize('elem_kwargs', ({'FamName': 'fam'},
                                         {'Length': 0.0, 'FamName': 'fam'}))
def test_find_Marker(elem_kwargs):
    assert find_class_name(elem_kwargs, True) == 'Marker'


@pytest.mark.parametrize('elem_kwargs', (
        {'Class': 'Marker', 'PassMethod': 'IdentityPass', 'Length': 1,
         'FamName': 'fam'},
        {'Class': 'Quadrupole', 'PassMethod': 'CavityPass', 'Length': 1,
         'FamName': 'fam'},
        {'Class': 'Marker', 'PassMethod': 'StrMPoleSymplectic4Pass',
         'FamName': 'fam'},
        {'Class': 'Monitor', 'PassMethod': 'StrMPoleSymplectic4Pass',
         'FamName': 'fam'},
        {'Class': 'RingParam', 'PassMethod': 'StrMPoleSymplectic4Pass',
         'Energy': 3E9, 'FamName': 'fam'},
        {'Class': 'Drift', 'PassMethod': 'StrMPoleSymplectic4Pass',
         'Length': 1.0, 'FamName': 'fam'},
        {'Class': 'Drift', 'PassMethod': 'InvalidPass'}))
def test_sanitise_class_error(elem_kwargs):
    with pytest.raises(AttributeError):
        elem = element_from_dict(elem_kwargs)
