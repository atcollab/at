import at
import numpy
import pytest
from at.load_mat import find_class_name, element_from_dict


def test_invalid_class_raises_AttributeError():
    elem_kwargs = {'Class': 'Invalid'}
    with pytest.raises(AttributeError):
        find_class_name(elem_kwargs)


def test_correct_class_names():
    elem_kwargs = {'FamName': 'fam'}
    for class_name in at.load_mat.CLASSES:
        elem_kwargs['Class'] = class_name
        assert find_class_name(elem_kwargs) == class_name


def test_class_mapping():
    elem_kwargs = {'FamName': 'fam'}
    for class_name in at.load_mat.CLASS_MAPPING.keys():
        elem_kwargs['Class'] = class_name
        assert find_class_name(elem_kwargs) == at.load_mat.CLASS_MAPPING[class_name]


def test_family_mapping():
    elem_kwargs = {}
    for family_name in at.load_mat.CLASSES:
        elem_kwargs['FamName'] = family_name
        assert find_class_name(elem_kwargs) == family_name
    for family_name in at.load_mat.FAMILY_MAPPING.keys():
        elem_kwargs['FamName'] = family_name
        assert find_class_name(elem_kwargs) == at.load_mat.FAMILY_MAPPING[family_name]
        elem_kwargs['FamName'] = family_name.upper()
        assert find_class_name(elem_kwargs) == at.load_mat.FAMILY_MAPPING[family_name]


def test_PassMethod_mapping():
    elem_kwargs = {'FamName': 'fam'}
    for pass_method in at.load_mat.PASSMETHOD_MAPPING.keys():
        elem_kwargs['PassMethod'] = pass_method
        assert find_class_name(elem_kwargs) == at.load_mat.PASSMETHOD_MAPPING[pass_method]


def test_find_Aperture():
    elem_kwargs = {'Limits': [-0.5, 0.5, -0.5, 0.5], 'FamName': 'fam'}
    assert find_class_name(elem_kwargs) == 'Aperture'


def test_find_RFCavity():
    elem_kwargs = {'Voltage': 2.5e+6, 'FamName': 'fam'}
    assert find_class_name(elem_kwargs) == 'RFCavity'
    elem_kwargs = {'Frequency': 5.e+8, 'FamName': 'fam'}
    assert find_class_name(elem_kwargs) == 'RFCavity'
    elem_kwargs = {'HarmNumber': 1000, 'FamName': 'fam'}
    assert find_class_name(elem_kwargs) == 'RFCavity'
    elem_kwargs = {'PhaseLag': 0, 'FamName': 'fam'}
    assert find_class_name(elem_kwargs) == 'RFCavity'
    elem_kwargs = {'TimeLag': 0.0, 'FamName': 'fam'}
    assert find_class_name(elem_kwargs) == 'RFCavity'


def test_find_Monitor():
    elem_kwargs = {'GCR': [1, 1, 0, 0], 'FamName': 'fam'}
    assert find_class_name(elem_kwargs) == 'Monitor'


@pytest.mark.parametrize('elem_kwargs', (
        {'FullGap': 0.05, 'FamName': 'fam'},
        {'FringeInt1': 0.5, 'FamName': 'fam'},
        {'FringeInt2': 0.5, 'FamName': 'fam'},
        {'EntranceAngle': 0.05, 'FamName': 'fam'},
        {'ExitAngle': 0.05, 'FamName': 'fam'},
        {'PassMethod': 'BndMPoleSymplectic4Pass', 'PolynomB': [0, 0, 0, 0], 'FamName': 'fam'}
))
def test_find_Dipole(elem_kwargs):
    assert find_class_name(elem_kwargs) == 'Dipole'


def test_find_Corrector():
    elem_kwargs = {'KickAngle': [0, 0], 'FamName': 'fam'}
    assert find_class_name(elem_kwargs) == 'Corrector'


def test_find_RingParam():
    elem_kwargs = {'Periodicity': 1, 'FamName': 'fam'}
    assert find_class_name(elem_kwargs) == 'RingParam'


def test_find_M66():
    elem_kwargs = {'M66': numpy.eye(6), 'FamName': 'fam'}
    assert find_class_name(elem_kwargs) == 'M66'


def test_find_Quadrupole():
    elem_kwargs = {'K': -0.5, 'FamName': 'fam'}
    assert find_class_name(elem_kwargs) == 'Quadrupole'


def test_find_Multipole():
    elem_kwargs = {'PolynomB': [0, 0, 0, 0],
                   'PassMethod': 'StrMPoleSymplectic4Pass', 'FamName': 'fam'}
    assert find_class_name(elem_kwargs) == 'Multipole'


def test_find_Drift():
    elem_kwargs = {'Length': 1.0, 'FamName': 'fam'}
    assert find_class_name(elem_kwargs) == 'Drift'


def test_find_Sextupole():
    elem_kwargs = {'PolynomB': [0, 0, 1, 0], 'FamName': 'fam'}
    assert find_class_name(elem_kwargs) == 'Sextupole'


def test_find_Octupole():
    elem_kwargs = {'PolynomB': [0, 0, 0, 1], 'PolynomA': [0, 0, 0, 0],
                   'FamName': 'fam'}
    assert find_class_name(elem_kwargs) == 'Octupole'


def test_find_ThinMultipole():
    elem_kwargs = {'PolynomB': [0, 0, 0, 0, 1], 'PolynomA': [0, 0, 0, 0],
                   'FamName': 'fam'}
    assert find_class_name(elem_kwargs) == 'ThinMultipole'
    elem_kwargs = {'PolynomB': [0, 0, 0, 0, 1], 'PolynomA': [0, 0, 0, 0],
                   'Length': 0, 'FamName': 'fam'}
    assert find_class_name(elem_kwargs) == 'ThinMultipole'


def test_find_Marker():
    elem_kwargs = {'Length': 0.0, 'FamName': 'fam'}
    assert find_class_name(elem_kwargs) == 'Marker'
    elem_kwargs = {'FamName': 'fam'}
    assert find_class_name(elem_kwargs) == 'Marker'


@pytest.mark.parametrize('elem_kwargs', (
        {'FamName': '', 'Class': 'Marker', 'PassMethod': 'IdentityPass', 'Length': 1.0},
        {'FamName': '', 'Class': 'Quadrupole', 'PassMethod': 'CavityPass', 'Length': 1.0},
        {'FamName': '', 'Class': 'Marker', 'PassMethod': 'Invalid'},
        {'FamName': '', 'Class': 'Monitor', 'PassMethod': 'Invalid'},
        {'FamName': '', 'Class': 'RingParam', 'PassMethod': 'Invalid', 'Energy': 3E9},
        {'FamName': '', 'Class': 'Drift', 'PassMethod': 'IdentityPass', 'Length': 1.0}
))
def test_sanitise_class_error(elem_kwargs):
    with pytest.raises(AttributeError):
        elem = element_from_dict(elem_kwargs)
