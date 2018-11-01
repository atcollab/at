import at
import numpy
import pytest
from at.load_mat import find_class_name, sanitise_class


class_list = list(at.load_mat.CLASSES)
class_mapped = at.load_mat.CLASS_MAPPING
famname_mapped = at.load_mat.FAMILY_MAPPING
passmethod_mapped = at.load_mat.PASSMETHOD_MAPPING


def test_invalid_class_raises_AttributeError():
    elem_kwargs = {'Class': 'Invalid'}
    with pytest.raises(AttributeError):
        find_class_name(elem_kwargs)


def test_correct_class_names():
    elem_kwargs = {'FamName': 'fam'}
    for class_name in class_list:
        elem_kwargs['Class'] = class_name
        assert find_class_name(elem_kwargs) == class_name


def test_class_mapping():
    elem_kwargs = {'FamName': 'fam'}
    for class_name in class_mapped.keys():
        elem_kwargs['Class'] = class_name
        assert find_class_name(elem_kwargs) == class_mapped[class_name]


def test_family_mapping():
    elem_kwargs = {}
    for family_name in class_list:
        elem_kwargs['FamName'] = family_name
        assert find_class_name(elem_kwargs) == family_name
    for family_name in famname_mapped.keys():
        elem_kwargs['FamName'] = family_name
        assert find_class_name(elem_kwargs) == famname_mapped[family_name]
        elem_kwargs['FamName'] = family_name.upper()
        assert find_class_name(elem_kwargs) == famname_mapped[family_name]


def test_PassMethod_mapping():
    elem_kwargs = {'FamName': 'fam'}
    for pass_method in passmethod_mapped.keys():
        elem_kwargs['PassMethod'] = pass_method
        assert find_class_name(elem_kwargs) == passmethod_mapped[pass_method]


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


def test_find_Bend():
    elem_kwargs = {'FullGap': 0.05, 'FamName': 'fam'}
    assert find_class_name(elem_kwargs) == 'Bend'
    elem_kwargs = {'FringeInt1': 0.5, 'FamName': 'fam'}
    assert find_class_name(elem_kwargs) == 'Bend'
    elem_kwargs = {'FringeInt2': 0.5, 'FamName': 'fam'}
    assert find_class_name(elem_kwargs) == 'Bend'
    elem_kwargs = {'gK': 0.05, 'FamName': 'fam'}
    assert find_class_name(elem_kwargs) == 'Bend'
    elem_kwargs = {'EntranceAngle': 0.05, 'FamName': 'fam'}
    assert find_class_name(elem_kwargs) == 'Bend'
    elem_kwargs = {'ExitAngle': 0.05, 'FamName': 'fam'}
    assert find_class_name(elem_kwargs) == 'Bend'
    elem_kwargs = {'BendingAngle': 0.1, 'PolynomB': [0, 0, 0, 0],
                   'FamName': 'fam'}
    assert find_class_name(elem_kwargs) == 'Bend'


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
    elem_kwargs = {'PolynomB': [0, -0.5, 0, 0], 'FamName': 'fam'}
    assert find_class_name(elem_kwargs) == 'Quadrupole'


def test_find_Multipole():
    elem_kwargs = {'PolynomB': [0, 0, 0, 0],
                   'PassMethod': 'StrMPoleSymplectic4Pass', 'FamName': 'fam'}
    assert find_class_name(elem_kwargs) == 'Multipole'
    elem_kwargs = {'PolynomB': [0, 0, 0, 0, 1], 'PolynomA': [0, 0, 0, 0],
                   'Length': 0, 'FamName': 'fam'}
    assert find_class_name(elem_kwargs) == 'Multipole'


def test_find_Drift():
    elem_kwargs = {'PolynomB': [0, 0, 0, 0], 'BendingAngle': 0.0,
                   'FamName': 'fam'}
    assert find_class_name(elem_kwargs) == 'Drift'
    elem_kwargs = {'PolynomB': [0, 0, 0, 0], 'FamName': 'fam'}
    assert find_class_name(elem_kwargs) == 'Drift'
    elem_kwargs = {'Length': 1.0, 'FamName': 'fam'}
    assert find_class_name(elem_kwargs) == 'Drift'


def test_find_Sextupole():
    elem_kwargs = {'PolynomB': [0, 0, 1, 0], 'FamName': 'fam'}
    assert find_class_name(elem_kwargs) == 'Sextupole'
    elem_kwargs = {'PolynomB': [0, 0, 0, 0, 1], 'PolynomA': [0, 1, 0, 0],
                   'FamName': 'fam'}
    assert find_class_name(elem_kwargs) == 'Sextupole'


def test_find_Octupole():
    elem_kwargs = {'PolynomB': [0, 0, 0, 1], 'PolynomA': [0, 0, 0, 0],
                   'FamName': 'fam'}
    assert find_class_name(elem_kwargs) == 'Octupole'
    elem_kwargs = {'PolynomB': [0, 0, 0, 0, 1], 'PolynomA': [0, 0, 0, 1],
                   'FamName': 'fam'}
    assert find_class_name(elem_kwargs) == 'Octupole'


def test_find_ThinMultipole():
    elem_kwargs = {'PolynomB': [0, 0, 0, 0, 1], 'PolynomA': [0, 0, 0, 0],
                   'FamName': 'fam'}
    assert find_class_name(elem_kwargs) == 'ThinMultipole'


def test_find_Dipole():
    elem_kwargs = {'BendingAngle': 0.1, 'FamName': 'fam'}
    assert find_class_name(elem_kwargs) == 'Dipole'


def test_find_Marker():
    elem_kwargs = {'Length': 0.0, 'FamName': 'fam'}
    assert find_class_name(elem_kwargs) == 'Marker'
    elem_kwargs = {'FamName': 'fam'}
    assert find_class_name(elem_kwargs) == 'Marker'


def test_sanitise_class():
    elem_kwargs = {'PassMethod': 'IdentityPass', 'Class': 'Drift'}
    sanitise_class(elem_kwargs)
    assert elem_kwargs['Class'] == 'Monitor'
    elem_kwargs = {'PassMethod': 'CavityPass', 'Class': 'Drift'}
    with pytest.raises(AttributeError):
        sanitise_class(elem_kwargs)
    elem_kwargs = {'PassMethod': 'Invalid', 'Class': 'Marker'}
    with pytest.raises(AttributeError):
        sanitise_class(elem_kwargs)
    elem_kwargs = {'PassMethod': 'Invalid', 'Class': 'Monitor'}
    with pytest.raises(AttributeError):
        sanitise_class(elem_kwargs)
    elem_kwargs = {'PassMethod': 'Invalid', 'Class': 'Drift'}
    with pytest.raises(AttributeError):
        sanitise_class(elem_kwargs)
    elem_kwargs = {'PassMethod': 'Invalid', 'Class': 'RingParam'}
    with pytest.raises(AttributeError):
        sanitise_class(elem_kwargs)
