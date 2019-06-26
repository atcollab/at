import sys
import numpy
import pytest
from at import elements
from at.lattice import Lattice, AtWarning, AtError


def test_lattice_creation_gets_attributes_from_arguments():
    lat = Lattice(name='lattice', energy=3.e+6, periodicity=32, an_attr=12)
    assert len(lat) == 0
    assert lat.name == 'lattice'
    assert lat.energy == 3.e+6
    assert lat.periodicity == 32
    assert lat._radiation is False
    assert lat.an_attr == 12


def test_lattice_creation_short_scan_reads_radiation_status_correctly():
    d = elements.Dipole('d1', 1, BendingAngle=numpy.pi, Energy=5.e+6,
                        PassMethod='BndMPoleSymplectic4RadPass')
    lat = Lattice([d], name='lattice', energy=3.e+6, periodicity=32)
    assert len(lat) == 1
    assert lat.name == 'lattice'
    assert lat.energy == 3.e+6
    assert lat.periodicity == 32
    assert lat._radiation is True


def test_lattice_creation_from_lattice_inherits_attributes():
    d = elements.Dipole('d1', 1, BendingAngle=numpy.pi, Energy=5.e+6,
                        PassMethod='BndMPoleSymplectic4RadPass')
    lat1 = Lattice([d], name='lattice', energy=3.e+6, periodicity=32,
                   an_attr=12)
    lat1.another_attr = 5
    lat2 = Lattice(lat1)
    assert id(lat1) != id(lat2)
    assert len(lat2) == 1
    assert lat2.name == 'lattice'
    assert lat2.energy == 3.e+6
    assert lat2.periodicity == 32
    assert lat2._radiation is True
    assert lat2.an_attr == 12
    assert lat2.another_attr == 5


def test_lattice_energy_is_not_defined_raises_AtError():
    with pytest.raises(AtError):
        Lattice()


def test_item_is_not_an_AT_element_warns_correctly():
    with pytest.warns(AtWarning):
        Lattice(['a'], energy=0, periodicity=1)


def test_lattice_string_ordering():
    lat = Lattice([elements.Drift('D0', 1.0, attr1=numpy.array(0))],
                  name='lat', energy=5, periodicity=1, attr2=3)
    # Default dictionary ordering is only in Python >= 3.6
    if sys.version_info < (3, 6):
        assert lat.__str__().startswith("Lattice(<1 elements>, ")
        assert lat.__str__().endswith(", attr2=3)")
        assert lat.__repr__().startswith("Lattice([Drift('D0', 1.0, "
                                         "attr1=array(0))], ")
        assert lat.__repr__().endswith(", attr2=3)")
    else:
        assert lat.__str__() == ("Lattice(<1 elements>, name='lat', energy=5,"
                                 " periodicity=1, attr2=3)")
        assert lat.__repr__() == ("Lattice([Drift('D0', 1.0, attr1=array(0))],"
                                  " name='lat', energy=5, periodicity=1,"
                                  " attr2=3)")


def test_getitem(simple_lattice, simple_ring):
    assert simple_lattice[-1] == simple_ring[-1]
    bool_refs = numpy.array([False, True, False, True, False, True])
    assert simple_lattice[bool_refs] == simple_ring[1::2]


def test_setitem(simple_lattice):
    new = elements.Monitor('M2')
    old = simple_lattice[5]
    simple_lattice[5] = new
    assert simple_lattice[5] != old
    assert simple_lattice[5] == new
    bool_refs = numpy.array([False, False, False, False, False, True])
    simple_lattice[bool_refs] = old
    assert simple_lattice[5] != new
    assert simple_lattice[5] == old


def test_delitem(simple_lattice, simple_ring):
    mon = elements.Monitor('M2')
    simple_lattice.append(mon)
    assert len(simple_lattice) == 7
    del simple_lattice[-1]
    assert simple_lattice[:] == simple_ring
    bool_refs = numpy.array([False, False, False, False, False, False, True])
    simple_lattice.append(mon)
    assert len(simple_lattice) == 7
    del simple_lattice[bool_refs]
    assert simple_lattice[:] == simple_ring


def test_copy(hmba_lattice):
    assert id(hmba_lattice.copy()) != id(hmba_lattice)
    assert id(hmba_lattice.copy()[0]) == id(hmba_lattice[0])


def test_deepcopy(hmba_lattice):
    assert id(hmba_lattice.deepcopy()) != id(hmba_lattice)
    assert id(hmba_lattice.deepcopy()[0]) != id(hmba_lattice[0])


def test_property_values_against_known(hmba_lattice):
    assert hmba_lattice.voltage == 6000000
    assert hmba_lattice.harmonic_number == 992
    assert hmba_lattice.radiation is False
    numpy.testing.assert_almost_equal(hmba_lattice.energy_loss,
                                      2526188.713461808)


def test_radiation_change(hmba_lattice):
    rfs = [elem for elem in hmba_lattice if isinstance(elem,
                                                       elements.RFCavity)]
    dipoles = [elem for elem in hmba_lattice if isinstance(elem,
                                                           elements.Dipole)]
    quads = [elem for elem in hmba_lattice if isinstance(elem,
                                                         elements.Quadrupole)]
    hmba_lattice.radiation_on(None, 'pass2', 'auto')
    assert hmba_lattice.radiation is True
    for elem in rfs:
        assert elem.PassMethod == 'IdentityPass'
    for elem in dipoles:
        assert elem.PassMethod == 'pass2'
    for elem in quads:
        assert elem.PassMethod == 'StrMPoleSymplectic4RadPass'
    hmba_lattice.radiation_off(None, 'BndMPoleSymplectic4Pass', 'auto')
    assert hmba_lattice.radiation is False
    for elem in rfs:
        assert elem.PassMethod == 'IdentityPass'
    for elem in dipoles:
        assert elem.PassMethod == 'BndMPoleSymplectic4Pass'
    for elem in quads:
        assert elem.PassMethod == 'StrMPoleSymplectic4Pass'


def test_radiation_state_errors(hmba_lattice):
    hmba_lattice.radiation_on()
    with pytest.raises(AtError):
        hmba_lattice.linopt()
    hmba_lattice.radiation_off()
    hmba_lattice.linopt()
    with pytest.raises(AtError):
        hmba_lattice.ohmi_envelope()
    hmba_lattice.radiation_on()
    hmba_lattice.ohmi_envelope()
    with pytest.raises(AtError):
        hmba_lattice.get_mcf()
    hmba_lattice.radiation_off()
    hmba_lattice.get_mcf()
