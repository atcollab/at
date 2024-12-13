import numpy
import pytest
from numpy.testing import assert_allclose, assert_equal

from at import elements
from at.lattice import Lattice, AtWarning, AtError


def test_lattice_creation_gets_attributes_from_arguments():
    lat = Lattice(name="lattice", energy=3.0e6, periodicity=32, an_attr=12)
    assert len(lat) == 0
    assert lat.name == "lattice"
    assert lat.energy == 3.0e6
    assert lat.periodicity == 32
    assert lat.radiation is False
    assert lat.particle.name == "relativistic"
    assert lat.an_attr == 12


def test_lattice_energy_radiation_periodicity():
    d = elements.Dipole(
        "d1",
        1,
        BendingAngle=numpy.pi / 16,
        Energy=5.0e6,
        PassMethod="BndMPoleSymplectic4RadPass",
    )
    lat = Lattice([d], name="lattice", periodicity=0, energy=3.0e6)
    assert lat.energy == 3.0e6
    assert lat.periodicity == 32
    assert lat.radiation is True
    assert lat.is_6d is True


def test_lattice_voltage_harmonic_number():
    rf = elements.RFCavity("rf", 0, 0.2e6, 0.5e9, 5, 3.0e6)
    d = elements.Dipole("d1", 2.99792458, BendingAngle=numpy.pi / 5)
    lat = Lattice([rf, d], name="lattice", periodicity=0)
    assert lat.energy == 3.0e6
    assert lat.periodicity == 10
    assert lat.rf_voltage == 2e6
    assert lat.revolution_frequency == 10.0e6
    assert lat.harmonic_number == 50
    assert lat.radiation is True


def test_lattice_creation_from_lattice_inherits_attributes():
    d = elements.Dipole(
        "d1",
        1,
        BendingAngle=numpy.pi,
        Energy=5.0e6,
        PassMethod="BndMPoleSymplectic4RadPass",
    )
    lat1 = Lattice([d], name="lattice", energy=3.0e6, periodicity=32, an_attr=12)
    lat2 = Lattice(lat1, another_attr=5)
    assert id(lat1) != id(lat2)
    assert len(lat2) == 1
    assert lat2.name == "lattice"
    assert lat2.energy == 3.0e6
    assert lat2.periodicity == 32
    assert lat2.is_6d is True
    assert lat2.another_attr == 5
    with pytest.raises(AttributeError):
        assert lat2.an_attr == 12


def test_lattice_energy_is_not_defined_raises_AtError():
    d = elements.Dipole("d1", 1, BendingAngle=numpy.pi)
    with pytest.raises(AtError):
        Lattice([d])


def test_item_is_not_an_AT_element_warns_correctly():
    with pytest.warns(AtWarning):
        Lattice(["a"], energy=0, periodicity=1)


def test_lattice_string_ordering():
    lat = Lattice(
        [elements.Drift("D0", 1.0, attr1=numpy.array(0))],
        name="lat",
        energy=5,
        periodicity=1,
        attr2=3,
    )
    latstr = str(lat)
    assert latstr.startswith(
        "Lattice(<1 elements>, name='lat', "
        "energy=5, particle=Particle('relativistic'), "
        "periodicity=1, beam_current=0.0, nbunch=1"
    )
    assert latstr.endswith("attr2=3)")

    latrepr = repr(lat)
    assert latrepr.startswith(
        "Lattice([Drift('D0', 1.0, attr1=array(0))], "
        "name='lat', "
        "energy=5, particle=Particle('relativistic'), "
        "periodicity=1, beam_current=0.0, nbunch=1"
    )
    assert latrepr.endswith("attr2=3)")


def test_getitem(simple_lattice, simple_ring):
    assert simple_lattice[-1] == simple_ring[-1]
    bool_refs = numpy.array([False, True, False, True, False, True])
    assert simple_lattice[bool_refs] == simple_ring[1::2]


def test_setitem(simple_lattice):
    new = elements.Monitor("M2")
    old = simple_lattice[5]
    simple_lattice[5] = new
    assert simple_lattice[5] != old
    assert simple_lattice[5] == new
    bool_refs = numpy.array([False, False, False, False, False, True])
    simple_lattice[bool_refs] = old
    assert simple_lattice[5] != new
    assert simple_lattice[5] == old


def test_delitem(simple_lattice, simple_ring):
    mon = elements.Monitor("M2")
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
    assert hmba_lattice.rf_voltage == 6000000
    assert hmba_lattice.harmonic_number == 992
    assert hmba_lattice.radiation is False
    numpy.testing.assert_almost_equal(
        hmba_lattice.energy_loss, 2526188.658758993, decimal=1
    )


def test_radiation_change(hmba_lattice):
    rfs = [elem for elem in hmba_lattice if isinstance(elem, elements.RFCavity)]
    dipoles = [elem for elem in hmba_lattice if isinstance(elem, elements.Dipole)]
    quads = [elem for elem in hmba_lattice if isinstance(elem, elements.Quadrupole)]
    hmba_lattice.radiation_on(None, "pass2", "auto")
    assert hmba_lattice.radiation is True
    assert hmba_lattice.has_cavity is False
    for elem in rfs:
        assert elem.PassMethod == "IdentityPass"
    for elem in dipoles:
        assert elem.PassMethod == "pass2"
    for elem in quads:
        assert elem.PassMethod == "StrMPoleSymplectic4RadPass"
    hmba_lattice.radiation_off(None, "BndMPoleSymplectic4Pass", "auto")
    assert hmba_lattice.radiation is False
    for elem in rfs:
        assert elem.PassMethod == "IdentityPass"
    for elem in dipoles:
        assert elem.PassMethod == "BndMPoleSymplectic4Pass"
    for elem in quads:
        assert elem.PassMethod == "StrMPoleSymplectic4Pass"


def test_radiation_state_errors(hmba_lattice):
    hmba_lattice.radiation_on()
    with pytest.raises(AtError):
        hmba_lattice.linopt()
    hmba_lattice.radiation_off()
    hmba_lattice.linopt()
    with pytest.raises(AtError):
        hmba_lattice.ohmi_envelope()
    hmba_lattice.enable_6d()
    hmba_lattice.ohmi_envelope()
    with pytest.raises(AtError):
        hmba_lattice.get_mcf()
    hmba_lattice.disable_6d()
    hmba_lattice.get_mcf()


@pytest.mark.parametrize("ring", ["hmba_lattice"])
def test_develop(request, ring):
    ring = request.getfixturevalue(ring)
    newring = ring.develop()
    assert_allclose(ring.circumference, newring.circumference)
    assert_allclose(ring.rf_voltage, newring.rf_voltage)
    assert_allclose(
        ring.revolution_frequency * ring.harmonic_number,
        newring.revolution_frequency * newring.harmonic_number,
    )
    rp1 = ring.radiation_parameters()
    rp2 = newring.radiation_parameters()
    assert_allclose(rp1.fulltunes, rp2.fulltunes)
    assert_allclose(rp1.U0, rp2.U0)


@pytest.mark.parametrize("ring", ["hmba_lattice"])
def test_operators(request, ring):
    ring = request.getfixturevalue(ring)
    newring1 = ring + ring
    newring2 = ring.concatenate(ring, copy=True)
    assert_equal(len(newring1), len(ring) * 2)
    assert_equal(len(newring2), len(ring) * 2)
    newring1 += ring
    newring2.concatenate(ring, copy=False)
    assert_equal(len(newring1), len(ring) * 3)
    assert_equal(len(newring2), len(ring) * 3)
    newring1 = ring * 2
    newring2 = ring.repeat(2)
    assert_equal(len(newring1), len(ring) * 2)
    assert_equal(len(newring2), len(ring) * 2)
    assert_equal(newring1.harmonic_number, ring.harmonic_number)
    newring1 = ring.reverse(copy=True)
    assert_equal(newring1[-1].FamName, ring[0].FamName)
