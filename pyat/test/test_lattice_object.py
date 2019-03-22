import pytest
import numpy
from at import elements
from at.lattice import Lattice, AtWarning, AtError


def test_lattice_creation():
    l = Lattice(name='lattice', energy=3.e+6, periodicity=32)
    assert len(l) == 0
    assert l.name == 'lattice'
    assert l.energy == 3.e+6
    assert l.periodicity == 32
    assert l._radiation is False
    d = elements.Dipole('d1', 1, BendingAngle=numpy.pi, Energy=5.e+6,
                        PassMethod='BndMPoleSymplectic4RadPass')
    l = Lattice([d], name='lattice', energy=3.e+6, periodicity=32)
    assert len(l) == 1
    assert l.name == 'lattice'
    assert l.energy == 3.e+6
    assert l.periodicity == 32
    assert l._radiation is True
    l.an_attr = 12
    lat = Lattice(l, another_attr=5)
    assert id(l) != id(lat)
    assert len(lat) == 1
    assert lat.name == 'lattice'
    assert lat.energy == 3.e+6
    assert lat.periodicity == 32
    assert lat._radiation is True
    assert lat.an_attr == 12
    assert lat.another_attr == 5


def test_lattice_gets_attributes_from_RingParam():
    rp1 = elements.RingParam('lattice_name', 3.e+6, Periodicity=32)
    rp2 = elements.RingParam('not_lattice_name', -12, Periodicity=0)
    with pytest.warns(AtWarning):
        l = Lattice([rp1, rp2])
    assert len(l) == 0
    assert l.name == 'lattice_name'
    assert l.energy == 3.e+6
    assert l.periodicity == 32
    assert l._radiation is False
    assert len(Lattice([rp1], keep_all=True)) == 1


def test_lattice_gets_attributes_from_elements():
    d = elements.Dipole('d1', 1, BendingAngle=numpy.pi, Energy=3.e+6,
                        PassMethod='BndMPoleSymplectic4RadPass')
    l = Lattice([d])
    assert len(l) == 1
    assert l.name == ''
    assert l.energy == 3.e+6
    assert l.periodicity == 2
    assert l._radiation is True


def test_lattice_creation_warnings_and_errors():
    # Lattice energy is not defined
    with pytest.raises(AtError):
        Lattice()
    # No bending in the cell, set "Periodicity" to 1
    with pytest.warns(AtWarning):
        Lattice([], energy=0)
    # item 0 (a) is not an AT element: ignored
    with pytest.warns(AtWarning):
        Lattice(['a'], energy=0, periodicity=1)
    # Non-integer number of cells: 12.5663706144 -> 13
    d = elements.Dipole('d1', 1, BendingAngle=0.5)
    with pytest.warns(AtWarning):
        Lattice([d], energy=0)
    # Inconsistent energy values, "energy" set to 5.0
    m1 = elements.Marker('m1', Energy=5)
    m2 = elements.Marker('m2', Energy=3)
    with pytest.warns(AtWarning):
        Lattice([m1, m2], periodicity=1)
    # >1 RingParams already tested
