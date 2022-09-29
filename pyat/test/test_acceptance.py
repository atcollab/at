import at
import numpy
from numpy.testing import assert_allclose


def test_1d_acceptance(hmba_lattice):
    hmba_lattice = hmba_lattice.radiation_off(copy=True)
    acceptance, _, _ = hmba_lattice.get_horizontal_acceptance(1e-3, 1.0e-3)
    expected = numpy.array([-0.011, 0.011])
    assert_allclose(acceptance, expected, atol=1e-6)


def test_1d_acceptance_refpts(hmba_lattice):
    hmba_lattice = hmba_lattice.radiation_off(copy=True)
    acceptance, _, _ = hmba_lattice.get_horizontal_acceptance(1e-3, 1.0e-3,
                                                              refpts=10)
    expected = numpy.array([-0.002, 0.002])
    assert_allclose(acceptance, expected, atol=1e-6)


def test_2d_acceptance(hmba_lattice):
    hmba_lattice = hmba_lattice.radiation_off(copy=True)
    acceptance, _, _ = hmba_lattice.get_acceptance(['x', 'y'], [10, 5],
                                                   [10.0e-3, 10.0e-3],
                                                   grid_mode=at.GridMode.RECURSIVE)
    expected = numpy.array([[-0.01125, -0.0053033, 0., 0.00574524, 0.010625],
                            [0., 0.0053033, 0.006875, 0.00574524, 0.]])
    assert_allclose(acceptance, expected, atol=1e-6, rtol=1.0e-4)
