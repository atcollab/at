from at import linopt2, linopt6
from numpy.testing import assert_allclose as assert_close
import pytest

@pytest.mark.parametrize('lattice',
                         [pytest.lazy_fixture('dba_lattice'),
                          pytest.lazy_fixture('hmba_lattice')])
def test_linopt6_norad(lattice):
    """Compare the results of linopt2 and linopt6 in 4d"""
    refpts = range(len(lattice) + 1)
    ld02, rd2, ld2 = linopt2(lattice, refpts, get_w=True)
    ld06, rd6, ld6 = linopt6(lattice, refpts, get_w=True)
    assert_close(rd2.tune, rd6.tune, atol=1e-12, rtol=0)
    assert_close(rd2.chromaticity, rd6.chromaticity, atol=1e-12, rtol=0)

    for field in ['s_pos', 'closed_orbit', 'dispersion', 'alpha', 'beta', 'mu']:
        assert_close(ld2[field], ld6[field], atol=1e-10, rtol=0)
    assert_close(ld2.W, ld6.W, atol=1e-6, rtol=0)


@pytest.mark.parametrize('lattice', [pytest.lazy_fixture('hmba_lattice')])
def test_linopt6_rad(lattice):
    """Compare the results with and without radiation"""
    refpts=range(len(lattice)+1)
    radlattice = lattice.radiation_on(dipole_pass=None, copy=True)
    ld04, rd4, ld4 = linopt6(lattice, refpts, get_w=True)
    ld06, rd6, ld6 = linopt6(radlattice, refpts, get_w=True)

    assert_close(rd4.tune, rd6.tune[:2], atol=1e-10, rtol=0)
    assert_close(rd4.chromaticity, rd6.chromaticity[:2], atol=1e-8, rtol=0)

    for field in ['s_pos', 'closed_orbit', 'dispersion', 'alpha', 'beta', 'W']:
        assert_close(ld4[field], ld6[field], atol=1.e-6, rtol=0)
    assert_close(ld4.mu, ld6.mu[:, :2], atol=1.e-8, rtol=0)
