import pytest
from numpy.testing import assert_allclose as assert_close

from at import linopt2, linopt4, linopt6, get_optics


@pytest.mark.parametrize(
    "lattice",
    ["dba_lattice", "hmba_lattice", "noringparam_lattice"],
)
def test_linopt6_norad(request, lattice):
    """Compare the results of linopt2 and linopt6 in 4d"""
    lattice = request.getfixturevalue(lattice)
    refpts = range(len(lattice) + 1)
    ld02, rd2, ld2 = linopt2(lattice, refpts, get_w=True)
    ld06, rd6, ld6 = linopt6(lattice, refpts, get_w=True)
    assert_close(rd2.tune, rd6.tune, atol=1e-10, rtol=0)
    assert_close(rd2.chromaticity, rd6.chromaticity, atol=1e-10, rtol=0)

    for field in ["s_pos", "closed_orbit", "dispersion", "alpha", "beta", "mu"]:
        assert_close(ld2[field], ld6[field], atol=1e-10, rtol=0, err_msg=field)
    assert_close(ld2.W, ld6.W, atol=1e-6, rtol=0)


@pytest.mark.parametrize(
    "lattice",
    ["hmba_lattice", "noringparam_lattice"],
)
def test_linopt6_rad(request, lattice):
    """Compare the results with and without radiation"""
    lattice = request.getfixturevalue(lattice)
    refpts = range(len(lattice) + 1)
    # Turn cavity ON, without radiation
    radlattice = lattice.enable_6d(dipole_pass=None, copy=True)
    ld04, rd4, ld4 = linopt6(lattice, refpts, get_w=True)
    ld06, rd6, ld6 = linopt6(radlattice, refpts, get_w=True)

    assert_close(rd4.tune, rd6.tune[:2], atol=1e-10, rtol=0)
    assert_close(rd4.chromaticity, rd6.chromaticity[:2], atol=1e-8, rtol=0)

    for field in ["s_pos", "closed_orbit", "dispersion", "alpha", "beta"]:
        assert_close(ld4[field], ld6[field], atol=1.0e-8, rtol=0, err_msg=field)
    assert_close(ld4.mu, ld6.mu[:, :2], atol=1.0e-8, rtol=0)
    assert_close(ld4.W, ld6.W, atol=1e-6, rtol=2e-5)


@pytest.mark.parametrize("dp", (-0.005, 0.0, 0.005))
@pytest.mark.parametrize("lattice", ["dba_lattice", "hmba_lattice"])
@pytest.mark.parametrize("method", [linopt2, linopt4, linopt6])
def test_linopt6_line(request, lattice, dp, method):
    """Run a ring as a transfer line and compare with the ring results"""
    lattice = request.getfixturevalue(lattice)
    refpts = lattice.uint32_refpts(range(len(lattice) + 1))
    ld04, rd4, ld4 = get_optics(lattice, refpts, dp=dp, method=method)

    twin = {
        "alpha": ld04.alpha,
        "beta": ld04.beta,
        "dispersion": ld04.dispersion,
        "closed_orbit": ld04.closed_orbit,
    }
    # twin = ld04
    tr04, bd4, tr4 = get_optics(lattice, refpts, dp=dp, twiss_in=twin, method=method)
    for field in ["s_pos", "closed_orbit", "dispersion", "alpha", "beta", "mu"]:
        assert_close(ld04[field], tr04[field], atol=1.0e-10, rtol=0, err_msg=field)
        assert_close(ld4[field], tr4[field], atol=1.0e-8, rtol=0, err_msg=field)
