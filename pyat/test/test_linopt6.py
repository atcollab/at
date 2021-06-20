from at import linopt, linopt6
from numpy.testing import assert_allclose as assert_close
import pytest

@pytest.mark.parametrize('dp', (0.0, 0.01))
@pytest.mark.parametrize('lattice',
                         [pytest.lazy_fixture('dba_lattice'),
                          pytest.lazy_fixture('hmba_lattice')])
def test_linopt6(lattice, dp):
    """Compare the results of linopt and linopt6 in 4d"""
    refpts=range(21)
    lindata0, tune, chrom, lindata = linopt(lattice, dp, refpts,
                                                    get_chrom=True)
    ld0, rd, ld = linopt6(lattice, refpts, dp=dp, get_chrom=True)
    assert_close(tune, rd.tune, atol=1e-12)
    assert_close(chrom, rd.chromaticity, atol=1e-8)

    assert_close(ld.beta, lindata.beta,  rtol=1e-6)
    assert_close(ld.alpha, lindata.alpha,  atol=1e-6, rtol=1e-6)
    assert_close(ld.dispersion, lindata.dispersion,  atol=1e-9)

