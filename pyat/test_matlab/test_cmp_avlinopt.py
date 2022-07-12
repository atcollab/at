import at
import matlab
import numpy
from numpy.testing import assert_allclose as assert_close
import pytest


def _ml_refs(refpts, nelems):
    """Convert refpoints to Matlab"""
    uintrefs = at.uint32_refpts(refpts, nelems)
    return matlab.double([float(ref+1) for ref in uintrefs])


@pytest.mark.parametrize('dp', (-0.01, 0.0, 0.01))
@pytest.mark.parametrize('lattices',
                         [pytest.lazy_fixture('dba'),
                          pytest.lazy_fixture('hmba')])
def test_avlinopt(engine, lattices, dp):
    py_lattice, ml_lattice, _ = lattices
    nelems = len(py_lattice)
    refpts = range(nelems)
    mlrefs = _ml_refs(refpts, nelems)

    # Python call
    lindata, avebeta, avemu, avedisp, aves, tune, chrom = \
        at.avlinopt(py_lattice, dp, refpts=refpts)
    # Matlab call
    ml_data, ml_avebeta, ml_avemu, ml_avedisp, ml_tune, ml_chrom = \
        engine.pyproxy('atavedata', ml_lattice, dp, mlrefs, nargout=6)
    # Comparison
    assert_close(avebeta, numpy.asarray(ml_avebeta), rtol=1e-7)
    assert_close(avemu, numpy.asarray(ml_avemu), rtol=0, atol=1e-7)
    assert_close(avedisp, numpy.asarray(ml_avedisp), rtol=0, atol=1e-8)
