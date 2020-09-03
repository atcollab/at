import at
import matlab
import numpy
import pytest


def _ml_refs(refpts, nelems):
    """Convert refpoints to Matlab"""
    uintrefs = at.uint32_refpts(refpts, nelems)
    # noinspection PyUnresolvedReferences
    return matlab.double([ref+1 for ref in uintrefs])


# noinspection PyUnresolvedReferences
@pytest.mark.parametrize('dp', (-0.01, 0.0, 0.01))
@pytest.mark.parametrize('ml_lattice, py_lattice',
                         [(pytest.lazy_fixture('ml_hmba'),
                           pytest.lazy_fixture('py_hmba'))])
def test_avlinopt(engine, ml_lattice, py_lattice, dp):
    """N.B. a 'mu' comparison is left out for twiss data as the values for 'mu'
        returned by 'twissring' in Matlab are inconsistent with those from
        'get_twiss' and 'linopt' in Python as well as those returned from
        'atlinopt' in Matlab.
    """
    nelems = len(py_lattice)
    refpts = range(nelems)
    mlrefs = _ml_refs(refpts, nelems)

    # Python call
    lindata, avebeta, avemu, avedisp, aves, tune, chrom = \
        at.avlinopt(py_lattice, dp, refpts=refpts, ddp=1.e-6)
    # Matlab call
    ml_data, ml_avebeta, ml_avemu, ml_avedisp, ml_tune, ml_chrom = \
        engine.pyproxy('atavedata', ml_lattice, dp, mlrefs, nargout=6)
    # Comparison
    numpy.testing.assert_allclose(avebeta, numpy.asarray(ml_avebeta), rtol=1.0e-11, atol=1.0e-20)
    numpy.testing.assert_allclose(avemu, numpy.asarray(ml_avemu), rtol=1.0e-11, atol=1.0e-20)
    numpy.testing.assert_allclose(avedisp, numpy.asarray(ml_avedisp), rtol=1.0e-11, atol=1.0e-10)
