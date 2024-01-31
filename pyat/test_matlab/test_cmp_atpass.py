import matlab
import numpy
import pytest
from numpy.testing import assert_allclose as assert_close

import at


@pytest.mark.parametrize(
    "lattices",
    [
        "dba",
        "hmba",
        "hmba_cav",
        "hmba_rad",
    ],
)
def test_atpass(engine, request, lattices):
    py_lattice, ml_lattice, _ = request.getfixturevalue(lattices)
    xy_step = 1.0e-8
    scaling = [xy_step, xy_step, xy_step, xy_step, xy_step, xy_step]

    ml_rin = engine.diag(matlab.double(scaling))
    ml_rout = engine.atpass(ml_lattice, ml_rin, 1, 1)

    py_rin = numpy.asfortranarray(numpy.diag(scaling))
    at.lattice_track(py_lattice, py_rin, 1, in_place=True)

    assert_close(py_rin, ml_rout, rtol=0, atol=1.0e-30)


@pytest.mark.parametrize(
    "lattices",
    [
        "dba",
        "hmba",
        "hmba_cav",
        "hmba_rad",
    ],
)
def test_linepass(engine, request, lattices):
    py_lattice, ml_lattice, _ = request.getfixturevalue(lattices)
    xy_step = 1.0e-8
    scaling = [xy_step, xy_step, xy_step, xy_step, xy_step, xy_step]

    ml_rin = engine.diag(matlab.double(scaling))
    ml_rout = numpy.array(engine.linepass(ml_lattice, ml_rin))

    py_rin = numpy.asfortranarray(numpy.diag(scaling))
    py_rout, *_ = at.lattice_track(py_lattice, py_rin, refpts=len(py_lattice))
    py_rout = numpy.squeeze(py_rout)

    assert_close(py_rout, ml_rout, rtol=0, atol=1.0e-30)
