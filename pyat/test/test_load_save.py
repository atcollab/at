"""pyat test."""

import os
import sys
from pathlib import Path
from tempfile import mkstemp

import machine_data
import pytest
from at.lattice import Lattice
from at.lattice import elements as elt
from at.lattice.elements.idtable_element import InsertionDeviceKickMap
from numpy.testing import assert_allclose, assert_equal

if sys.version_info < (3, 9):
    from importlib_resources import files
else:
    from importlib.resources import files


@pytest.fixture()
def simple_hmba(hmba_lattice: Lattice) -> None:
    """Modify hmba_lattice to make it compatible with MAD-X and Elegant.

    Arguments:
        hmba_lattice: AT Lattice.

    Returns:
        Lattice: AT lattice
    """
    ring = hmba_lattice.deepcopy()
    # Set NumIntSteps to default, remove FringeQuadEntrance and FringeQuadExit
    for mult in ring.select(elt.Multipole):
        mult.NumIntSteps = 10
        if isinstance(mult, (elt.Dipole, elt.Quadrupole)):
            del mult.FringeQuadEntrance
            del mult.FringeQuadExit
    # Replace Multipoles by Octupoles
    for idx in ring.get_uint32_index("OF*"):
        _q = ring[idx]
        ring[idx] = elt.Octupole(_q.FamName, _q.Length, _q.PolynomA, _q.PolynomB)
    # Disable useless e;ements
    for idx in ring.get_uint32_index("SH*"):
        _q = ring[idx]
        ring[idx] = elt.Drift(_q.FamName, _q.Length)
    return ring


@pytest.mark.parametrize("lattice", ["dba_lattice", "simple_hmba"])
@pytest.mark.parametrize(
    "suffix, options",
    (
        (".m", {}),
        (".repr", {}),
        (".mat", {"use": "abcd"}),
        (".mat", {"mat_key": "efgh"}),
        (".json", {}),
        (".seq", {"use": "ring"}),
        (".lte", {"use": "ring"}),
    ),
)
def test_m(request, lattice, suffix, options) -> None:
    """Test m."""
    ring0 = request.getfixturevalue(lattice)
    fhandle, fname = mkstemp(suffix=suffix)
    os.close(fhandle)

    # Create a new file
    ring0.save(fname, **options)

    # load the new file
    ring1 = Lattice.load(fname, **options)

    # Check that we get back the original lattice
    _, rg1, _ = ring0.linopt6()
    _, rg2, _ = ring1.linopt6()
    assert_allclose(rg1.tune, rg2.tune, atol=1.0e-12)
    assert_allclose(rg1.chromaticity, rg2.chromaticity, atol=1.0e-12)

    # delete temporary file
    temp_file = Path(fname)
    temp_file.unlink()


def test_long_arrays_in_m_file() -> None:
    """Test long array saving in .m files."""
    # create an element with long arrays
    thepath = files(machine_data).as_posix()
    theidfname = "/kickmap_w150_20mm.txt"
    elem = InsertionDeviceKickMap("idmap", 10, thepath + theidfname, 6.04)

    # save a ring with one element into a temporary file
    ring0 = Lattice([elem], energy=6.04e9)
    fhandle, fname = mkstemp(suffix=".m")
    os.close(fhandle)
    ring0.save(fname)

    # retrieve the ring
    ring1 = Lattice.load(fname)

    assert_equal(ring1[0].xkick.shape, ring0[0].xkick.shape)  # act

    # delete temporary file
    temp_file = Path(fname)
    temp_file.unlink()
