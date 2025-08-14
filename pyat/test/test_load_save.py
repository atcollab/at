import os
from pathlib import Path
from tempfile import mkstemp

import pytest
from at.lattice import Lattice
from at.lattice import elements as elt
from at.lattice.elements.idtable_element import InsertionDeviceKickMap
from numpy.testing import assert_allclose, assert_equal


@pytest.fixture
def simple_hmba(hmba_lattice):
    """Modify hmba_lattice to make it compatible with MAD-X and Elegant"""
    ring = hmba_lattice.deepcopy()
    # Set NumIntSteps to default, remove FringeQuadEntrance and FringeQuadExit
    for mult in ring.select(elt.Multipole):
        mult.NumIntSteps = 10
        if isinstance(mult, (elt.Dipole, elt.Quadrupole)):
            del mult.FringeQuadEntrance
            del mult.FringeQuadExit
    # Replace Multipoles by Octupoles
    for id in ring.get_uint32_index("OF*"):
        q = ring[id]
        ring[id] = elt.Octupole(q.FamName, q.Length, q.PolynomA, q.PolynomB)
    # Disable useless e;ements
    for id in ring.get_uint32_index("SH*"):
        q = ring[id]
        ring[id] = elt.Drift(q.FamName, q.Length)
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
def test_m(request, lattice, suffix, options):
    ring0 = request.getfixturevalue(lattice)
    fname = mkstemp(suffix=suffix)[1]

    # Create a new file
    ring0.save(fname, **options)

    # load the new file
    ring1 = Lattice.load(fname, **options)

    # Check that we get back the original lattice
    el1, rg1, _ = ring0.linopt6()
    el2, rg2, _ = ring1.linopt6()
    assert_allclose(rg1.tune, rg2.tune, atol=1.0e-12)
    assert_allclose(rg1.chromaticity, rg2.chromaticity, atol=1.0e-12)

    os.unlink(fname)


def test_long_arrays_in_m_file() -> None:
    """Test long array saving in .m files."""
    # create an element with long arrays
    elem = InsertionDeviceKickMap(
        "idmap", 10, "../../machine_data/IDs/kickmap_w150_20mm.txt", 6.04
    )

    # save a ring with one element into a temporary file
    ring0 = Lattice([elem], energy=6.04e9)
    fname = mkstemp(suffix=".m")[1]
    ring0.save(fname)

    # retrieve the ring
    ring1 = Lattice.load(fname)

    assert_equal(ring1[0].xkick.shape, ring0[0].xkick.shape)  # act

    # delete temporary file
    temp_file = Path(fname)
    temp_file.unlink()
