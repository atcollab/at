import os
from tempfile import mktemp

import pytest
from numpy.testing import assert_equal

from at.lattice import Lattice


@pytest.mark.parametrize("lattice", ["dba_lattice", "hmba_lattice"])
@pytest.mark.parametrize(
    "suffix, options",
    [
        (".m", {}),
        (".repr", {}),
        (".mat", {"use": "abcd"}),
        (".json", {}),
    ],
)
def test_m(request, lattice, suffix, options):
    ring0 = request.getfixturevalue(lattice)
    fname = mktemp(suffix=suffix)

    # Create a new .m or .repr file
    ring0.save(fname, **options)

    # load the new file
    ring1 = Lattice.load(fname, **options)

    # Check that we get back the original lattice
    el1, rg1, _ = ring0.linopt6()
    el2, rg2, _ = ring1.linopt6()
    assert_equal(rg1.tune, rg2.tune)
    assert_equal(rg1.chromaticity, rg2.chromaticity)

    os.unlink(fname)
