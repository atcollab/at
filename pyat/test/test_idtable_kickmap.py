"""Tests for InsertionDeviceKickMap multi-kickmap store interface."""

from __future__ import annotations

from importlib.resources import files

import machine_data
import pytest
from at.lattice.elements.idtable_element import InsertionDeviceKickMap
from numpy.testing import assert_array_equal


@pytest.fixture()
def idkm_file() -> str:
    """Path to the bundled test kickmap text file."""
    return files(machine_data).as_posix() + "/kickmap_w150_20mm.txt"


@pytest.fixture()
def idkm_elem(idkm_file: str) -> InsertionDeviceKickMap:
    """A base InsertionDeviceKickMap element created from the test file."""
    return InsertionDeviceKickMap("idmap", 10, idkm_file, 6.04)


# ---------------------------------------------------------------------------
# Multi-kickmap store
# ---------------------------------------------------------------------------

class TestKickmapStore:
    """Tests for add_kickmap / use_kickmap / list_kickmaps / active_kickmap."""

    def test_list_kickmaps_empty_on_new_element(self, idkm_elem):
        """A freshly constructed element has only the 'default' entry in the store."""
        assert idkm_elem.list_kickmaps() == ["default"]

    def test_active_kickmap_none_on_new_element(self, idkm_elem):
        """active_kickmap is 'default' immediately after construction."""
        assert idkm_elem.active_kickmap == "default"

    def test_default_kickmap_matches_initial_fields(self, idkm_elem):
        """The 'default' kickmap must contain the element's construction-time arrays."""
        assert_array_equal(
            idkm_elem._kickmap_store["default"]["xkick"], idkm_elem.xkick
        )
        assert_array_equal(
            idkm_elem._kickmap_store["default"]["ykick"], idkm_elem.ykick
        )

    def test_use_kickmap_default_restores_initial_fields(self, idkm_elem, idkm_file):
        """use_kickmap('default') restores the original kick arrays after a swap."""
        xkick_initial = idkm_elem.xkick.copy()
        idkm_elem.add_kickmap("half_e", 5, idkm_file, 3.0)
        idkm_elem.use_kickmap("half_e")
        idkm_elem.use_kickmap("default")
        assert_array_equal(idkm_elem.xkick, xkick_initial)
        assert idkm_elem.active_kickmap == "default"


        """add_kickmap makes the key visible in list_kickmaps."""
        idkm_elem.add_kickmap("mode_a", 10, idkm_file, 6.04)
        assert "mode_a" in idkm_elem.list_kickmaps()

    def test_add_multiple_kickmaps_independent(self, idkm_elem, idkm_file):
        """Multiple keys are stored independently alongside 'default'."""
        idkm_elem.add_kickmap("mode_a", 10, idkm_file, 6.04)
        idkm_elem.add_kickmap("mode_b", 5, idkm_file, 3.0)
        assert set(idkm_elem.list_kickmaps()) == {"default", "mode_a", "mode_b"}

    def test_use_kickmap_sets_active_key(self, idkm_elem, idkm_file):
        """use_kickmap updates active_kickmap to the chosen key."""
        idkm_elem.add_kickmap("mode_a", 10, idkm_file, 6.04)
        idkm_elem.use_kickmap("mode_a")
        assert idkm_elem.active_kickmap == "mode_a"

    def test_use_kickmap_switches_active_key(self, idkm_elem, idkm_file):
        """active_kickmap reflects the most recently activated key."""
        idkm_elem.add_kickmap("mode_a", 10, idkm_file, 6.04)
        idkm_elem.add_kickmap("mode_b", 5, idkm_file, 3.0)
        idkm_elem.use_kickmap("mode_a")
        idkm_elem.use_kickmap("mode_b")
        assert idkm_elem.active_kickmap == "mode_b"

    def test_use_kickmap_changes_tracking_arrays(self, idkm_elem, idkm_file):
        """use_kickmap replaces xkick/ykick/xkick1/ykick1/xtable/ytable."""
        idkm_elem.add_kickmap("norm", 10, idkm_file, 6.04)
        idkm_elem.add_kickmap("half_e", 5, idkm_file, 3.0)

        idkm_elem.use_kickmap("norm")
        xkick_norm = idkm_elem.xkick.copy()

        idkm_elem.use_kickmap("half_e")
        xkick_half = idkm_elem.xkick.copy()

        # kick values are 1/E² normalised; different energies must yield different tables
        assert not (xkick_norm == xkick_half).all()

    def test_use_kickmap_changes_nslice(self, idkm_elem, idkm_file):
        """use_kickmap also updates Nslice."""
        idkm_elem.add_kickmap("n10", 10, idkm_file, 6.04)
        idkm_elem.add_kickmap("n5", 5, idkm_file, 6.04)

        idkm_elem.use_kickmap("n10")
        assert int(idkm_elem.Nslice) == 10

        idkm_elem.use_kickmap("n5")
        assert int(idkm_elem.Nslice) == 5

    def test_use_kickmap_roundtrip_restores_exact_tables(self, idkm_elem, idkm_file):
        """Switching A→B→A restores the exact kick arrays from mode A."""
        idkm_elem.add_kickmap("norm", 10, idkm_file, 6.04)
        idkm_elem.add_kickmap("half_e", 5, idkm_file, 3.0)

        idkm_elem.use_kickmap("norm")
        xkick_before = idkm_elem.xkick.copy()
        ykick_before = idkm_elem.ykick.copy()

        idkm_elem.use_kickmap("half_e")
        idkm_elem.use_kickmap("norm")

        assert_array_equal(idkm_elem.xkick, xkick_before)
        assert_array_equal(idkm_elem.ykick, ykick_before)

    def test_use_kickmap_unknown_key_raises_keyerror(self, idkm_elem, idkm_file):
        """use_kickmap raises KeyError when the key is not in the store."""
        idkm_elem.add_kickmap("mode_a", 10, idkm_file, 6.04)
        with pytest.raises(KeyError, match="bad_key"):
            idkm_elem.use_kickmap("bad_key")

    def test_use_kickmap_on_empty_store_raises_keyerror(self, idkm_elem):
        """use_kickmap raises KeyError when the store is empty."""
        with pytest.raises(KeyError):
            idkm_elem.use_kickmap("anything")

    def test_list_kickmaps_excludes_internal_active_key(self, idkm_elem, idkm_file):
        """The internal '_active' sentinel must never appear in list_kickmaps."""
        idkm_elem.add_kickmap("m", 10, idkm_file, 6.04)
        idkm_elem.use_kickmap("m")
        assert "_active" not in idkm_elem.list_kickmaps()

    def test_add_kickmap_does_not_activate(self, idkm_elem, idkm_file):
        """add_kickmap alone must not change the active tracking fields."""
        xkick_orig = idkm_elem.xkick.copy()
        idkm_elem.add_kickmap("other", 5, idkm_file, 3.0)
        assert_array_equal(idkm_elem.xkick, xkick_orig)
        assert idkm_elem.active_kickmap == "default"


# ---------------------------------------------------------------------------
# PassMethod ↔ kickmap swap interaction
# ---------------------------------------------------------------------------

class TestKickmapPassMethodInteraction:
    """use_kickmap must not interfere with set_DriftPass / set_IdTablePass."""

    def test_idtablepass_preserved_after_use_kickmap(self, idkm_elem, idkm_file):
        """use_kickmap must not change an existing IdTablePass setting."""
        idkm_elem.set_IdTablePass()
        idkm_elem.add_kickmap("m", 10, idkm_file, 6.04)
        idkm_elem.use_kickmap("m")
        assert idkm_elem.PassMethod == "IdTablePass"

    def test_driftpass_preserved_after_use_kickmap(self, idkm_elem, idkm_file):
        """use_kickmap must not revert a DriftPass setting."""
        idkm_elem.set_DriftPass()
        idkm_elem.add_kickmap("m", 10, idkm_file, 6.04)
        idkm_elem.use_kickmap("m")
        assert idkm_elem.PassMethod == "DriftPass"

    def test_set_driftpass_after_use_kickmap(self, idkm_elem, idkm_file):
        """set_DriftPass works correctly after a kickmap swap."""
        idkm_elem.add_kickmap("m", 10, idkm_file, 6.04)
        idkm_elem.use_kickmap("m")
        idkm_elem.set_DriftPass()
        assert idkm_elem.PassMethod == "DriftPass"

    def test_set_idtablepass_after_driftpass(self, idkm_elem, idkm_file):
        """set_IdTablePass restores tracking pass after set_DriftPass."""
        idkm_elem.add_kickmap("m", 10, idkm_file, 6.04)
        idkm_elem.use_kickmap("m")
        idkm_elem.set_DriftPass()
        idkm_elem.set_IdTablePass()
        assert idkm_elem.PassMethod == "IdTablePass"

    def test_kick_tables_intact_after_passmethod_round_trip(
        self, idkm_elem, idkm_file
    ):
        """DriftPass → IdTablePass round-trip must leave kick arrays unchanged."""
        idkm_elem.add_kickmap("m", 10, idkm_file, 6.04)
        idkm_elem.use_kickmap("m")
        xkick_before = idkm_elem.xkick.copy()
        ykick_before = idkm_elem.ykick.copy()

        idkm_elem.set_DriftPass()
        idkm_elem.set_IdTablePass()

        assert_array_equal(idkm_elem.xkick, xkick_before)
        assert_array_equal(idkm_elem.ykick, ykick_before)

    def test_use_kickmap_then_passmethod_cycle_then_swap_again(
        self, idkm_elem, idkm_file
    ):
        """Full cycle: swap→Drift→IdTable→swap works with correct final state."""
        idkm_elem.add_kickmap("norm", 10, idkm_file, 6.04)
        idkm_elem.add_kickmap("half_e", 5, idkm_file, 3.0)

        idkm_elem.use_kickmap("norm")
        idkm_elem.set_DriftPass()
        idkm_elem.set_IdTablePass()

        # now switch to a different kickmap
        idkm_elem.use_kickmap("half_e")
        assert idkm_elem.PassMethod == "IdTablePass"
        assert idkm_elem.active_kickmap == "half_e"
        assert int(idkm_elem.Nslice) == 5
