"""Tests for InsertionDeviceKickMap multi-kickmap store interface."""

from __future__ import annotations

from importlib.resources import files
from inspect import Parameter, signature

import pytest
from at.lattice import Lattice
from at.lattice.elements.idtable_element import InsertionDeviceKickMap
from numpy.testing import assert_array_equal

import machine_data


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

    def test_constructor_has_explicit_arguments(self):
        """Required kickmap inputs are visible in the public signature."""
        parameters = signature(InsertionDeviceKickMap).parameters
        assert {"family_name", "nslice", "fname", "norm_energy", "kickmaps"} <= set(
            parameters
        )
        assert all(
            parameter.kind is not Parameter.VAR_POSITIONAL
            for parameter in parameters.values()
        )

    def test_empty_element(self):
        """An element may be created before any kickmaps are available."""
        elem = InsertionDeviceKickMap("idmap")
        assert elem.PassMethod == "DriftPass"
        assert elem.list_kickmaps() == []
        assert elem.active_kickmap is None

    def test_first_added_kickmap_activates_empty_element(self, idkm_file):
        """Adding the first map makes an empty element trackable."""
        elem = InsertionDeviceKickMap("idmap")
        elem.add_kickmap("mode_a", 10, idkm_file, 6.04)
        assert elem.PassMethod == "IdTablePass"
        assert elem.list_kickmaps() == ["mode_a"]
        assert elem.active_kickmap == "mode_a"

    def test_constructor_accepts_named_kickmap_mapping(self, idkm_file):
        """Several named kickmaps may be supplied at construction."""
        elem = InsertionDeviceKickMap(
            "idmap",
            kickmaps={
                "mode_a": (10, idkm_file, 6.04),
                "mode_b": (5, idkm_file, 3.0),
            },
        )
        assert elem.list_kickmaps() == ["mode_a", "mode_b"]
        assert elem.active_kickmap == "mode_a"
        elem.use_kickmap("mode_b")
        assert int(elem.Nslice) == 5

    def test_constructor_rejects_single_map_and_mapping(self, idkm_file):
        """Single-map arguments cannot be mixed with a kickmap mapping."""
        with pytest.raises(TypeError, match="cannot be combined"):
            InsertionDeviceKickMap(
                "idmap",
                10,
                idkm_file,
                6.04,
                kickmaps={"mode_a": (10, idkm_file, 6.04)},
            )

    @pytest.mark.parametrize(
        "kwargs",
        [
            {"nslice": 10},
            {"fname": "kickmap.txt"},
            {"norm_energy": 6.04},
            {"nslice": 10, "fname": "kickmap.txt"},
        ],
    )
    def test_constructor_rejects_incomplete_single_map(self, kwargs):
        """Single-map construction requires all three map arguments."""
        with pytest.raises(TypeError, match="must be supplied together"):
            InsertionDeviceKickMap("idmap", **kwargs)

    @pytest.mark.parametrize("suffix", [".json", ".mat", ".m"])
    def test_empty_element_round_trip(self, tmp_path, suffix):
        """An empty element remains empty when saved and loaded."""
        fname = tmp_path / f"empty_idmap{suffix}"
        Lattice([InsertionDeviceKickMap("idmap")], energy=6.04e9).save(fname)

        loaded = Lattice.load(fname, energy=6.04e9, periodicity=1)[0]

        assert loaded.PassMethod == "DriftPass"
        assert loaded.list_kickmaps() == []
        assert loaded.active_kickmap is None

    def test_legacy_m_file_loads(self, tmp_path):
        """The previous positional MATLAB constructor syntax remains readable."""
        fname = tmp_path / "legacy_idmap.m"
        fname.write_text(
            """function ring = legacy_idmap()
ring = {...
atinsertiondevicekickmap('idmap','IdTablePass','',6.04,10,1,[1],[2],[0],[0],[0],[0]);...
};
end
""",
            encoding="utf-8",
        )

        loaded = Lattice.load(fname, energy=6.04e9, periodicity=1)[0]

        assert loaded.PassMethod == "IdTablePass"
        assert loaded.list_kickmaps() == ["default"]
        assert loaded.active_kickmap == "default"

    def test_list_kickmaps_on_new_element(self, idkm_elem):
        """A freshly constructed element has only the 'default' entry in the store."""
        assert idkm_elem.list_kickmaps() == ["default"]

    def test_active_kickmap_on_new_element(self, idkm_elem):
        """active_kickmap is 'default' immediately after construction."""
        assert idkm_elem.active_kickmap == "default"

    def test_default_kickmap_matches_initial_fields(self, idkm_elem):
        """The 'default' kickmap must contain the element's construction-time arrays."""
        assert_array_equal(
            idkm_elem.KickmapStore["default"]["xkick"], idkm_elem.xkick
        )
        assert_array_equal(
            idkm_elem.KickmapStore["default"]["ykick"], idkm_elem.ykick
        )

    def test_use_kickmap_default_restores_initial_fields(self, idkm_elem, idkm_file):
        """use_kickmap('default') restores the original kick arrays after a swap."""
        xkick_initial = idkm_elem.xkick.copy()
        idkm_elem.add_kickmap("half_e", 5, idkm_file, 3.0)
        idkm_elem.use_kickmap("half_e")
        idkm_elem.use_kickmap("default")
        assert_array_equal(idkm_elem.xkick, xkick_initial)
        assert idkm_elem.active_kickmap == "default"

    def test_add_kickmap_lists_new_key(self, idkm_elem, idkm_file):
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

        # Kicks are normalized by 1/E², so different energies change the tables.
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
        idkm_elem.KickmapStore = {}
        with pytest.raises(KeyError):
            idkm_elem.use_kickmap("anything")

    def test_active_is_available_as_a_kickmap_key(self, idkm_elem, idkm_file):
        """The active marker does not reserve a key in the kickmap store."""
        idkm_elem.add_kickmap("active", 10, idkm_file, 6.04)
        idkm_elem.use_kickmap("active")
        assert idkm_elem.active_kickmap == "active"
        assert "active" in idkm_elem.list_kickmaps()

    def test_add_kickmap_does_not_activate(self, idkm_elem, idkm_file):
        """add_kickmap alone must not change the active tracking fields."""
        xkick_orig = idkm_elem.xkick.copy()
        idkm_elem.add_kickmap("other", 5, idkm_file, 3.0)
        assert_array_equal(idkm_elem.xkick, xkick_orig)
        assert idkm_elem.active_kickmap == "default"

    @pytest.mark.parametrize("suffix", [".json", ".mat", ".m"])
    def test_kickmap_store_round_trip(
        self, idkm_elem, idkm_file, tmp_path, suffix
    ):
        """All native lattice formats preserve every stored kickmap."""
        idkm_elem.add_kickmap("mode_a", 10, idkm_file, 6.04)
        smaller_map = {
            "Length": idkm_elem.Length,
            "xkick": idkm_elem.xkick[:7, :5],
            "ykick": idkm_elem.ykick[:7, :5],
            "xkick1": idkm_elem.xkick1[:7, :5],
            "ykick1": idkm_elem.ykick1[:7, :5],
            "xtable": idkm_elem.xtable[:5],
            "ytable": idkm_elem.ytable[:7],
        }
        idkm_elem.add_kickmap("mode_b", 5, smaller_map, 3.0)
        idkm_elem.use_kickmap("mode_b")
        fname = tmp_path / f"idmap{suffix}"

        Lattice([idkm_elem], energy=6.04e9).save(fname)
        loaded = Lattice.load(fname)[0]

        assert loaded.list_kickmaps() == ["default", "mode_a", "mode_b"]
        assert loaded.active_kickmap == "mode_b"
        assert loaded.xkick.shape == (7, 5)
        loaded.use_kickmap("mode_a")
        assert_array_equal(
            loaded.xkick, idkm_elem.KickmapStore["mode_a"]["xkick"]
        )


# ---------------------------------------------------------------------------
# PassMethod ↔ kickmap swap interaction
# ---------------------------------------------------------------------------

class TestKickmapPassMethodInteraction:
    """use_kickmap must not interfere with enable / disable."""

    def test_enabled_state_preserved_after_use_kickmap(
        self, idkm_elem, idkm_file
    ):
        """use_kickmap must not change the enabled pass method."""
        idkm_elem.add_kickmap("m", 10, idkm_file, 6.04)
        idkm_elem.use_kickmap("m")
        assert idkm_elem.PassMethod == "IdTablePass"

    def test_disabled_state_preserved_after_use_kickmap(
        self, idkm_elem, idkm_file
    ):
        """use_kickmap must not re-enable a disabled element."""
        idkm_elem.disable()
        idkm_elem.add_kickmap("m", 10, idkm_file, 6.04)
        idkm_elem.use_kickmap("m")
        assert idkm_elem.PassMethod == "DriftPass"

    def test_disable_after_use_kickmap(self, idkm_elem, idkm_file):
        """disable works correctly after a kickmap swap."""
        idkm_elem.add_kickmap("m", 10, idkm_file, 6.04)
        idkm_elem.use_kickmap("m")
        idkm_elem.disable()
        assert idkm_elem.PassMethod == "DriftPass"

    def test_enable_after_disable(self, idkm_elem, idkm_file):
        """enable restores the tracking pass method after disable."""
        idkm_elem.add_kickmap("m", 10, idkm_file, 6.04)
        idkm_elem.use_kickmap("m")
        idkm_elem.disable()
        idkm_elem.enable()
        assert idkm_elem.PassMethod == "IdTablePass"

    def test_kick_tables_intact_after_passmethod_round_trip(
        self, idkm_elem, idkm_file
    ):
        """disable / enable round-trip must leave kick arrays unchanged."""
        idkm_elem.add_kickmap("m", 10, idkm_file, 6.04)
        idkm_elem.use_kickmap("m")
        xkick_before = idkm_elem.xkick.copy()
        ykick_before = idkm_elem.ykick.copy()

        idkm_elem.disable()
        idkm_elem.enable()

        assert_array_equal(idkm_elem.xkick, xkick_before)
        assert_array_equal(idkm_elem.ykick, ykick_before)

    def test_use_kickmap_then_passmethod_cycle_then_swap_again(
        self, idkm_elem, idkm_file
    ):
        """Full cycle: swap→Drift→IdTable→swap works with correct final state."""
        idkm_elem.add_kickmap("norm", 10, idkm_file, 6.04)
        idkm_elem.add_kickmap("half_e", 5, idkm_file, 3.0)

        idkm_elem.use_kickmap("norm")
        idkm_elem.disable()
        idkm_elem.enable()

        # now switch to a different kickmap
        idkm_elem.use_kickmap("half_e")
        assert idkm_elem.PassMethod == "IdTablePass"
        assert idkm_elem.active_kickmap == "half_e"
        assert int(idkm_elem.Nslice) == 5
