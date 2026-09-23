"""ID table :py:class:`.Element`."""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any
from warnings import warn

import numpy as np

from ...constants import clight, e_mass
from .element_object import Element


def _anyarray(value: np.ndarray) -> np.ndarray:
    # Ensure proper ordering(F) and alignment(A) for "C" access in integrators
    return np.require(value, dtype=np.float64, requirements=["F", "A"])


class InsertionDeviceKickMap(Element):
    """Insertion device kick-map element for a parallel electron beam.

    The element implements tracking through integrated first- and second-order
    magnetic field maps normalized to a reference energy. It may be created
    empty, from one kickmap, or from a mapping of named kickmaps.

    Args:
        family_name: Element family name.
        length: Insertion device length in m. A zero value uses the length
            read from a supplied kickmap.
        norm_energy: Normalization energy in GeV.

    Keyword Args:
        nslice: Number of integration slices for a single kickmap.
        fname: Radia text-file path or dictionary containing a single kickmap.
        kickmaps: Mapping of kickmap names to ``(nslice, source, energy)``
            tuples. The first entry is initially active.
        **kwargs: Additional element attributes. This is also used internally
            when restoring an element from a lattice file.

    If no kickmap arguments are supplied, an empty element with
    ``PassMethod="DriftPass"`` is created. Adding the first kickmap activates
    it and changes the pass method to ``IdTablePass``. Supplying ``nslice``,
    ``fname`` and ``norm_energy`` creates one kickmap named ``"default"``.

    First-order maps are optional. Positive and negative signs are not applied
    by this implementation, so the input data must already contain the desired
    signs.

    See P. Elleaume, "A New Approach to the Electron Beam Dynamics in
    Undulators and Wigglers", EPAC 1992, 0661.

    Examples:
        Create an empty element and add a kickmap later:

        >>> elem = InsertionDeviceKickMap("ID")
        >>> elem.add_kickmap("LH", 25, "lh_kickmap.txt", 2.75)

        Create an element with one kickmap:

        >>> elem = InsertionDeviceKickMap(
        ...     "ID",
        ...     norm_energy=2.75,
        ...     nslice=25,
        ...     fname="lh_kickmap.txt",
        ... )

        Create an element with several named kickmaps:

        >>> elem = InsertionDeviceKickMap(
        ...     "ID",
        ...     kickmaps={
        ...         "LH": (25, "lh_kickmap.txt", 2.75),
        ...         "LV": (25, "lv_kickmap.txt", 2.75),
        ...     },
        ... )
    """

    _BUILD_ATTRIBUTES = [
        *Element._BUILD_ATTRIBUTES,
        "Length",
        "Normalization_energy",
    ]
    _LEGACY_M_BUILD_ATTRIBUTES = [
        "FamName",
        "PassMethod",
        "Filename_in",
        "Normalization_energy",
        "Nslice",
        "Length",
        "xkick",
        "ykick",
        "xkick1",
        "ykick1",
        "xtable",
        "ytable",
    ]

    _conversions = dict(
        Element._conversions,
        Normalization_energy=float,
        Nslice=int,
        xkick=_anyarray,
        ykick=_anyarray,
        xkick1=_anyarray,
        ykick1=_anyarray,
        xtable=_anyarray,
        ytable=_anyarray,
    )

    def __init__(
        self: InsertionDeviceKickMap,
        family_name: str,
        length: float = 0.0,
        norm_energy: float = 0.0,
        *,
        nslice: int | None = None,
        fname: str | Path | dict[str, Any] | None = None,
        kickmaps: Mapping[
            str, tuple[int, str | Path | dict[str, Any], float]
        ]
        | None = None,
        **kwargs: Any,
    ) -> None:
        """Initialize an insertion device kick-map element."""
        file_data = "xkick" in kwargs or "KickmapStore" in kwargs
        source_values = (nslice, fname)

        if file_data:
            if kickmaps is not None or any(
                value is not None for value in source_values
            ):
                msg = "Kickmap arguments cannot be combined with serialized data"
                raise TypeError(msg)
            kwargs.setdefault("Length", length)
            kwargs.setdefault("Normalization_energy", norm_energy)
            super().__init__(family_name, **kwargs)
            if hasattr(self, "KickmapStore"):
                self._normalise_kickmap_store()
            else:
                self.KickmapStore = {"default": self._snapshot()}
                self.ActiveKickmap = "default"
            return

        if kickmaps is not None and (
            length != 0.0
            or norm_energy != 0.0
            or any(value is not None for value in source_values)
        ):
            msg = (
                "kickmaps cannot be combined with length, norm_energy, "
                "nslice or fname"
            )
            raise TypeError(msg)
        if any(value is not None for value in source_values) and (
            not all(value is not None for value in source_values)
            or norm_energy == 0.0
        ):
            msg = "nslice, fname and norm_energy must be supplied together"
            raise TypeError(msg)

        if kickmaps:
            entries = list(kickmaps.items())
            first_key, first_spec = entries[0]
            first_data = self._load_kickmap(
                *self._validate_kickmap_spec(first_key, first_spec)
            )
            first_data.update(kwargs)
            super().__init__(family_name, **first_data)
            self.KickmapStore = {first_key: self._snapshot()}
            self.ActiveKickmap = first_key
            for key, spec in entries[1:]:
                self._store_kickmap(key, *self._validate_kickmap_spec(key, spec))
        elif all(value is not None for value in source_values):
            elemargs = self._load_kickmap(nslice, fname, norm_energy)
            if length != 0.0:
                elemargs["Length"] = length
            elemargs.update(kwargs)
            super().__init__(family_name, **elemargs)
            self.KickmapStore = {"default": self._snapshot()}
            self.ActiveKickmap = "default"
        else:
            kwargs.setdefault("PassMethod", "DriftPass")
            kwargs.setdefault("Length", length)
            kwargs.setdefault("Normalization_energy", norm_energy)
            super().__init__(family_name, **kwargs)
            self.KickmapStore = {}
            self.ActiveKickmap = ""
            self._disable_passmethod = "DriftPass"

    def to_dict(self) -> dict[str, Any]:
        """Return serializable attributes, omitting the default empty store."""
        attributes = super().to_dict()
        if not self.KickmapStore:
            attributes.pop("KickmapStore")
            attributes.pop("ActiveKickmap")
        return attributes

    @staticmethod
    def _validate_kickmap_spec(
        key: str,
        spec: tuple[int, str | Path | dict[str, Any], float],
    ) -> tuple[int, str | Path | dict[str, Any], float]:
        if not isinstance(key, str):
            msg = "Kickmap names must be strings"
            raise TypeError(msg)
        if not isinstance(spec, (tuple, list)) or len(spec) != 3:
            msg = (
                f"Kickmap {key!r} must be defined as "
                "(nslice, source, norm_energy)"
            )
            raise TypeError(msg)
        map_nslice, map_source, map_energy = spec
        return map_nslice, map_source, map_energy

    def _store_kickmap(
        self,
        key: str,
        nslice: int,
        fname: str | Path | dict[str, Any],
        norm_energy: float,
    ) -> dict:
        if not isinstance(key, str):
            msg = "Kickmap names must be strings"
            raise TypeError(msg)
        data = self._load_kickmap(nslice, fname, norm_energy)
        self.KickmapStore[key] = {
            field: data[field]
            for field in (
                "Filename_in",
                "Normalization_energy",
                "Nslice",
                "Length",
                "xkick",
                "ykick",
                "xkick1",
                "ykick1",
                "xtable",
                "ytable",
            )
        }
        return data

    def _snapshot(self: InsertionDeviceKickMap) -> dict:
        """Return a dict of the element's current tracking-field values."""
        return {
            "Filename_in": getattr(self, "Filename_in", ""),
            "Normalization_energy": float(self.Normalization_energy),
            "Nslice": int(self.Nslice),
            "Length": float(self.Length),
            "xkick":  self.xkick.copy(),
            "ykick":  self.ykick.copy(),
            "xkick1": self.xkick1.copy(),
            "ykick1": self.ykick1.copy(),
            "xtable": self.xtable.copy(),
            "ytable": self.ytable.copy(),
        }

    def _normalise_kickmap_store(self: InsertionDeviceKickMap) -> None:
        """Restore kickmap types after loading from a lattice file."""
        for data in self.KickmapStore.values():
            data["Filename_in"] = str(data["Filename_in"])
            data["Normalization_energy"] = float(data["Normalization_energy"])
            data["Nslice"] = int(data["Nslice"])
            data["Length"] = float(data["Length"])
            for field in ("xkick", "ykick", "xkick1", "ykick1", "xtable", "ytable"):
                data[field] = _anyarray(data[field])

    def _load_kickmap(
        self: InsertionDeviceKickMap,
        nslice: int,
        fname: str | Path | dict[str, Any],
        norm_energy: float,
    ) -> dict:
        """Load and normalize one Radia field map."""

        def sorted_table(
            table_in: np.ndarray, sorted_index: np.ndarray, order_axis: str
        ) -> np.ndarray:
            """Return ordered table.

            Arguments:
                table_in: input table..
                sorted_index: index to sort the table.
                order_axis: sort on 'col' or 'row'.

            Returns:
                Sorted table as fortran array.
            """
            # np.asfortranarray makes a copy of contiguous memory positions
            table_out = np.copy(table_in)
            for i, iis in zip(range(len(sorted_index)), sorted_index):
                if order_axis == "col":
                    table_out[:, i] = table_in[:, iis]
                if order_axis == "row":
                    table_out[i, :] = table_in[iis, :]
            return np.asfortranarray(table_out)

        if isinstance(fname, dict):
            thefields = self.read_dict_radia_field_map(fname)
            fname = ""
        else:
            # assume text file
            fname = str(fname)
            thefields = self.read_text_radia_field_map(fname)

        (
            el_length,
            hkickmap2,
            vkickmap2,
            table_colshkick,
            table_rowshkick,
            table_colsvkick,
            table_rowsvkick,
            hkickmap1,
            vkickmap1,
            _,
            _,
        ) = thefields

        # set to float
        table_colshkickarray = np.array(table_colshkick, dtype="float64")
        table_rowshkickarray = np.array(table_rowshkick, dtype="float64")
        table_colsvkickarray = np.array(table_colsvkick, dtype="float64")
        table_rowsvkickarray = np.array(table_rowsvkick, dtype="float64")

        # Reorder table_axes
        cols1sorted_index = np.argsort(table_colshkickarray)
        table_colshkickarray.sort()
        rows1sorted_index = np.argsort(table_rowshkickarray)
        table_rowshkickarray.sort()
        cols2sorted_index = np.argsort(table_colsvkickarray)
        table_colsvkickarray.sort()
        rows2sorted_index = np.argsort(table_rowsvkickarray)
        table_rowsvkickarray.sort()
        # Reorder kickmap2
        hkickmap2_a = sorted_table(hkickmap2, cols1sorted_index, "col")
        hkickmap2 = sorted_table(hkickmap2_a, rows1sorted_index, "row")
        vkickmap2_a = sorted_table(vkickmap2, cols2sorted_index, "col")
        vkickmap2 = sorted_table(vkickmap2_a, rows2sorted_index, "row")
        # Reorder kickmap1
        hkickmap1_a = sorted_table(hkickmap1, cols1sorted_index, "col")
        hkickmap1 = sorted_table(hkickmap1_a, rows1sorted_index, "row")
        vkickmap1_a = sorted_table(vkickmap1, cols2sorted_index, "col")
        vkickmap1 = sorted_table(vkickmap1_a, rows2sorted_index, "row")

        # Field to kick factors
        e_mass_gev = e_mass * 1e-9
        brho = 1e9 * np.sqrt(norm_energy**2 - e_mass_gev**2) / clight
        # kick2 vars
        factor2 = 1.0 / (brho**2)
        xkick = factor2 * hkickmap2
        ykick = factor2 * vkickmap2
        # kick1 vars
        factor1 = 1.0 / (brho)
        xkick1 = factor1 * hkickmap1
        ykick1 = factor1 * vkickmap1
        # axes
        xtable = table_colshkickarray.T
        ytable = table_rowshkickarray.T

        return {
            "PassMethod": "IdTablePass",
            "Filename_in": fname,
            "Normalization_energy": norm_energy,
            "Nslice": np.uint8(nslice),
            "Length": el_length,
            "xkick": xkick,
            "ykick": ykick,
            "xkick1": xkick1,
            "ykick1": ykick1,
            "xtable": xtable,
            "ytable": ytable,
        }

    def read_text_radia_field_map(
        self: InsertionDeviceKickMap, file_in_name: str
    ) -> tuple:
        """
        Read a RadiaField map in text format and return.

        A File, where :
        - comments start with #.
        - the first data line is the length in meters.
        - the second data line is the number of points in the h. plane.
        - the third data line is the number of points in the v. plane.
        - each block is a table with axes.
        - each data block comes after a START.
        - first the horizontal data block, and second the vertical data
        block with the second order kicks.
        There might be two other blocks with the horizontal and
        vertical first order kicks.

        File example (ignore the !SPACE):
        ! #comment in line 1
        ! #comment in line 2
        ! Length_in_m
        ! #comment in line 4
        ! Number of points in horizontal plane :nh
        ! #comment in line 6
        ! Number of points in vertical plane :nv
        ! #comment in line 8
        ! START
        !             pos_point1h pos_point2h ... pos_pointnh
        ! pos_point1v
        ! ...                    horizontal kick_map(nv,nh)
        ! pos_pointnv
        ! START
        !             pos_point1h pos_point2h ... pos_pointnh
        ! pos_point1v
        ! ...                    vertical kick_map(nv,nh)
        ! pos_pointnv
        ! (EOL)

        Arguments:
            file_in_name: the file name.

        Returns:
            Tuple with file tables and axes.

        Raises:
            ValueError: if the number of blocks in less than 2 or equal to 3.
        """
        thepath = Path(file_in_name)
        with thepath.open(encoding="utf-8") as thefile:
            lines = thefile.readlines()
        thefile.close()
        data_lines = 0  # line not starting with '#'
        header_lines = 0  # line starting with '#'
        block_counter = 0  # START of the h.map, START of the v.map
        kick_block_list = []
        kick_haxes_list = []
        kick_vaxes_list = []
        for line in lines:
            sline = line.split()
            if sline[0] == "#":  # line is comment
                header_lines += 1
            else:
                data_lines += 1
                if data_lines == 1:  # get the element length
                    el_length = float(sline[0])
                elif data_lines == 2:  # get the number of hor. points
                    h_points = int(sline[0])
                elif data_lines == 3:  # get the number of ver. points
                    v_points = int(sline[0])
                    # initialize element kicks and table_axes
                    kick_block = np.zeros((v_points, h_points))
                    haxis = np.zeros(h_points)
                    vaxis = np.zeros(v_points)
                else:
                    # read block of data
                    if sline[0] == "START" or sline[0] == "START\n":
                        block_counter += 1
                        block_lines = 0
                    if block_lines == 1:
                        haxis = sline
                    if block_lines > 1:
                        # minus one due to python index starting at 0
                        # and minus another one due
                        # to the column labels in first line
                        vaxis[block_lines - 2] = float(sline[0])
                        kick_block[block_lines - 2][:] = sline[1:]
                    if block_lines > v_points:
                        block_lines = 0
                        kick_block_list.append(np.copy(kick_block))
                        kick_haxes_list.append(np.copy(haxis))
                        kick_vaxes_list.append(np.copy(vaxis))
                    block_lines += 1
        # checking how many kick blocks were added
        lenkick_block_list = len(kick_block_list)
        if lenkick_block_list < 2 or lenkick_block_list == 3:
            _minimumblocknumbererrormsg = (
                "Input file contains only " f"{len(kick_block_list)} block"
            )
            raise ValueError(_minimumblocknumbererrormsg)
        if lenkick_block_list == 2:
            # first order kick not in file
            kick_block_list.append(0.0 * np.copy(kick_block))
            kick_block_list.append(0.0 * np.copy(kick_block))
        elif lenkick_block_list > 4:
            # file contains more blocks that required
            _warn4kickblocks = (
                "Input file contains more than 4 blocks. Additional blocks ignored"
            )
            warn(_warn4kickblocks)

        return (
            el_length,
            kick_block_list[0],
            kick_block_list[1],
            kick_haxes_list[0],
            kick_vaxes_list[0],
            kick_haxes_list[1],
            kick_vaxes_list[1],
            kick_block_list[2],
            kick_block_list[3],
            h_points,
            v_points,
        )

    def read_dict_radia_field_map(
        self: InsertionDeviceKickMap, id_input: dict
    ) -> tuple:
        """Read a dictionary with Radia field map tables.

        The required keys are "Length", "xkick" and "ykick"
        for the second order maps, "xtable" and "ytable" for
        the grid, and "xkick1" and "ykick1" for the first order
        maps.

        Arguments:
            id_input: Radia field map input.

        Returns:
            Tuple with Insertion Device parameters.
        """
        (v_points, h_points) = id_input["xkick"].shape
        return (
            id_input["Length"],
            id_input["xkick"],
            id_input["ykick"],
            id_input["xtable"],
            id_input["ytable"],
            id_input["xtable"],
            id_input["ytable"],
            id_input["xkick1"],
            id_input["ykick1"],
            h_points,
            v_points,
        )
    def _apply_kickmap_data(self: InsertionDeviceKickMap, data: dict) -> None:
        """Apply kickmap data dict to the active tracking fields."""
        self.Filename_in = data["Filename_in"]
        self.Normalization_energy = data["Normalization_energy"]
        self.Nslice = data["Nslice"]
        self.Length = data["Length"]
        self.xkick = data["xkick"]
        self.ykick = data["ykick"]
        self.xkick1 = data["xkick1"]
        self.ykick1 = data["ykick1"]
        self.xtable = data["xtable"]
        self.ytable = data["ytable"]

    def add_kickmap(
        self: InsertionDeviceKickMap,
        key: str,
        nslice: int,
        fname: str | Path | dict[str, Any],
        norm_energy: float,
    ) -> None:
        """Store a kickmap under a string key without activating it.

        Use :meth:`use_kickmap` to make a stored kickmap active for tracking.

        Arguments:
            key: identifier string for this kickmap.
            nslice: number of slices in integrator.
            fname: input filename (text file or dict).
            norm_energy: normalization energy in GeV.
        """
        was_empty = not self.KickmapStore
        data = self._store_kickmap(key, nslice, fname, norm_energy)
        if was_empty:
            self.Filename_in = data["Filename_in"]
            self.Normalization_energy = data["Normalization_energy"]
            self.ActiveKickmap = key
            self._apply_kickmap_data(self.KickmapStore[key])
            self.PassMethod = data["PassMethod"]
            self._enable_passmethod = data["PassMethod"]
            self._disable_passmethod = "DriftPass"

    def use_kickmap(self: InsertionDeviceKickMap, key: str) -> None:
        """Activate a stored kickmap by key for tracking.

        Swaps xtable, ytable, xkick, ykick, xkick1, ykick1 and Nslice/Length
        to the data stored under *key*.

        Arguments:
            key: identifier of a kickmap previously added with
                :meth:`add_kickmap`.

        Raises:
            KeyError: if *key* is not found in the store.
        """
        if not hasattr(self, "KickmapStore") or key not in self.KickmapStore:
            available = self.list_kickmaps()
            raise KeyError(
                f"Kickmap '{key}' not found. Available keys: {available}"
            )
        self.ActiveKickmap = key
        self._apply_kickmap_data(self.KickmapStore[key])

    @property
    def active_kickmap(self: InsertionDeviceKickMap) -> str | None:
        """The key of the currently active kickmap, or None if not set."""
        return getattr(self, "ActiveKickmap", "") or None

    def list_kickmaps(self: InsertionDeviceKickMap) -> list[str]:
        """Return the list of stored kickmap keys.

        Returns:
            List of key strings.
        """
        return list(getattr(self, "KickmapStore", {}))

# EOF
