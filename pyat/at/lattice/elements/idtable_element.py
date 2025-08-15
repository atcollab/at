"""ID table :py:class:`.Element`."""

from __future__ import annotations

from pathlib import Path
from warnings import warn

import numpy as np

from ...constants import clight, e_mass
from .element_object import Element

# 2023jan15 orblancog first release
# 2023jan18 orblancog fix bug with element print
# 2023apr30 orblancog redefinition to function
# 2023jul04 orblancog changing class to save in .mat and .m
# 2025ago15 orblancog include first order kick maps from text files


def _anyarray(value: np.ndarray) -> None:
    # Ensure proper ordering(F) and alignment(A) for "C" access in integrators
    return np.require(value, dtype=np.float64, requirements=["F", "A"])


class InsertionDeviceKickMap(Element):
    """
    Insertion Device Element. Valid for a parallel electron beam.

      Pascale ELLEAUME, "A New Approach to the Electron Beam Dynamics in
        Undulators and  Wigglers". EPAC1992 0661.
        European Synchrotron Radiation Facility.
        BP 220, F-38043 Grenoble, France
    """

    _BUILD_ATTRIBUTES = Element._BUILD_ATTRIBUTES + [
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
        Nslice=int,
        xkick=_anyarray,
        ykick=_anyarray,
        xkick1=_anyarray,
        ykick1=_anyarray,
    )

    def __init__(self, family_name: str, *args, **kwargs) -> None:
        """
        Init IdTable.

        This __init__ takes the input to initialize an InsertionDeviceKickMap
        from an input file with arguments, for example at the moment of the element creation,
        or from a dictionary with all parameters, for example when reading a Lattice.

        Arguments:
            family_name: the family name
        """
        _argnames = [
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
        if len(args) < 11:
            # get data from text file
            elemargs = self.from_text_file(*args)
        else:
            # get data from arguments
            elemargs = dict(zip(_argnames, args))
        elemargs.update(kwargs)
        super().__init__(family_name, **elemargs)

    def set_DriftPass(self) -> None:
        """Set DriftPass tracking pass method."""
        self.PassMethod = "DriftPass"

    def set_IdTablePass(self) -> str:
        """Set IdTablePass tracking pass method."""
        self.PassMethod = "IdTablePass"

    def get_PassMethod(self) -> str:
        """Get the current tracking pass method.

        Returns:
            String with the current tracking pass method.
        """
        warn(UserWarning("get_PassMethod is deprecated; do not use"), stacklevel=2)
        return self.PassMethod

    def from_text_file(self, nslice: int, fname: str, norm_energy: float) -> tuple:
        """
        Create an Insertion Device Kick Map from a Radia field map file in text format.

        FamilyName is part of the base class, and all other arguments are described below.

        Arguments:
            nslice: number of slices in integrator.
            fname: input filename.
            norm_energy: particle energy in GeV.

        Returns:
            A tuple.
            Element parameters, default PassMethod: ``IdTablePass``.
        """

        def read_radia_field_map(file_in_name: str) -> tuple:
            """
            Read a RadiaField map in text format and return.

            A File, where :
            - comments start with #.
            - the first data line is the length in meters.
            - the second data line is the number of points in the h. plane.
            - the third data line is the number of points in the v. plane.
            - each block is a table with axes.
            - each data block comes after a START.
            - first the horizontal data block, and second the
              vertical data block with the second order kicks.
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
                            kick_block_list.append(kick_block)
                            kick_haxes_list.append(haxis)
                            kick_vaxes_list.append(vaxis)
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
                    "Input file contains more than 4 blocks. "
                    "Additional blocks ignored"
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

        def sorted_table(
            table_in: np.ndarray, sorted_index: np.ndarray, order_axis: str
        ) -> np.ndarray:
            # np.asfortranarray makes a copy of contiguous memory positions
            table_out = np.copy(table_in)
            for i, iis in zip(range(len(sorted_index)), sorted_index):
                if order_axis == "col":
                    table_out[:, i] = table_in[:, iis]
                if order_axis == "row":
                    table_out[i, :] = table_in[iis, :]
            return np.asfortranarray(table_out)

        # read the input data
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
        ) = read_radia_field_map(fname)

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


# EOF
