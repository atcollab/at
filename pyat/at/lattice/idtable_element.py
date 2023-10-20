from .elements import Element
import numpy
import io
from ..constants import clight
from warnings import warn


def _anyarray(value):
    # Ensure proper ordering(F) and alignment(A) for "C" access in integrators
    return numpy.require(value, dtype=numpy.float64, requirements=['F', 'A'])


class InsertionDeviceKickMap(Element):
    """
    Insertion Device Element. Valid for a parallel electron beam.

      Pascale ELLEAUME, "A New Approach to the Electron Beam Dynamics in
        Undulators and  Wigglers". EPAC1992 0661.
        European Synchrotron Radiation Facility.
        BP 220, F-38043 Grenoble, France
    """
    _BUILD_ATTRIBUTES = Element._BUILD_ATTRIBUTES + ['PassMethod',
                                                     'Filename_in',
                                                     'Normalization_energy',
                                                     'Nslice',
                                                     'Length',
                                                     'xkick',
                                                     'ykick',
                                                     'xkick1',
                                                     'ykick1',
                                                     'xtable',
                                                     'ytable']

    _conversions = dict(Element._conversions,
                        Nslice=int,
                        xkick=_anyarray,
                        ykick=_anyarray,
                        xkick1=_anyarray,
                        ykick1=_anyarray
                        )

    def set_DriftPass(self):
        setattr(self, 'PassMethod', 'DriftPass')

    def set_IdTablePass(self):
        setattr(self, 'PassMethod', 'IdTablePass')

    def get_PassMethod(self):
        warn(UserWarning('get_PassMethod is deprecated; do not use'))
        return getattr(self, 'PassMethod')

    def from_text_file(self, Nslice, Filename_in, Energy):
        """
        This function creates an Insertion Device Kick Map
        from a Radia field map file.

        Args:
            family_name:    family name
            Nslice:         number of slices in integrator
            Filename_in:    input filename
            Energy:         particle energy in GeV

        Returns:
            KickMap element
            Default PassMethod: ``IdTablePass``
        """
        # 2023jul04 changing class to save in .mat and .m
        # 2023apr30 redefinition to function
        # 2023jan18 fix bug with element print
        # 2023jan15 first release
        # orblancog
        def readRadiaFieldMap(file_in_name):
            """
            Read a RadiaField map and return
            """
            with io.open(file_in_name, mode="r", encoding="utf-8") as f:
                """
                File, where :
                - the first data line is the length in meters
                - the second data line is the number of points in the h. plane
                - the third data line is the number of points in the v. plane
                - each data block comes after a START
                - first the horizontal data block, and second the
                      vertical data block
                - each block is a table with axes
                - comments start with #
                """
                # File example:
                # #comment in line 1
                # #comment in line 2
                # Length_in_m
                # #comment in line 4
                # Number of points in horizontal plane :nh
                # #comment in line 6
                # Number of points in vertical plane :nv
                # #comment in line 8
                # START
                #             pos_point1h pos_point2h ... pos_pointnh
                # pos_point1v
                # ...                    horizontal kick_map(nv,nh)
                # pos_pointnv
                # START
                #             pos_point1h pos_point2h ... pos_pointnh
                # pos_point1v
                # ...                    vertical kick_map(nv,nh)
                # pos_pointnv
                # (EOL)

                data_lines = 0     # line not starting with '#'
                header_lines = 0   # line starting with '#'
                block_counter = 0  # START of the h.map, START of the v.map
                for line in f:
                    sline = line.split()
                    if sline[0] == '#':  # line is comment
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
                            kick_map = numpy.zeros((v_points, h_points))
                            haxis = numpy.zeros(h_points)
                            vaxis = numpy.zeros(v_points)
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
                                kick_map[block_lines - 2][:] = sline[1:]
                            if block_lines > v_points:
                                block_lines = 0
                                if block_counter == 1:
                                    hkickmap = numpy.copy(kick_map)
                                    table_cols1 = haxis
                                    table_rows1 = vaxis
                                if block_counter == 2:
                                    vkickmap = numpy.copy(kick_map)
                                    table_cols2 = haxis
                                    table_rows2 = vaxis
                                if block_counter > 2:
                                    print('atWarning: only two tables read')
                            block_lines += 1
            # dummy variables not implemented in the reading function
            # but required
            hkickmap1 = 0.0 * numpy.copy(hkickmap)
            vkickmap1 = 0.0 * numpy.copy(vkickmap)

            return el_length, hkickmap, vkickmap, table_cols1, table_rows1, \
                table_cols2, table_rows2, h_points, v_points, \
                hkickmap1, vkickmap1

        def sorted_table(table_in, sorted_index, order_axis):
            # numpy.asfortranarray makes a copy of contiguous memory positions
            table_out = numpy.copy(table_in)
            for i, iis in zip(range(len(sorted_index)), sorted_index):
                if order_axis == 'col':
                    table_out[:, i] = table_in[:, iis]
                if order_axis == 'row':
                    table_out[i, :] = table_in[iis, :]
            table_out2 = numpy.asfortranarray(table_out)
            return table_out2

        # read the input data
        el_length, hkickmap, vkickmap, \
            table_cols1, table_rows1, \
            table_cols2, table_rows2, \
            NumX, NumY, \
            hkickmap1, vkickmap1 \
            = readRadiaFieldMap(Filename_in)

        # set to float
        table_cols1array = numpy.array(table_cols1, dtype='float64')
        table_rows1array = numpy.array(table_rows1, dtype='float64')
        table_cols2array = numpy.array(table_cols2, dtype='float64')
        table_rows2array = numpy.array(table_rows2, dtype='float64')

        # Reorder table_axes
        cols1sorted_index = numpy.argsort(table_cols1array)
        table_cols1array.sort()
        rows1sorted_index = numpy.argsort(table_rows1array)
        table_rows1array.sort()
        cols2sorted_index = numpy.argsort(table_cols2array)
        table_cols2array.sort()
        rows2sorted_index = numpy.argsort(table_rows2array)
        table_rows2array.sort()
        # Reorder kickmap
        hkickmap_a = sorted_table(hkickmap, cols1sorted_index, 'col')
        hkickmap = sorted_table(hkickmap_a, rows1sorted_index, 'row')
        vkickmap_a = sorted_table(vkickmap, cols2sorted_index, 'col')
        vkickmap = sorted_table(vkickmap_a, rows2sorted_index, 'row')
        # Reorder kickmap1
        hkickmap1_a = sorted_table(hkickmap1, cols1sorted_index, 'col')
        hkickmap1 = sorted_table(hkickmap1_a, rows1sorted_index, 'row')
        vkickmap1_a = sorted_table(vkickmap1, cols2sorted_index, 'col')
        vkickmap1 = sorted_table(vkickmap1_a, rows2sorted_index, 'row')

        # Field to kick factors
        Brho = 1e9 * Energy/clight
        factor = 1.0/(Brho**2)
        xkick = factor * hkickmap
        ykick = factor * vkickmap
        # kick1 vars set to zero, not yet implemented
        factor1 = -1.0/(Brho)
        xkick1 = factor1 * hkickmap1
        ykick1 = factor1 * vkickmap1
        xtable = table_cols1array.T
        ytable = table_rows1array.T

        args_dict = {'PassMethod': 'IdTablePass',
                     'Filename_in': Filename_in,
                     'Normalization_energy': Energy,
                     'Nslice': numpy.uint8(Nslice),
                     'Length': el_length,
                     'xkick': xkick,
                     'ykick': ykick,
                     'xkick1': xkick1,
                     'ykick1': ykick1,
                     'xtable': xtable,
                     'ytable': ytable,
                     }
        return args_dict

    def __init__(self, family_name: str, *args, **kwargs):
        _argnames = ['PassMethod',
                     'Filename_in',
                     'Normalization_energy',
                     'Nslice',
                     'Length',
                     'xkick',
                     'ykick',
                     'xkick1',
                     'ykick1',
                     'xtable',
                     'ytable']
        if len(args) < 11:
            # get data from text file
            elemargs = self.from_text_file(*args)
        else:
            # get data from arguments
            elemargs = dict(zip(_argnames, args))
        elemargs.update(kwargs)
        super().__init__(family_name, **elemargs)
# EOF
