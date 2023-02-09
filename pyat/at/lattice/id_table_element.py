import numpy
import io
from scipy.constants import c as clight

class InsertionDevice(Element):
    """
    Insertion Device Element. Valid for a parallel electron beam.

      Pascale ELLEAUME, "A New Approach to the Electron Beam Dynamics in
        Undulators and  Wigglers". EPAC1992 0661.
        European Synchrotron Radiation Facility.
        BP 220, F-38043 Grenoble, France
    """
    _BUILD_ATTRIBUTES = Element._BUILD_ATTRIBUTES + ['Nslice',
                                                     'Filename_in',
                                                     'Energy']

    def readRadiaFieldMap(self, file_in_name):
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
            - first the horizontal data block, and second the vertical data block
            - each block is a table with axes
            - comments start with #
            """
            ### File example:
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

            data_lines = 0;     # line not starting with '#'
            header_lines = 0;   # line starting with '#'
            block_counter  = 0; # START of the hor. map and START of the vert. map
            for line in f:
                sline = line.split()
                if sline[0] == '#': # line is comment
                    header_lines += 1;
                else:
                    data_lines += 1;
                    if   data_lines == 1: el_length = float(sline[0]); # get the element length
                    elif data_lines == 2: h_points  = int(sline[0]);   # get the number of horizontal points
                    elif data_lines == 3:
                        # get the number of vertical points
                        v_points  = int(sline[0]);
                        # initialize element kicks and table_axes
                        kick_map = numpy.zeros((v_points,h_points));
                        haxis = numpy.zeros(h_points);
                        vaxis = numpy.zeros(v_points);
                    else:
                        # read block of data
                        if sline[0] == "START" or sline[0] == "START\n" :
                            block_counter += 1;
                            block_lines = 0;
                        if block_lines == 1: haxis = sline;
                        if block_lines > 1:
                            # minus one due to python index starting at zero
                            # and minus another one due to the column labels in first line
                            vaxis[block_lines - 2]       = float(sline[0]);
                            kick_map[block_lines - 2][:] = sline[1:];
                        if block_lines > v_points:
                            block_lines = 0;
                            if block_counter == 1:
                                hkickmap    = numpy.copy(kick_map);
                                table_cols1 = haxis;
                                table_rows1 = vaxis;
                            if block_counter == 2:
                                vkickmap    = numpy.copy(kick_map);
                                table_cols2 = haxis;
                                table_rows2 = vaxis;
                        block_lines += 1;


        return el_length, hkickmap, vkickmap, table_cols1, table_rows1, table_cols2, table_rows2

    def sorted_table(self, table_in, sorted_index, order_axis):
        # numpy.asfortranarray makes a copy of contiguous memory positions
        table_out = numpy.copy(table_in);
        for i,iis in zip(range(len(sorted_index)), sorted_index):
            if order_axis == 'col': table_out[:,i] = table_in[ : ,iis]
            if order_axis == 'row': table_out[i,:] = table_in[iis, : ]
        table_out2 = numpy.asfortranarray(table_out)
        return table_out2

    def __init__(self, family_name: str, Nslice: float, Filename_in: str, Energy: float, **kwargs):
        """
        Args:
            family_name:    family name
            Nslice:         number of slices in ... tracking??? integrator ???
            Filename_in:    input filename
            Energy:         particle energy in GeV

        Default PassMethod: ``IdTablePass``
        """
        # 2023jan18 fix bug with element print
        # 2023jan15 first release
        # orblancog
        ## read the input data
        el_length, hkickmap, vkickmap, table_cols1, table_rows1, table_cols2, table_rows2 \
                = self.readRadiaFieldMap(Filename_in);

        ## set to float
        table_cols1array  = numpy.array(table_cols1, dtype='float64');
        table_rows1array  = numpy.array(table_rows1, dtype='float64');
        table_cols2array  = numpy.array(table_cols2, dtype='float64');
        table_rows2array  = numpy.array(table_rows2, dtype='float64');

        ## Reorder table_axes and kick maps
        cols1sorted_index = numpy.argsort(table_cols1array); table_cols1array.sort();
        rows1sorted_index = numpy.argsort(table_rows1array); table_rows1array.sort();
        cols2sorted_index = numpy.argsort(table_cols2array); table_cols2array.sort();
        rows2sorted_index = numpy.argsort(table_rows2array); table_rows2array.sort();
        hkickmap_a  = self.sorted_table(hkickmap,   cols1sorted_index, 'col');
        hkickmap    = self.sorted_table(hkickmap_a, rows1sorted_index, 'row');
        vkickmap_a  = self.sorted_table(vkickmap,   cols2sorted_index, 'col');
        vkickmap    = self.sorted_table(vkickmap_a, rows2sorted_index, 'row');

        ## Field to kick factors
        Brho     = 1e9 * Energy/clight
        factor   =  1/(Brho**2)
        xkick    = factor * hkickmap
        ykick    = factor * vkickmap
        factor1  = -1/(Brho)
        xtable   = table_cols1array
        ytable   = table_rows1array

        ## Create element properties
        kwargs.setdefault('PassMethod', 'IdTablePass')
        kwargs.setdefault('Filename_in', Filename_in)
        kwargs.setdefault('Energy',      Energy)
        kwargs.setdefault('Nslice',      Nslice)
        kwargs.setdefault('Length',      el_length)
        kwargs.setdefault('xkick',       xkick)
        kwargs.setdefault('ykick',       ykick)
        kwargs.setdefault('xtable',      xtable)
        kwargs.setdefault('ytable',      ytable)

        ## Set Element inherited properties
        # pyat issue #522
        #elem = Element('name', PassMethod='IdTablePass', xkick=x0, ykick=x1, \
        #                       xtable=x2, ytable=x3, Nslice=x4)
        super(InsertionDevice, self).__init__(family_name, **kwargs)

