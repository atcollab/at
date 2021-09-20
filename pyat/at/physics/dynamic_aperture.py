import at
import copy
import at.plot
from itertools import compress
import numpy as np
import math
import pickle
from scipy.io import savemat
import time

__all__ = ['Acceptance6D', 'dynamic_aperture', 'off_energy_dynamic_aperture', 'momentum_acceptance']


class Acceptance6D(object):
    """
    Compute acceptance
    """

    verbose = True

    # possible acceptance planes
    planes = ('x', 'xp', 'y', 'yp', 'delta', 'ct')

    # 2D combination of planes
    modes = ('x-y',
             'delta-x',
             'xp-yp',
             'x-xp',
             'y-yp',
             'ct-delta',
             '6D')

    dict_units = {'x': [1e3, 'mm'],
                  'xp': [1e3, 'mrad'],
                  'y': [1e3, 'mm'],
                  'yp': [1e3, 'mrad'],
                  'delta': [1e2, '%'],
                  'ct': [1e3, 'mm'],
                  }

    dict_def_range = {'x': [-3e-2, 3e-2],
                      'xp': [-1e-3, 1e-3],
                      'y': [-1e-2, 1e-2],
                      'yp': [-1e-3, 1e-3],
                      'delta': [-2e-1, 2e-1],
                      'ct': [-1e-1, 1e-1],
                      }
    """
    @property
    def dict_def_range(self):
        return self._dict_def_range

    @dict_def_range.setter
    def dict_def_range(self, new_value):
        self._dict_def_range = new_value
        self.test_points_grid()  # recompute grid points if ranges are changed
    """

    def __init__(self,
                 ring,
                 dp=0.0,
                 mode='x-y',
                 start_index=0,
                 grid_mode='grid',
                 compute_range=False,
                 n_turns=2**10,
                 search_divider=10,
                 use_radiation=True,
                 parallel = False
                ):
        """
        Acceptance6D computes the 6D acceptance for a pyAT lattice

        :param ring: pyAT lattice
        :param dp: momentum deviation
        :param mode: mode for computation. None = 6D, '6D', 'x-y' (default), 'delta-x'. 'x-xp',...
        :param start_index: index in ring where to start the computation
        :param grid_mode: 'grid' or 'radial' set of test points to compute Acceptance
        :param compute_range: compute range for each plane with more than 1 point to scan
        :param nturns: number of turns that a particle must survive
        :param search_divider: division of range used by recursive range
        :param use_radiation: default true (turn on radiation, in future version, use radiation state of input lattice)
        :param parallel: default False if True, use patpass when usefull
        : n_point: dictionary to determine the # points to scan in each dimension
        : dict_def_range: default range for each dimension
        : dict_unit:  default units for each dimension

        """

        self.parallel_computation = parallel

        # rotate lattice to given element
        self.ring = at.Lattice(ring[start_index:] + ring[:start_index])

        self.dp = dp
        self.mode = mode

        # if explicit input, turn off radiation
        if not use_radiation:
            self.ring.radiation_off()
        else:
            # if no input use state of radiation of the input lattice
            if self.verbose:
                print('using input lattice radiation state: {}'.format(ring.radiation))

        use_radiation = ring.radiation

        # define orbit about wich to compute DA
        self.compute_orbit()

        # define coordinates dictionary structure
        self.coordinates = {'x': [],
                            'xp': [],
                            'y': [],
                            'yp': [],
                            'delta': [],
                            'ct': []
                            }
        self.survived = []

        self.number_of_turns = n_turns
        self.search_divider = search_divider

        # define grid modes for later use
        self.grid_modes = ['grid', 'radial']
        self.grid_mode = grid_mode
        if not (self.grid_mode in self.grid_modes):
            raise('grid mode must be: grid or radial')

        # point to scan in each dimension
        if mode == 'x-y':
            self.n_points = {'x': 13,
                             'xp': 1,
                             'y': 13,
                             'yp': 1,
                             'delta': 1,
                             'ct': 1
                             }
        elif mode == 'delta-x':
            self.n_points = {'x': 13,
                             'xp': 1,
                             'y': 1,
                             'yp': 1,
                             'delta': 13,
                             'ct': 1
                             }
        elif mode == 'ct-delta':
            self.n_points = {'x': 1,
                             'xp': 1,
                             'y': 1,
                             'yp': 1,
                             'delta': 13,
                             'ct': 13
                             }
        elif mode == 'x-xp':
            self.n_points = {'x': 13,
                             'xp': 13,
                             'y': 1,
                             'yp': 1,
                             'delta': 1,
                             'ct': 1
                             }
        elif mode == 'y-yp':
            self.n_points = {'x': 1,
                             'xp': 1,
                             'y': 13,
                             'yp': 13,
                             'delta': 1,
                             'ct': 1
                             }
        elif mode == 'xp-yp':
            self.n_points = {'x': 1,
                             'xp': 13,
                             'y': 1,
                             'yp': 13,
                             'delta': 1,
                             'ct': 1
                             }
        elif mode == '6D':
            self.n_points = {'x': 5,
                             'xp': 5,
                             'y': 5,
                             'yp': 5,
                             'delta': 5,
                             'ct': 5
                             }
        else:  # default for no mode = x-y
            self.n_points = {'x': 13,
                             'xp': 1,
                             'y': 13,
                             'yp': 1,
                             'delta': 1,
                             'ct': 1
                             }

        self.test_points_grid(compute_limits=compute_range)

        pass

    def compute_orbit(self):
        """
        computes orbit for given dp
        """

        """
        L.F.
        you have to change the RF frequency by alpha * f0 * dpp, 
        compute the 6D closed orbit
        """
        elems = self.ring.get_elements('CAV*')
        ind_rf = [self.ring.index(elem) for elem in elems]
        # set RF frequency
        self.ligth_speed = 2.99792458e08
        self.Circumference = self.ring.s_range[-1]
        self.harm = self.ring[self.ind_rf[0]].HarmNumber
        self.rf_frequency = self.ring[self.ind_rf[0]].Frequency


        rad = self.ring.radiation
        self.ring.radiation_off()

        self.orbit = self.ring.find_orbit4(dp=self.dp)  # dpp is added to orbit here

        if rad:
            self.ring.radiation_on()

        return self.orbit

    def recursive_search_of_directional_limit(self, direction=(1, 0, 0, 0, 0, 0),
                                              number_of_recursions=3,
                                              back_step=2):
        """

        recursively search for maximum coordinate along a given direction where particle survive

        :param direction: 6D array to deterine a given direction
        :param number_of_recursions:  number of recursion. at each recursion, the last step in the search is replaced by a finer mesh
        :param back_step:  number of steps to decrement to start with next search.
        :return: 6D coordinates array of last surviving particle
        """

        # disable parallel computation
        init_parcomp = self.parallel_computation
        self.parallel_computation = False

        # define initial coordinate
        coord0 = {'x': 1e-6, 'xp': 0.0, 'y': 1e-6, 'yp': 0.0, 'delta': 0.0, 'ct': 0.0}

        # define step in each plane
        step = {pl: direction[ip]*(self.dict_def_range[pl][1] - self.dict_def_range[pl][0])
                for ip, pl in enumerate(self.planes)}

        limit = [v for _, v in coord0.items()]

        for ir in range(number_of_recursions):

            if self.verbose:
                print('search {d} recursively step {s}/{ts}'.format(d=direction, s=ir, ts=number_of_recursions))

            # define initial coordinates for search
            #if ir > 0:
            #    for pl in self.planes:
            #        coord[pl] -= step[pl]  # go back by one of the previous size steps.
            #else:
            #    coord = copy.deepcopy(coord0)
            if ir == 0:
                coord = {'x': 1e-6, 'xp': 0.0, 'y': 1e-6, 'yp': 0.0, 'delta': 0.0, 'ct': 0.0}

            if self.verbose:
                print(coord)

            # reduce step (initial = full range)
            step = {k: v / self.search_divider for k, v in step.items()}

            # do not recompute already computed point
            if ir > 0:
                for pl in self.planes:
                    coord[pl] += step[pl]  # go back by one of the previous size steps.

            # search limit
            while self.test_survived(coord)[0]:
                for pl in self.planes:
                    coord[pl] += step[pl]
                    # update tested points
                    self.coordinates[pl].append(coord[pl])
                if self.verbose:
                    print('step forward')
                    print(coord)

                self.survived.append(True)

            # last survived is previous point
            for pl in self.planes:
                _c = back_step * step[pl]
                # limit to 0.0
                if coord[pl] > 0.0:
                    if coord[pl]-_c < 0.0:
                        _c = coord[pl]
                if coord[pl] < 0.0:
                    if coord[pl]-_c > 0.0:
                        _c = coord[pl]
                # assign new start coord
                coord[pl] -= _c

            if self.verbose:
                print('step back')
                print(coord)

            # last coordinates are lost and add an additional test point lost (for display and contour computation
            if self.survived:
                self.survived[-1] = False

            coord_lost = copy.deepcopy(coord)
            for pl in self.planes:
                coord_lost[pl] += step[pl]
                # update tested points
                self.coordinates[pl].append(coord_lost[pl])
            self.survived.append(False)

            # define limit 6 element array
            limit = [copy.deepcopy(v) for _, v in coord.items()]

            if self.verbose:
                print('present limit: {l}'.format(l=limit))

        # restore initial parallel computation value
        self.parallel_computation = init_parcomp

        return limit

    def compute_range(self):
        """
        computes maximum +/- range for all dimensions with more than 1 point to scan
        result is stored in dict_def_range
        ranges are limited to [-1 1] in any case without notice.
        :return:
        """
        # disable parallel computation
        init_parcomp = self.parallel_computation
        self.parallel_computation = False

        # relevant_planes = []
        # if self.mode:
        #     relevant_planes = self.mode.split('-', 2)
        coord0 = {'x': 1e-6, 'xp': 0.0, 'y': 1e-6, 'yp': 0.0, 'delta': 0.0, 'ct': 0.0}

        for p in self.planes:
            if self.n_points[p] > 1:

                step = (self.dict_def_range[p][1] - self.dict_def_range[p][0]) / self.search_divider

                coord = {'x': 1e-6, 'xp': 0.0, 'y': 1e-6, 'yp': 0.0, 'delta': 0.0, 'ct': 0.0} # copy.deepcopy(coord0)
                while self.test_survived(coord)[0] and coord[p]< 1.0:
                    coord[p] += step
                self.dict_def_range[p][1] = coord[p]+step

                coord = {'x': 1e-6, 'xp': 0.0, 'y': 1e-6, 'yp': 0.0, 'delta': 0.0, 'ct': 0.0} # copy.deepcopy(coord0)
                while self.test_survived(coord)[0] and coord[p]> -1.0:
                    coord[p] -= step
                self.dict_def_range[p][0] = coord[p]-step
                if self.verbose:
                    print('Computed range for {pl}: [{mi}, {ma}]'.format(pl=p,
                                                                         mi=self.dict_def_range[p][0],
                                                                         ma=self.dict_def_range[p][1]))

        # restore initial parallel computation value
        self.parallel_computation = init_parcomp

        pass

    def test_points_grid(self, compute_limits=False):
        """
        define fixed grid of test points
        the grid of points is stored in a dictionary of coordinates, one array for each dimension

        :param compute_limits calls compute_range to deterine the maximum ranges before defining the grid

        :return:
        """

        if compute_limits:
            self.compute_range()

        # define grid
        if self.grid_mode == 'radial':

            d_ = {'x': [], 'xp': [], 'y': [], 'yp': [], 'delta': [], 'ct': []}

            # find n_points more than 1
            two_planes = [p for p in self.planes if self.n_points[p] > 1]

            if len(two_planes) != 2:
                raise IndexError('There are less or more than 2 planes with more than 1 point to scan')

            ellipse_axes = [self.dict_def_range[two_planes[0]][1],
                            self.dict_def_range[two_planes[1]][1]]

            if self.dict_def_range[two_planes[1]][0] == 0:
                theta = np.linspace(0, math.pi, self.n_points[two_planes[0]])
            else:
                theta = np.linspace(0, 2*math.pi, self.n_points[two_planes[0]])

            ip = 0

            for ind_p, p in enumerate(self.planes):  # must loop all plane to define all 6D ranges (1 point)

                if self.n_points[p] > 1:  # if here, there are only 2 planes available

                    # scan predefined range
                    for ea in np.sqrt(np.linspace(ellipse_axes[ip]**2/self.n_points[two_planes[1]],
                                          ellipse_axes[ip]**2, self.n_points[two_planes[1]])):
                        if ip == 0:
                            d_[p].append([ea * math.cos(t) for t in theta])
                        elif ip == 1:
                            d_[p].append([ea * math.sin(t) for t in theta])

                    ip += 1

                else:  # this dimension is not part of the scan
                    d_[p].append([(self.dict_def_range[p][0] + self.dict_def_range[p][1]) / 2])

            flat_dd = []
            for k, val in d_.items():
                flat_dd.append([item for sublist in val for item in sublist])
                print('{} has {} elements'.format(k, len(flat_dd[-1])))

            num_points = max([len(dd) for dd in flat_dd])

            # define lists for DA scan
            for ip, p in enumerate(self.planes):
                if p in two_planes:
                     self.coordinates[p] = [item for sublist in d_[p] for item in sublist]
                else:
                    self.coordinates[p] = [item for sublist in d_[p] for item in sublist]*num_points

                # print(len(self.coordinates[p]))

        else:

            d_ = []

            for p in self.planes:
                if self.n_points[p] > 1:
                    # scan predefined range
                    d_.append(np.linspace(self.dict_def_range[p][0],
                                          self.dict_def_range[p][1],
                                          self.n_points[p])
                              )
                else:  # this dimension is not part of the scan
                    d_.append((self.dict_def_range[p][0] +
                               self.dict_def_range[p][1])/2
                              )

            # define mesh of points
            xx, xpxp, yy, ypyp, deltadelta, ctct = np.meshgrid(d_[0],
                                                               d_[1],
                                                               d_[2],
                                                               d_[3],
                                                               d_[4],
                                                               d_[5],)

            self.coordinates['x'] = xx.flatten()
            self.coordinates['xp'] = xpxp.flatten()
            self.coordinates['y'] = yy.flatten()
            self.coordinates['yp'] = ypyp.flatten()
            self.coordinates['delta'] = deltadelta.flatten()
            self.coordinates['ct'] = ctct.flatten()

        if len(self.coordinates['x']) == 0:
            raise IOError('grid mode must be grid or radial')

        self.survived = [False]*len(self.coordinates['x'])

        pass

    def test_survived(self, coordinates):
        """
        test if a set of particle coordinates survived
        returns a boolean if the coordinates survived
        """

        """
        if delta != 0.0:
            # recompute orbit
            self.ring.radiation_off()
            orbit_delta = self.ring.find_sync_orbit()
            self.ring.radiation_on()

            test_coord = list(orbit_delta[0])
        else:
            test_coord = list(self.orbit[0])
        """

        # nothing to do if no coordinate to test
        if len(coordinates) == 0:
            return []

        # if coordinates to test:

        # create 6xN matrix add orbit (with dpp) to each coordinate
        rin = np.concatenate(([coordinates['x'] + self.orbit[0][0]],
                              [coordinates['xp'] + self.orbit[0][1]],
                              [coordinates['y'] + self.orbit[0][2]],
                              [coordinates['yp'] + self.orbit[0][3]],
                              [coordinates['delta'] + self.orbit[0][4]],
                              [coordinates['ct'] + self.orbit[0][5]]),
                              axis=0)

        # track coordinates for N turns
        if not self.parallel_computation:

            t, losses = at.atpass(self.ring,
                                  copy.deepcopy(np.asfortranarray(rin)),
                                  self.number_of_turns,
                                  refpts=np.array(np.uint32(0)),  # np.array(np.uint32([len(self.ring)])))
                                  losses=True)

        else:
            print('parallel computation')

            # track coordinates
            t, losses = at.patpass(self.ring,
                                   copy.deepcopy(rin),
                                   self.number_of_turns,
                                   refpts=np.array(np.uint32(0)),  # np.array(np.uint32([len(self.ring)])))
                                   losses=True)

        survived = [not s for s in losses['islost']]

        # print if survived for each test particle
        if self.verbose:
            if type(coordinates['x']) is float:
                print('[{x:2.2g}, {xp:2.2g}, {y:2.2g}, {yp:2.2g}, {delta:2.2g}, {ct:2.2g}] '
                      '| {tot:02d} coord: {surv} for {nt} turns'.format(
                    tot=1, surv=survived, nt=self.number_of_turns,
                    x=coordinates['x'], y=coordinates['y'], xp=coordinates['xp'],
                    ct=coordinates['ct'], yp=coordinates['yp'], delta=coordinates['delta']))
            else:
                for x, xp, y, yp, delta, ct, s in zip(coordinates['x'],
                                                      coordinates['xp'],
                                                      coordinates['y'],
                                                      coordinates['yp'],
                                                      coordinates['delta'],
                                                      coordinates['ct'], survived):
                    print('[{x:2.2g}, {xp:2.2g}, {y:2.2g}, {yp:2.2g}, {delta:2.2g}, {ct:2.2g}] '
                          '| {tot:02d} coord: {surv} for {nt} turns'.format(
                        tot=len(coordinates['x']), surv=s, nt=self.number_of_turns,
                        x=x, y=y, xp=xp, ct=ct, yp=yp, delta=delta))

        """
        test_coord = copy.deepcopy(list(self.orbit[0]))
        # print(test_coord[4])
        # print(self.dp)
        # test_coord[4] = test_coord[4] + self.dp

        rin = np.asfortranarray(np.zeros((6, 1)))

        for ip, p in enumerate(self.planes):

            if p == 'x':
                test_coord[ip] += coordinates['x']
            if p == 'xp':
                test_coord[ip] += coordinates['xp']
            if p == 'y':
                test_coord[ip] += coordinates['y']
            if p == 'yp':
                test_coord[ip] += coordinates['yp']
            if p == 'delta':
                test_coord[ip] += coordinates['delta']
            if p == 'ct':
                test_coord[ip] += coordinates['ct']

        rin[0][0] = test_coord[0]
        rin[1][0] = test_coord[1]
        rin[2][0] = test_coord[2]
        rin[3][0] = test_coord[3]
        rin[4][0] = test_coord[4]
        rin[5][0] = test_coord[5]

        # track
        t = at.atpass(self.ring,
                      copy.deepcopy(rin),
                      self.number_of_turns,
                      refpts=np.array(np.uint32(0)))  #np.array(np.uint32([len(self.ring)])))
        # atpass help:
        # rout: 6 x n_particles x n_refpts x n_turns Fortran - ordered numpy array of particle coordinates

        if self.verbose:
            print('test if {r} survives {n} turns: {s}'.format(
                r=rin.T, n=self.number_of_turns, s=not(np.isnan(t[0][0][0][-1]))))
        """

        return survived

    def compute(self):
        """
        compute 2D DA limits

        fills the self.survived property with booleans corresponding to each coordinate in coordinates

        :return: h_s, v_s lists of dim1 and dim2 coordinates defining the 2D limit of Acceptance
        """

        # display summary of work to do
        print('Computing acceptance {m}'.format(m=self.mode))

        [print('test {np:02d} points in {pl}'.format(np=self.n_points[p], pl=p))
         for ip, p in enumerate(self.planes)]

        # test all coordinates at once
        self.survived = self.test_survived(self.coordinates)

        h_s = []
        v_s = []

        if self.mode in self.modes[0:-1]:
            # simplify 6D structure to 2D
            h, v, sel = self.select_test_points_based_on_mode()

            # find maximum of each column and return as border
            h_s, v_s = self.get_border(h, v, sel)

            try:
                h_s, v_s = self.get_border(h, v, sel)
            except Exception:
                print('DA limit could not be computed (probably no closed contour)')
        else:
            print('DA limit could not be computed for this mode')

        return h_s, v_s

    def get_border(self, h, v, sel):
        """
        find border of list of points.
        works only for grid mode.
        """
        h_border = []
        v_border = []

        # [print(h_, v_, s_) for h_, v_, s_ in zip(h, v, sel)]

        # loop columns of grid
        if self.grid_mode == 'grid':

            for hc in np.unique(h):
                col=[]
                for i, h_ in enumerate(h):
                    if h_ == hc and sel[i]:
                        col.append(v[i])
                if len(col) > 0:
                    h_border.append(hc)
                    v_border.append(np.max(col))
                    h_border.insert(0,hc)
                    v_border.insert(0,np.min(col))
        elif self.grid_mode == 'radial':
            # find radial grid extremes.  not implemented
            print('radial mode, no border computed')
            pass
        else:
            print('grid_mode must be grid or radial')

        # remove bottom line if any. not implemented

        return h_border, v_border

    def select_test_points_based_on_mode(self):
        """
        reduces the 6D coordinated according to the mode to 3 lists for plotting:
        horizontal, vertical and boolean array of survival.
        """

        if self.mode == None or self.mode=='6D':
            if self.verbose:
                print('no selection of test points for modes None and 6D')
            return [], [], []

        if self.mode == 'x-y':
            h = self.coordinates['x']
            v = self.coordinates['y']

            sel = [a and b == 0 and c == 0 and d == 0 and e == 0 for a, b, c, d, e in zip(
                self.survived,
                self.coordinates['delta'],
                self.coordinates['xp'],
                self.coordinates['yp'],
                self.coordinates['ct'])]

        if self.mode == 'delta-x':
            h = self.coordinates['delta']
            v = self.coordinates['x']

            sel = [a and b == 0 and c == 0 and d == 0 and e == 0 for a, b, c, d, e in zip(
                self.survived,
                self.coordinates['y'],
                self.coordinates['ct'],
                self.coordinates['yp'],
                self.coordinates['xp'])]

        if self.mode == 'xp-yp':
            h = self.coordinates['xp']
            v = self.coordinates['yp']

            sel = [a and b == 0 and c == 0 and d == 0 and e == 0 for a, b, c, d, e in zip(
                self.survived,
                self.coordinates['y'],
                self.coordinates['ct'],
                self.coordinates['delta'],
                self.coordinates['x'])]

        if self.mode == 'ct-delta':
            h = self.coordinates['ct']
            v = self.coordinates['delta']

            sel = [a and b == 0 and c == 0 and d == 0 and e == 0 for a, b, c, d, e in zip(
                self.survived,
                self.coordinates['y'],
                self.coordinates['x'],
                self.coordinates['xp'],
                self.coordinates['yp'])]

        if self.mode == 'x-xp':
            h = self.coordinates['x']
            v = self.coordinates['xp']

            sel = [a and b == 0 and c == 0 and d == 0 and e == 0 for a, b, c, d, e in zip(
                self.survived,
                self.coordinates['y'],
                self.coordinates['ct'],
                self.coordinates['delta'],
                self.coordinates['yp'])]

        if self.mode == 'y-yp':
            h = self.coordinates['y']
            v = self.coordinates['yp']

            sel = [a and b == 0 and c == 0 and d == 0 and e == 0 for a, b, c, d, e in zip(
                self.survived,
                self.coordinates['x'],
                self.coordinates['ct'],
                self.coordinates['delta'],
                self.coordinates['xp'])]

        return h, v, sel

    pass



def off_energy_dynamic_aperture(sr_ring,
                                deltaps=np.linspace(-0.1, 0.1, 11),
                                step_size=1e-3,
                                search=True,
                                inject_from_inside=True,
                                n_turns=2**10,
                                start_index=0,
                                file_name_save=None,
                                num_recursions=5):
    """
    maximum negative horizontal DA for sereral energy deviations

    :param sr_ring: pyAT lattice
    :param deltaps: list of deltap/p to test
    :param step_size: step for search (not used if search = True)
    :param search: if True recursively search for limit
    :param inject_from_inside: if True look fo maximum negative x, if False look for maximum positive x
    :param num_recursions: number of recursion to refine limit search
    :param n_turns: number of turns for survival
    :param start_index: index where to compute off-energy DA
    :param file_name_save: if given, save to file
    :return: max_neg_x list of maximum negative da for each deltap/p required
    :return: deltaps, da: for use with plotting function in at.plot.dynamic_aperture.plot_off_energy_dynamic_aperture
    """
    da = Acceptance6D(copy.deepcopy(sr_ring), start_index=start_index, n_turns=n_turns)
    da.verbose = False

    max_neg_x = []

    coord0 = {'x': 1e-6, 'xp': 0.0, 'y': 1e-6, 'yp': 0.0, 'delta': 0.0, 'ct': 0.0}

    if search:
        da.coordinates = {'x': [], 'xp': [], 'y': [], 'yp': [], 'delta': [], 'ct': []}
        da.survived = []

    ring_rad = da.ring.radiation

    for deltap in deltaps:
        print('{d:2.1f}%'.format(d=deltap*100))

        # change dp
        da.dp = deltap

        # refresh reference orbit with new dp
        da.compute_orbit()

        if search:
            if inject_from_inside:
                lim = da.recursive_search_of_directional_limit(
                    direction=(-1, 0, 0, 0, 0, 0),
                    number_of_recursions=num_recursions)
            else:
                lim = da.recursive_search_of_directional_limit(
                    direction=(1, 0, 0, 0, 0, 0),
                    number_of_recursions=num_recursions)
            max_neg_x.append(lim[0])
        else:
            coord = copy.deepcopy(coord0)
            while da.test_survived(coord)[0]:
                if inject_from_inside:
                    coord['x'] -= step_size
                else:
                    coord['x'] += step_size
            max_neg_x.append(coord['x'])

    if file_name_save:
        # save data
        file_to_store = open(file_name_save + '.pickle', "wb")
        pickle.dump(max_neg_x, file_to_store)
        pickle.dump(deltaps, file_to_store)
        pickle.dump(n_turns, file_to_store)
        pickle.dump(start_index, file_to_store)
        file_to_store.close()

        m_dict = {'max_neg_x': max_neg_x, 'deltaps': deltaps, 'n_turns': n_turns, 'start_index': start_index}
        savemat(file_name_save + '.mat', m_dict)


    return max_neg_x, deltaps, da


def momentum_acceptance(sr_ring, n_turns=2**10, ref_pts=None, file_name_save=None, num_recursions=5):
    """
    compute momentum acceptance at given reference positions along the lattice
    :param sr_ring: pyAT lattice
    :param n_turns: number of turns to survive
    :param ref_pts: list of positions at which to compute momentum acceptance
    :param num_recursions: number of recursion to refine search
    :param file_name_save: if given, save figure

    :return: [[max deltap/p], [min deltap/p]], s locations
    :return: da, for use with plotting function in at.plot.dynamic_aperture.plot_momentum_acceptance
    """

    mom_acc = [[], []]

    for el in ref_pts:

        print('ind {d}'.format(d=el))

        da = Acceptance6D(copy.deepcopy(sr_ring), start_index=el)
        da.number_of_turns = n_turns
        da.verbose = False

        da.coordinates = {'x': [], 'xp': [], 'y': [], 'yp': [], 'delta': [], 'ct': []}
        da.survived = []

        # pos acceptance
        lim = da.recursive_search_of_directional_limit(
            direction=(0, 0, 0, 0, 1, 0),
            number_of_recursions=num_recursions)
        mom_acc[0].append(lim[4])
        # neg acceptance
        lim = da.recursive_search_of_directional_limit(
            direction=(0, 0, 0, 0, -1, 0),
            number_of_recursions=num_recursions)
        mom_acc[1].append(lim[4])

    s = sr_ring.get_s_pos(ref_pts)

    if file_name_save:
        # save data
        file_to_store = open(file_name_save + '.pickle', "wb")
        pickle.dump(mom_acc, file_to_store)
        pickle.dump(s, file_to_store)
        pickle.dump(ref_pts, file_to_store)
        pickle.dump(n_turns, file_to_store)
        file_to_store.close()

        m_dict = {'mom_ac': mom_acc, 's': s, 'n_turns': n_turns, 'ref_pts': ref_pts}
        savemat(file_name_save + '.mat', m_dict)

    return mom_acc, s, da


def dynamic_aperture(sr_ring,
                     dp=0.0,
                     n_turns=2**10,
                     start_index=0,
                     n_radii=13,
                     n_theta=13,
                     grid_mode='radial',
                     parallel=False,
                     search=True,
                     file_name_save=None,
                     num_recursions=5):
    """
    compute Dynamic aperture x-y searching recursively

    :param sr_ring: pyAT lattice
    :param dp: momentum deviation
    :param n_turns: number of turns to survive
    :param start_index: index in sr_ring where to start the tracking
    :param n_radii: number of radial points to test  (not used if search = True)
    :param n_theta: number of angular divisions
    :param grid_mode: 'grid' or 'radial' fixed test points grid
    :param search: True (default), ignore grid and search recursively for maximum along given radial direction
    :param num_recursions: number of recursion for search
    :param parallel : uses patpass rather than atpass
    :param file_name_save: if given, save to file

    :return: h,v horizontal and vertical coordinates of DA limit
    :return: da, search : for use with plotting function in at.plot.dynamic_aperture.plot_dynamic_aperture
    """
    h = []
    v = []

    da = Acceptance6D(copy.deepcopy(sr_ring), mode='x-y', grid_mode=grid_mode, start_index=start_index)
    da.number_of_turns = n_turns
    da.dp = dp
    da.verbose = False
    da.parallel_computation = parallel
    
    if grid_mode == 'radial':

        da.compute_range()  # implement an init at change of mode, npoint, or range.

        da.n_points['x'] = n_theta  # used as theta, if grid radial
        da.n_points['y'] = n_radii  # used as radii, if grid radial
        da.dict_def_range['y'][0] = 0.0  # make y single sided

        if search:
            theta = np.linspace(0, math.pi, da.n_points['x'])
            lim = []

            # list of coordinates tested is emptied (will be filled by recursive search)
            da.coordinates = {'x': [], 'xp': [], 'y': [], 'yp': [], 'delta': [], 'ct': []}
            da.survived = []

            for t in theta:
                lim = da.recursive_search_of_directional_limit(
                           direction=(math.cos(t), 0, math.sin(t), 0, 0, 0),
                           number_of_recursions=num_recursions)

                h.append(lim[0])
                v.append(lim[2])

            # # convert coordinates to np array
            # for k, v in da.coordinates.items():
            #     v = np.array(v)

        else:
            da.compute_range()  # implement an init at change of mode, npoint, or range.

            da.n_points['x'] = n_theta  # used as theta, if grid radial
            da.n_points['y'] = n_radii  # used as radii, if grid radial
            da.dict_def_range['y'][0] = 0.0  # make y single sided

            da.test_points_grid()
            h, v = da.compute()

    if grid_mode == 'grid':
        da.compute_range()  # implement an init at change of mode, npoint, or range.
        da.n_points['x'] = n_theta  # used as theta, if grid radial
        da.n_points['y'] = n_radii  # used as radii, if grid radial
        da.dict_def_range['y'][0] = 0.0  # make y single sided
        da.test_points_grid()
        h, v = da.compute()

    # plot
    if file_name_save:

        # save data
        file_to_store = open(file_name_save + '.pickle', "wb")
        pickle.dump(h, file_to_store)
        pickle.dump(v, file_to_store)
        pickle.dump(n_turns, file_to_store)
        pickle.dump(start_index, file_to_store)
        pickle.dump(dp, file_to_store)
        file_to_store.close()

        m_dict = {'h': h, 'v': v, 'n_turns': n_turns, 'start_index': start_index, 'dp': dp}
        savemat(file_name_save + '.mat', m_dict)

    return h, v, da, search
