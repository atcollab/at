import at
from at.physics.dynamic_aperture import Acceptance6D, dynamic_aperture, off_energy_dynamic_aperture, momentum_acceptance
import time
import copy
import pytest


def test_acceptance_x_xp(hmba_lattice):

    sr_ring = at.Lattice(hmba_lattice*32)

    sr_ring.radiation_off()

    # x-xp acceptance at index 120 with dpp -0.5% for 32 turns
    print('x-xp acceptance, using class directly')
    da = Acceptance6D(copy.deepcopy(sr_ring),
                      mode='x-xp', start_index=120, n_turns=2**5, dp=-0.005)

    da.verbose = True
    da.parallel_computation = True

    da.dict_def_range = {'x': [-2e-3, 2e-3],
                      'xp': [-1e-4, 1e-4],
                      'y': [-5e-3, 5e-3],
                      'yp': [-1e-3, 1e-3],
                      'delta': [-2e-1, 2e-1],
                      'ct': [-1e-1, 1e-1],
                      }
    da.n_points = {'x': 7,
                     'xp': 7,
                     'y': 1,
                     'yp': 1,
                     'delta': 1,
                     'ct': 1
                     }
    # search range for scan
    da.compute_range()
    # define grid of test particles in the range
    da.test_points_grid()
    # test if particles of the grid survive or not and produce a curve arround the survived ones
    h, v = da.compute()
    [print(h_, v_) for h_, v_ in zip(h, v)]
    print(da.survived)

    pass

def test_dynamic_aperture(hmba_lattice):

    sr_ring = at.Lattice(hmba_lattice*32)

    sr_ring.radiation_off()

    # on-energy DA
    t = time.time()
    print('On-energy DA')
    h, v, _, _= dynamic_aperture(sr_ring, n_turns=2**5, n_radii=15, num_recursions=4,
                            file_name_save='') # 'on_en_da')
    [print(h_,v_) for h_, v_ in zip(h,v)]
    elapsed = time.time() - t
    print('On-energy DA took {s:2.1f}s'.format(s=elapsed))
    pass

def test_dynamic_aperture_grid(hmba_lattice):

    sr_ring = at.Lattice(hmba_lattice*32)

    sr_ring.radiation_off()

    # on-energy DA
    t = time.time()
    print('On-energy DA grid')
    h, v, da, _= dynamic_aperture(sr_ring, grid_mode='grid',
                                 n_turns=2**7, n_radii=25, n_theta=25, num_recursions=4,
                            file_name_save='') # 'on_en_da')
    [print(h_,v_) for h_, v_ in zip(h,v)]
    print(da.survived)
    print(len(da.survived))
    print(len([a for a in da.survived if a]))
    elapsed = time.time() - t
    print('On-energy DA grid took {s:2.1f}s'.format(s=elapsed))
    pass


def test_off_en_dynamic_aperture(hmba_lattice):
    sr_ring = at.Lattice(hmba_lattice * 32)

    sr_ring.radiation_off()
    # off-en DA
    t = time.time()
    print('Off-energy DA')
    off_energy_dynamic_aperture(sr_ring, n_turns=2**4,
                                file_name_save='', # 'off_en_da_inside',
                                inject_from_inside=True,
                                num_recursions=3)
    elapsed = time.time() - t
    print('Off-energy DA took {s:2.1f}s'.format(s=elapsed))
    pass

def test_mom_acc(hmba_lattice):
    sr_ring = at.Lattice(hmba_lattice * 32)

    sr_ring.radiation_off()
    # momentum acceptance
    t = time.time()
    print('Momentum acceptance')
    momentum_acceptance(sr_ring, n_turns=2**4,
                        file_name_save='', # 'mom_acc',
                        ref_pts=range(0, 100, 10), num_recursions=3)
    elapsed = time.time() - t
    print('Momentum acceptance took {s:2.1f}s'.format(s=elapsed))

    # plt.show()

    pass


def test_search(hmba_lattice):

    n_turns = 100
    dpp = 0.0
    start_index = 0
    num_recursions = 5

    sr_ring = hmba_lattice

    # sr_lattice_file = '../test_matlab/hmba.mat'
    # sr_lattice_variable = 'ring'
    # sr_ring = at.load_mat(sr_lattice_file, mat_key=sr_lattice_variable)

    sr_ring.radiation_off()

    h = []
    v = []

    da = Acceptance6D(copy.deepcopy(sr_ring), mode='x-y', start_index=start_index)
    da.number_of_turns = n_turns
    da.dpp = dpp
    da.verbose = True
    da.compute_range()  # implement an init at change of mode, npoint, or range.

    da.n_points['x'] = 13  # used as theta, if grid radial
    da.n_points['y'] = 5  # used as radii, if grid radial
    da.dict_def_range['y'][0] = 0.0  # make y single sided

    # list of coordinates tested is emptied (will be filled by recursive search)
    da.coordinates = {'x': [], 'xp': [], 'y': [], 'yp': [], 'delta': [], 'ct': []}
    da.survived = []

    lim = da.recursive_search_of_directional_limit(
        direction=(1, 0, 0, 0, 0, 0),
        number_of_recursions=num_recursions)

    print(lim)

    return lim
