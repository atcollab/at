import at
from at.physics.dynamic_aperture import Acceptance6D, dynamic_aperture, off_energy_dynamic_aperture, momentum_acceptance
import time
import copy
import pytest


def test_dynamic_aperture():

    # load and display optics
    sr_lattice_file = '../test_matlab/hmba.mat'
    sr_lattice_variable = 'ring'
    sr_ring = at.load_mat(sr_lattice_file, mat_key=sr_lattice_variable)

    sr_ring.radiation_off()
    # sr_ring.plot_beta()
    # plt.savefig('./optics.png', dpi=600)

    # x-xp acceptance at index 120 with dpp -0.5% for 32 turns
    print('x-xp acceptance, using class directly')
    da = Acceptance6D(copy.deepcopy(sr_ring),
                      mode='x-xp', start_index=120, n_turns=2**5, dpp=-0.005)
    da.verbose = False
    da.compute_range()
    h, v = da.compute()
    da.plot(file_name_save='test')

    # on-energy DA
    t = time.time()
    print('On-energy DA')
    dynamic_aperture(sr_ring, n_turns=2**5, n_radii=14, file_name_save='on_en_da', num_recursions=5)
    elapsed = time.time() - t
    print('On-energy DA took {s:2.1f}s'.format(s=elapsed))

    # off-en DA
    t = time.time()
    print('Off-energy DA')
    off_energy_dynamic_aperture(sr_ring, n_turns=2**5,
                                file_name_save='off_en_da_inside',
                                inject_from_inside=True,
                                num_recursions=3)
    elapsed = time.time() - t
    print('Off-energy DA took {s:2.1f}s'.format(s=elapsed))

    # momentum acceptance
    t = time.time()
    print('Momentum acceptance')
    momentum_acceptance(sr_ring, n_turns=2**5,
                        file_name_save='mom_acc',
                        ref_pts=range(0, 100, 10), num_recursions=3)
    elapsed = time.time() - t
    print('Momentum acceptance took {s:2.1f}s'.format(s=elapsed))

    # plt.show()

    pass


def test_search(n_turns=100, dpp=0.0, start_index=0, num_recursions=5):

    sr_lattice_file = '../test_matlab/hmba.mat'
    sr_lattice_variable = 'ring'
    sr_ring = at.load_mat(sr_lattice_file, mat_key=sr_lattice_variable)

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
