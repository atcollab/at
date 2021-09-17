import at
import at.physics.dynamic_aperture as da
import at.plot.dynamic_aperture as daplot
import time
import numpy as np

folder_data = '.'  # to store files and figures

sr_lattice_file = '../test_matlab/hmba.mat'
sr_lattice_variable = 'RING'
N = 2**5 # number of turns for tracking

sr_ring = at.load_mat(sr_lattice_file, mat_key=sr_lattice_variable)

sr_ring.radiation_on()
t = time.time()
print('On-energy DA')
h0, v0, da_, search = da.dynamic_aperture(sr_ring,
                    n_turns=N,
                    n_radii=14,
                    file_name_save=folder_data + '/on_en_da',
                    num_recursions=5)
daplot.plot_dynamic_aperture(h0, v0, da_, search, file_name_save=folder_data + '/on_en_da')

elapsed = time.time() - t
print('On-energy DA took {s:2.1f}s'.format(s=elapsed))

# off-en DA
t = time.time()
print('Off-energy DA')
maxneg0, dpp, da_ = da.off_energy_dynamic_aperture(sr_ring,
                               n_turns=N,
                               deltaps=np.linspace(-0.1, 0.1, 21),
                               file_name_save=folder_data + '/off_en_da',
                               inject_from_inside=True,
                               num_recursions=5)

daplot.plot_off_energy_dynamic_aperture(maxneg0, dpp, da_, file_name_save=folder_data + '/off_en_da')

elapsed = time.time() - t
print('Off-energy DA took {s:2.1f}s'.format(s=elapsed))

# MA
t = time.time()
print('Momentum acceptance')
mom, s, da_ = da.momentum_acceptance(sr_ring, n_turns=N,
                       file_name_save=folder_data + '/mom_acc',
                       ref_pts=range(0, 125, 5), num_recursions=2)
daplot.plot_momentum_acceptance(mom, s, da_, file_name_save=folder_data + '/mom_acc')

elapsed = time.time() - t
print('Momentum acceptance took {s:2.1f}s'.format(s=elapsed))
