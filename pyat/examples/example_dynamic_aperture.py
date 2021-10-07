import at
import at.physics.dynamic_aperture as da
import at.plot.dynamic_aperture as daplot
import at.lattice.cavity_access
import at.plot
import time
import numpy as np

folder_data = '.'  # to store files and figures

sr_lattice_file = '../test_matlab/hmba.mat'
sr_lattice_variable = 'RING'
N = 2**2 # number of turns for tracking

sr_arc = at.load_mat(sr_lattice_file, mat_key=sr_lattice_variable)

sr_ring = at.Lattice(sr_arc*32)

sr_ring.periodicity =1
sr_ring.set_rf_harmonic_number(992)
sr_ring.set_rf_voltage(6.0e6)
sr_ring.radiation_on()

Nex = 7 # number of examples to run
i = 1

# x - xp Acceptance
t = time.time()
print('{}/{}  Computing: x-xp acceptance, using class directly ...'.format(i,Nex))
i +=1
DA = da.Acceptance6D(sr_ring,
                  mode='x-xp', start_index=120, n_turns=N, dp=-0.0)
DA.verbose = False
DA.compute_range()
# define grid of test particles in the range
DA.test_points_grid()
# test if particles of the grid survive or not and produce a curve arround the survived ones
_, _, h, v, s = DA.compute()
if DA.verbose:
    [print(h_, v_, s_) for h_, v_, s_ in zip(h, v, s)]
elapsed = time.time() - t
daplot.plot(DA, file_name_save=folder_data + '/Acceptance')
print('x-xp Acceptance took {s:2.1f}s'.format(s=elapsed))

# Dynamic aperture on-energy
t = time.time()
print('{}/{}  Computing: On-energy DA ...'.format(i,Nex))
i +=1

h0, v0, da_, search = da.dynamic_aperture(
                    sr_ring,
                    n_turns=N,
                    n_radii=15,
                    file_name_save=folder_data + '/on_en_da',
                    num_recursions=3)
daplot.plot_dynamic_aperture(h0, v0, da_, search, file_name_save=folder_data + '/on_en_da')

elapsed = time.time() - t
print('On-energy DA took {s:2.1f}s'.format(s=elapsed))

print('{}/{}  Computing: On-energy DA grid ...'.format(i,Nex))
i +=1
h0, v0, da_, search = da.dynamic_aperture(
                    sr_ring,
                    n_turns=N,
                    n_radii=21, # here = N Horizontal grid points
                    n_theta=21, # here = N Vertical grid points
                    grid_mode='cartesian', # search input will be ignored
                    parallel=False,  # can be set to True on Unix to use patpass, see next example
                    file_name_save=folder_data + '/on_en_da',
                    num_recursions=5)
daplot.plot_dynamic_aperture(h0, v0, da_, search, file_name_save=folder_data + '/on_en_da_grid')

elapsed = time.time() - t
print('On-energy DA grid took {s:2.1f}s'.format(s=elapsed))


"""  FOLLOWING WORKS ONLY IN UNIX.
print('{}/{}  Computing: On-energy DA grid (parallel computation)  ...'.format(i,Nex))
i +=1
h0, v0, da_, _ = da.dynamic_aperture(sr_ring,
                    parallel=True,
                    n_turns=N,
                    n_radii=21,
                    n_theta=21,
                    grid_mode='grid',  # grid mode also allows parallel computation
                    file_name_save=folder_data + '/on_en_da_grid',
                    num_recursions=5)
# daplot.plot_dynamic_aperture(h0, v0, da_, search, file_name_save=folder_data + '/on_en_da_grid')

elapsed = time.time() - t
print('On-energy DA grid (parallel) took {s:2.1f}s'.format(s=elapsed))

"""


# max hor. position vs momenutm deviation (off-energy DA)
t = time.time()
print('{}/{}  Computing: Off-energy DA ...'.format(i,Nex))
i +=1
maxneg0, dp, da_ = da.off_energy_dynamic_aperture(sr_ring,
                               n_turns=N,
                               deltaps=np.linspace(-0.1, 0.1, 21),
                               file_name_save=folder_data + '/off_en_da',
                               inject_from_inside=True,
                               num_recursions=5)

daplot.plot_off_energy_dynamic_aperture(maxneg0, dp, da_, file_name_save=folder_data + '/off_en_da')

elapsed = time.time() - t
print('Off-energy DA took {s:2.1f}s'.format(s=elapsed))


# max hor. position vs momenutm deviation (off-energy DA)
t = time.time()
print('{}/{}  Computing: Off-energy DA at a different location ...'.format(i,Nex))
i +=1
maxneg0, dp, da_ = da.off_energy_dynamic_aperture(sr_ring,
                               n_turns=N,
                               start_index=56,
                               deltaps=np.linspace(-0.1, 0.1, 21),
                               file_name_save=folder_data + '/off_en_da',
                               inject_from_inside=True,
                               num_recursions=5)

daplot.plot_off_energy_dynamic_aperture(maxneg0, dp, da_, file_name_save=folder_data + '/off_en_da_somewhere')

elapsed = time.time() - t
print('Off-energy DA took {s:2.1f}s'.format(s=elapsed))



#  Local Momentum Acceptance
t = time.time()
print('{}/{} Computing: Momentum acceptance ...'.format(i,Nex))
i +=1
mom, s, da_ = da.momentum_acceptance(sr_ring, n_turns=N,
                       file_name_save=folder_data + '/mom_acc',
                       ref_pts=range(0, 125, 5), num_recursions=2)
daplot.plot_momentum_acceptance(mom, s, da_, file_name_save=folder_data + '/mom_acc')

elapsed = time.time() - t
print('Momentum acceptance took {s:2.1f}s'.format(s=elapsed))



# generic multidimentsional Acceptance
print('{}/{} Computing: 6D acceptance, using class directly ...'.format(i,Nex))
i +=1
DA = da.Acceptance6D(sr_ring,
                  mode='6D', start_index=1, n_turns=N, dp=-0.0)
DA.n_points['ct'] = 3 # do not scan ct coordinate
DA.n_points['delta'] = 1 # do not scan ct coordinate
DA.n_points['x'] = 9 # do not scan ct coordinate
DA.n_points['y'] = 7 # do not scan ct coordinate
DA.n_points['xp'] = 11 # do not scan ct coordinate
DA.n_points['yp'] = 3 # do not scan ct coordinate
DA.grid_mode='grid'
DA.verbose = False
DA.parallel_computation = False
# define range in each dimentsion
DA.compute_range()
# define grid of test particles in the range
DA.test_points_grid()
# test if particles of the grid survive or not and produce a curve arround the survived ones
h, v, h_all, v_all, survived_all = DA.compute()
[print(h_, v_) for h_, v_ in zip(h, v)]
print(DA.survived)
elapsed = time.time() - t
daplot.plot6d(DA,file_name_save=folder_data + '/Acceptance_5D')
print('6D Acceptance took {s:2.1f}s'.format(s=elapsed))

