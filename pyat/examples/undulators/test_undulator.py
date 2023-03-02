import at
# import at.plot
import matplotlib.pyplot as plt
from at.lattice.undulator_dipoles import undulator_dipoles, insert_undulator

lattice_file_for_test = '/machfs/liuzzo/EBS/beamdyn/matlab/optics/sr/theory/betamodel.mat'
lattice_variable_name = 'betamodel'

ring = at.load_mat(lattice_file_for_test, mat_key=lattice_variable_name)


und = undulator_dipoles(0.150*10, 10,
              BendAngle=None,
              B0=1.6497127666295313,
              Energy=6e9,
              magnetmodel='rectangularbend',
              PoleDistance=0.075/2)  # from G.Le Bec (approximately)

und = at.Lattice(und, energy=6e9)


def split_element(ring, idx, insert_elem=None):
    frac = 0.5
    enew = ring[idx].divide([frac, 1 - frac])
    ring = ring[:idx] + enew + ring[idx + 1:]
    if insert_elem is not None:
        ring.insert(idx + 1, insert_elem)
    return ring

insertIDindex = ring.get_refpts('ID17')[0]
ring = split_element(ring, insertIDindex+1, insert_elem=at.Marker('ID17downstream'))
insertIDindex = ring.get_refpts('ID17downstream')[0]

ring_w150, ring_w150_unmatched = insert_undulator(ring, und, insertIDindex)

ring_w150.plot_beta()

plt.show()

ring_w150.radiation_on(dipole_pass='BndMPoleSymplectic4RadPass')
emit0, bbb, eee = ring_w150.ohmi_envelope()
params = ring_w150.radiation_parameters()

# env = self.ring.envelope_parameters()
# env.emittances[0]

# get global lattice parameters
emittance_h_u = emit0['emitXY'][0]
emittance_v_u = emit0['emitXY'][1]
U0_u = params.U0

ring.radiation_on(dipole_pass='BndMPoleSymplectic4RadPass')
emit0, bbb, eee = ring.ohmi_envelope()
params = ring.radiation_parameters()

# env = self.ring.envelope_parameters()
# env.emittances[0]

# get global lattice parameters
emittance_h_0 = emit0['emitXY'][0]
emittance_v_0 = emit0['emitXY'][1]
U0_0 = params.U0

print(f'emittance Hor w/o undulator: {emittance_h_0}')
print(f'emittance Hor w/  undulator: {emittance_h_u}')


print(f'U0 w/o undulator: {U0_0}')
print(f'U0 w/  undulator: {U0_u}')
print(f'Delta U0 = {U0_u-U0_0} eV')

if __name__=='__main__':
    pass