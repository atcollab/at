import numpy as np
from at.constants import clight, e_mass, qe
import matplotlib.pyplot as plt
import at
from at.collective import Wake, LongResonatorElement

# Set up the ring
nturns = 10000
Npart = 10000
Nbunches = 1
current = 0.1   # A

ring = at.load_m('../../../machine_data/esrf.m')
ring.radiation_off()

ring.set_rf_frequency()
freqres = ring.rf_frequency

ring.set_fillpattern(Nbunches)
ring.beam_current = current

fring, _ = at.fast_ring(ring)
# Switch on RF cavity but have all other radiation sources off
fring.radiation_on(cavity_pass='RFCavityPass',
                   dipole_pass=None, quadrupole_pass=None)

# Define the resonator parameters and current
wturns = 50
srange = Wake.build_srange(0., 0.3, 1.0e-5, 1.0e-2,
                           ring.circumference, ring.circumference * wturns)

detuneHz = -5e4
fr = freqres + detuneHz
qfactor = 4500
rshunt = 6e6
bucket_size = clight/freqres


welem = LongResonatorElement('wake', ring, srange, fr,
                             qfactor, rshunt, Nturns=wturns, Nslice=1)
fring.append(welem)

# Define beam monitor
bmon = at.BeamMoments('mon')
fring.append(bmon)

# Particle generation and tracking
sigm = at.sigma_matrix(ring.radiation_on(copy=True))
part = at.beam(Npart, sigm)

fring.track(part, nturns=nturns, refpts=None, in_place=True)

dp_all = bmon.means[4, 0, :]

# Fit the results to obtain the simulated growth rate
width = 1000
step = 600
xr = (2 * np.arange(int(nturns/step)) * step + width) / 2
amp1 = np.array([np.amax(dp_all[ii * step: ii * step + width])
                 for ii in np.arange(int(nturns / step))])
fit1 = np.polyfit(xr, np.log(amp1), 1)
plotrange = np.arange(nturns)


# Define the parameters we need for the analytical value
gamma = ring.gamma
t_rev = 1.0/ring.revolution_frequency
eta = ring.slip_factor
params = ring.radiation_parameters()
fs = params.f_s
intensity = current * t_rev / qe / Nbunches


# Compute the growth rate
Rp = rshunt / (1 + 1j * qfactor * ((freqres + fs) / fr - fr / (freqres + fs)))
Rm = rshunt / (1 + 1j * qfactor * ((freqres - fs) / fr - fr / (freqres - fs)))

numer = e_mass * 2 * gamma * t_rev**2 * 2 * np.pi * fs
denom = -qe * intensity * eta * 2 * np.pi * fr * np.real(Rp - Rm)
tau_RS = numer / denom
gr = tau_RS/t_rev


lab = (r'$\tau_{{fit}}={:.1f}\ \mathrm{{turns}}$' '\n'
       r'$\tau_{{Chao}}={:.1f}\ \mathrm{{turns}}$').format(1/fit1[0], gr)
# Plot all together
fig = plt.figure(figsize=(9, 5))
ax1 = fig.add_subplot(111)
ax1.plot(dp_all,
         label=r'$\Delta f = {:.1f}\ \mathrm{{kHz}}$'.format(detuneHz/1e3),
         color='b')

ax1.plot(plotrange, np.exp(fit1[1])*np.exp(fit1[0]*plotrange),
         color='k', linestyle='solid', marker='None',
         label=lab)
ax1.plot(xr, amp1, linestyle='None', marker='x', color='r')

ax1.legend(bbox_to_anchor=(1.01, 1))
ax1.set_xlabel('Turn')
ax1.set_ylabel('dp/p')
fig.subplots_adjust(right=0.7)
plt.show()
