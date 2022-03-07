import numpy
import matplotlib.pyplot as plt
import at
from at.collective import BeamLoadingElement, get_anal_qs
from at.collective import build_srange
from at.physics import harmonic_analysis


ring = at.load_lattice('../../../machine_data/esrf.m')
ring.radiation_off(cavity_pass='CavityPass')
ring.set_rf_frequency()

wturns = 50
beta = 2.8
qfactor = 37500/(1+beta) 
rshunt = 145*qfactor*11
phil = 0.0

fring, _ = at.fast_ring(ring)
qs0 = at.get_tune(ring)[2]
u0 = ring.get_energy_loss(method=at.ELossMethod.TRACKING)
frev = ring.revolution_frequency
volt = ring.rf_voltage
rffreq = ring.rf_frequency
harm = ring.harmonic_number 
phis = numpy.arcsin(u0/volt)


srange = build_srange(0, 0.3, 1e-5, 1e-2, fring.circumference, fring.circumference*wturns)
bl_elem = BeamLoadingElement('bl', fring, srange, qfactor, rshunt, Nturns=wturns, Nslice=1)
fring.append(bl_elem)

print('u0', u0)
print('RF Freq.', rffreq)
print('Harm.', harm)
print('Phi_l angle: ', phil/(2*numpy.pi)*360)
print('Synch. phase: ', phis)
print('Synch. freq.: ', qs0*frev)
print('Q', qfactor)
print('Rs', rshunt)

Npart = 500
part = at.beam(Npart, at.sigma_matrix(ring.radiation_on(copy=True)))
Nturns = 1000
allCurrents = numpy.linspace(0.0, 0.5, 11)
allQs = numpy.zeros(len(allCurrents))
allQsa = numpy.zeros(len(allCurrents))

for p, Ib in enumerate(allCurrents):
    bl_elem.Current = Ib
    qsa = get_anal_qs(qs0,Ib,volt,rshunt,phil,u0)
    print('\n')
    print('Current : ',Ib)
    print('Generator voltage: ', fring.rf_voltage)
    print('Resonator frequency: ', fring[-1].ResFrequency)
    print('Detuning angle: ', fring.get_elements(at.RFCavity)[0].PhaseLag/numpy.pi*180)
    print('Synch. tune: ', qsa)   
    part = numpy.squeeze(at.lattice_pass(fring, part, nturns=Nturns))
    tbt_dp = numpy.mean(part[4,:,:], axis=0)
    tbt_s = numpy.mean(part[5,:,:], axis=0)
    try:
        qs_dp = harmonic_analysis.get_tunes_harmonic(tbt_dp,'laskar', num_harmonics=20, fmin=1e-5, fmax=0.1)
        qs_s = harmonic_analysis.get_tunes_harmonic(tbt_s,'laskar', num_harmonics=20, fmin=1e-5, fmax=0.1)
    except ValueError:
        qs_dp = 0
        qs_s = 0
    allQs[p] = (qs_dp+qs_s)/2
    allQsa[p] = qsa
    part = part[:,:,-1]

fig = plt.figure()
ax1 = fig.add_subplot(211)
ax1.plot(allCurrents, allQs, marker='x', linestyle='dashed', color='k', label='Tracking')
ax1.plot(allCurrents, allQsa, marker='x', linestyle='dashed', color='r', label='Analytic')
ax1.set_xlabel('Ib [mA]')
ax1.set_ylabel('Qs')
ax1.legend()
ax2 = fig.add_subplot(212)
ax2.plot(tbt_dp, color='k', label='dp')
ax2.plot(tbt_s, color='r', label='ct')
ax2.set_xlabel('Turn #')
ax2.set_ylabel('dp/ ct [m]')
ax2.legend()
plt.show()







    


