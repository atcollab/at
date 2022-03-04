import numpy
import matplotlib.pyplot as plt
import at
from at.collective.beam_loading import BeamLoadingElement
from at.collective.wake_object import build_srange
from at.physics import harmonic_analysis

#ring = at.load_m('../../../machine_data/esrf.m')
ring = at.load_mat('/machfs/carver/Repo/EBSlattices/AT/S28F.mat', key='LOW_EMIT_RING_INJ')

ring.radiation_on()
ring.set_rf_voltage(6e6)
ring.set_cavity_phase()
_,fring = at.fast_ring(ring)

cavpts = at.get_refpts(fring, at.RFCavity)
wturns=250
srange = build_srange(-0.36, 0.36, 1e-5, 1e-2, fring.circumference, fring.circumference*wturns)

beta = 2.8
qfactor = 37500/(1+beta)
rshunt = 145*qfactor*11

bl_elem = BeamLoadingElement('bl', fring, srange, qfactor, rshunt, cavpts=cavpts)
fring.append(bl_elem)
print('U0=',bl_elem.u0)

#Npart = 100
#part = at.beam(Npart, at.sigma_matrix(ring))
part = numpy.array([[0.,0.,0.,0.,0.,0.]]).T
Nturns = 10000
allCurrents = numpy.arange(0,201,50)
allQs = numpy.zeros(len(allCurrents))

for p, Ib in enumerate(allCurrents):
    bl_elem.Current = -Ib*1e-3
    for cav in cavpts:
        print(fring[cav].Voltage, fring[cav].PhaseLag)
    
    tbt_dp = numpy.zeros(Nturns)
    for i in numpy.arange(Nturns):
        part = at.lattice_pass(fring, part)[:,:,0,0]
        tbt_dp[i] = numpy.mean(part[4,:])

    qs = harmonic_analysis.get_tunes_harmonic(tbt_dp,'laskar', num_harmonics=20, fmin=1e-4, fmax=0.1)
    allQs[p] = qs

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(allCurrents, allQs, marker='x', linestyle='dashed', color='k')
ax1.set_xlabel('Ib [mA]')
ax1.set_ylabel('Qs')
plt.show()






    


