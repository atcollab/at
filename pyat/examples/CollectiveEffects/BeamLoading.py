import numpy
import matplotlib.pyplot as plt
import at
from at.collective import BeamLoadingElement, add_beamloading
from at.collective import get_qs_beamloading, get_params_beamloading, BLMethod
from at.physics import harmonic_analysis
from at.constants import clight, qe

ring = at.load_lattice('../../../machine_data/esrf.m')
ring.radiation_on(cavity_pass='RFCavityPass')
ring.set_rf_frequency()


# In fact, the following resonator parameters are for the EBS machine
# Not the old ESRF SR, however for benchmarking this doesn't matter
beta = 2.8
qfactor = 37500/(1+beta)
rshunt = 145*qfactor*13
phil = 0.0
o6, _ = ring.find_orbit6()


# Generate the fast ring, but remove the quantum diffusion element in
# order to get an easy and clean Qs measurement
_, fring = at.fast_ring(ring)
fring.pop(-1)
cavpts = at.get_refpts(fring, at.RFCavity)


# Covnert the cavity element into a beam loading element
Nslice = 1
Npart = 1
bl_elem = add_beamloading(fring, cavpts[0], qfactor, rshunt,
                          Nslice=Nslice, Nturns=50, mode=BLMethod.AUTO)


qs0 = at.get_tune(ring)[2]
u0 = ring.get_energy_loss(method=at.ELossMethod.TRACKING)
frev = ring.revolution_frequency
volt = ring.rf_voltage
rffreq = ring.rf_frequency
harm = ring.harmonic_number
phis = numpy.pi - numpy.arcsin(u0/volt)


print('u0', u0)
print('RF Freq.', rffreq)
print('Harm.', harm)
print('Phi_l angle: ', phil/(2*numpy.pi)*360)
print('Synch. phase: ', phis)
print('Synch. freq.: ', qs0*frev)
print('Q', qfactor)
print('Rs', rshunt)

# Number of turns for the tracking
Nturns = 2**14

allCurrents = numpy.linspace(0., 0.30, 31)
allQs = numpy.zeros(len(allCurrents))

for p, Ib in enumerate(allCurrents):

    # Initialise Particles
    part = numpy.zeros((6, Npart))
    part = (part.T + o6).T

    # Set the current within the beam loading element
    # This wil automatically set the initial values assuming
    # Vb of 2*Ib*Rshunt
    bl_elem.Current = Ib
    psi0 = bl_elem.Psi
    vgen0 = bl_elem.Vgen
    fres0 = bl_elem.ResFrequency
    Vcav0 = bl_elem.Vcav

    print('\n')
    print('Current : ', Ib)
    allPart = numpy.zeros((6, Npart, Nturns))
    for i in numpy.arange(Nturns):
        _ = at.lattice_pass(fring, part, nturns=1)
        allPart[:, :, i] = part

    # Print the initial assumptions and the final converged values
    print('Psi:', psi0, bl_elem.Psi)
    print('Amplitude:', 2*Ib*rshunt*numpy.cos(psi0),
          bl_elem.Vbeam*numpy.cos(bl_elem.Psi))
    print('Phase:', phis,
          2*numpy.pi*numpy.mean(part, axis=1)[5]*rffreq/clight + phis)
    print('Vcav:', Vcav0, bl_elem.Vcav)
    print('Vgen:', vgen0, bl_elem.Vgen)
    print('Vb:', 2*Ib*rshunt*numpy.cos(psi0)**2,
          bl_elem.Vbeam*numpy.cos(bl_elem.Psi)**2)

    # Perform a harmonic analysis and save the average of
    # the synchrotron tune coming from dp and s
    tbt_dp = numpy.mean(allPart[4, :, :], axis=0)
    tbt_s = numpy.mean(allPart[5, :, :], axis=0)
    try:
        qs_dp = harmonic_analysis.get_tunes_harmonic(
                    tbt_dp, 'laskar', num_harmonics=20, fmin=1e-5, fmax=0.1)
        qs_s = harmonic_analysis.get_tunes_harmonic(
                    tbt_s, 'laskar', num_harmonics=20, fmin=1e-5, fmax=0.1)
    except ValueError:
        qs_dp = 0
        qs_s = 0
    allQs[p] = (qs_dp + qs_dp)/2
    print('Qs:', allQs[p])

# Now we compute the analytical values as a function of current
# to make the comparison
analCurrents = numpy.arange(0, numpy.amax(allCurrents)+0.001, 0.001)
analVbeam = 2*analCurrents*rshunt
analVgen, analFres, analPsi, analVcav = get_params_beamloading(
        rffreq, analCurrents, volt, qfactor, rshunt, phil, phis, vb=None)
analQs = get_qs_beamloading(qs0, analVbeam, volt, numpy.pi - phis, -analPsi)

# The two lines start to diverge at high current due to the fact
# that we are lumping all the current into one single bunch.
# This will be improved when simulating with multiple bunches
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(analCurrents*1e3, analQs, 'r', label='Analytical')
ax1.plot(allCurrents*1e3, allQs,
         color='b', marker='x', label='PyAT', linestyle='None')
ax1.set_ylabel('Qs')
ax1.set_xlabel('I [mA]')
ax1.legend()
plt.show()
