import numpy
import matplotlib.pyplot as plt
import matplotlib.cm as cmap
import at
from at.collective import BeamLoadingElement, add_beamloading, BLMode
from mpi4py import MPI
from at.tracking.utils import get_bunches_std_mean
from at.constants import qe
from at.physics import harmonic_analysis


def analytical_qs(ring, I):

    E0 = ring.energy * qe
    alpha_c = at.get_mcf(ring.radiation_off(copy=True))
    omega_rf = 2 * numpy.pi * ring.get_rf_frequency()
    Vc = ring.get_rf_voltage()

    envel = at.envelope_parameters(ring.radiation_on(copy=True))
    U0 = envel.U0 * qe

    Q0 = bl_elem.Qfactor * (1 + beta)
    Rsh = bl_elem.Rshunt * (1 + beta)  # Difference in definition

    phi_s = numpy.pi - numpy.arcsin(U0 / (qe * Vc))  # sync. phase in radian
    psi = numpy.arctan(-2 * Rsh * I * numpy.cos(numpy.pi - phi_s) /
                         (Vc * (1 + beta)))  # tuning angle
                         
    omega_res = (omega_rf /
                 (1 - 2 * Rsh * I * (numpy.cos(numpy.pi - phi_s)) /
                  (2 * Q0 * Vc)))  # resonant freq.

    h = ring.harmonic_number  # harmonic number
    omega0 = omega_rf / h  # rev. freq.
    T0 = 2 * numpy.pi / omega0  # rev. time

    omega_s0 = numpy.sqrt(qe * Vc * omega_rf * alpha_c * numpy.cos(numpy.pi - phi_s) /
                          (E0 * T0))  # Synch. freq. for a single particle

    K = alpha_c * qe * I / (E0 * T0) * Rsh / (1 + beta)

    # Synch. freq. with beam loading
    tilde_omega_s = numpy.sqrt(omega_s0**2 +
                               K * omega_rf * numpy.sin(2 * psi) + 1j * 0)
    return tilde_omega_s / omega0, omega_s0/omega0


comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

print(size, rank)

ring = at.load_lattice('../../../machine_data/esrf.m')
ring.enable_6d(cavity_pass='RFCavityPass')
ring.set_rf_frequency()


# In fact, the following resonator parameters are for the EBS machine
# Not the old ESRF SR, however for benchmarking this doesn't matter
beta = 2.8
q0 = 37500  # unloaded
qfactor = q0 / (1 + beta)  # loaded
rshunt = 145 * qfactor * 13  # loaded
o6, _ = ring.find_orbit6()

ring.set_rf_voltage(8.e6)

Nbunches = 16
ring.set_fillpattern(Nbunches)

_, fring = at.fast_ring(ring)
fring.pop(-1)  # drop diffusion element

# Here we specify whether we want to use PHASOR or WAKE
# beam loading models.
mode = 'WAKE'
if mode == 'WAKE':
    blm = BLMode.WAKE
else:
    blm = BLMode.PHASOR

# Npart must be at least Nbunches per core
Npart = Nbunches

# Now we give the fring and convert it
# into a beam loaded cavity.
add_beamloading(fring, qfactor, rshunt, Nslice=1,
                Nturns=50, blmode=blm,
                VoltGain=0.1, PhaseGain=0.1)

bl_elem = fring[at.get_refpts(fring, BeamLoadingElement)[0]]

# Specify some simulation parameters
kickTurn = 500
Nturns = 2**14 + kickTurn
current = 150e-3
fring.beam_current = current

part = numpy.zeros((6, Npart))
part = (part.T + o6).T

if rank == 0:

    print('\n')
    print('Current:', current)

    z_all = numpy.zeros((Nturns, Nbunches))
    dp_all = numpy.zeros((Nturns, Nbunches))


# Here it should be considered that there are 2 ways to run
# simulations in pyat. Either with ring.track(part, nturns=Nturns)
# or with for i in np.arange(Nturns): ring.track(part, nturns=1).
# With the former, you should use the BeamMoments element to acquire
# the means and stds of each bunch turn by turn when using MPI.
# However, when you want to kick after a certain turn, you have 1
# of 2 choices. Either you use a BeamMoments element, and split the 
# tracking into 2 (before and after). You can then concatenate the
# the results. Or you can use the code below to gather all particles
# yourself with the MPICOMM after each turn.

for i in numpy.arange(Nturns):
    # Apply a kick to ensure the coherent tune is large
    if i == kickTurn:
        part[4, :] += 3e-3

    fring.track(part, nturns=1, refpts=None, in_place=True)

    # Gather particles over all cores (compatible with MPI on or off)
    allPartsg = comm.gather(part)
    if rank == 0:
        allPartsg = numpy.concatenate(allPartsg, axis=1)
        stds, means = get_bunches_std_mean(allPartsg, Nbunches)

        dp_all[i] = numpy.array([x[4] for x in means])
        z_all[i] = numpy.array([x[5] for x in means])
        if i % 1000 == 0:
            # Print the turn number and number of lost particles
            print(i, numpy.sum(numpy.isnan(allPartsg)))


# Now we compute the tune of each bunch, and compare with the analytical
if rank == 0:
    qscoh = numpy.zeros(Nbunches)
    for ib in numpy.arange(Nbunches):
        dp_dat = dp_all[kickTurn:, ib] - numpy.mean(dp_all[kickTurn:, ib])
        qs = harmonic_analysis.get_tunes_harmonic(dp_dat,
                                                  num_harmonics=20,
                                                  fmin=1e-5, fmax=0.1)
        qscoh[ib] = qs

    qs_mn, qs_std = numpy.array([numpy.mean(qscoh), numpy.std(qscoh)])

    qs_theory, qs_zerocurrent = analytical_qs(ring, current)

    print('Analytical:', numpy.real(qs_theory))
    print('Simulated:', qs_mn, 'pm', qs_std)

    freq = numpy.fft.rfftfreq(Nturns - kickTurn)
    for i in numpy.arange(Nbunches):
        fftdat = numpy.abs(numpy.fft.rfft(dp_all[kickTurn:, i]))
        if i == 0:
            plt.semilogy(freq, fftdat,
                         color=cmap.jet(float(i + 1) / Nbunches),
                         label='Bunch by Bunch FFT')
        else:
            plt.semilogy(freq, fftdat,
                         color=cmap.jet(float(i + 1) / Nbunches))

    plt.xlim([0, 8e-3])
    plt.xlabel('Tune')
    plt.ylabel('FFT Amp [A.U.]')
    plt.axvline(numpy.mean(qscoh), label='Coherent Tune',
                linestyle='dashed', color='k')
    plt.axvline(numpy.real(qs_theory), label='Analytical Beam Loaded Tune',
                linestyle='dashed', color='r')
    plt.axvline(qs_zerocurrent, label='Zero Current Tune',
                linestyle='dashed', color='b')
    plt.legend()
    plt.show()
