import numpy as np
import matplotlib.pyplot as plt
from at.constants import clight
import at


# In this example, we will add a passive cavity which has the same frequency
# as the main cavity. This is similar to having a ' parked' cavity, except that
# here we will have the voltage feedback on. This means the frequency will be
# modified to fulfil a voltage setpoint.

# set point for the passive cavity
volt_set = 1e5

parking_detune = 50e3  # initial detuning for the passive cavity

I0 = 400e-3  # Total Current [A]

Vc = 1e6  # main RF Voltage [V]

# MAXIV
ring_dict = {
    "circumference": 528,
    "harmonic_number": 176,
    "ac": 0.000306,
    "energy": 3e9,
    "sigma_e": 0.000769,
    "U0": 363.8e3,
}

# main cavity resonator parameters
beta_main = 1.9
QL = 3688
Rsl = 0.32e6
Nmain = 4


Nbunches = ring_dict["harmonic_number"]
Nslice = 1

nparts_per_bunch = 1
Nparts = int(Nbunches * nparts_per_bunch)

Nturns = 50000

sigma_matrix = at.sigma_matrix(
    betax=1,
    betay=1,
    alphax=0,
    alphay=0,
    emitx=100e-12,
    emity=10e-12,
    espread=ring_dict["sigma_e"],
    blength=20e-3,
)

t0 = ring_dict["circumference"] / clight
tauz = 25.194e-3 / t0  # convert damping time into turns

simple_ring = at.simple_ring(
    ring_dict["energy"],
    ring_dict["circumference"],
    ring_dict["harmonic_number"],
    0.1,
    0.1,
    Vc,
    ring_dict["ac"],
    U0=ring_dict["U0"],
    tauz=tauz,
    espread=ring_dict["sigma_e"],
)

simple_ring.set_fillpattern(Nbunches)
simple_ring.set_beam_current(I0)

# The following line defines the element 0 as the main cavity. This takes the rf
# frequency of this element as the nominal, and defines things like bunch_spos
# from it. Very very important, as if you define a cavity with a frequency lower
# than the main, it will take the lowest frequency as default!!
simple_ring.cavpts = [0]

# Setup main cavity
simple_ring.set_cavity_phase()

## setup parked cavity
psi_parked = np.arctan(
    2
    * QL
    * (1 - simple_ring.rf_frequency / (simple_ring.rf_frequency + parking_detune))
)
parked_cav = at.RFCavity(
    "parkyboi",
    0.0,
    0.0,
    simple_ring.rf_frequency,
    simple_ring.harmonic_number,
    simple_ring.energy,
)
simple_ring.insert(1, parked_cav)


# add beamloading
at.add_beamloading(
    simple_ring,
    QL,
    Nmain * Rsl,
    cavpts=[0],
    Nslice=Nslice,
    VoltGain=1.0,
    PhaseGain=0.001,
    cavitymode=at.CavityMode.ACTIVE,
    fbmode=at.FeedbackMode.ONETURN,
)

at.add_beamloading(
    simple_ring,
    QL,
    Rsl,
    detune=parking_detune,
    cavpts=[1],
    Nslice=Nslice,
    cavitymode=at.CavityMode.PASSIVE_SETVOLTAGE,
    PhaseGain=0.00,
    passive_voltage=volt_set,
    buffersize=1000,
    windowlength=1000,
    fbmode=at.FeedbackMode.WINDOW,
)

# at.add_beamloading(simple_ring, QL, Rsl, detune=parking_detune, cavpts=[1],
#                   Nslice=Nslice, cavitymode=at.CavityMode.PASSIVE)

# Here we initialise the parked cavity with NO feedback on, we track and allow
# the beam to stabilise, then we switch it on
turn_passive_fb_on = 5000

bmon = at.BeamMoments("bmon")
simple_ring.append(bmon)

all_z_means = np.zeros((Nbunches, Nturns))
all_dp_means = np.zeros((Nbunches, Nturns))
all_z_stds = np.zeros((Nbunches, Nturns))
all_dp_stds = np.zeros((Nbunches, Nturns))

all_Vg_main = np.zeros((4, Nturns))
all_Vb_main = np.zeros((2, Nturns))

all_Vb_park = np.zeros((2, Nturns))

all_resfreq_park = np.zeros((Nturns))

parts = at.beam(Nparts, sigma_matrix)

for iturn in np.arange(Nturns):
    if iturn % 100 == 0:
        print(iturn)

    _ = simple_ring.track(parts, in_place=True, refpts=None, nturns=1)
    if iturn == turn_passive_fb_on:
        simple_ring[1].PhaseGain = 0.001

    all_z_means[:, iturn] = bmon.means[5, :, 0]
    all_dp_means[:, iturn] = bmon.means[4, :, 0]
    all_z_stds[:, iturn] = bmon.stds[5, :, 0]
    all_dp_stds[:, iturn] = bmon.stds[4, :, 0]

    all_Vg_main[:, iturn] = simple_ring[0].Vgen
    all_Vb_main[:, iturn] = simple_ring[0].Vbeam

    all_Vb_park[:, iturn] = simple_ring[1].Vbeam
    all_resfreq_park[iturn] = simple_ring[1].ResFrequency

fig, (ax1, ax2) = plt.subplots(2, 1)
ax1.plot(all_Vb_park[0, :] / 1e3, "r", label="Parked Cavity Voltage")

ax1.axhline(volt_set / 1e3, color="k", linestyle="dashed")
ax2.plot((all_resfreq_park - simple_ring.rf_frequency) / 1e3, "r")
ax2.set_ylabel(r"$\Delta$(Res Freq) [kHz]")
ax1.set_title("Vset={:.1f} kV".format(volt_set / 1e3))
ax2.set_xlabel("Turn")
ax1.set_ylabel("Parked Cavity Volt [kV]")
ax1.legend()
plt.show()
