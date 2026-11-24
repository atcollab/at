import numpy as np
import matplotlib.pyplot as plt
import at

CUR, NPART, NTURNS = 1.0, 2000, 2000

ring = at.load_mat("../../../machine_data/hmba.mat").repeat(32)
ring.disable_6d()
at.fit_tune(ring, "QF1*", "QD2*", [0.24, 0.29])
ring[ring.get_uint32_index("QF1*")[0]].PolynomA[1] = 0.05
ring.enable_6d()
ring.set_fillpattern(1)
ring.set_beam_current(CUR)
frev = ring.revolution_frequency

res = at.ibs_equilibrium(ring, history=True)
ep = ring.envelope_parameters()
emit = res["history"]
scale = np.sqrt(emit[:, 2] / emit[0, 2])
ode = np.column_stack((emit, ep.sigma_l * scale, ep.sigma_e * scale))
ode_turns = res["time"] * frev

fr = at.fast_ring(ring)
fr.insert(0, at.IBSElement("IBS", ring, update_turns=10))
fr.set_fillpattern(1)
fr.set_beam_current(CUR)
orbit0 = fr.find_orbit()[0]
beam = at.beam(NPART, at.sigma_matrix(ring), orbit=orbit0)


def stats(b, il=5, ie=4):
    sig = at.sigma_matrix(beam=b)
    emits, _ = at.emittances_from_beam(sigma_mat=sig)
    return np.concatenate((emits, np.sqrt([sig[il, il], sig[ie, ie]])))


rout, *_ = fr.track(beam, nturns=NTURNS)
turns = np.arange(NTURNS + 1)
hist = np.array([stats(beam)] + [stats(rout[:, :, 0, t]) for t in range(NTURNS)])


labels = [r"$\epsilon_1$ [pm]", r"$\epsilon_2$ [pm]", r"$\epsilon_3$ [$\mu$m]",
          r"$\sigma_l$ [mm]", r"$\sigma_\delta$ [$10^{-3}$]"]
units = [1e12, 1e12, 1e6, 1e3, 1e3]
fig, axes = plt.subplots(1, 5, figsize=(18, 3.8))
for k, (ax, label, u) in enumerate(zip(axes, labels, units)):
    ax.plot(turns, u * hist[:, k], label="AT tracking")
    ax.plot(ode_turns, u * ode[:, k], "k--", label="ibs_equilibrium history")
    ax.set_xlim(0, NTURNS)
    ax.set_xlabel("turn")
    ax.set_ylabel(label)
axes[0].legend()
fig.tight_layout()
plt.show()
