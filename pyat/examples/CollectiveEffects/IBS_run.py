"""
Example python script to test the IBS effect
"""

import at
import numpy as np
import matplotlib.pyplot as plt
from at.constants import qe, clight
from at import Wake, WakeType, WakeComponent, WakeElement
import pickle


class IBS_run:

    def __init__(self, bunch_current_mA, emity_pm):

        self.bunch_current = float(bunch_current_mA) * 1e-3  # Bunch current [A]
        self.emity = float(emity_pm) * 1e-12  # Vertical emittance [m]
        self.emitx = 141e-12
        self.model = "CIMP"  # model to compute the growth rates.
        self.get_opt = "markers"  # method to compute the lattice optics.
        self.update_turns = (
            500  # number of turns after which the grwoth rates for IBS are computed.
        )
        self.nturns = 10000  # total number of turns for tracking

        self.ring, self.mcf = self._define_ring()

    def _define_ring(self):
        """
        This method is used to load the lattice and set it by defining the number of bunches , bunch current and including
        radiation damping and quantum diffusion element.
        """

        ring_hmba = at.load_mat('../../../machine_data/hmba.mat')  # single cell for the EBS lattice
        ring = ring_hmba.repeat(32)  # cell repeated 32 times to get the complete ring
        ring.disable_6d()
        mcf = at.get_mcf(ring)
        ring.enable_6d()
        ring.set_rf_voltage(5.5e6)  # voltage of the RF cavity (5.5 MV for EBS ring)
        ring.set_fillpattern(1)
        ring.set_beam_current(self.bunch_current)
        ring.set_cavity_phase()
        qdelem = at.gen_quantdiff_elem(ring)  # adds the quantum diffusion element
        ring.append(qdelem)
        return ring, mcf

    def _build_simple_ring(self, ibs_elem, params, optics, chromaticity):
        """
        This method is to define a simple ring (with an RF cavity element ,a 6x6 linear transfer map with no radiation damping,
        a detuning and chromaticity element, a simple radiation damping element, and a simplified quantum diffusion element
         which contains equilibrium emittance) for fast tracking.

        Further, an ibs element is also added to model the intrabeam scattering effect.

        """
        sring = at.simple_ring(
            energy=self.ring.energy,
            circumference=self.ring.circumference,
            harmonic_number=self.ring.harmonic_number,
            Qx=params.tunes6[0],
            Qy=params.tunes6[1],
            Vrf=self.ring.rf_voltage,
            alpha=self.mcf,
            betax=optics[0].beta[0],
            betay=optics[0].beta[1],
            alphax=optics[0].alpha[0],
            alphay=optics[0].alpha[1],
            dispx=optics[0].dispersion[0],
            dispxp=optics[0].dispersion[1],
            dispy=optics[0].dispersion[2],
            dispyp=optics[0].dispersion[3],
            Qpx=chromaticity[0],
            Qpy=chromaticity[1],
            emitx=params.emittances[0],
            emity=self.emity,
            espread=params.sigma_e,
            taux=params.Tau[0] * self.ring.revolution_frequency,
            tauy=params.Tau[1] * self.ring.revolution_frequency,
            tauz=params.Tau[2] * self.ring.revolution_frequency,
            U0=self.ring.energy_loss,
            name="esrf",
            particle="electron",
            TimeLag=0.0,
        )

        sring.enable_6d()
        sring.set_fillpattern(1)
        sring.set_beam_current(self.bunch_current)
        sring.set_cavity_phase()

        ring_ibs = sring.deepcopy()
        ring_ibs.append(ibs_elem)
        return ring_ibs

    def _generate_bunch(self, params, optics):
        """
        Generate a particle bunch using the sigma matrix based on given optics and beam parameters.

        """
        sigma_matrix = at.sigma_matrix(
            betax=optics[0].beta[0],
            betay=optics[0].beta[1],
            alphax=optics[0].alpha[0],
            alphay=optics[0].alpha[1],
            emitx=self.emitx,
            emity=self.emity,
            blength=params.sigma_l,
            espread=params.sigma_e,
        )
        return at.beam(
            1000, sigma_matrix
        )  # 400000 is the number of particles in the bunch.

    def _track(self, ring_ibs, bunch):
        """
        Track the bunch through the ring including radiation damping, quantum diffusion,
        and intrabeam scattering effects.

        Calculates bunch length, energy spread, and emittances every turn and store the values in lists.

        """

        emit_ibs = []
        blength_ibs = []
        espread_ibs = []

        for _ in range(self.nturns):
            ring_ibs.track(bunch, refpts=None, in_place=True)
            sigma = at.sigma_matrix(beam=bunch)
            eigs, _ = np.linalg.eig(sigma @ at.jmat(3))
            emits = np.sort(np.abs(eigs)[::2])
            emit_ibs.append([emits[1], emits[0], emits[2]])
            blength_ibs.append(np.std(bunch[5, :]))
            espread_ibs.append(np.std(bunch[4, :]))

        return blength_ibs, espread_ibs, emit_ibs

    def run_simulation(self):
        """
        This method defines the IBS element with the choice of model ('CIMP', 'Bane', 'PM', 'PS'),
        the method to calculate optics ('mbtrack_eqv', 'markers', 'average'), and update_turns
        (i.e., how many turns between IBS growth rate computations) etc. These parameters are defined in the __init__ method.

        It then calls the previously defined functions for ring and bunch definition, and performs tracking.

        """

        ibs_elem = at.IBSElement(
            "Ibs",
            0,
            self.ring,
            model=self.model,
            get_opt=self.get_opt,
            update_turns=self.update_turns,
            PassMethod="pyIBSRadPass",
        )
        ibs_elem.compute_optics_params()
        params = self.ring.envelope_parameters()
        _, rd, _ = at.get_optics(self.ring, get_chrom=True)
        chromaticity = rd.chromaticity
        optics = at.get_optics(self.ring.disable_6d(copy=True))

        ring_ibs = self._build_simple_ring(ibs_elem, params, optics, chromaticity)
        bunch = self._generate_bunch(params, optics)

        return self._track(ring_ibs, bunch)


############################################################################################################################

## Run simulation with given bunch current and horizontal emittance values.
## Adjust these parameters to study how IBS effects vary.
##
## Outputs bunch length, energy spread, and emittance, then saves results to a pickle file.


def launch():
    bunch_current = 10.0  # mA
    emity = 1.0  # pm

    sim = IBS_run(bunch_current, emity)
    blength, espread, emit = sim.run_simulation()

    outdict = {
        "bunch_current_mA": bunch_current,
        "emity_pm": emity,
        "bunch_length": blength,
        "energy_spread": espread,
        "emittances": emit,
    }

    with open("./IBS_run_output.pkl", "wb") as f:
        pickle.dump(outdict, f)

    print("Simulation complete. Output saved to IBS_run_output.pkl")


if __name__ == "__main__":
    launch()
