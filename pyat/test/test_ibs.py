import numpy
import at
import pytest


def test_two_particles_for_one_turn(hmba_lattice):

    ring = hmba_lattice.enable_6d(copy=True)
    bunch_current = 10e-3
    ring.set_fillpattern(1)
    ring.set_beam_current(bunch_current)
    optics = at.get_optics(ring.disable_6d(copy=True), refpts=0)
    ibs_element = at.IBSElement(
        "Ibs_element",
        0,
        ring,
        model="CIMP",
        get_opt="markers",
        n_points=100,
        n_bin=1,
        PassMethod="pyIBSRadPass",
    )
    ibs_element.compute_optics_params()
    ring.append(ibs_element)

    numpy.random.seed(42)
    sigma_matrix = at.sigma_matrix(
        betax=optics[0].beta[0],
        betay=optics[0].beta[1],
        alphax=optics[0].alpha[0],
        alphay=optics[0].alpha[1],
        emitx=141e-12,
        emity=1e-12,
        blength=0.003,
        espread=0.001,
    )
    rin = at.beam(2, sigma_matrix)

    rout_particle1_expected = numpy.array(
        [
            -7.25596327e-07,
            -7.64355188e-06,
            1.14092197e-06,
            7.33024711e-07,
            3.94815955e-04,
            6.27494565e-04,
        ]
    )
    rout_particle2_expected = numpy.array(
        [
            -2.80118387e-05,
            -3.36321436e-06,
            -1.19125693e-06,
            -2.81718065e-07,
            6.31241787e-04,
            -5.73637434e-03,
        ]
    )

    ring.track(rin)

    numpy.testing.assert_almost_equal(rin[:, 0], rout_particle1_expected, decimal=6)
    numpy.testing.assert_almost_equal(rin[:, 1], rout_particle2_expected, decimal=6)
