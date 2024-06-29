import numpy as np
from numpy.testing import assert_allclose as assert_close

import at
from at import (
    RingObservable,
    ObservableList,
    OrbitObservable,
    GlobalOpticsObservable,
    LocalOpticsObservable,
    MatrixObservable,
    TrajectoryObservable,
    EmittanceObservable,
    LatticeObservable,
    GeometryObservable,
)


def test_observables(hmba_lattice):
    # noinspection PyUnusedLocal
    def phase_advance(elemdata):
        mu = elemdata.mu
        return mu[-1] - mu[0]

    def circumference(ring):
        return ring.cell_length

    ring = hmba_lattice.enable_6d(copy=True)
    ring.set_rf_frequency(dp=0.0)

    # Create an empty ObervableList
    allobs = ObservableList()
    # Populate with all kinds of Observables
    allobs.append(OrbitObservable(at.Monitor, axis="x"))
    allobs.append(
        LocalOpticsObservable(
            at.Monitor, "beta", plane=1, target=7.0, bounds=(-np.inf, 0.0)
        )
    )
    allobs.append(MatrixObservable("BPM_02"))
    allobs.append(LocalOpticsObservable(at.Monitor, "beta", plane="v", statfun=np.amax))
    allobs.append(
        LocalOpticsObservable(
            at.Quadrupole, "closed_orbit", plane=slice(4), target=0.0, weight=1.0e-6
        )
    )
    allobs.append(LocalOpticsObservable(at.Quadrupole, "s_pos"))
    allobs.append(
        LocalOpticsObservable([33, 101], phase_advance, all_points=True, summary=True)
    )
    allobs.append(GlobalOpticsObservable("tune", plane=0, use_integer=True))
    allobs.append(LocalOpticsObservable(at.End, "mu"))
    allobs.append(GlobalOpticsObservable("chromaticity"))
    allobs.append(LatticeObservable(at.Sextupole, "H", statfun=np.mean))
    allobs.append(LatticeObservable(at.Sextupole, "PolynomB", index=2))
    allobs.append(EmittanceObservable("emittances", plane="x"))
    allobs.append(RingObservable(circumference))
    allobs.append(TrajectoryObservable(at.Monitor, axis="px"))
    allobs.append(GeometryObservable(at.Monitor, "x"))

    # Evaluate the Observables
    r_in = np.zeros(6)
    r_in[0] = 0.001
    r_in[2] = 0.001
    allobs.evaluate(ring, r_in=r_in, initial=True)

    # Get the expected values
    o0, o = ring.find_orbit(refpts=at.All)
    el0, rg, el = ring.get_optics(refpts=at.All, orbit=o0, get_chrom=True)
    m66, ms = ring.find_m66(refpts="BPM_02", orbit=o0)
    prms = ring.envelope_parameters(orbit=o0)
    rout, _, _ = ring.track(r_in, refpts=at.Monitor)
    geodata, _ = ring.get_geometry(refpts=at.Monitor)

    monitors = ring.get_bool_index(at.Monitor)
    quadrupoles = ring.get_bool_index(at.Quadrupole)

    # Compare the results
    assert_close(allobs.values[0], o[monitors, 0])
    assert_close(allobs.values[1], el.beta[monitors, 1])
    assert_close(allobs.values[2], ms)
    assert_close(allobs.values[3], np.amax(el.beta[monitors, 1]))
    assert_close(allobs.values[4], el.closed_orbit[quadrupoles, :4])
    assert_close(allobs.values[5], el.s_pos[quadrupoles])
    assert_close(allobs.values[6], el.mu[101] - el.mu[33])
    assert_close(allobs.values[7], el.mu[-1, 0] / 2.0 / np.pi)
    assert_close(allobs.values[8], [el.mu[-1]])
    assert_close(allobs.values[9], rg.chromaticity)
    assert_close(
        allobs.values[10], np.mean([elem.H for elem in ring.select(at.Sextupole)])
    )
    assert_close(
        allobs.values[11], [elem.PolynomB[2] for elem in ring.select(at.Sextupole)]
    )
    assert_close(allobs.values[12], prms.emittances[0])
    assert_close(allobs.values[13], ring.cell_length)
    assert_close(allobs.values[14], rout[1, 0, :, 0])
    assert_close(allobs.values[15], geodata.x)
