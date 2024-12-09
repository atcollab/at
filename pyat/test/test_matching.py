"""Test the mew matching."""

from __future__ import annotations

import numpy as np
import pytest

# import at
from at.lattice import Lattice, End
from at.latticetools import LocalOpticsObservable, GlobalOpticsObservable
from at.latticetools import EmittanceObservable, ObservableList
from at.future import RefptsVariable, VariableList, match


@pytest.fixture()
def mring(hmba_lattice: Lattice):
    """Return hmba lattice with focusing sextupoles split."""
    ring = hmba_lattice.deepcopy()
    sf = ring.get_uint32_index("SF*")
    sf1 = ring[sf[0]].divide([0.5, 0.5])
    sf2 = ring[sf[1]].divide([0.5, 0.5])
    ring.pop(sf[1])
    ring.insert(sf[1], sf2[1])
    ring.insert(sf[1], sf2[0])
    ring.pop(sf[0])
    ring.insert(sf[0], sf1[1])
    ring.insert(sf[0], sf1[0])
    ring.periodicity = 1
    return ring


@pytest.fixture()
def mline(hmba_lattice: Lattice):
    """Return hmba lattice with focusing sextupoles split."""
    ring = hmba_lattice.deepcopy()
    twiss_in, _, _ = ring.linopt6()
    sf = ring.get_uint32_index("SF*")
    sf1 = ring[sf[0]].divide([0.5, 0.5])
    ring.pop(sf[0])
    ring.insert(sf[0], sf1[0])
    return ring[:sf[0]+1], twiss_in


def test_linopt_matching(mring: Lattice):
    """Test linear optics matching."""

    # noinspection PyUnusedLocal
    def phase_advance(elemdata):
        # Evaluation function for the phase advance (both planes)
        mu = elemdata.mu
        return (mu[-1] - mu[0]) / 2.0 / np.pi

    # Define the location of constraint
    sf = mring.get_uint32_index("SF*")[1::2]
    center = mring.get_uint32_index("CellCenter")

    # Define the variables
    names = ["QF1*", "QD2*", "QD3*", "QF4*", "QD5*", "QF6*", "QF8*"]
    bounds = [[0, 5], [-5, 0], [-5, 0], [0, 5], [-5, 0], [0, 5], [0, 5]]
    variables = VariableList(
        RefptsVariable(nm, "PolynomB", index=1, bounds=bnd, name=nm, ring=mring)
        for nm, bnd in zip(names, bounds)
    )

    # Define the observables
    obs = ObservableList()
    # Tunes
    obs.append(GlobalOpticsObservable("tune", target=[0.38, 0.85]))
    # Start of the ring
    obs.append(LocalOpticsObservable(0, "alpha", target=[0.0, 0.0]))
    obs.append(LocalOpticsObservable(0, "dispersion", plane=1, target=0.0))
    # Center of the cell
    obs.append(LocalOpticsObservable(center, "beta", plane=1, target=4.69999257))
    obs.append(LocalOpticsObservable(center, "alpha", target=[0.0, 0.0]))
    obs.append(LocalOpticsObservable(center, "dispersion", plane=1, target=0.0))
    # Focusing sextupoles
    obs.append(LocalOpticsObservable(sf, "beta", plane=1, target=[5.4, 5.4]))
    obs.append(
        LocalOpticsObservable(sf, "alpha", plane=1, target=[0.68000392, -0.67999686])
    )
    obs.append(
        LocalOpticsObservable(
            sf, "dispersion", plane=0, target=[0.08820002, 0.08819999]
        )
    )
    # Phase advance
    obs.append(
        LocalOpticsObservable(sf, phase_advance, target=[1.49721338, 0.48228011])
    )

    # Perform the fit
    newring = match(mring, variables, obs, copy=True)

    # check the residuals
    obs.evaluate(newring)
    assert obs.sum_residuals < 1.0e-08


def test_envelope_matching(mring: Lattice):
    """Test linear and envelope."""
    # Define the variables: QF1; QF2, cavity voltage
    names = ["QF1*", "QD2*"]
    variables = VariableList(
        RefptsVariable(nm, "PolynomB", index=1, name=nm, ring=mring) for nm in names
    )
    variables.append(RefptsVariable(0, "Voltage", name="voltage"))

    # Define the constraints: fit the tunes from "tunes6"
    obs = ObservableList([EmittanceObservable("tunes6", target=[0.38, 0.85, 1.1e-4])])

    # Perform the fit
    ring = mring.enable_6d(copy=True)
    newring = match(ring, variables, obs, copy=True)

    # check the residuals
    obs.evaluate(newring)
    assert obs.sum_residuals < 1.0e-12

    # Define the constraints: fit transv. tunes form "tune" and long. tune from tunes6
    obs = ObservableList()
    obs.append(GlobalOpticsObservable("tune", plane=slice(2), target=[0.38, 0.85]))
    obs.append(EmittanceObservable("tunes6", plane="z", target=1.1e-04))

    # Perform the fit
    ring = mring.enable_6d(copy=True)
    newring = match(ring, variables, obs, copy=True)

    # check the residuals
    obs.evaluate(newring)
    assert obs.sum_residuals < 1.0e-12


def test_line_matching(mline):
    line, twiss_in = mline
    names = ["QF1*", "QD2*"]
    variables = VariableList(
        RefptsVariable(nm, "PolynomB", index=1, name=nm, ring=line) for nm in names
    )

    # Define the constraints: set alpha = [0.0, 0.0] at the end of the line
    obs = ObservableList(twiss_in=twiss_in)
    obs.append(LocalOpticsObservable(End, "alpha", target=[0.0, 0.0]))

    # Perform the fit
    newline = match(line, variables, obs, copy=True)

    # check the residuals
    obs.evaluate(newline)
    assert obs.sum_residuals < 1.0e-08
