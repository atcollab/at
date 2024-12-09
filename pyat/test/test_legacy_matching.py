"""Test the legacy matching."""

from __future__ import annotations

import numpy as np
from numpy.testing import assert_allclose as assert_close
import pytest

from at.lattice import Lattice
from at.physics import linopt2
from at.deprecated import LinoptConstraints, EnvelopeConstraints
from at.deprecated import ElementVariable, match


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


def test_linopt_matching(mring: Lattice):
    """Test linear optics matching."""
    # Define the location of constraint
    sf = mring.get_uint32_index("SF*")[1::2]
    center = mring.get_uint32_index("CellCenter")

    # Define the variables
    names = ["QF1*", "QD2*", "QD3*", "QF4*", "QD5*", "QF6*", "QF8*"]
    bounds = [[0, 5], [-5, 0], [-5, 0], [0, 5], [-5, 0], [0, 5], [0, 5]]
    variables = [
        ElementVariable(mring.get_uint32_index(nm), "PolynomB", index=1, name=nm)
        for nm, bnd in zip(names, bounds)
    ]

    # Define an evaluation function for the phase advance (both planes)
    # noinspection PyUnusedLocal
    def mu_diff(lindata, tune, chrom):
        delta_mu = (lindata[1].mu - lindata[0].mu) / (2 * np.pi)
        return delta_mu % 1.0

    # Define the constraints
    cst1 = LinoptConstraints(mring, method=linopt2)
    cst1.add("tunes", [0.38, 0.85], name="tunes")
    # Start of the ring
    cst1.add("alpha", [0.0, 0.0], refpts=0, name="alpha_0")
    cst1.add("dispersion", 0.0, refpts=0, index=1, name="eta'_0")
    # Center of the cell
    cst1.add("beta", 4.69999257, refpts=center, index=1, name="betay_c")
    cst1.add("alpha", [0.0, 0.0], refpts=center, name="alpha_c")
    cst1.add("dispersion", 0.0, refpts=center, index=1, name="eta'_c")
    # Focusing sextupoles
    cst1.add("beta", [5.40000677, 5.39998124], refpts=sf, index=1, name="betay_sf")
    cst1.add("alpha", [0.68000392, -0.67999686], refpts=sf, index=1, name="alphay_sf")
    cst1.add("dispersion", [0.08820002, 0.08819999], refpts=sf, index=0, name="eta_sf")
    # Phase advance
    cst1.add(mu_diff, [0.49721338, 0.48228011], refpts=sf, name="delta phi")

    # Perform the fit
    newring = match(mring, variables, (cst1,), method="lm", verbose=1)

    # check the residuals
    assert_close(cst1.evaluate(newring), 0, rtol=0.0, atol=1e-4)


def test_envelope_matching(mring: Lattice):
    """Test linear and envelope."""
    # Define the variables
    ring = mring.radiation_on(copy=True)
    names = ["QF1*", "QD2*"]
    variables = [
        ElementVariable(ring.get_uint32_index(nm), "PolynomB", index=1, name=nm)
        for nm in names
    ]
    variables += [ElementVariable([0], "Voltage", name="voltage")]
    # [0] instead of 0 because of bug in utils.py

    # Define the constraints
    lopcst = EnvelopeConstraints(ring)
    lopcst.add("tunes", [0.38, 0.85, 1.1e-4])

    # Perform the fit
    newring = match(ring, variables, (lopcst,), verbose=1)

    # check the residuals
    residual = lopcst.evaluate(newring.radiation_on(copy=True))
    assert_close(residual, 0, rtol=0.0, atol=6e-9)

    # Define the constraints
    lincst = LinoptConstraints(ring)
    lincst.add("tunes", [0.38, 0.85], name="tunes")
    lopcst = EnvelopeConstraints(ring)
    lopcst.add("tunes", 1.1e-4, index=2, name="sync. tune")

    # Perform the fit
    newring = match(ring, variables, (lincst, lopcst), verbose=1)

    # check the residuals
    linresidual = lincst.evaluate(newring)
    lopresidual = lopcst.evaluate(newring.radiation_on(copy=True))
    assert_close(linresidual, 0, rtol=0.0, atol=6e-9)
    assert_close(lopresidual, 0, rtol=0.0, atol=3e-8)
