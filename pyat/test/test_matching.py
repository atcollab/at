import at
from numpy.testing import assert_allclose as assert_close
import pytest
import numpy as np


@pytest.fixture(scope='session')
def test_ring(hmba_lattice):
    ring = hmba_lattice.deepcopy()
    ring.periodicity = 1
    sf = ring.get_uint32_index('SF*')
    ring[sf[1]:sf[1] + 1] = ring[sf[1]].divide([0.5, 0.5])
    ring[sf[0]:sf[0] + 1] = ring[sf[0]].divide([0.5, 0.5])
    return ring


def test_linopt_matching(test_ring):

    def varqp(refs):
        name = refs
#       bnds = [0., 5.] if refs[1] == 'F' else [-5., 0]
#       return at.latticetools.ElementVariable(ik, 'K', name=name, bounds=bnds)
        return at.latticetools.ElementVariable(refs, 'K', name=name)

    # Define the location of constraint
    sf = test_ring.get_uint32_index('SF*')[1::2]
    center = test_ring.get_uint32_index('CellCenter')
    # Define the variables
    names = ['QF1*', 'QD2*', 'QD3*', 'QF4*', 'QD5*', 'QF6*', 'QF8*']
    variables = at.VariableList([varqp(refs) for refs in names])

    # Define an evaluation function for the phase advance (both planes)
    # noinspection PyUnusedLocal
    def phase_advance(ring, elemdata):
        mu = elemdata.mu
        return ((mu[-1] - mu[0]) / 2.0 / np.pi) % 1.0

    # Define the constraints
    constraints = at.ObservableList(test_ring)
    constraints.append(at.GlobalOpticsObservable('tune', target=[0.38, 0.85]))
    # Start of the ring
    constraints.append(
        at.LocalOpticsObservable(0, 'alpha', target=[0., 0.]))
    constraints.append(
        at.LocalOpticsObservable(0, 'dispersion', plane='px', target=0.))
    # Center of the cell
    constraints.append(
        at.LocalOpticsObservable(center, 'alpha', target=[0., 0.]))
    constraints.append(
        at.LocalOpticsObservable(center, 'beta', plane='v', target=4.67))
    constraints.append(
        at.LocalOpticsObservable(center, 'dispersion', plane='px', target=0.))
    # Focusing sextupoles
    constraints.append(
        at.LocalOpticsObservable(sf, 'alpha', plane='v', target=[0.68, -0.68]))
    constraints.append(
        at.LocalOpticsObservable(sf, 'beta', plane='v', target=[5.4, 5.4]))
    constraints.append(
        at.LocalOpticsObservable(sf, 'dispersion', plane='x', target=0.0882))
    # Phase advance
    constraints.append(at.LocalOpticsObservable(sf, phase_advance,
                                                target=[0.49721338, 0.48228011],
                                                summary=True))
    # Perform the fit
    newring = at.match(test_ring, variables, constraints, method='lm')

    # check the residuals
    constraints.evaluate(newring)
    assert_close(constraints.sum_residuals, 0.0, rtol=0.0, atol=5e-8)


def test_envelope_matching(test_ring):

    def varqp(refs):
        name = refs
        return at.latticetools.ElementVariable(refs, 'K', name=name)

    # Define the variables
    test_ring = test_ring.radiation_on(copy=True)
    names = ['QF1*', 'QD2*']
    variables = at.VariableList([varqp(refs) for refs in names])
    variables.append(at.ElementVariable(0, 'Voltage', name='voltage'))

    # Define the constraints
    constraints = at.ObservableList(test_ring)
    constraints.append(
        at.EmittanceObservable('tunes6', target=[0.38, 0.85, 1.1e-4],
                               weight=[1., 1., 0.0001]))

    # Perform the fit
    newring = at.match(test_ring, variables, constraints)

    # check the residuals
    constraints.evaluate(newring)
    assert_close(constraints.sum_residuals, 0, rtol=0.0, atol=1.e-10)

    # Define the constraints
    constraints = at.ObservableList(test_ring)
    constraints.append(
        at.GlobalOpticsObservable('tune', plane=slice(2), target=[0.38, 0.85]))
    constraints.append(
        at.EmittanceObservable('f_s', target=1300.0, weight=1000.0))

    # Perform the fit
    newring = at.match(test_ring, variables, constraints)

    # check the residuals
    constraints.evaluate(newring)
    assert_close(constraints.sum_residuals, 0, rtol=0.0, atol=1.e-10)
