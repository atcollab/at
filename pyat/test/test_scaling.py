import at
from at.lattice.elements import Corrector, Dipole, Quadrupole, ThinMultipole
from math import pi
from numpy.testing import assert_allclose as assert_close
import pytest


def simple_fodo(ncells: int = 64) -> at.Lattice:
    """Return a simple F0D0 lattice"""
    # dipole
    dipole_length = 1.0
    bending_angle = pi / ncells
    # Focusing quad
    kf = 0.5
    lf = 1.0
    # Defocusing quad
    kd = -0.5
    ld = 1.0
    # Drifts
    drift_f = 0.5
    drift_d = 0.5 + 0.5 * ld

    drf = at.Drift('drf', drift_f)
    drd = at.Drift('drd', drift_d)
    qf0 = Quadrupole('qf0', lf, kf)
    qd0 = ThinMultipole('qd0', [0, 0], [0, kd])
    dip0 = Dipole('dip0', dipole_length, bending_angle, 0)
    kick0 = Corrector('kick0', 0.0, [0.001, -0.001])
    # build the nominal lattice
    fodo0 = at.Lattice([drd, kick0, dip0, drf, qf0, drf, dip0, drd, qd0],
                       energy=6.0E9, name='fodo0')
    return fodo0


@pytest.mark.parametrize('dp', [0.01, -0.01])
@pytest.mark.parametrize('ncells', [64, 6])
def test_scaling(dp, ncells):
    """Check the FieldScaling attribute"""
    scaling = 1.0 + dp
    fodo0 = simple_fodo(ncells)
    magnets = fodo0.get_bool_index(
        at.checktype((Corrector, Dipole, Quadrupole, ThinMultipole)))

    # Nominal optics:
    _, rd0, el0 = at.get_optics(fodo0, refpts=at.All)

    # Off-momentum optics with scaled magnets
    fodoa = fodo0.replace(magnets)
    for elem in fodoa[magnets]:
        elem.FieldScaling = scaling
    _, rda, ela = at.get_optics(fodoa, dp=dp, refpts=at.All)
    assert_close(rd0.tune, rda.tune, atol=1.e-10, rtol=0)
    assert_close(el0.closed_orbit[:, [0, 2]],
                 ela.closed_orbit[:, [0, 2]], atol=1.e-10, rtol=0)

    # Off-momentum optics with scaling the whole lattice
    scin = at.Element('scin', PassMethod='ChangePRefPass',
                      FieldScaling=scaling)
    scout = at.Element('scout', PassMethod='ChangePRefPass',
                       FieldScaling=1.0 / scaling)
    fodob = fodo0.copy()
    fodob.insert(0, scin)
    fodob.append(scout)
    _, rdb, elb = at.get_optics(fodob, dp=dp, refpts=at.All)
    assert_close(rd0.tune, rdb.tune, atol=1.e-10, rtol=0)
    assert_close(el0.closed_orbit[:,    [0, 2]],
                 elb.closed_orbit[1:-1, [0, 2]], atol=1.e-10, rtol=0)
