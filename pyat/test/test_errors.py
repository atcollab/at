import pytest
from at import checktype, Monitor, Quadrupole
import numpy as np
from numpy.testing import assert_allclose as assert_close


@pytest.mark.parametrize('offset', ((0.001, 0.0), ([0.001, -0.001], 0.0)))
@pytest.mark.parametrize('gain', (([1.1, 0.9], 0.0), (1.05, 0.0)))
@pytest.mark.parametrize('tilt', ((-0.002, 0.0),))
def test_systematic_bpm_errors(hmba_lattice, offset, gain, tilt):
    ring = hmba_lattice.copy()
    bpms = ring.get_cells(checktype(Monitor))
    ring.assign_errors(bpms, BPMOffset=offset, BPMGain=gain, BPMTilt=tilt)
    ring = ring.enable_errors()
    bpmoff = np.vstack([el._BPMOffset for el in ring.select(bpms)])
    bpmgain = np.vstack([el._BPMGain for el in ring.select(bpms)])
    bpmtilt = np.array([el._BPMTilt for el in ring.select(bpms)])

    # Check that all values are correctly assigned
    assert np.all(bpmoff == offset[0])
    assert np.all(bpmgain == gain[0])
    assert np.all(bpmtilt == tilt[0])


def test_random_bpm_errors(hmba_lattice):

    def _rotmat(theta):
        cs = np.cos(theta)
        sn = np.sin(theta)
        return np.array([[cs, sn], [-sn, cs]])

    ring = hmba_lattice.copy()
    bpms = ring.get_cells(checktype(Monitor))

    offset = (0.001, [0.001, 0.002])
    gain = ([1.1, 0.9], 0.01)
    tilt = 0.001

    ring.assign_errors(bpms, BPMOffset=offset, BPMGain=gain, BPMTilt=tilt)

    xyorbit0 = ring.find_orbit(bpms)[1][:, [0, 2]]
    ring = ring.enable_errors()

    bpmoff = np.vstack([el._BPMOffset for el in ring.select(bpms)])
    bpmgain = np.vstack([el._BPMGain for el in ring.select(bpms)])
    bpmtilt = np.array([el._BPMTilt for el in ring.select(bpms)])

    # Build all rotation matrices
    rotm = np.stack([_rotmat(v).T for v in bpmtilt], axis=0)
    # reshape offset,
    of = np.reshape(xyorbit0 - bpmoff, (-1, 1, 2))
    expected = np.squeeze(of @ rotm) * bpmgain

    xyorbit = ring.find_orbit(bpms, monitor_errors=True)[1][:, [0, 2]]
    assert_close(xyorbit, expected, rtol=0, atol=0)

    _, _, el = ring.get_optics(bpms, monitor_errors=True)
    xyorbit = el.closed_orbit[:, [0, 2]]
    assert_close(xyorbit, expected, rtol=0, atol=0)


def test_field_errors(hmba_lattice):
    ring = hmba_lattice.copy()
    qpoles = ring.get_cells(checktype(Quadrupole))
    # Nominal quadrupoles strengths
    k = np.array([el.PolynomB[1] for el in ring.select(qpoles)])
    shifterr = 0.0001            # Random shift error both planes
    rotationerr = 0.0002         # Random tilt
    octu = 100.0
    PAErr = 0.0001               # random vertical dipole
    systPBErr = [0, 0, 0, octu]  # Systematic scaling octupole
    randPBErr = [0, 0.001]       # Random strength error
    # Assign errors on quadrupoles
    ring.assign_errors(qpoles, ShiftErr=shifterr, RotationErr=rotationerr,
                       ScalingPolynomAErr=PAErr,
                       ScalingPolynomBErr=(systPBErr, randPBErr))
    # Enable the errors
    errlat = hmba_lattice.enable_errors()
    # Get the PolynomB including errors
    pb = np.vstack([el.PolynomB[:4] for el in errlat.select(qpoles)])
    # Check the systematic octupole
    assert_close(pb[:, 3], k*octu, rtol=0, atol=0)
