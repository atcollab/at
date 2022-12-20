import pytest
from at import checktype, Monitor
import numpy as np
from numpy.testing import assert_allclose as assert_close


@pytest.mark.parametrize('offset', ((0.001, 0.0), ([0.001, -0.001], 0.0)))
@pytest.mark.parametrize('gain', (([1.1, 0.9], 0.0), (1.05, 0.0)))
@pytest.mark.parametrize('tilt', ((-0.002, 0.0),))
def test_systematic_bpm_errors(hmba_lattice, offset, gain, tilt):
    ring = hmba_lattice.copy()
    bpms = ring.get_cells(checktype(Monitor))
    ring.assign_errors(bpms, BPMOffset=offset, BPMGain=gain, BPMTilt=tilt)
    bpmoff = np.vstack([el.BPMOffset for el in ring.select(bpms)])
    bpmgain = np.vstack([el.BPMGain for el in ring.select(bpms)])
    bpmtilt = np.array([el.BPMTilt for el in ring.select(bpms)])

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
    xyorbit = ring.find_orbit_err(bpms, all=False)[1][:, [0, 2]]

    bpmoff = np.vstack([el.BPMOffset for el in ring.select(bpms)])
    bpmgain = np.vstack([el.BPMGain for el in ring.select(bpms)])
    bpmtilt = np.array([el.BPMTilt for el in ring.select(bpms)])

    # Build all rotation matrices
    rotm = np.stack([_rotmat(v).T for v in bpmtilt], axis=0)
    # reshape offset,
    of = np.reshape(xyorbit0+bpmoff, (-1, 1, 2))
    expected = np.squeeze(of @ rotm) * bpmgain
    assert_close(xyorbit, expected, rtol=0, atol=0)
