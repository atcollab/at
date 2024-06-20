import at
import numpy
import pytest
import warnings
from numpy.testing import assert_allclose as assert_close
from at.collective import Wake, WakeElement, ResonatorElement
from at.collective import WakeComponent, ResWallElement
from at.collective import add_beamloading, remove_beamloading, BLMode
from at import lattice_track
from at import lattice_pass, internal_lpass


_issorted = lambda a: numpy.all(a[:-1] <= a[1:])


def test_fillpattern(hmba_lattice):
    ring = hmba_lattice.enable_6d(copy=True)
    nbunch = 16
    current = 0.2
    ring.set_fillpattern(nbunch)
    ring.beam_current = current
    assert ring.fillpattern.shape == (ring.harmonic_number, )
    assert numpy.sum(ring.fillpattern) == 1
    assert ring.bunch_list.shape == (nbunch, )
    assert ring.nbunch == nbunch
    assert ring.bunch_currents.shape == (nbunch, )
    assert numpy.sum(ring.bunch_currents) == current
    

def test_build_srange():
    bspacing = 10
    wlength = 20
    bsize = 0.1
    steps = 1.0e-3
    stepl = 1.0e-2
    srange = Wake.build_srange(-bsize, bsize, steps, stepl, bspacing, wlength)
    assert len(numpy.unique(srange)) == len(srange)
    assert _issorted(srange)
    assert_close(max(srange), wlength + bsize - bspacing - stepl)
    assert min(srange) == -bsize
    assert_close(min(numpy.diff(srange)) ,steps)
    assert_close(min(numpy.diff(srange[srange>bsize])) ,stepl)
    
                  
def test_wake_object():
    srange = [0, 1, 1]
    with pytest.raises(Exception):
        Wake(srange)
    srange = Wake.build_srange(-0.36, 0.36, 1.0e-5, 1.0e-2, 844, 844)
    with pytest.raises(Exception):
        Wake.resistive_wall(srange, Wake.WakeComponent.X, 1.0, 1.0e-3, 1.0)
    srange = Wake.build_srange(0.0, 0.36, 1.0e-5, 1.0e-2, 844, 844) 
    long_res =  Wake.long_resonator(srange, 1.0e9, 1.0, 1.0e3, 1.0)
    assert 1.0e-24 in long_res.srange
    assert long_res.Z.shape == long_res.srange.shape
    
    
def test_wake_element(hmba_lattice):
    ring = hmba_lattice.enable_6d(copy=True)
    srange = Wake.build_srange(0.0, 0.36, 1.0e-5, 1.0e-2, 844, 844) 
    long_res =  Wake.long_resonator(srange, 1.0e9, 1.0, 1.0e3, 1.0) 
    welem = WakeElement('WELEM', ring, long_res)
    assert welem.TurnHistory.shape == (ring.nbunch*welem.Nslice*welem.Nturns, 4)
    ring.append(welem)
    ring.set_fillpattern(16)
    assert welem.TurnHistory.shape == (ring.nbunch*welem.Nslice*welem.Nturns, 4)
    ring.disable_6d(at.Collective)
    assert welem.PassMethod == 'IdentityPass'
    
    
def test_resonator_element(hmba_lattice):
    ring = hmba_lattice.enable_6d(copy=True)
    srange = Wake.build_srange(0.0, 0.36, 1.0e-5, 1.0e-2, 844, 844)
    welem = ResonatorElement('WELEM', ring, srange, WakeComponent.Z, 1.0e9, 1.0, 1.0e3)  
    assert welem.ResFrequency == 1.0e9
    assert welem.Qfactor == 1
    assert welem.Rshunt == 1.0e3
    assert_close(welem.WakeZ[[0,1]], [3.14159265e+12, 6.28318531e+12], rtol = 1e-8)
    welem.ResFrequency += 100
    assert welem.ResFrequency == 1.0e9 + 100
    
    
def test_resistive_wall_element(hmba_lattice):
    ring = hmba_lattice.enable_6d(copy=True)
    srange = Wake.build_srange(0.0, 0.36, 1.0e-5, 1.0e-2, 844, 844)
    welem = ResWallElement('WELEM', ring, srange, WakeComponent.DX, 1.0, 1.0e-2, 1.0e6)  
    assert welem.RWLength == 1.0
    assert welem.Conductivity == 1.0e6
    assert_close(welem.WakeDX[[0,1]], [ 0.0000000e+00, -1.0449877e+24], rtol = 1e-8)
    welem.Conductivity += 100
    assert welem.Conductivity == 1.0e6 + 100
    
    
def test_beamloading(hmba_lattice):
    ring = hmba_lattice.enable_6d(copy=True)
    with pytest.raises(Exception):
        add_beamloading(ring, 44e3, 400, cavpts=range(len(ring)))
    add_beamloading(ring, 44e3, 400)
    cavs = ring.get_elements(at.RFCavity)  
    for cav in cavs:
        assert cav.PassMethod == 'BeamLoadingCavityPass'
        assert hasattr(cav, 'Vbeam') 
        assert hasattr(cav, 'Vgen')
        assert hasattr(cav, 'Vcav') 
    ring.disable_6d(at.RFCavity)
    for cav in cavs:
        assert cav.PassMethod == 'IdentityPass'  
    ring.enable_6d(at.RFCavity)
    for cav in cavs:
        assert cav.PassMethod == 'BeamLoadingCavityPass'          
    remove_beamloading(ring)
    cavs = ring.get_elements(at.RFCavity)  
    for cav in cavs:
        assert cav.PassMethod == 'RFCavityPass' 
    

@pytest.mark.parametrize('func', (lattice_track, lattice_pass))
def test_track_beamloading(hmba_lattice, func):
    ring = hmba_lattice.enable_6d(copy=True)
    rin0 = numpy.zeros(6)
    func(ring, rin0, refpts=[])
    add_beamloading(ring, 44e3, 400, blmode=BLMode.WAKE)
    rin1 = numpy.zeros(6)
    func(ring, rin1, refpts=[])
    assert_close(rin0, rin1, atol=1e-21)
    ring.set_fillpattern(2)
    ring.beam_current = 0.2
    rin = numpy.zeros((6, 1))
    with pytest.raises(Exception):
        func(ring, rin, refpts=[])
    rin = numpy.zeros((6, 2))
    if func == lattice_track:
        func(ring, rin, refpts=[], in_place=True)
    else:
        func(ring, rin, refpts=[])
    assert_close(rin[:, 0], numpy.array([-2.318948e-08, -1.599715e-09,
                                        0.000000e+00,  0.000000e+00,
                                        -1.313306e-05, -1.443748e-08]
                                        ), atol=5e-10)
 

def test_buffers(hmba_lattice):
    ring = hmba_lattice.enable_6d(copy=True)
    ring.periodicity = 1
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ring.harmonic_number = 32
    nturns = 11
    nbunch = 4
    nslice = 51
    ns = nbunch*nslice
    ls = ns*ring.circumference/ring.periodicity
    add_beamloading(ring, 44e3, 400, Nturns=nturns, Nslice=nslice,
                    buffersize=nturns, blmode=BLMode.WAKE)
    ring.set_fillpattern(nbunch)
    ring.beam_current = 0.2
    rin = numpy.zeros((6, nbunch)) + 1.0e-6
    bl_elem = ring.get_elements('*_BL')[0]
    th = numpy.zeros((nturns, ) + bl_elem.TurnHistory.shape)
    vbh = numpy.zeros((nturns, ) + bl_elem.Vbeam_buffer.shape)
    vgh = numpy.zeros((nturns, ) + bl_elem.Vgen_buffer.shape)
    vbbh = numpy.zeros((nturns, ) + bl_elem.Vbunch_buffer.shape)
    for i in numpy.arange(nturns):
        ring.track(rin, nturns=1, refpts=None, in_place=True)
        th[i] = bl_elem.TurnHistory
        vbh[i] = bl_elem.Vbeam_buffer
        vgh[i] = bl_elem.Vgen_buffer
        vbbh[i] = bl_elem.Vbunch_buffer
    for i in numpy.arange(1, nturns):
        dth = numpy.sum(th[i-1, (nturns-i)*ns:]
                        - th[i, (nturns-i-1)*ns:(nturns-1)*ns]) - i*ls
        dvbh = numpy.sum(vbh[i-1, :, (nturns-i):]
                         - vbh[i, :, (nturns-i-1):(nturns-1)])
        dvgh = numpy.sum(vgh[i-1, :, (nturns-i):]
                         - vgh[i, :, (nturns-i-1):(nturns-1)])
        dvbbh = numpy.sum(vbbh[i-1, :, (nturns-i):]
                          - vbbh[i, :, (nturns-i-1):(nturns-1)])
    assert_close([dth, dvbh, dvgh, dvbbh], numpy.zeros(4), atol=1e-9)
