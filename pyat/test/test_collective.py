import at
import numpy
import pytest
from numpy.testing import assert_allclose as assert_close
from numpy.testing import assert_equal
from at.collective import Wake, WakeElement, ResonatorElement
from at.collective import WakeComponent, ResWallElement
from at.collective import add_beamloading, remove_beamloading, BLMode


_issorted = lambda a: numpy.all(a[:-1] <= a[1:])


def test_fillpattern(hmba_lattice):
    ring = hmba_lattice.radiation_on(copy=True)
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
    ring = hmba_lattice.radiation_on(copy=True)
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
    ring = hmba_lattice.radiation_on(copy=True)
    srange = Wake.build_srange(0.0, 0.36, 1.0e-5, 1.0e-2, 844, 844)
    welem = ResonatorElement('WELEM', ring, srange, WakeComponent.Z, 1.0e9, 1.0, 1.0e3)  
    assert welem.ResFrequency == 1.0e9
    assert welem.Qfactor == 1
    assert welem.Rshunt == 1.0e3
    print(welem.WakeZ[[0,1]])
    assert_close(welem.WakeZ[[0,1]], [3.14159265e+12, 6.28318531e+12], rtol = 1e-8)
    welem.ResFrequency += 100
    assert welem.ResFrequency == 1.0e9 + 100
    
    
def test_resistive_wall_element(hmba_lattice):
    ring = hmba_lattice.radiation_on(copy=True)
    srange = Wake.build_srange(0.0, 0.36, 1.0e-5, 1.0e-2, 844, 844)
    welem = ResWallElement('WELEM', ring, srange, WakeComponent.DX, 1.0, 1.0e-2, 1.0e6)  
    assert welem.RWLength == 1.0
    assert welem.Conductivity == 1.0e6
    print(welem.WakeDX[[0,1]])
    assert_close(welem.WakeDX[[0,1]], [ 0.0000000e+00, -1.0449877e+24], rtol = 1e-8)
    welem.Conductivity += 100
    assert welem.Conductivity == 1.0e6 + 100
    
    
def test_beamloading(hmba_lattice):
    ring = hmba_lattice.radiation_on(copy=True)
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
    
        
def test_track_beamloading(hmba_lattice):
    ring = hmba_lattice.radiation_on(copy=True)
    rin0 = numpy.zeros(6)
    at.lattice_pass(ring, rin0, refpts=[])
    add_beamloading(ring, 44e3, 400, mode=BLMode.WAKE)
    rin1 = numpy.zeros(6)
    at.lattice_pass(ring, rin1, refpts=[])    
    assert_close(rin0, rin1, atol=1e-21)
    ring.set_fillpattern(2)
    ring.beam_current = 0.2
    rin = numpy.zeros((6,1))
    with pytest.raises(Exception):
        at.lattice_pass(ring, rin, refpts=[])
    rin = numpy.zeros((6,2))
    at.lattice_pass(ring, rin, refpts=[])
    assert_close(rin[:,0], numpy.array([-2.318948e-08, -1.599715e-09,
                                        0.000000e+00,  0.000000e+00,
                                        -1.313306e-05, -1.443748e-08]
                                        ), atol=1e-10)
    
