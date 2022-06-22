import numpy
import matplotlib.pyplot as plt
import at
from at.collective import BeamLoadingElement, add_beamloading, remove_beamloading
from at.collective import get_params_beamloading, get_qs_beamloading
from at.collective import long_resonator_wf
from at.physics import harmonic_analysis
from at.constants import clight, qe


ring = at.load_lattice('../../../machine_data/esrf.m')
ring.radiation_on(cavity_pass='RFCavityPass')
ring.set_rf_frequency()
ring.set_rf_voltage(5.5e6)

o6,_ = ring.find_orbit6()
u0 = ring.get_energy_loss()
qs0 = at.get_tune(ring)[2]
rffreq = ring.get_rf_frequency()
timelag = ring.get_rf_timelag()
volt = ring.get_rf_voltage()
circ = ring.circumference
_, phis = at.get_timelag_fromU0(ring)
phis *= 2*numpy.pi*rffreq/clight

nturns = 1000
nturnsw = 50
nslice = 1
npart = 1

beta = 2.8
qfactor = 37500/(1+beta)
rshunt = 145*qfactor*11
current = 0.2

fring, fringr = at.fast_ring(ring)
cavpts = at.get_refpts(fringr, at.RFCavity)
bl_elem = add_beamloading(fringr, cavpts[0], qfactor, rshunt, Nslice=nslice, Nturns=nturnsw)
bl_elem.Current = current

qs_anal = get_qs_beamloading(qs0, 2*current*rshunt, volt, phis, bl_elem.Psi)
vgen0, fres0, psi0, vcav0 = get_params_beamloading(rffreq, current, volt, qfactor, rshunt, 0.0, phis)

dp_dat = numpy.zeros(nturns)
dz_dat = numpy.zeros(nturns)

part = numpy.zeros((6,npart))
part = (part.T + o6).T
#part = at.beam(npart,at.sigma_matrix(ring), o6)


for t in range(nturns):                                   
    _ = at.lattice_pass(fringr, part, nturns=1)        
    mns = numpy.mean(part,axis=1)
    dp_dat[t] = mns[4]
    dz_dat[t] = mns[5]

qs_anal_ss = get_qs_beamloading(qs0, bl_elem.Vbeam, volt, phis, bl_elem.Psi)

try:
    qs_dp = harmonic_analysis.get_tunes_harmonic(dp_dat,'laskar', num_harmonics=20, fmin=1e-5, fmax=0.1)
    qs_dz = harmonic_analysis.get_tunes_harmonic(dz_dat,'laskar', num_harmonics=20, fmin=1e-5, fmax=0.1)
except ValueError:
    qs_dp = 0
    qs_dz = 0  
    
    
print('Vbeam*cos(psi):',2*current*rshunt*numpy.cos(psi0), bl_elem.Vbeam*numpy.cos(bl_elem.Psi)) 
print('Psi:', psi0, bl_elem.Psi)
print('ResFreq', fres0, bl_elem.ResFrequency)   
print('Vcav:', vcav0, bl_elem.Vcav)
print('Vgen:', vgen0, bl_elem.Vgen) 
print('Qs:', qs_anal, qs_anal_ss, qs_dp, qs_dz)








    


