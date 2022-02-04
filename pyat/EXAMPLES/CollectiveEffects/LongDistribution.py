import numpy
import at
from at.constants import clight, e_mass, qe
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from at.collective.wake_object import build_srange
from at.collective.wake_elements import LongResonatorElement
from at.collective.Haissinski import Haissinski
import time


# First we define the ring, the BB resonator, the current and the wake element
ring = at.load_m('../../../machine_data/esrf.m')

freq = 10e9
qfactor = 1
Rs = 1e4
current = 5e-4
m = 50 #30 is quite coarse, 70 or 80 is very fine. 50 is middle
kmax = 8

srange = build_srange(-0.36, 0.36, 1.0e-5, 1.0e-2, ring.circumference, ring.circumference)
welem = LongResonatorElement('wake', ring, srange, freq, qfactor, Rs, Nslice=300)
welem.Current = current

# Now we initialise the Haissinski class, and solve, then we normalise the distribution and shift the charge center to be at 0
ha = Haissinski(welem, ring, m=m, kmax=kmax, current=current, numIters = 30, eps=1e-13)
ha.solve()

ha_prof = ha.phi_1/ha.I
ha_prof /= numpy.trapz(ha.phi_1/ha.I, x=ha.q_array*ha.sigma_l)
ha_cc = numpy.average(ha.q_array*ha.sigma_l, weights=ha_prof)
ha_x = (ha.q_array*ha.sigma_l - ha_cc)[::-1] #This should be integrated into code but the sign reversal is needed


# Now we set up and run the tracking. The final distribution is an average of the last numAve turns 
ring.radiation_on()
ring.set_cavity_phase()
_,fring = at.fast_ring(ring)
fring.append(welem)


sigm = at.sigma_matrix(ring.radiation_on(copy=True))
Nparts = 40000
Nturns = 20000
nbins = 30
part = at.beam(Nparts, sigm)
part[:,2] += 10e-6

histAveFlag = False
numAve = 5000
id0 = 0
for t in numpy.arange(Nturns):
    if t%1000==0:
        print('Tracking turn ', t+1, ' of ', Nturns)
    part = at.lattice_pass(fring, part)[:,:,0,0]
    
    if t > Nturns-numAve:  
        if not histAveFlag:
            histAveFlag = True

            nt1, et1 = numpy.histogram(part[5,:], bins=nbins)

            etc1 = numpy.array([(et1[i] + et1[i+1])/2 for i in numpy.arange(len(et1)-1)]) 

            allData = numpy.zeros((len(etc1), numAve + 1))
            allData[:, 0] = etc1
            allData[:, id0 + 1] =  nt1

        else:
            n1, e1 = numpy.histogram(part[5,:], bins=et1)
            ec1 = numpy.array([(e1[i] + e1[i+1])/2 for i in numpy.arange(len(e1)-1)]) 
            allData[:, id0 + 1] =  n1
        id0 += 1   
        

# Identical post processing of tracking distribution (to the Haissinski)
zr = allData[:,0]
prof = numpy.mean(allData[:,1:],axis=1) 

cc = numpy.average(zr, weights=prof)
norm = numpy.trapz(prof,x=zr)
prof = nt1/norm

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(1e3*ha_x, ha_prof, color='r', linestyle='solid', label='Haissinski Solution')
ax1.plot(1e3*(zr+cc), prof, color='k', linestyle='dashed', label='Tracking')
ax1.set_xlabel('z [mm]')
ax1.set_ylabel(r'$\rho(z)$')
ax1.legend()
plt.show()



