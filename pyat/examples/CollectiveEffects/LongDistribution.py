import numpy
import at
import matplotlib.pyplot as plt
from at.collective import Wake, LongResonatorElement
from at.collective.haissinski import Haissinski

# First we define the ring, the BB resonator, the current and the wake element
ring = at.load_m('../../../machine_data/esrf.m')



freq = 10e9
qfactor = 1
Rs = 1e4
current = 5e-4
m = 60  # 30 is quite coarse, 70 or 80 is very fine. 50 is middle
kmax = 8


srange = Wake.build_srange(0., 0.36, 1.0e-5, 1.0e-2,
                           ring.circumference, ring.circumference)

# Now we initialise the Haissinski class, and solve
# then we normalise the distribution and shift the charge center to be at 0
wobj = Wake.long_resonator(srange, freq, qfactor, Rs, ring.beta)
ha = Haissinski(wobj, ring, m=m, kmax=kmax,
                current=current, numIters=30, eps=1e-13)
ha.solve()

ha_x_tmp = ha.q_array*ha.sigma_l
ha_prof = ha.res/ha.Ic
ha_prof /= numpy.trapz(ha_prof, x=ha_x_tmp)
ha_cc = numpy.average(ha_x_tmp, weights=ha_prof)
ha_x = (ha_x_tmp - ha_cc)

'''
currents = numpy.arange(0, 1.1e-3, 2e-4)
ha.solve_steps(currents)

fig = plt.figure()
ax1 = fig.add_subplot(111)
for i in numpy.arange(len(currents)):
    ax1.plot(1e3*ha.q_array*ha.sigma_l, ha.res_steps[i,:]/ha.I_steps[i],
             label='Ib={:f}mA'.format(currents[i]*1e3),
             color=cm.jet(float(i)/len(currents)))
ax1.legend()
ax1.set_xlabel('z [mm]')
ax1.set_ylabel(r'$\rho(z)$')
plt.show()
'''


# Now we set up and run the tracking.
# The final distribution is an average of the last numAve turns
welem = LongResonatorElement('wake', ring,
                             srange, freq, qfactor, Rs, Nslice=300)


Nbunches = 1
ring.beam_current = current
ring.set_fillpattern(Nbunches)

ring.enable_6d()
ring.set_cavity_phase()
_, fring = at.fast_ring(ring)
fring.append(welem)


sigm = at.sigma_matrix(ring.radiation_on(copy=True))
Nparts = 10000
Nturns = 20000
nbins = 60
part = at.beam(Nparts, sigm)

histAveFlag = False
numAve = 5000
id0 = 0
for t in numpy.arange(Nturns):
    if t % 1000 == 0:
        print('Tracking turn ', t, ' of ', Nturns)
    fring.track(part, in_place=True, refpts=None)

    if t > Nturns-numAve:

        if not histAveFlag:
            histAveFlag = True

            nt1, et1 = numpy.histogram(part[5, :], bins=nbins)

            etc1 = numpy.array([(et1[i] + et1[i + 1]) / 2
                                for i in numpy.arange(len(et1) - 1)])

            allData = numpy.zeros((len(etc1), numAve + 1))
            allData[:, 0] = etc1
            allData[:, id0 + 1] = nt1

        else:
            n1, e1 = numpy.histogram(part[5, :], bins=et1)
            ec1 = numpy.array([(e1[i] + e1[i + 1]) / 2
                               for i in numpy.arange(len(e1) - 1)])
            allData[:, id0 + 1] = n1
        id0 += 1


# Identical post processing of tracking distribution (to the Haissinski)
zr_tmp = allData[:, 0]
prof = numpy.mean(allData[:, 1:], axis=1)
prof /= numpy.trapz(prof, x=zr_tmp)

cc = numpy.average(zr_tmp, weights=prof)
zr = zr_tmp - cc

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(1e3*ha_x, ha_prof,
         color='r', linestyle='solid', label='Haissinski Solution')
ax1.plot(1e3*zr, prof, color='k', linestyle='dashed', label='Tracking')
ax1.set_xlabel('z [mm]')
ax1.set_ylabel(r'$\rho(z)$')
ax1.legend()
plt.savefig('./haissinski_dist.png')
plt.show()
