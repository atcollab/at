import numpy as np
import at
from at.constants import clight
import matplotlib.pyplot as plt
import pickle as pkl
from scipy.optimize import curve_fit


def get_growth_rates():
    '''
    Function to compute mode 0 TRW mode growth rates
    R. Nagaoka, K. Bane, Collective Effects in Diffraction
    limited storage rings, eq31

    '''

    omega0 = 2*np.pi*ring.revolution_frequency
    energy = ring.energy
    R = ring.circumference/(2*np.pi)

    gr = (betax*omega0*current*R) / (4*np.pi*energy*beff**3) \
        * ((2*clight*Z0*rho_material) / ((1-qx)*omega0))**0.5

    return gr


def sin_fit(x, amp, freq, phase, off):
    return amp * np.sin(freq * x + phase) + off


def getTimeVectors(tbtMatrix):

    print('Do SVD')
    u, s, v = np.linalg.svd(tbtMatrix)

    # time vector of strongest mode
    vector = np.array([v[0][k] for k in np.arange(v.shape[1])])

    # SVD does not give mode number, so we fit a sin to the spatial pattern
    fit, _ = curve_fit(sin_fit, np.arange(Nbunches), np.real(u[:, 0]),
                       p0=[np.amax(np.real(u[:, 0])),
                       -2*np.pi/Nbunches,
                       -np.pi/4, 0])
    mode = fit[1]/2/np.pi*Nbunches

    bunchPlot = np.arange(0, Nbunches, 0.01)

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(np.real(u[:, 0]),
             color='r', marker='x', linestyle='None', label='Spatial Pattern')
    ax1.plot(bunchPlot, sin_fit(bunchPlot, fit[0], fit[1], fit[2], fit[3]),
             color='b', linestyle='solid', label='Sinusoidal Fit')
    ax1.set_xlabel('Bunch')
    ax1.legend()
    plt.show()

    return vector, mode, fit


def fitTimeVectors(tv):

    width = 400
    step = 200

    numsteps = int((len(tv) - width)/step) + 1
    starts = np.arange(0, len(tv) - width + 1, step)
    ends = np.arange(width, len(tv) + 1, step)

    dat = np.real(np.abs(tv))

    newdat = np.array([np.amax(dat[starts[i]:ends[i]])
                       for i in np.arange(numsteps)])

    turns = np.arange(len(dat))
    fitVal = np.polyfit((starts+ends) / 2, np.log(newdat), 1)
    chi_squared = np.sum((
                  np.polyval(fitVal, (starts+ends)/2) - np.log(newdat))**2)

    if 1/fitVal[0] > 1e6:
        fitVal[0] = np.inf

    return fitVal, newdat, starts, ends, chi_squared


# First we define everything we need and get the theoretical growth rates
Z0 = 377.7  # Vacuum impedance [Ohm]

resDict = pkl.load(open('./TRW_output.pkl', 'rb'))
ring = resDict['ring']
Nbunches = resDict['M']
current = resDict['current']
beff = resDict['rvac']
rho_material = 1/resDict['sigmac']


ld0, bd, ld = ring.get_optics()
betax = ld0.beta[0]
alphax = ld0.alpha[0]

f0 = ring.revolution_frequency
alpha = at.get_mcf(ring.radiation_off(copy=True))
eta = at.get_slip_factor(ring.radiation_off(copy=True))
envel = at.envelope_parameters(ring.radiation_on(copy=True))
qx, qy, qs = ring.get_tune()
fs0 = qs*f0

gr_theory = get_growth_rates()

# Find the last turn in the data where there were no lost particles
tstart = 4000
try:
    tend = np.where(not resDict['lostturn'])[0][0]
except:
    tend = len(resDict['xp'])
tend = 7000

data_xp = resDict['xp'].T
data_x = resDict['x'].T

# Plot all data dp, just to see how it looks
bbs = np.arange(0, Nbunches)
for bb in bbs:
    plt.plot(data_x[bb, tstart:tend])
plt.show()

# Here we use xp and x and combine with betax into one complex number
# this allows us to distinguish from + and - modes
inv = np.array([[1/np.sqrt(betax), 0],
                [alphax/np.sqrt(betax), np.sqrt(betax)]])

data = np.zeros(data_xp[:, tstart:tend].shape) + 1j * 0
for i in np.arange(data_xp.shape[0]):
    xp = data_xp[i, tstart:tend]
    x = data_x[i, tstart:tend]
    datpos, datmom = np.matmul(inv, [x, xp])
    data[i, :] = datpos + 1j * datmom

# Do the SVD and get all the time vectors
tv, mode, _ = getTimeVectors(data)
mode = int(np.round(mode))

# Fit the time vector
gr, newdat, starts, ends, chis = fitTimeVectors(tv)

turns = np.arange(0, tv.shape[0], 50)

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(np.real(tv), color='r', linestyle='solid', marker='None')
ax1.plot(turns, np.exp(gr[1])*np.exp(turns*gr[0]),
         linestyle='solid', marker='None',
         label=r'$\tau_{{\mathrm{{sim}}}}={:.2f}\ \mathrm{{turns}}$'
               .format(1/gr[0]) + '\n' +
               r'$\tau_{{\mathrm{{analytic}}}}={:.2f}\ \mathrm{{turns}}$'
               .format(f0/gr_theory))
ax1.plot((starts + ends)/2, newdat, color='k', marker='x', linestyle='None')
ax1.legend()
ax1.set_title('Mode {:d}, mu={:.3f}'.format(mode, mode/Nbunches))
ax1.set_xlabel('Turn')
plt.show()
