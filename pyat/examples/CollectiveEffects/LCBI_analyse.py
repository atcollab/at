import numpy as np
import at
import matplotlib.pyplot as plt
import pickle as pkl
from scipy.optimize import curve_fit


def get_growth_rates():
    '''
    Function to compute all LCBI mode growth rates
    Largely stolen from mbtrack2 (Alexis Gamelin)
    '''

    def Zr(omega):
        df = (omega_res/omega - omega/omega_res)
        return np.real(rshunt / (1 + 1j * qfactor * df))

    factor = eta * current / (4 * np.pi * ring.energy * qs)

    growth_rates = np.zeros(Nbunches)
    for mu in np.arange(Nbunches):
        omega0 = 2 * np.pi * f0
        omega_res = 2 * np.pi * fres
        n_max = int(10 * omega_res / (omega0 * Nbunches))

        n0 = np.arange(n_max)
        n1 = np.arange(1, n_max)
        omega_p = omega0 * (n0 * Nbunches + mu + qs)
        omega_m = omega0 * (n1 * Nbunches - mu - qs)

        sum_val = np.sum(omega_p * Zr(omega_p)) - np.sum(omega_m * Zr(omega_m))
        growth_rates[mu] = factor * sum_val
    return growth_rates


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
                       2*np.pi*expected_mode/Nbunches,
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
resDict = pkl.load(open('./LCBI_output.pkl', 'rb'))
ring = resDict['ring']
Nbunches = resDict['M']
current = resDict['current']
fres = resDict['fres']
qfactor = resDict['qfactor']
rshunt = resDict['rshunt']


f0 = ring.revolution_frequency
alpha = at.get_mcf(ring.radiation_off(copy=True))
eta = at.get_slip_factor(ring.radiation_off(copy=True))
envel = at.envelope_parameters(ring.radiation_on(copy=True))
qx, qy, qs = ring.get_tune()
betaz = ring.circumference * alpha / (2 * np.pi * qs)
fs0 = qs*f0

gr_theory = get_growth_rates()
expected_mode = np.argmax(gr_theory)
gr_peak = gr_theory[expected_mode]


data_dp = resDict['dp'].T
data_z = resDict['z'].T

# Plot all data dp, just to see how it looks
bbs = np.arange(0, Nbunches)
for bb in bbs:
    plt.plot(data_dp[bb, :])
plt.xlabel('Turn')
plt.ylabel('dp/p')
plt.show()

# Here we use dp and z and combine with betaz into one complex number
# this allows us to distinguish from + and - modes
inv = np.array([[1/np.sqrt(betaz), 0],
                [0, np.sqrt(betaz)]])

data = np.zeros(data_dp[:, :].shape) + 1j * 0
for i in np.arange(data_dp.shape[0]):
    dp = data_dp[i, :]
    z = data_z[i, :]
    datpos, datmom = np.matmul(inv, [z, dp])
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
         label=r'tau={:.2f} turns'.format(1/gr[0]))
ax1.plot((starts + ends)/2, newdat, color='k', marker='x', linestyle='None')
ax1.legend()
ax1.set_title('Mode {:d}, mu={:.3f}'.format(mode, mode/Nbunches))
plt.show()


gr_pyat = gr[0] * f0  # growth rate in seconds

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot([mode], [gr_pyat],
         color='r', marker='x', linestyle='None', label='PyAT')
ax1.plot(gr_theory, color='b', marker='o', linestyle='None', label='Theory')
ax1.set_xlabel('Mode')
ax1.set_ylabel('Growth Rate [1/s]')
ax1.legend()
plt.show()
