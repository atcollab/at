import numpy
from at import radiation_parameters
from at.constants import clight, qe
from scipy.interpolate import interp1d
import time
from at.collective import Wake
from at.lattice import Lattice


class Haissinski(object):
    r"""
    Class to find the longitudinal distribution in the presence
    of a short range wakefield.

    Parameters: 
    
        wake_object:    class initialized with Wake that contains a 
          longitudinal wake Z and equivalent srange.
        ring:           a lattice object that is used to extract 
         machine parameters

    Keyword Args:
        m (int):        the number of points in the full
          distribution
        kmax:           the min and max of the range of the distribution.
          See equation 27. In units of [:math:`\sigma_z0`]
        current:        bunch current.
        numIters (int): the number of iterations
        eps:            the convergence criteria.

    This class is a direct implementation of the following paper:
    "Numerical solution of the Haïssinski equation for the equilibrium
    state of  a stored electron beam", R. Warnock, K.Bane, Phys. Rev. Acc.
    and Beams 21, 124401 (2018)

    The reader is referred to this paper for a description of the methods.
    The equation number of key formula are written next to the relevant
    function.

    The functions solve or solve_steps can be used after initialisation
    An example usage can be found in:

    at/pyat/examples/Collective/LongDistribution.py

    Future developments of this class:
        Adding LR wake or harmonic cavity as done at SOLEIL. Needs
        to be added WITH this class which is just for short range wake.
    """

    def __init__(self, wake_object: Wake, ring: Lattice, m: int = 12, kmax: float = 1,
                 current: float = 1e-4, numIters: int = 10, eps: float = 1e-10):

        self.circumference = ring.circumference
        self.energy = ring.energy

        radpars = radiation_parameters(ring)
        self.f_s = radpars.f_s
        self.nu_s = self.f_s / (clight/self.circumference)
        self.sigma_e = radpars.sigma_e * self.energy
        self.sigma_l = radpars.sigma_l

        #  The paper uses alpha and eta with same sign.
        #  AT uses opposite signs. So negative sign is needed for eta.
        self.eta = -ring.radiation_off(copy=True).slip_factor
        self.ga = ring.gamma
        self.betrel = ring.beta

        self.numIters = numIters
        self.eps = eps

        #  negative s to be consistent with paper and negative Wz
        s = wake_object.srange/self.sigma_l
        self.ds = numpy.diff(s)[-1]
        self.s = -s
        self.wtot_fun = interp1d(self.s, -wake_object.Z,
                                 bounds_error=False, fill_value = 0)

        if m % 2 != 0:
            raise AttributeError('m must be even and int')
        self.m = m
        self.npoints = 2*self.m + 1

        self.kmax = kmax
        self.dq = self.kmax/self.m

        self.q_array = -self.kmax + numpy.arange(self.npoints)*self.dq
        self.set_weights()
        self.dtFlag = False

        self.set_I(current)
        self.initial_phi()

        try:
            self.Sfun
        except AttributeError:
            self.precompute_S()

        self.compute_Smat()

    def set_weights(self):
        '''
        Page 7 second paragraph, in the text
        '''
        self.weights = numpy.ones(self.npoints)*2
        self.weights[1::2] = self.weights[1::2]**2
        self.weights[0] = 1
        self.weights[-1] = 1
        self.weights *= self.dq/3

    def precompute_S(self):
        '''
        Equation 16
        '''
        print('Computing integrated wake potential')
        sr = numpy.arange(2*numpy.amin(self.q_array),
                          numpy.abs(2*numpy.amin(self.q_array))
                          + self.ds, self.ds)
        res = numpy.zeros(len(sr))

        topend = numpy.trapz(self.wtot_fun(numpy.arange(numpy.amax(sr),
                             numpy.amax(self.s), self.ds)),
                             x=numpy.arange(numpy.amax(sr),
                             numpy.amax(self.s), self.ds))
        for p, low in enumerate(sr):
            if p % 10000 == 0:
                print(p, ' out of ', len(sr))
            rr = numpy.arange(low, numpy.amax(sr), self.ds)
            val = numpy.trapz(self.wtot_fun(rr), x=rr) + topend
            res[p] = val
        #  Not used except for plotting
        self.Sfun_range = sr
        self.Sfun = interp1d(sr, res)

    def compute_Smat(self):
        '''
        The sampling of the integrated wake potential S
        is only made at certain places. So all possibilities
        are loaded into a matrix for speed.
        '''
        self.Smat = numpy.zeros((self.npoints, self.npoints))
        for iq1, q1 in enumerate(self.q_array):
            for iq2, q2 in enumerate(self.q_array):
                self.Smat[iq1, iq2] = self.Sfun(q1-q2)

    def set_I(self, current):
        '''
        Equation 11
        '''
        self.N = current*self.circumference/(clight*self.betrel*qe)
        self.Ic = numpy.sign(self.eta)*qe*self.N/(2*numpy.pi *
                                                  self.nu_s*self.sigma_e)
        #  removed one qe for eV to joule
        print('Normalised current: ', self.Ic*1e12, ' pC/V')

    def initial_phi(self):
        '''
        Simply a gaussian but using the normalised units.
        Page 5 top right, in the text.
        '''
        self.phi = numpy.exp(-self.q_array**2/2)*self.Ic/numpy.sqrt(2*numpy.pi)
        self.phi_0 = self.phi.copy()

    def Fi(self):
        '''
        Equation 28
        '''
        self.allFi = numpy.zeros(self.npoints)
        for i in numpy.arange(self.npoints):
            sum1 = 0
            for j in numpy.arange(self.npoints):
                sum1b = 0
                for k in numpy.arange(self.npoints):
                    sum1b += self.weights[k] * self.Smat[j, k] * self.phi[k]
                sum1 += self.weights[j] * \
                    numpy.exp(-self.q_array[j]**2/2 + sum1b)
            sum2 = 0
            for j in numpy.arange(self.npoints):
                sum2 += self.weights[j] * self.Smat[i, j] * self.phi[j]
            self.allFi[i] = self.phi[i] * sum1 - self.Ic * \
                numpy.exp(-self.q_array[i]**2/2 + sum2)

    def dFi_ij(self, i, j):

        sum1 = 0
        for k in numpy.arange(self.npoints):

            sum1b = 0
            for li in numpy.arange(self.npoints):

                sum1b += self.weights[li] * self.Smat[k, li] * self.phi[li]

            #  in paper i starts from 1, but kron(0,0) is 0 in python.
            #  No good for python, so i need this to be able to loop from 0
            kron = 0 if i != j else 1

            sum1 += self.weights[k] * \
                (kron + self.phi[i] * self.weights[j] * self.Smat[k, j]) * \
                numpy.exp(-self.q_array[k]**2/2 + sum1b)

        sum2 = 0
        for k in numpy.arange(self.npoints):
            sum2 += self.weights[k] * self.Smat[i, k] * self.phi[k]

        val = sum1 - self.Ic*self.weights[j] * self.Smat[i, j] \
            * numpy.exp(-self.q_array[i]**2/2 + sum2)

        return val

    def dFi_dphij(self):
        '''
        Equation 30
        '''
        if not self.dtFlag:
            t0 = time.time()
            self.dFi_ij(0, 0)
            dt = time.time() - t0
            print('Computing dFi/dphij, will take approximately ',
                  numpy.round(dt*self.npoints**2/60, 3),
                  ' minutes per iteration')
            self.dtFlag = True
        allDat = [self.dFi_ij(i, j) for i in numpy.arange(self.npoints)
                  for j in numpy.arange(self.npoints)]
        self.alldFi_dphij = numpy.reshape(allDat, (self.npoints, self.npoints))

    def compute_new_phi(self):
        ainv = numpy.linalg.inv(self.alldFi_dphij)
        self.pseudo_inv = numpy.dot(ainv, -self.allFi)
        self.phi_1 = self.pseudo_inv + self.phi

    def update(self):
        self.phi = self.phi_1.copy()

    def convergence(self):
        '''
        Equation 32
        '''
        self.conv = (numpy.sum(numpy.abs(self.phi_1 - self.phi))) \
            / numpy.sum(numpy.abs(self.phi))

    def set_output(self):
        self.res = (self.phi_1.copy())[::-1]

    def solve(self):
        for it in numpy.arange(self.numIters):
            self.Fi()
            self.dFi_dphij()
            self.compute_new_phi()
            self.convergence()
            self.set_output()
            print('Iteration: ', it, ', Delta: ', self.conv)
            if self.conv < self.eps:
                print('Converged')
                break
            if it != self.numIters-1:
                self.update()
        if self.conv > self.eps:
            print('Did not converge, maybe solve for an intermediate '
                  'current then use the solution as a starting point')

    def solve_steps(self, currents):
        '''
        INPUT:
            currents  an array of currents to solve. If 0 is given,
                      a current of 10uA is used to prevent failure.
        '''
        self.I_steps = numpy.zeros(len(currents))
        self.res_steps = numpy.zeros((len(currents), len(self.q_array)))
        for ii, Ib in enumerate(currents):
            print('Running step ', ii+1, ' out of ', len(currents))
            #  If Ib = 0, the gaussian is zero. Small epsilon is given
            if Ib == 0:
                self.set_I(1e-5)
            else:
                self.set_I(Ib)
            self.I_steps[ii] = self.Ic
            self.solve()
            self.res_steps[ii, :] = self.res
