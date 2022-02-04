import numpy
from at import envelope_parameters, get_mcf
from at.constants import clight, e_mass, qe
from scipy.interpolate import interp1d
import time

class Haissinski(object):

    '''
    Class to find the longitudinal distribution in the presence 
    of a short range wakefield.

    This class is a direct implementation of the following paper:
    "Numerical solution of the Haïssinski equation for the equilibrium 
    state of  a stored electron beam", R. Warnock, K.Bane, Phys. Rev. Acc.
    and Beams 21, 124401 (2018)

    The reader is referred to this paper for a description of the methods.
    The equation number of key formula are written next to the relevant function.

    As input:
    wake_element is a wake_object that contains a WakeT and WakeZ array. 
    ring is a ring instance which is needed for machine parameters (sigma_l, sigma_e, etc)
    
    m is the number of points in the full distribution that you want
    kmax is the min and max of the range of the distribution. See equation 27. In units of sigma_z0
    current is the bunch current.
    numIters is the number of iterations
    eps is the convergence criteria.

    Future developments of this class:
        When a high current is provided, it is better to solve for intermediate currents and use this solution as start for higher current otherwise it may not converge
        Adding LR wake or harmonic cavity as done at SOLEIL. Needs to be added WITH this class which is just for short range wake.
    '''

    def __init__(self, wake_element, ring, m=12, kmax=1, current=1e-4, numIters = 10, eps = 1e-10):

        self.circumference = ring.circumference
        self.energy = ring.energy

        envelpars = envelope_parameters(ring.radiation_on(copy=True))
        self.f_s = envelpars.f_s
        self.nu_s = self.f_s / (clight/self.circumference)
        self.sigma_e = envelpars.sigma_e * self.energy
        self.sigma_l = envelpars.sigma_l

        self.eta = get_mcf(ring.radiation_off(copy=True))
        self.ga = self.energy/e_mass
        self.betrel = numpy.sqrt(1.0-1.0/self.ga/self.ga)


        self.numIters = numIters
        self.eps = eps    


        #negative s to be consistent with paper and negative Wz
        s = wake_element.WakeT
        self.ds = numpy.diff(s)[0]/self.sigma_l
        self.wtot_fun = interp1d(-s/self.sigma_l, -wake_element.WakeZ) 

        self.s = -s/self.sigma_l

        if m%2!=0:
            raise AttributeError('m must be even and int')
        self.m = m
        self.npoints = 2*self.m + 1
        
        self.kmax = kmax
        self.dq = self.kmax/self.m

        self.q_array = -self.kmax + numpy.arange(self.npoints)*self.dq
        self.set_weights()


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
        sr = numpy.arange(2*numpy.amin(self.q_array), numpy.abs(2*numpy.amin(self.q_array)) + self.ds, self.ds)
        res = numpy.zeros(len(sr))

        topend = numpy.trapz(self.wtot_fun(numpy.arange(numpy.amax(sr), numpy.amax(self.s), self.ds)), x=numpy.arange(numpy.amax(sr), numpy.amax(self.s), self.ds))
        for p, low in enumerate(sr):
            if p % 10000 == 0:
                print(p, ' out of ', len(sr))
            rr = numpy.arange(low, numpy.amax(sr), self.ds)
            val = numpy.trapz(self.wtot_fun(rr), x=rr) + topend
            res[p] = val
    
        self.Sfun_range = sr #Not used except for plotting
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
        self.I = numpy.sign(self.eta)*qe*self.N/(2*numpy.pi*self.nu_s*self.sigma_e) #removed one qe for eV to joule
        print('Normalised current: ', self.I*1e12, ' pC/V')

    def initial_phi(self):
        '''
        Simply a gaussian but using the normalised units.
        Page 5 top right, in the text.
        '''
        self.phi = numpy.exp(-self.q_array**2/2)*self.I/numpy.sqrt(2*numpy.pi)
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
                    sum1b += self.weights[k] * self.Smat[j,k] * self.phi[k]    
                sum1 += self.weights[j] * numpy.exp(-self.q_array[j]**2/2 + sum1b)
            sum2 = 0
            for j in numpy.arange(self.npoints):
                sum2 += self.weights[j] * self.Smat[i,j] * self.phi[j]
            self.allFi[i] = self.phi[i] * sum1 - self.I * numpy.exp(-self.q_array[i]**2/2 + sum2)
   
    def dFi_ij(self, i, j):

        sum1 = 0
        for k in numpy.arange(self.npoints):

            sum1b = 0
            for l in numpy.arange(self.npoints):

                sum1b += self.weights[l] * self.Smat[k,l] * self.phi[l]

            #in paper i starts from 1, but kron(0,0) is 0 in python. 
            #No good for python, so i need this to be able to loop from 0
            kron = 0 if i != j else 1 

            sum1 += self.weights[k] * (kron + self.phi[i] * self.weights[j] * self.Smat[k,j]) * numpy.exp(-self.q_array[k]**2/2 + sum1b)

        sum2 = 0
        for k in numpy.arange(self.npoints):
            sum2 += self.weights[k] * self.Smat[i, k] * self.phi[k]

        val = sum1 - self.I*self.weights[j] * self.Smat[i, j] * numpy.exp(-self.q_array[i]**2/2 + sum2)    
       
        return val

    def dFi_dphij(self):
        '''
        Equation 30
        '''
        try:
            self.dtFlag
        except:
            t0 = time.time()
            self.dFi_ij(0,0)
            dt = time.time() - t0
            print('Computing dFi/dphij, will take approximately ', numpy.round(dt*self.npoints**2/60, 3), ' minutes per iteration') 
            self.dtFlag = True
        allDat = [self.dFi_ij(i,j) for i in numpy.arange(self.npoints) for j in numpy.arange(self.npoints)]
        self.alldFi_dphij = numpy.reshape(allDat, (self.npoints, self.npoints))

    def compute_new_phi(self):
        self.pseudo_inv = numpy.linalg.inv(self.alldFi_dphij).dot(-self.allFi) 
        self.phi_1 = self.pseudo_inv + self.phi

    def update(self):
        self.phi = self.phi_1.copy()
        #self.phi_1 = None         
   
    
    def convergence(self):
        '''
        Equation 32
        '''
        self.conv = (numpy.sum(numpy.abs(self.phi_1 - self.phi)))/numpy.sum(numpy.abs(self.phi))


    def solve(self):
        for it in numpy.arange(self.numIters):
            self.Fi()
            self.dFi_dphij()
            self.compute_new_phi()
            self.convergence()
            print('Iteration: ', it, ', Delta: ', self.conv)
            if self.conv < self.eps:
                print('Converged')
                break
            if it != self.numIters-1:
                self.update()
        if self.conv > self.eps:
            print('Did not converge, maybe solve for an intermediate current then use the solution as a starting point')
          







