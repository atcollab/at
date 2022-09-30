import numpy as np
from at.constants import clight, e_mass, qe
import matplotlib.pyplot as plt
import at
from at.collective import Wake, LongResonatorElement
from at.tracking.utils import get_bunches_std_mean
from mpi4py import MPI
import time
import sys
import pickle as pkl


class Logger(object):
    def __init__(self, fName):

        self.terminal = sys.stdout
        self.log = open(fName, "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        # this flush method is needed for python 3 compatibility.
        # this handles the flush command by doing nothing.
        # you might want to specify some extra behavior here.
        pass


def launch():

    sys.stdout = Logger('./logfiler.txt')
    sys.stderr = Logger('./logfiler.txt')

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    print(size, rank)

    # Set up the ring
    ring = at.load_m('../../../machine_data/esrf.m')

    nturns = 5000
    Nbunches = 992

    # Number of particles PER CORE must be a multiple of Nbunches
    Npart = int(Nbunches)

    current = 200e-3

    bunchSpacing = ring.circumference / Nbunches
    ring.beam_current = current
    ring.set_rf_voltage(6e6)

    ring.radiation_off(cavity_pass='RFCavityPass')
    ring.set_fillpattern(Nbunches)

    ring.set_rf_frequency()
    freqres = ring.rf_frequency

    qx, qy, qs = ring.get_tune()
    fs0 = qs * ring.revolution_frequency

    ld0, ld = ring.find_orbit()
    fring, _ = at.fast_ring(ring)

    # Define the resonator parameters and current
    wturns = 40
    srange = Wake.build_srange(0., 0.41, 1.0e-3, 1.0e-1,
                               bunchSpacing, ring.circumference*wturns)

    # This is needed to fix an interpolation issue with LongResElem
    srange = np.sort(np.concatenate((srange, [1e-24])))

    fres = ring.get_rf_frequency() - ring.revolution_frequency - fs0
    qfactor = 1e4
    rshunt = 2e6

    welem = LongResonatorElement('wake', ring, srange, fres, qfactor,
                                 rshunt, Nturns=wturns, Nslice=1)
    fring.append(welem)

    # Particle generation and tracking
    sigm = at.sigma_matrix(ring.radiation_on(copy=True))
    part = at.beam(Npart, sigm, orbit=ld0)
    part[4, :] = 0
    part[5, :] = 0

    if rank == 0:
        z_all = np.zeros((nturns, Nbunches))
        dp_all = np.zeros((nturns, Nbunches))
        turnMsk = np.zeros(nturns, dtype=bool)

    for i in np.arange(nturns):
        _ = at.lattice_pass(fring, part)
        allPartsg = comm.gather(part)
        if rank == 0:
            allPartsg = np.concatenate(allPartsg, axis=1)
            stds, means = get_bunches_std_mean(allPartsg, Nbunches)
            dp_all[i] = np.array([x[4] for x in means])
            z_all[i] = np.array([x[5] for x in means])

            sumNan = np.sum(np.isnan(allPartsg))
            if sumNan > 0:
                turnMsk[i] = True
            print(i, sumNan, np.mean(np.abs(dp_all[i])))

    if rank == 0:
        outDict = pkl.dump({'dp': dp_all, 'z': z_all, 'lostturn': ~turnMsk,
                            'ring': ring, 'fres': fres, 'qfactor': qfactor,
                            'rshunt': rshunt, 'M': Nbunches,
                            'current': current},
                           open('./LCBI_output.pkl', 'wb'))

if __name__ == '__main__':
    launch()
