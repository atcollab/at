import sys
import numpy as np
import at
from at.collective import Wake, ResWallElement, WakeComponent
from mpi4py import MPI
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

    sys.stdout = Logger('./logfile.txt')
    sys.stderr = Logger('./logfile.txt')

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    print(size, rank)

    # Set up the ring
    ring = at.load_m('../../../machine_data/esrf.m')
    ring.disable_6d(cavity_pass='RFCavityPass')

    current = 200e-3
    nturns = 10000
    Nbunches = 992
    Npart = int(1 * Nbunches)

    bunchSpacing = ring.circumference/Nbunches
    ring.beam_current = current

    ring.set_fillpattern(Nbunches)

    ring.set_rf_frequency()
    freqres = ring.rf_frequency

    allSextInds = at.get_refpts(ring, at.elements.Sextupole)
    sdInds = []
    sfInds = []
    for i in allSextInds:
        if ring[i].PolynomB[2] > 0:
            sfInds.append(i)
        else:
            sdInds.append(i)
    sdInds = np.array(sdInds)
    sfInds = np.array(sfInds)
      
    at.fit_chrom(ring, sfInds, sdInds, [0, 0])

    l0, _, _ = ring.get_optics()

    fring, _ = at.fast_ring(ring)

    # Define the resonator parameters and current
    wturns = 30
    srange = Wake.build_srange(1e-4, 0.3, 1.0e-4, 1.0e-1, bunchSpacing,
                               ring.circumference * wturns)

    sigmac = 1e8
    rvac = 15e-3
    beta = 1
    welem = ResWallElement('wake', ring, srange, WakeComponent.DX,
                           ring.circumference, rvac, sigmac, beta,
                           Nslice=1, Nturns=wturns)
    fring.append(welem)

    bmon = at.BeamMoments('mon')
    fring.append(bmon)

    # Particle generation and tracking
    sigm = at.sigma_matrix(ring.enable_6d(copy=True))
    part = at.beam(Npart, sigm)
    part[0, :] = 1e-6
    part[1, :] = 0

    fring.track(part, nturns=nturns, refpts=None, in_place=True)

    x_all = bmon.means[0, :, :]
    xp_all = bmon.means[1, :, :]
    if rank == 0:
        outdict = {'x': x_all.T, 'xp': xp_all.T,
                   'ring': ring, 'sigmac': sigmac, 'rvac': rvac,
                   'M': Nbunches, 'alphax': l0['alpha'][0],
                   'betax': l0['beta'][0], 'current': current}
        pkl.dump(outdict, open('./TRW_output.pkl', 'wb'))


if __name__ == '__main__':
    launch()
