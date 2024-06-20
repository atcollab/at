import at
from at.collective import Wake, LongResonatorElement
from mpi4py import MPI
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

    ring.disable_6d(cavity_pass='RFCavityPass')
    ring.set_fillpattern(Nbunches)

    ring.set_rf_frequency()

    qx, qy, qs = ring.get_tune()
    fs0 = qs * ring.revolution_frequency

    ld0, ld = ring.find_orbit()
    fring, _ = at.fast_ring(ring)

    # Define the resonator parameters and current
    wturns = 40
    srange = Wake.build_srange(0., 0.41, 1.0e-3, 1.0e-1,
                               bunchSpacing, ring.circumference*wturns)

    fres = ring.get_rf_frequency() - ring.revolution_frequency - fs0
    qfactor = 1e4
    rshunt = 2e6

    welem = LongResonatorElement('wake', ring, srange, fres, qfactor,
                                 rshunt, Nturns=wturns, Nslice=1)
    fring.append(welem)

    bmon = at.BeamMoments('mon')
    fring.append(bmon)

    # Particle generation and tracking
    sigm = at.sigma_matrix(ring.enable_6d(copy=True))
    part = at.beam(Npart, sigm, orbit=ld0)
    part[4, :] = 0
    part[5, :] = 0

    fring.track(part, nturns=nturns, refpts=None, in_place=True)

    dp_all = bmon.means[4, :, :]
    z_all = bmon.means[5, :, :]
    if rank == 0:
        outdict = {'dp': dp_all.T, 'z': z_all.T,
                   'ring': ring, 'fres': fres, 'qfactor': qfactor,
                   'rshunt': rshunt, 'M': Nbunches,
                   'current': current}
        pkl.dump(outdict, open('./LCBI_output.pkl', 'wb'))


if __name__ == '__main__':
    launch()
