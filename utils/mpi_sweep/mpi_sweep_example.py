#!/usr/bin/env python3

import at
import numpy
from mpi_sweep import mpi_sweep

# this code will run on all MPI nodes
d1 = at.elements.Drift('DR01', 3)
qf = at.elements.Quadrupole('QF', 1, 0.2)
d2 = at.elements.Drift('DR02', 3)
qd = at.elements.Quadrupole('QD', 1, -0.2)
fodocell = numpy.array([d1, qf, d2, qd])
thering = at.lattice.Lattice(fodocell, energy=1e6)

def generator():
    """
    Generate input parameter list for worker function.
    It will only be executed on node 0 to save on allocations.
    """
    return [(0.1,-0.1), (0.2,-0.2), (0.3,-0.3)]

def worker(arg, **kwargs):
    """
    Function to run on multiple MPI nodes with argument arg provided
    by the generator function.
    This function may run multiple times on one node.
    """
    # additional arguments may be used
    # kwargs['name']  # name of the processor
    # kwargs['rank']  # rank (id) of the node
    # kwargs['size']  # total number of nodes
    qf.K = arg[0]
    qd.K = arg[1]
    return at.find_m44(thering)[0]

def collector(args):
    """
    Collect computed data.
    This function will only be executed on node 0.
    """
    with open('output.csv', 'w') as f:
        f.write('Kqf,Kqd')
        [f.write(f',m{i}{j}') for i in range(0, 4) for j in range(0, 4)]
        f.write('\n')
        for case in args:
            params, output = case
            f.write(f'{params[0]},{params[1]}')
            [f.write(f',{output[i,j]}') for i in range(0, 4) for j in range(0, 4)]
            f.write('\n')

mpi_sweep(generator, worker, collector)
