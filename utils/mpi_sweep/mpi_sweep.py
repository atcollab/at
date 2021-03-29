from mpi4py import MPI
import numpy

def mpi_sweep(generator, worker, collector):
    """
    Evaluate function on multiple MPI nodes and collect the result.
    generator: function returning list of worker arguments
    worker: worker function
    collector: function collecting calculation results
    """
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    name = MPI.Get_processor_name()

    data = None
    if rank == 0:
        data = numpy.array_split(generator(), size)

    chunk = comm.scatter(data)

    output = [worker(par, name=name, rank=rank, size=size)
              for par in chunk]

    data = comm.gather(list(zip(chunk, output)))

    if rank == 0:
        return collector([d for sublist in data for d in sublist])
