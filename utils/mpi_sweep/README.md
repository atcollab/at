`mpi_sweep.py` script allows parameter sweeping on computational clusters via MPI.

You may run example script locally `mpi_sweep_example.py`:
```bash
mpirun --host localhost:$(nproc) -n 3 ./mpi_sweep_example.py
```
