.. _parallel:

Parallel processing
===================

In most passmethods, particles are independent. So multi-particle tracking is
well adapted for parallel processing. But a few elements, like impedance elements,
need the knowledge of the whole particle distribution. These are subclasses of
the :py:class:`.Collective` class. Depending on the multiprocessing option
chosen, they must be excluded from parallel tracking or be protected by a fence.

PyAT has several options for parallel computing. Though nothing wrong should
happen when simultaneously activating several options, it is usually advised
not to do so.

Python :py:obj:`multiprocessing` package
----------------------------------------
Some PyAT functions are based on the standard Python :obj:`multiprocessing`
package.

* Always available, on all platforms, with the default installation,
* Must be explicitly called, either by calling :py:func:`.patpass` instead
  of :py:func:`.lattice_pass` for tracking, or by specifying a ``use_mp``
  boolean keyword argument in some functions (see :py:func:`.get_acceptance` for
  example).
* :py:class:`.Collective` element classes are not allowed.

`OpenMP <https://www.openmp.org>`_
----------------------------------
PyAT C integrators can be compiled with `OpenMP <https://www.openmp.org>`_.
Parallelism is activated inside each integrator, when possible.

* Available on all platforms with a special installation procedure,
* Fully transparent once installed: `OpenMP <https://www.openmp.org>`_ is
  automatically activated when tracking more than ``OMP_PARTICLE_THRESHOLD``
  particles (10 by default). As a consequence, many internal PyAT functions
  (closed orbit, transfer matrices…) benefit from parallelism,
* No restriction on element classes,
* Uses shared-memory multiprocessing: small overhead, but limited to
  single-node multicore machines,

Installation
............
`OpenMP <https://www.openmp.org>`_ is integrated in most modern compilers and
does not require any preliminary installation, except for macOS:

.. admonition:: macOS preliminary OpenMP installation

   While the ``clang`` compiler integrates `OpenMP <https://www.openmp.org>`_,
   ``libomp.dylib`` must be installed with your favorite package manager,
   `Homebrew <https://brew.sh>`_ or `MacPorts <https://www.macports.org>`_:

   | Use ``$ brew install libomp``
   | or  ``$ sudo port install libomp``

From the PyPI repository, you must install from source (disable the binary
distribution)::

    pip install --no-binary accelerator-toolbox --config-settings openmp=1 accelerator-toolbox


From a local clone: when switching from a standard install
to an `OpenMP <https://www.openmp.org>`_ one or vice-versa, remove the
``build`` directory to force a re-compilation::

    cd <at>
    rm -rf build
    pip install --config-settings openmp=1 .

Configuration
.............
The minimum number of particles triggering parallelism can be configured at
compile time either by setting a ``OMP_PARTICLE_THRESHOLD`` environment variable
or by specifying the ``omp_particle_threshold`` configuration option::

    OPENMP=1 OMP_PARTICLE_THRESHOLD=4 pip install . # or
    pip install --config-settings openmp=1 --config-settings omp_particle_threshold=4 .

The number of threads is by default automatically selected by
`OpenMP <https://www.openmp.org>`_. It can be forced to any value with the
``omp_num_threads`` keyword argument of :py:func:`.lattice_pass`

`MPI <https://www.mpi-forum.org/docs/>`_
----------------------------------------
PyAT can be installed with `MPI <https://www.mpi-forum.org/docs/>`_
compatibility.

* Available on all platforms with a special installation procedure, but mostly
  interesting on computer clusters,
* Message-passing interface implies a significant overhead for data transfer,
  but allows large scale and heterogeneous computer infrastructures,
* No restriction on element classes,
* The whole python session needs to be started with ``mpirun``: better suited
  for massive multi-particle batch jobs than for interactive sessions,

Installation
............
`MPI <https://www.mpi-forum.org/docs/>`_ must be preliminary installed on the
computer.

From the PyPI repository, you must install from source (disable the binary
distribution)::

    pip install --no-binary accelerator-toolbox --config-settings mpi=1 "accelerator-toolbox[mpi]"

Installing PyAT from a local clone: when switching from a standard install
to a `MPI <https://www.mpi-forum.org/docs/>`_ one or vice-versa, remove the
``build`` directory to force a re-compilation::

    cd <at>
    rm -rf build
    pip install --config-settings mpi=1 ".[mpi]"

`GPU Tracking`_
---------------
PyAT can be installed with GPU support, either `OpenCL <https://github.com/KhronosGroup/OpenCL-Guide/tree/main>`_ or
`CUDA <https://developer.nvidia.com/cuda-toolkit>`_, compatibility. GPU are especially interesting for tracking large
number of particle. The performance of the tracking is mainly related to the GPU double precision arithmetic performance.

OpenCL Installation
...................
`OpenCL`_ and installable client driver (ICD) must be preliminary installed on the system using either linux packages or by
building the `OpenCL SDK <https://github.com/KhronosGroup/OpenCL-SDK>`_::

    sudo apt install opencl-headers ocl-icd-opencl-dev -y

Then you need to set the environment variable ``OCL_PATH`` to the SDK install path if you don't use a standard install::

    export OCL_PATH=<sdk_intall_path>

or (on Windows)::

    set OCL_PATH=C:\clpeak\build\sdk_instal

Note: `clpeak <https://github.com/krrishnarraj/clpeak>`_ is an OpenCL benchmarking tool that can be used to check system
performance (especially double precision floating point arithmetic) and to build the OpenCL SDK.

Install PyAT using the ``opencl`` flag::

    cd <at>
    rm -rf build
    pip install --config-settings opencl=1 .

You can check the install as bellow, without GPU support, the method ``at.tracking.gpu_info()`` will return an empty list::

    Z:\at>python
    Python 3.11.5 (tags/v3.11.5:cce6ba9, Aug 24 2023, 14:38:34) [MSC v.1936 64 bit (AMD64)] on win32
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import at
    >>> at.tracking.gpu_info()
    [['NVIDIA GeForce GTX 1050 Ti', '6.1', 6, 'NVIDIA CUDA OpenCL 3.0 CUDA 12.3.68'], ['Intel(R) UHD Graphics 630', '0.0', 24, 'Intel(R) OpenCL HD Graphics OpenCL 3.0 ']]

CUDA Installation
.................
NVidia `CUDA`_ toolkit must be preliminary installed on the system from `NVidia <https://developer.nvidia.com/cuda-downloads>`_.
Set the environment variable ``CUDA_PATH``::

    export CUDA_PATH=/cvmfs/hpc.esrf.fr/software/packages/ubuntu20.04/x86_64/cuda/12.3.1

or on Windows::

    set CUDA_PATH=C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.3

Install PyAT using the ``cuda`` flag::

    pip install --config-settings cuda=1 .

You can check the install using the method ``at.tracking.gpu_info()`` as described above.
