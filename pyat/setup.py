import glob
import os
from os.path import abspath, basename, dirname, exists, join, splitext
import sys
import shutil
from setuptools import setup, Extension

# Numpy build dependency defined in pyproject.toml.
import numpy


def select_omp():
    if exists('/usr/local/include/omp.h'):              # Homebrew
        return '-I/usr/local/include', '/usr/local/lib'
    elif exists('/opt/local/include/libomp/omp.h'):     # MacPorts
        return '-I/opt/local/include/libomp', '/opt/local/lib/libomp'
    else:
        raise FileNotFoundError('\n'.join(('',
          'libomp.dylib must be installed with your favourite package manager:',
          '',
          'Use "$ brew install libomp"',
          'Or  "$ sudo port install libomp"',
          ''
        )))


here = abspath(dirname(__file__))
macros = [('PYAT', None)]
with_openMP = False

cflags = []


mpi = os.environ.get('MPI', None)
if mpi is None:
    mpi_macros = []
    mpi_includes = []
else:
    mpi_macros = [('MPI', None)]
    os.environ["CC"] = 'mpicc'
    try:
        import mpi4py
    except ImportError:
        print('\npyAT with MPI requires mpi4py. '
              'Please install mpi4py: "pip install mpi4py"\n')
        sys.exit()
    mpi_includes = mpi4py.get_include()


omp = os.environ.get('OPENMP', None)
if omp is None:
    omp_cflags = []
    omp_lflags = []
    omp_macros = []
else:
    # Get the location of an alternate OpenMP library
    # Example: OMP_MATLAB=$MATLABROOT/sys/os/glnx64
    omp_path = os.environ.get('OMP_MATLAB', None)
    # Get the threshold on the number of particles
    omp_threshold = int(os.environ.get('OMP_PARTICLE_THRESHOLD', 10))
    omp_macros = [('OMP_PARTICLE_THRESHOLD', omp_threshold)]
    if sys.platform.startswith('win'):
        omp_cflags = ['/openmp']
        omp_lflags = []
    elif sys.platform.startswith('darwin'):
        omp_inc, omp_lib = select_omp()
        omp_cflags = ['-Xpreprocessor', '-fopenmp', omp_inc]
        if omp_path is None:
            omp_lflags = ['-L' + omp_lib, '-lomp']
        else:
            omp_lflags = ['-L' + omp_path, '-Wl,-rpath,' + omp_path, '-liomp5']
    else:
        omp_cflags = ['-fopenmp']
        if omp_path is None:
            omp_lflags = ['-lgomp']
        else:
            omp_lflags = ['-L' + omp_path, '-Wl,-rpath,' + omp_path, '-liomp5']

if not sys.platform.startswith('win32'):
    cflags += ['-Wno-unused-function']

# It is easier to copy the integrator files into a directory inside pyat
# for packaging. However, we cannot always rely on this directory being
# inside a clone of at and having the integrator sources available, and
# this file is executed each time any setup.py command is run.
# It appears that only copying the files when they are available is
# sufficient.
integrator_src_orig = abspath(join(here, '..', 'atintegrators'))
integrator_src = 'integrator-src'
diffmatrix_orig = abspath(join(here, '..', 'atmat', 'atphysics', 'Radiation'))

if exists(integrator_src_orig):
    # Copy files into pyat for distribution.
    source_files = glob.glob(join(integrator_src_orig, '*.[ch]'))
    source_files.extend(glob.glob(join(integrator_src_orig, '*.cc')))
    source_files.append(join(diffmatrix_orig, 'findmpoleraddiffmatrix.c'))
    if not exists(integrator_src):
        os.makedirs(integrator_src)
    for f in source_files:
        shutil.copy2(f, integrator_src)

pass_methods = glob.glob(join(integrator_src, '*Pass.c'))
pass_methods.extend(glob.glob(join(integrator_src, '*Pass.cc')))
diffmatrix_source = join(integrator_src, 'findmpoleraddiffmatrix.c')
at_source = 'at.c'


def integrator_ext(pass_method):
    name, _ = splitext(basename(pass_method))
    name = ".".join(('at', 'integrators', name))
    return Extension(
        name=name,
        sources=[pass_method],
        include_dirs=[numpy.get_include(), mpi_includes, integrator_src],
        define_macros=macros + omp_macros + mpi_macros,
        extra_compile_args=cflags + omp_cflags,
        extra_link_args=omp_lflags
    )


at = Extension(
    'at.tracking.atpass',
    sources=[at_source],
    define_macros=macros + omp_macros + mpi_macros,
    include_dirs=[numpy.get_include(), integrator_src],
    extra_compile_args=cflags + omp_cflags,
    extra_link_args=omp_lflags
)

diffmatrix = Extension(
    name='at.physics.diffmatrix',
    sources=[diffmatrix_source],
    include_dirs=[numpy.get_include(), integrator_src],
    define_macros=macros,
    extra_compile_args=cflags
)

setup(
    ext_modules=[at, diffmatrix] + [integrator_ext(pm) for pm in pass_methods],
)
