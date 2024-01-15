import glob
import os
from os.path import basename, exists, join, splitext
import sys
from setuptools import setup, Extension

# Numpy build dependency defined in pyproject.toml.
import numpy


def select_omp():
    if exists('/usr/local/include/omp.h'):              # Homebrew
        return '-I/usr/local/include', '/usr/local/lib'
    elif exists('/opt/local/include/libomp/omp.h'):     # MacPorts
        return '-I/opt/local/include/libomp', '/opt/local/lib/libomp'
    else:
        raise FileNotFoundError('\n'.join((
            '',
            'libomp.dylib must be installed with your package manager:',
            '',
            'Use "$ brew install libomp"',
            'Or  "$ sudo port install libomp"',
            ''
        )))


print("** Entering setup.py:", str(sys.argv))
print("** MPI:", os.environ.get('MPI', None))
print("** OPENMP:", os.environ.get('OPENMP', None))
print("** CUDA:", os.environ.get('CUDA', None))
print("** OPENCL:", os.environ.get('OPENCL', None))
macros = [('PYAT', None)]
with_openMP = False

cflags = ['-std=c99']
cppflags = []

mpi = eval(os.environ.get('MPI', 'None'))
if not mpi or (len(sys.argv) >= 2 and
               any(sys.argv[1] == arg for arg in ('egg_info', 'sdist'))):
    mpi_macros = []
    mpi_includes = []
else:
    mpi_macros = [('MPI', None)]
    os.environ["CC"] = 'mpicc'
    os.environ["CXX"] = 'mpicxx'
    try:
        import mpi4py
    except ImportError:
        print('\npyAT with MPI requires mpi4py. '
              'Please install mpi4py: "pip install mpi4py"\n')
        sys.exit()
    mpi_includes = mpi4py.get_include()


omp = eval(os.environ.get('OPENMP', 'None'))
if not omp:
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


cuda = eval(os.environ.get('CUDA', 'None'))
if not cuda:
    cuda_cppflags = []
    cuda_lflags = []
    cuda_macros = []
else:
    # Generate the shared include file for the GPU kernel
    exec(open('atgpu/genheader.py').read())
    cuda_path = os.environ.get('CUDA_PATH', None)
    cuda_macros = [('CUDA', None)]
    if cuda_path is None:
        raise RuntimeError('CUDA_PATH environment variable not defined')
    if sys.platform.startswith('win'):
        cuda_cppflags = ['-I' + cuda_path + '\\include']
        cuda_lflags = ['/LIBPATH:'+cuda_path+'\\lib\\x64', "cuda.lib", "nvrtc.lib"]
    else:
        cuda_cppflags = ['-I' + cuda_path + '/include']
        cuda_lflags = ['-L' + cuda_path + '/lib64', '-Wl,-rpath,' + cuda_path + '/lib64', '-lcuda', '-lnvrtc']

opencl = eval(os.environ.get('OPENCL', 'None'))
if not opencl:
    opencl_cppflags = []
    opencl_lflags = []
    opencl_macros = []
else:
    # Generate the shared include file for the GPU kernel
    exec(open('atgpu/genheader.py').read())
    opencl_ocl_path = os.environ.get('OCL_PATH', None)
    opencl_sdk_path = os.environ.get('OPENCL_SDK_PATH', None)
    opencl_macros = [('OPENCL', None)]
    if opencl_ocl_path is None:
        raise RuntimeError('OCL_PATH environment variable not defined')
    if opencl_sdk_path is None:
        raise RuntimeError('OPENCL_SDK_PATH environment variable not defined')
    if sys.platform.startswith('win'):
        opencl_cppflags = ['-I' + opencl_sdk_path + '\\include','-I' + opencl_ocl_path + '\\include']
        opencl_lflags = ['/LIBPATH:'+opencl_ocl_path+'\\lib\\x64', "OpenCL.lib"]
    else:
        opencl_cppflags = ['-I' + opencl_ocl_path + '/include','-I' + opencl_sdk_path + '/include']
        opencl_lflags = ['-L' + opencl_ocl_path + '/lib64', '-Wl,-rpath,' + opencl_ocl_path + '/lib64', '-lOpenCL']

if not sys.platform.startswith('win32'):
    cflags += ['-Wno-unused-function']

# It is easier to copy the integrator files into a directory inside pyat
# for packaging. However, we cannot always rely on this directory being
# inside a clone of at and having the integrator sources available, and
# this file is executed each time any setup.py command is run.
# It appears that only copying the files when they are available is
# sufficient.
integrator_src_orig = 'atintegrators'
diffmatrix_orig = join('atmat', 'atphysics', 'Radiation')

c_pass_methods = glob.glob(join(integrator_src_orig, '*Pass.c'))
cpp_pass_methods = glob.glob(join(integrator_src_orig, '*Pass.cc'))
diffmatrix_source = join(diffmatrix_orig, 'findmpoleraddiffmatrix.c')
at_source = join('pyat', 'at.c')


def c_integrator_ext(pass_method):
    name, _ = splitext(basename(pass_method))
    name = ".".join(('at', 'integrators', name))
    return Extension(
        name=name,
        sources=[pass_method],
        include_dirs=[numpy.get_include(), mpi_includes, integrator_src_orig],
        define_macros=macros + omp_macros + mpi_macros,
        extra_compile_args=cflags + omp_cflags,
        extra_link_args=omp_lflags
    )


def cpp_integrator_ext(pass_method):
    name, _ = splitext(basename(pass_method))
    name = ".".join(('at', 'integrators', name))
    return Extension(
        name=name,
        sources=[pass_method],
        include_dirs=[numpy.get_include(), mpi_includes, integrator_src_orig],
        define_macros=macros + omp_macros + mpi_macros,
        extra_compile_args=cppflags + omp_cflags,
        extra_link_args=omp_lflags
    )


at = Extension(
    'at.tracking.atpass',
    sources=[at_source],
    define_macros=macros + omp_macros + mpi_macros,
    include_dirs=[numpy.get_include(), integrator_src_orig],
    extra_compile_args=cflags + omp_cflags,
    extra_link_args=omp_lflags
)

cconfig = Extension(
    'at.cconfig',
    sources=[join('pyat', 'at', 'cconfig.c')],
    define_macros=macros + omp_macros + mpi_macros + cuda_macros + opencl_macros,
    extra_compile_args=cflags + omp_cflags,
)

diffmatrix = Extension(
    name='at.physics.diffmatrix',
    sources=[diffmatrix_source],
    include_dirs=[numpy.get_include(), integrator_src_orig],
    define_macros=macros,
    extra_compile_args=cflags
)

gpusource = [join('atgpu', 'AbstractGPU.cpp'),
             join('atgpu', 'AbstractInterface.cpp'),
             join('atgpu', 'PyInterface.cpp'),
             join('atgpu', 'Lattice.cpp'),
             join('atgpu', 'PassMethodFactory.cpp'),
             join('atgpu', 'SymplecticIntegrator.cpp'),
             join('atgpu', 'IdentityPass.cpp'),
             join('atgpu', 'DriftPass.cpp'),
             join('atgpu', 'StrMPoleSymplectic4Pass.cpp'),
             join('atgpu', 'BndMPoleSymplectic4Pass.cpp'),
             join('atgpu', 'StrMPoleSymplectic4RadPass.cpp'),
             join('atgpu', 'BndMPoleSymplectic4RadPass.cpp'),
             join('atgpu', 'CavityPass.cpp'),
             join('atgpu', 'RFCavityPass.cpp')]

cudaext = Extension(
    name='at.tracking.gpu',
    sources=gpusource + [join('atgpu', 'CudaGPU.cpp')],
    define_macros=macros + cuda_macros,
    extra_compile_args=cppflags + cuda_cppflags,
    extra_link_args=cuda_lflags
)

openclext = Extension(
    name='at.tracking.gpu',
    sources=gpusource + [join('atgpu', 'OpenCLGPU.cpp')],
    define_macros=macros + opencl_macros,
    extra_compile_args=cppflags + opencl_cppflags,
    extra_link_args=opencl_lflags
)

setup(
    ext_modules=[at, cconfig, diffmatrix] +
                ([cudaext] if cuda else []) +
                ([openclext] if opencl else []) +
                [c_integrator_ext(pm) for pm in c_pass_methods] +
                [cpp_integrator_ext(pm) for pm in cpp_pass_methods],
)

print("** Leaving setup.py")
