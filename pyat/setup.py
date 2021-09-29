import glob
from io import open
import os
from os.path import abspath, basename, dirname, exists, join, splitext
import sys
import shutil
from setuptools import setup, Extension, find_packages

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


# Get the long description from the README file
with open(join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()


# It is easier to copy the integrator files into a directory inside pyat
# for packaging. However, we cannot always rely on this directory being
# inside a clone of at and having the integrator sources available, and
# this file is executed each time any setup.py command is run.
# It appears that only copying the files when they are available is
# sufficient.
at_source = abspath(join(here, 'at.c'))
integrator_src_orig = abspath(join(here, '..', 'atintegrators'))
integrator_src = abspath(join(here, 'integrator-src'))
diffmatrix_source = abspath(
    join(here, '..', 'atmat', 'atphysics', 'Radiation')
)

if exists(integrator_src_orig):
    # Copy files into pyat for distribution.
    source_files = glob.glob(join(integrator_src_orig, '*.[ch]'))
    source_files.extend(glob.glob(join(integrator_src_orig, '*.cc')))
    source_files.extend(
        glob.glob(join(diffmatrix_source, 'findmpoleraddiffmatrix.c'))
    )
    if not exists(integrator_src):
        os.makedirs(integrator_src)
    for f in source_files:
        shutil.copy2(f, integrator_src)

pass_methods = glob.glob(join(integrator_src, '*Pass.c'))
pass_methods.extend(glob.glob(join(integrator_src, '*Pass.cc')))
diffmatrix_method = join(integrator_src, 'findmpoleraddiffmatrix.c')


def integrator_ext(pass_method):
    name, _ = splitext(basename(pass_method))
    name = ".".join(('at', 'integrators', name))
    return Extension(
        name=name,
        sources=[pass_method],
        include_dirs=[numpy.get_include(), integrator_src, diffmatrix_source],
        define_macros=macros + omp_macros,
        extra_compile_args=cflags + omp_cflags,
        extra_link_args=omp_lflags
    )


at = Extension(
    'at.tracking.atpass',
    sources=[at_source],
    define_macros=macros + omp_macros,
    include_dirs=[numpy.get_include(), integrator_src, diffmatrix_source],
    extra_compile_args=cflags + omp_cflags,
    extra_link_args=omp_lflags
)

diffmatrix = Extension(
    name='at.physics.diffmatrix',
    sources=[diffmatrix_method],
    include_dirs=[numpy.get_include(), integrator_src, diffmatrix_source],
    define_macros=macros,
    extra_compile_args=cflags
)

setup(
    name='accelerator-toolbox',
    version='0.2.1',
    description='Accelerator Toolbox',
    long_description=long_description,
    author='The AT collaboration',
    author_email='atcollab-general@lists.sourceforge.net',
    url='https://github.com/atcollab/at',
    # Numpy 1.16.6 is the oldest version that builds with Python 3.9.
    install_requires=['numpy>=1.16.6', 'scipy>=0.16'],
    packages=find_packages(),
    ext_modules=[at, diffmatrix] + [integrator_ext(pm) for pm in pass_methods],
    zip_safe=False,
    python_requires='>=3.6.0',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ]
)
