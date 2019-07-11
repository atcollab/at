import glob
from io import open
import os
from os.path import abspath, basename, dirname, exists, join, splitext
import sys
import shutil
try:
    import numpy
except ImportError:
    print('\npyAT requires numpy. '
          'Please install numpy: "pip install numpy"\n')
    sys.exit()
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext


here = abspath(dirname(__file__))
macros = [('PYAT', None)]

cflags = []

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
at_source = abspath(join(here,'at.c'))
integrator_src_orig = abspath(join(here, '..', 'atintegrators'))
integrator_src = abspath(join(here, 'integrator-src'))
diffmatrix_source = abspath(
    join(here, '..', 'atmat', 'atphysics', 'Radiation')
)

if exists(integrator_src_orig):
    diffmatrix_source = abspath(join(here, '../atmat/atphysics/Radiation'))
    # Copy files into pyat for distribution.
    source_files = glob.glob(join(integrator_src_orig, '*.[ch]'))
    source_files.extend(
        glob.glob(join(diffmatrix_source, 'findmpoleraddiffmatrix.c'))
    )
    if not exists(integrator_src):
        os.makedirs(integrator_src)
    for f in source_files:
        shutil.copy2(f, integrator_src)

pass_methods = glob.glob(join(integrator_src, '*Pass.c'))
diffmatrix_method = join(integrator_src, 'findmpoleraddiffmatrix.c')


def integrator_ext(pass_method):
    name, _ = splitext(basename(pass_method))
    name = ".".join(('at', 'integrators', name))
    return Extension(
        name=name,
        sources=[pass_method],
        include_dirs=[numpy.get_include(), integrator_src, diffmatrix_source],
        define_macros=macros,
        extra_compile_args=cflags
    )


at = Extension(
    'at.tracking.atpass',
    sources=[at_source],
    define_macros=macros,
    include_dirs=[numpy.get_include(), integrator_src, diffmatrix_source],
    extra_compile_args=cflags
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
    version='0.0.2',
    description='Accelerator Toolbox',
    long_description=long_description,
    author='The AT collaboration',
    author_email='atcollab-general@lists.sourceforge.net',
    url='https://pypi.org/project/accelerator-toolbox/',
    install_requires=['numpy>=1.10', 'scipy>=0.16'],
    packages=find_packages(),
    ext_modules=[at, diffmatrix] + [integrator_ext(pm) for pm in pass_methods],
    zip_safe=False,
    python_requires='>=2.7.4'
)
