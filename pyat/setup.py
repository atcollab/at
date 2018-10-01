import sys
try:
    import numpy
except ImportError:
    print('\npyAT requires numpy. '
          'Please install numpy: "pip install numpy"\n')
    sys.exit()
from setuptools import setup, Extension
import os
import glob
import shutil


macros = [('PYAT', None)]

integrator_src_orig = os.path.abspath('../atintegrators')
integrator_src = './integrator-src'

# Copy files into pyat for distribution.
source_files = glob.glob(os.path.join(integrator_src_orig, '*.[ch]'))
if not os.path.exists(integrator_src):
    os.makedirs(integrator_src)
for f in source_files:
    shutil.copy2(f, integrator_src)

pass_methods = glob.glob(os.path.join(integrator_src, '*Pass.c'))

cflags = []

if not sys.platform.startswith('win32'):
    cflags += ['-Wno-unused-function']


def integrator_extension(pass_method):
    name, _ = os.path.splitext(os.path.basename(pass_method))
    name = ".".join(('at', 'integrators', name))
    return Extension(name=name,
                     sources=[pass_method],
                     include_dirs=[numpy.get_include(), integrator_src],
                     define_macros=macros,
                     extra_compile_args=cflags)


at = Extension('at.atpass',
               sources=['at.c'],
               define_macros=macros,
               include_dirs=[numpy.get_include(), integrator_src_orig],
               extra_compile_args=cflags)

diffmatrix = Extension(name='at.physics.diffmatrix',
                       sources=['../atmat/atphysics/Radiation/findmpoleraddiffmatrix.c'],
                       include_dirs=[numpy.get_include(), integrator_src_orig],
                       define_macros=macros,
                       extra_compile_args=cflags)

setup(name='at-python',
      version='0.0.1',
      description='Accelerator Toolbox',
      author='The AT collaboration',
      author_email='atcollab-general@lists.sourceforge.net',
      install_requires=['numpy>=1.10'],
      packages=['at', 'at.integrators', 'at.lattice', 'at.physics'],
      ext_modules=[at, diffmatrix] + [integrator_extension(pm) for pm in pass_methods],
      zip_safe=False)
