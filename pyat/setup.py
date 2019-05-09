import sys
try:
    import numpy
except ImportError:
    print('\npyAT requires numpy. '
          'Please install numpy: "pip install numpy"\n')
    sys.exit()
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
import os
import glob
import shutil


here = os.path.abspath(os.path.dirname(__file__))
macros = [('PYAT', None)]


class CopyDuringBuild(build_ext):
    """Only copy integrators and diffmatrix during build. This avoids filepath
    issues during wheel building where a temporary copy of pyAT is made.
    """
    def run(self):
        here = os.path.abspath(os.path.dirname(__file__))
        integrator_src_orig = os.path.abspath(os.path.join(here, '../atintegrators'))
        integrator_src = os.path.abspath(os.path.join(here, 'integrator-src'))
        diffmatrix_source = os.path.abspath(os.path.join(here, '../atmat/atphysics/Radiation'))
        # Copy files into pyat for distribution.
        source_files = glob.glob(os.path.join(integrator_src_orig, '*.[ch]'))
        source_files.extend(glob.glob(os.path.join(diffmatrix_source, 'findmpoleraddiffmatrix.c')))
        if not os.path.exists(integrator_src):
            os.makedirs(integrator_src)
        for f in source_files:
            shutil.copy2(f, integrator_src)
        build_ext.run(self)


at_source = os.path.abspath(os.path.join(here,'at.c'))
integrator_src_orig = os.path.abspath(os.path.join(here, '../atintegrators'))
integrator_src = os.path.abspath(os.path.join(here, 'integrator-src'))
diffmatrix_source = os.path.abspath(os.path.join(here, '../atmat/atphysics/Radiation'))
pass_methods = glob.glob(os.path.join(integrator_src, '*Pass.c'))
diffmatrix_method = os.path.join(integrator_src, 'findmpoleraddiffmatrix.c')

cflags = []

if not sys.platform.startswith('win32'):
    cflags += ['-Wno-unused-function']


def integrator_extension(pass_method):
    name, _ = os.path.splitext(os.path.basename(pass_method))
    name = ".".join(('at', 'integrators', name))
    return Extension(name=name,
                     sources=[pass_method],
                     include_dirs=[numpy.get_include(), integrator_src, diffmatrix_source],
                     define_macros=macros,
                     extra_compile_args=cflags)


at = Extension('at.tracking.atpass',
               sources=[at_source],
               define_macros=macros,
               include_dirs=[numpy.get_include(), integrator_src, diffmatrix_source],
               extra_compile_args=cflags)

diffmatrix = Extension(name='at.physics.diffmatrix',
                       sources=[diffmatrix_method],
                       include_dirs=[numpy.get_include(), integrator_src, diffmatrix_source],
                       define_macros=macros,
                       extra_compile_args=cflags)

setup(cmdclass={'build_ext': CopyDuringBuild},
      name='accelerator-toolbox',
      version='0.0.1',
      description='Accelerator Toolbox',
      author='The AT collaboration',
      author_email='atcollab-general@lists.sourceforge.net',
      install_requires=['numpy>=1.10', 'scipy>=0.16'],
      packages=find_packages(),
      ext_modules=[at, diffmatrix] + [integrator_extension(pm) for pm in pass_methods],
      zip_safe=False,
      python_requires='>=2.7.4')
