from distutils.core import setup, Extension
from distutils import sysconfig
import numpy
import sys
import os
import glob

macros = [('PYAT', None)]

integrator_src = os.path.abspath('../atintegrators')
integrator_build = None

cflags = []

suffix = sysconfig.get_config_var('EXT_SUFFIX')
if suffix is None:
    if sys.platform.startswith('win32'):
        suffix = '.pyd'
    else:
        suffix = '.so'

if not sys.platform.startswith('win32'):
    cflags += ['-Wno-unused-function']


def integrator_extension(pass_method):
    name = ".".join(('at', 'integrators', os.path.basename(pass_method)[:-2]))
    return Extension(name=name,
                     sources=[pass_method],
                     include_dirs=[numpy.get_include(), integrator_src],
                     define_macros=macros,
                     extra_compile_args=cflags)


integ_list = glob.glob(os.path.join(integrator_src, '*Pass.c'))
dist = setup(name='at.integrators', package_dir={'at': ''}, packages=['at.integrators'],
             ext_modules=[integrator_extension(pm) for pm in integ_list])
try:
    install_location = dist.command_obj['install'].install_platlib
    if integrator_build is None:
        integrator_build = '"{}"'.format(os.path.join(install_location, 'at', 'integrators', '%s{}'.format(suffix)))
        macros.append(('INTEGRATOR_PATH', integrator_build))
    at = Extension('at.atpass', sources=['at.c'],
                   define_macros=macros,
                   include_dirs=[numpy.get_include(), integrator_src],
                   extra_compile_args=cflags)
    setup(name='at', package_dir={'at': ''}, packages=['at'], ext_modules=[at])

except KeyError:
    print('\npyat should be built in one step by calling "setup.py install"\n')
