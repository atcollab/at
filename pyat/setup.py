from distutils.core import setup, Extension
import numpy
import sys
import os
import glob

macros = [('PYAT', None)]

integrator_src = os.path.abspath('../atintegrators')
integrator_build = None

cflags = []
if sys.platform.startswith('win32'):
    suffix = 'pyd'
else:
    cflags += ['-Wno-unused-function']
    suffix = 'so'


def integrator_extension(pass_method):
    name = ".".join(('pyat', 'atintegrators', os.path.basename(pass_method)[:-2]))
    return Extension(name=name,
                     sources=[pass_method],
                     include_dirs=[numpy.get_include(), integrator_src],
                     define_macros=macros,
                     extra_compile_args=cflags)


integ_list = glob.glob(os.path.join(integrator_src, '*Pass.c'))
dist = setup(name='pyat.atintegrators', package_dir={'pyat': ''}, packages=['pyat.atintegrators'],
             ext_modules=[integrator_extension(pm) for pm in integ_list])
try:
    install_location = dist.command_obj['install'].install_platlib
    if integrator_build is None:
        integrator_build = os.path.join(install_location, 'pyat', 'atintegrators', '.'.join(('%s', suffix)))
        macros.append(('INTEGRATOR_PATH', integrator_build))
    print(integrator_build)
    at = Extension('pyat.at', sources=['at.c'],
                   define_macros=macros,
                   include_dirs=[numpy.get_include(), integrator_src],
                   extra_compile_args=cflags)
    setup(name='pyat', package_dir={'pyat': ''}, packages=['pyat'], py_modules=['elements', 'load_mat'],
          ext_modules=[at])

except KeyError:
    print('\npyat should be built in one step by calling "setup.py install"\n')
    pass
