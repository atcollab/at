from distutils.core import setup, Extension
import numpy
import sys
import os
import glob

macros=[('PYAT', None)]

mach_home = os.getenv('MACH_HOME')
mach_arch = os.getenv('MACHARCH')
pyversion = "python{0.major}.{0.minor}".format(sys.version_info)

integrator_src = os.path.abspath('../atintegrators')
integrator_build = None


for pass_method in glob.glob(os.path.join(integrator_src, '*Pass.c')):
    name = os.path.basename(pass_method)[:-2]
    ext = Extension(name=name,
          sources=[pass_method],
          define_macros=macros,
          include_dirs=[numpy.get_include(),
                        integrator_src])
    dist = setup(name=name, ext_modules=[ext])
    try:
        # if installing, fetch the Python path
        install_location = dist.command_obj['install'].install_platlib
        if integrator_build is None:
            integrator_build = install_location.replace("\\", "/")
            macros.append(('INTEGRATOR_PATH', integrator_build))
    except KeyError:
        pass

print(integrator_build)
at = Extension('at', sources=['at.c'],
               define_macros=macros,
               include_dirs=[numpy.get_include(), integrator_src])

setup(name='at', py_modules=['elements','load_mat'], ext_modules=[at])

