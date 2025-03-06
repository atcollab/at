import sys
import os
from os.path import basename, exists, join, splitext
import glob

# Numpy build dependency defined in pyproject.toml.
import numpy as np
from setuptools import setup, Extension


def select_omp():
    if os.uname().machine.startswith("arm"):
        homeb = "/opt/homebrew/opt/libomp"
    else:
        homeb = "/usr/local"

    if exists(os.path.join(homeb, "include/omp.h")):  # Homebrew
        return "-I" + os.path.join(homeb, "include"), os.path.join(homeb, "lib")
    elif exists("/opt/local/include/libomp/omp.h"):  # MacPorts
        return "-I/opt/local/include/libomp", "/opt/local/lib/libomp"
    else:
        raise FileNotFoundError(
            "\n".join(
                (
                    "",
                    "libomp.dylib must be installed with your package manager:",
                    "",
                    'Use "$ brew install libomp"',
                    'Or  "$ sudo port install libomp"',
                    "",
                )
            )
        )


print("** Entering setup.py:", str(sys.argv))
print("** MPI:", os.environ.get("MPI", None))
print("** OPENMP:", os.environ.get("OPENMP", None))

platform = sys.platform
macros = [("PYAT", None)]

cflags = []
cppflags = []

if not platform.startswith("win32"):
    cflags += ["-std=c99", "-Wno-unused-function", "-Wno-unknown-pragmas"]
    cppflags += ["-Wno-unused-function", "-Wno-unknown-pragmas"]

if platform.startswith("darwin"):
    cflags += ["-ffp-contract=off"]
    cppflags += ["-ffp-contract=off"]

mpi = eval(os.environ.get("MPI", "None"))
if not mpi or (len(sys.argv) >= 2 and sys.argv[1] in {"egg_info", "sdist"}):
    mpi_macros = []
    mpi_includes = []
else:
    mpi_macros = [("MPI", None)]
    os.environ["CC"] = "mpicc"
    os.environ["CXX"] = "mpicxx"
    try:
        import mpi4py
    except ImportError:
        print(
            "\npyAT with MPI requires mpi4py. "
            'Please install mpi4py: "pip install mpi4py"\n'
        )
        sys.exit()
    mpi_includes = mpi4py.get_include()


omp = eval(os.environ.get("OPENMP", "None"))
if not omp:
    omp_cflags = []
    omp_lflags = []
    omp_macros = []
else:
    # Get the location of an alternate OpenMP library
    # Example: OMP_MATLAB=$MATLABROOT/sys/os/glnx64
    omp_path = os.environ.get("OMP_MATLAB", None)
    # Get the threshold on the number of particles
    omp_threshold = int(os.environ.get("OMP_PARTICLE_THRESHOLD", 10))
    omp_macros = [("OMP_PARTICLE_THRESHOLD", omp_threshold)]
    if platform.startswith("win"):
        omp_cflags = ["/openmp"]
        omp_lflags = []
    elif platform.startswith("darwin"):
        omp_inc, omp_lib = select_omp()
        omp_cflags = ["-Xpreprocessor", "-fopenmp", omp_inc]
        if omp_path is None:
            omp_lflags = ["-L" + omp_lib, "-lomp"]
        else:
            omp_lflags = ["-L" + omp_path, "-Wl,-rpath," + omp_path, "-liomp5"]
    else:
        omp_cflags = ["-fopenmp"]
        if omp_path is None:
            omp_lflags = ["-lgomp"]
        else:
            omp_lflags = ["-L" + omp_path, "-Wl,-rpath," + omp_path, "-liomp5"]

# It is easier to copy the integrator files into a directory inside pyat
# for packaging. However, we cannot always rely on this directory being
# inside a clone of at and having the integrator sources available, and
# this file is executed each time any setup.py command is run.
# It appears that only copying the files when they are available is
# sufficient.
integrator_src_orig = "atintegrators"
diffmatrix_orig = join("atmat", "atphysics", "Radiation")

c_pass_methods = glob.glob(join(integrator_src_orig, "*Pass.c"))
cpp_pass_methods = glob.glob(join(integrator_src_orig, "*Pass.cc"))
diffmatrix_source = join(diffmatrix_orig, "findmpoleraddiffmatrix.c")
at_source = join("pyat", "at.c")


def integrator_ext(pass_method, flags):
    name, _ = splitext(basename(pass_method))
    name = ".".join(("at", "integrators", name))
    return Extension(
        name=name,
        sources=[pass_method],
        include_dirs=[np.get_include(), mpi_includes, integrator_src_orig],
        define_macros=macros + omp_macros + mpi_macros,
        extra_compile_args=flags + omp_cflags,
        extra_link_args=omp_lflags,
    )


at = Extension(
    "at.tracking.atpass",
    sources=[at_source],
    define_macros=macros + omp_macros + mpi_macros,
    include_dirs=[np.get_include(), integrator_src_orig],
    extra_compile_args=cflags + omp_cflags,
    extra_link_args=omp_lflags,
)

cconfig = Extension(
    "at.cconfig",
    sources=[join("pyat", "at", "cconfig.c")],
    define_macros=macros + omp_macros + mpi_macros,
    extra_compile_args=cflags + omp_cflags,
)

diffmatrix = Extension(
    name="at.physics.diffmatrix",
    sources=[join(diffmatrix_orig, "findmpoleraddiffmatrix.c")],
    include_dirs=[np.get_include(), integrator_src_orig],
    define_macros=macros,
    extra_compile_args=cflags,
)

wiggdiffmatrix = Extension(
    name="at.physics.wiggdiffmatrix",
    sources=[join(diffmatrix_orig, "FDW.c")],
    include_dirs=[np.get_include(), integrator_src_orig],
    define_macros=macros,
    extra_compile_args=cflags,
)

setup(
    ext_modules=[at, cconfig, diffmatrix, wiggdiffmatrix]
    + [integrator_ext(pm, cflags) for pm in c_pass_methods]
    + [integrator_ext(pm, cppflags) for pm in cpp_pass_methods],
)

print("** Leaving setup.py")
