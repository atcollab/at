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
print("** CUDA:", os.environ.get("CUDA", None))
print("** OPENCL:", os.environ.get("OPENCL", None))

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


cuda = eval(os.environ.get("CUDA", "None"))
if not cuda:
    cuda_cppflags = []
    cuda_lflags = []
    cuda_macros = []
else:
    # Generate the shared include file for the GPU kernel
    exec(open("atgpu/genheader.py").read())
    cuda_path = os.environ.get("CUDA_PATH", None)
    cuda_macros = [("CUDA", None)]
    if cuda_path is None:
        raise RuntimeError("CUDA_PATH environment variable not defined")
    if sys.platform.startswith("win"):
        cuda_cppflags = ["-I" + cuda_path + "\\include"]
        cuda_lflags = ["/LIBPATH:" + cuda_path + "\\lib\\x64", "cuda.lib", "nvrtc.lib"]
    else:
        cuda_cppflags = ["-I" + cuda_path + "/include"]
        cuda_lflags = [
            "-L" + cuda_path + "/lib64",
            "-Wl,-rpath," + cuda_path + "/lib64",
            "-lcuda",
            "-lnvrtc",
        ]

opencl = eval(os.environ.get("OPENCL", "None"))
if not opencl:
    opencl_cppflags = []
    opencl_lflags = []
    opencl_macros = []
else:
    # Generate the shared include file for the GPU kernel
    exec(open("atgpu/genheader.py").read())
    opencl_ocl_path = os.environ.get("OCL_PATH", None)
    opencl_macros = [("OPENCL", None)]
    if sys.platform.startswith("win"):
        # Static link
        if opencl_ocl_path is None:
            raise RuntimeError("OCL_PATH environment variable not defined")
        opencl_cppflags = ["-I" + opencl_ocl_path + "\\include"]
        opencl_lflags = [
            "/LIBPATH:" + opencl_ocl_path + "\\lib",
            "OpenCL.lib",
            "cfgmgr32.lib",
            "runtimeobject.lib",
            "Advapi32.lib",
            "ole32.lib",
        ]
    elif sys.platform.startswith("darwin"):
        opencl_cppflags = ["-std=c++11"]
        opencl_lflags = ["-Wl,-framework,OpenCL"]
    else:
        if opencl_ocl_path is not None:
            # Private install
            opencl_cppflags = ["-I" + opencl_ocl_path + "/include"]
            opencl_lflags = [
                "-L" + opencl_ocl_path + "/lib",
                "-Wl,-rpath," + opencl_ocl_path + "/lib",
                "-lOpenCL",
            ]
        elif exists("/usr/include/CL/opencl.h"):
            # Standard install
            opencl_cppflags = []
            opencl_lflags = ["-lOpenCL"]
        else:
            raise RuntimeError(
                "Install OpenCL include and driver (ICD) in standard path or set OCL_PATH environment variable"
            )

if not sys.platform.startswith("win32"):
    cflags += ["-Wno-unused-function"]

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
gpu_pass_methods = glob.glob(join("atgpu", "*Pass.cpp"))
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
    define_macros=macros + omp_macros + mpi_macros + cuda_macros + opencl_macros,
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

gpusource = gpu_pass_methods + [
    join("atgpu", "AbstractGPU.cpp"),
    join("atgpu", "AbstractInterface.cpp"),
    join("atgpu", "PyInterface.cpp"),
    join("atgpu", "Lattice.cpp"),
    join("atgpu", "PassMethodFactory.cpp"),
    join("atgpu", "SymplecticIntegrator.cpp"),
]

cudaext = Extension(
    name="at.tracking.gpu",
    sources=gpusource + [join("atgpu", "CudaGPU.cpp")],
    define_macros=macros + cuda_macros,
    include_dirs=[np.get_include()],
    extra_compile_args=cppflags + cuda_cppflags,
    extra_link_args=cuda_lflags,
)

openclext = Extension(
    name="at.tracking.gpu",
    sources=gpusource + [join("atgpu", "OpenCLGPU.cpp")],
    define_macros=macros + opencl_macros,
    include_dirs=[np.get_include()],
    extra_compile_args=cppflags + opencl_cppflags,
    extra_link_args=opencl_lflags,
)

setup(
    ext_modules=[at, cconfig, diffmatrix, wiggdiffmatrix]
    + [integrator_ext(pm, cflags) for pm in c_pass_methods]
    + [integrator_ext(pm, cppflags) for pm in cpp_pass_methods]
    + ([cudaext] if cuda else [])
    + ([openclext] if opencl else []),
)

print("** Leaving setup.py")
