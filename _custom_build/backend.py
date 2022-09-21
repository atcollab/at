from setuptools import build_meta as _orig
import os

# Mandatory hooks
build_sdist = _orig.build_sdist


def _get_option(config, key, default_value):
    val = os.environ.get(key.upper(), default_value)
    if config is not None:
        val = config.get(key.lower(), val)
    return eval(val)


def build_wheel(wheel_dir, config_settings=None, metadata_dir=None):
    print("** Entering build_wheel. "
          "Config settings:", config_settings)
    os.environ["MPI"] = str(_get_option(config_settings, 'MPI', '0'))
    os.environ["OPENMP"] = str(_get_option(config_settings, 'OPENMP', '0'))
    omp_threshold = _get_option(config_settings,
                                'OMP_PARTICLE_THRESHOLD', 'None')
    if omp_threshold is not None:
        os.environ["OMP_PARTICLE_THRESHOLD"] = str(omp_threshold)
    print("** MPI:", os.environ.get('MPI', 'None'))
    print("** OPENMP:", os.environ.get('OPENMP', 'None'))
    ret = _orig.build_wheel(wheel_dir, config_settings, metadata_dir)
    print("** Leaving build_wheel")
    return ret


# Optional hook
def get_requires_for_build_wheel(config_settings=None):
    print("** Entering get_requires_for_build_wheel. "
          "Config settings:", config_settings)
    mpi = _get_option(config_settings, 'MPI', '0')
    print("** MPI:", mpi)
    addlist = ["mpi4py"] if mpi else []
    print("** Additional modules:", addlist)
    ret = _orig.get_requires_for_build_wheel(config_settings) + addlist
    print("** Leaving get_requires_for_build_wheel")
    return ret
