"""Implements PEP 517 and PEP 660 hooks for custom build"""

import os

from setuptools import build_meta as _orig
from setuptools.build_meta import *  # noqa: F401, F403


def _get_option(config, key, default_value):
    val = os.environ.get(key.upper(), default_value)
    if config is not None:
        val = config.get(key.lower(), val)
    return eval(val)


def _set_environ(func):
    """Set environment for building according to config"""

    def build(wheel_dir, config_settings=None, metadata_dir=None):
        print("** Entering set_environ. " "Config settings:", config_settings)
        os.environ["MPI"] = str(_get_option(config_settings, "MPI", "0"))
        os.environ["OPENMP"] = str(_get_option(config_settings, "OPENMP", "0"))
        os.environ["CUDA"] = str(_get_option(config_settings, "CUDA", "0"))
        os.environ["OPENCL"] = str(_get_option(config_settings, "OPENCL", "0"))
        omp_threshold = _get_option(config_settings, "OMP_PARTICLE_THRESHOLD", "None")
        if omp_threshold is not None:
            os.environ["OMP_PARTICLE_THRESHOLD"] = str(omp_threshold)
        print("** MPI:", os.environ.get("MPI", "None"))
        print("** OPENMP:", os.environ.get("OPENMP", "None"))
        print("** CUDA:", os.environ.get("CUDA", "None"))
        print("** OPENCL:", os.environ.get("OPENCL", "None"))
        ret = func(wheel_dir, config_settings, metadata_dir)
        print("** Leaving set_environ")
        return ret

    return build


def _requirements(func):
    """Add build requirements according to config"""

    def get_requires(config_settings=None):
        print("** Entering get_requires. Config settings:", config_settings)
        mpi = _get_option(config_settings, "MPI", "0")
        print("** MPI:", mpi)
        addlist = ["mpi4py"] if mpi else []
        print("** Additional modules:", addlist)
        ret = func(config_settings) + addlist
        print("** Leaving get_requires")
        return ret

    return get_requires


build_wheel = _set_environ(_orig.build_wheel)
build_editable = _set_environ(_orig.build_editable)

get_requires_for_build_wheel = _requirements(_orig.get_requires_for_build_wheel)
get_requires_for_build_editable = _requirements(_orig.get_requires_for_build_editable)
