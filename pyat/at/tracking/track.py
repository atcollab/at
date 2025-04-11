from __future__ import annotations

__all__ = [
    "lattice_track",
    "element_track",
    "internal_lpass",
    "internal_epass",
    "internal_plpass",
    "gpu_info",
    "gpu_core_count",
    "MPMode",
]

import multiprocessing
from collections.abc import Iterable
from functools import partial
from warnings import warn
from enum import Enum

import numpy

from .atpass import atpass as _atpass, elempass as _elempass, reset_rng
from .utils import fortran_align, has_collective, format_results
from .utils import initialize_lpass, disable_varelem, variable_refs
from ..lattice import AtError, AtWarning, DConstant, random
from ..lattice import Lattice, Element, Refpts, End
from ..lattice import get_uint32_index
from ..cconfig import iscuda
from ..cconfig import isopencl


class MPMode(Enum):
    """
    Multi Processing mode
    """

    CPU = 1  #: CPU multiprocessing
    GPU = 2  #: GPU multiprocessing


# Possible values for the use_mp flag, converts to use_mp,use_gpu bool pair
_mp_table = {
    False: (False, False),
    True: (True, False),
    MPMode.CPU: (True, False),
    MPMode.GPU: (False, True),
}

if iscuda() or isopencl():
    from .gpu import gpupass as _gpupass
    from .gpu import gpuinfo as _gpuinfo

_imax = numpy.iinfo(int).max
_globring: list[Element] | None = None


def _atpass_fork(seed, rank, rin, **kwargs):
    """Single forked job"""
    reset_rng(rank=rank, seed=seed)
    result = _atpass(_globring, rin, **kwargs)
    return rin, result


def _atpass_spawn(ring, seed, rank, rin, **kwargs):
    """Single spawned job"""
    reset_rng(rank=rank, seed=seed)
    result = _atpass(ring, rin, **kwargs)
    return rin, result


def _pass(ring, r_in, pool_size, start_method, seed, **kwargs):
    ctx = multiprocessing.get_context(start_method)
    # Split input in as many slices as processes
    args = enumerate(numpy.array_split(r_in, pool_size, axis=1))
    # Generate a new starting point for C RNGs
    global _globring
    _globring = ring
    if ctx.get_start_method() == "fork":
        passfunc = partial(_atpass_fork, seed, **kwargs)
    else:
        passfunc = partial(_atpass_spawn, ring, seed, **kwargs)
    # Start the parallel jobs
    with ctx.Pool(pool_size) as pool:
        results = pool.starmap(passfunc, args)
    _globring = None
    # Gather the results
    losses = kwargs.pop("losses", False)
    return format_results(results, r_in, losses)


@fortran_align
def _element_pass(element: Element, r_in, **kwargs):
    return _elempass(element, r_in, **kwargs)


@fortran_align
def _lattice_pass(
    lattice: list[Element],
    r_in,
    nturns: int = 1,
    refpts: Refpts = End,
    no_varelem=True,
    seed: int | None = None,
    use_gpu: bool = False,
    **kwargs,
):
    kwargs["reuse"] = kwargs.pop("keep_lattice", False)
    if no_varelem:
        lattice = disable_varelem(lattice)
    else:
        if sum(variable_refs(lattice)) > 0:
            kwargs["reuse"] = False
    refs = get_uint32_index(lattice, refpts)
    if seed is not None:
        reset_rng(seed=seed)
    if use_gpu:
        if not (iscuda() or isopencl()):
            raise AtError("No GPU support enabled")
        else:
            return _gpupass(lattice, r_in, nturns, refpts=refs, **kwargs)
    else:
        return _atpass(lattice, r_in, nturns, refpts=refs, **kwargs)


@fortran_align
def _plattice_pass(
    lattice: list[Element],
    r_in,
    nturns: int = 1,
    refpts: Refpts = End,
    seed: int | None = None,
    pool_size: int = None,
    start_method: str = None,
    **kwargs,
):
    verbose = kwargs.pop("verbose", False)
    refpts = get_uint32_index(lattice, refpts)
    any_collective = has_collective(lattice)
    kwargs["reuse"] = kwargs.pop("keep_lattice", False)
    rshape = r_in.shape
    if len(rshape) >= 2 and rshape[1] > 1 and not any_collective:
        if seed is None:
            seed = random.common.integers(0, high=_imax, dtype=int)
        if pool_size is None:
            pool_size = min(
                len(r_in[0]), multiprocessing.cpu_count(), DConstant.patpass_poolsize
            )
        if start_method is None:
            start_method = DConstant.patpass_startmethod
        return _pass(
            lattice,
            r_in,
            pool_size,
            start_method,
            seed=seed,
            nturns=nturns,
            refpts=refpts,
            **kwargs,
        )
    else:
        if seed is not None:
            reset_rng(seed=seed)
        if any_collective:
            warn(
                AtWarning("Collective PassMethod found: use single process"),
                stacklevel=2,
            )
        else:
            warn(
                AtWarning("no parallel computation for a single particle"), stacklevel=2
            )
        return _atpass(lattice, r_in, nturns=nturns, refpts=refpts, **kwargs)


def lattice_track(
    lattice: Iterable[Element],
    r_in,
    nturns: int = 1,
    refpts: Refpts = End,
    in_place: bool = False,
    **kwargs,
):
    """
    :py:func:`track_function` tracks particles through each element of a
    lattice or throught a single Element calling the element-specific
    tracking function specified in the Element's *PassMethod* field.

    Usage:
      >>> lattice_track(lattice, r_in)
      >>> lattice.track(r_in)

    Parameters:
        lattice: list of elements
        r_in: (6, N) array: input coordinates of N particles.
          *r_in* is modified in-place only if *in_place* is
          :py:obj:`True` and reports the coordinates at
          the end of the element. For the best efficiency, *r_in*
          should be given as F_CONTIGUOUS numpy array.

    Keyword arguments:
        nturns: number of turns to be tracked
        refpts: Selects the location of coordinates output.
          See ":ref:`Selecting elements in a lattice <refpts>`"
        in_place (bool): If True *r_in* is modified in-place and
          reports the coordinates at the end of the element.
          (default: False)
        seed (int | None): Seed for the random generators. If None (default)
          continue the sequence
        keep_lattice (bool):    Use elements persisted from a previous
          call. If :py:obj:`True`, assume that the lattice has not changed
          since the previous call.
        keep_counter (bool):    Keep the turn number from the previous
          call.
        turn (int):             Starting turn number. Ignored if
          *keep_counter* is :py:obj:`True`. The turn number is necessary to
          compute the absolute path length used in RFCavityPass.
        losses (bool):          Boolean to activate loss maps output
        omp_num_threads (int):  Number of OpenMP threads
          (default: automatic)
        use_mp (bool | MPMode):  Flag to activate multiprocessing
          (default: False)
        pool_size:              number of processes used when
          *use_mp* is :py:obj:`True` or :py:attr:`MPMode.CPU`. If None,
          ``min(npart,nproc)`` is used. It can be globally set using the
          variable *at.lattice.DConstant.patpass_poolsize*
        gpu_pool: List of GPU to use when *use_mp* is :py:attr:`MPMode.GPU`.
          If none, first GPU is selected.
        start_method:           python multiprocessing start method.
          :py:obj:`None` uses the python default that is considered safe.
          Available values: ``'fork'``, ``'spawn'``, ``'forkserver'``.
          Default for linux is ``'fork'``, default for macOS and  Windows is
          ``'spawn'``. ``'fork'`` may be used on macOS to speed up the
          calculation or to solve Runtime Errors, however it is considered
          unsafe. Used only when *use_mp* is :py:obj:`True`. It can be globally
          set using the variable *at.lattice.DConstant.patpass_startmethod*

    The following keyword arguments overload the lattice values

    Keyword arguments:

        particle (Optional[Particle]): circulating particle.
          Default: :code:`lattice.particle` if existing,
          otherwise :code:`Particle('relativistic')`
        energy (Optiona[float]): lattice energy. Default 0.
        unfold_beam (bool): Internal beam folding activate, this
          assumes the input particles are in bucket 0, works only
          if all bucket see the same RF Voltage.
          Default: :py:obj:`True`

    If *energy* is not available, relativistic tracking if forced,
    *rest_energy* is ignored.

    Returns:
        r_out: (6, N, R, T) array containing output coordinates of N particles
          at R reference points for T turns
        trackparam: A dictionary containing tracking input parameters with the
          following keys:

          ==============    ===================================================
          **npart**         number of particles
          **rout**          final particle coordinates
          **turn**          starting turn
          **refpts**        array of index where particle coordinate are saved
                            (only for lattice tracking)
          **nturns**        number of turn
          ==============    ===================================================

        trackdata: A dictionary containing tracking data with the following
          keys:

          ==============    ===================================================
          **loss_map**:     recarray containing the loss_map (only for lattice
                            tracking)
          ==============    ===================================================


          The **loss_map** is filled only if *losses* is :py:obj:`True`,
          it contains the following keys:

          ==============    ===================================================
          **islost**        (npart,) bool array indicating lost particles
          **turn**          (npart,) int array indicating the turn at
                            which the particle is lost
          **element**       (npart,) int array indicating the element at
                            which the particle is lost
          **coord**         (npart, 6) float array giving the coordinates at
                            which the particle is lost (zero for surviving
                            particles)
          ==============    ===================================================


    .. note::

       * :pycode:`track_function(lattice, r_in, refpts=len(line))` is the same
         as :pycode:`track_function(lattice, r_in)` since the reference point
         len(line) is the exit of the last element.
       * :pycode:`track_function(lattice, r_in, refpts=0)` is a copy of *r_in*
         since the reference point 0 is the entrance of the first element.
       * To resume an interrupted tracking (for instance to get intermediate
         results), one must use one of the *turn* or *keep_counter*
         keywords to ensure the continuity of the turn number.
       * For multiparticle tracking with large number of turn the size of
         *r_out* may increase excessively. To avoid memory issues
         :pycode:`track_function(lattice, r_in, refpts=None, in_place=True)`
         can be used. An empty list is returned and the tracking results of
         the last turn are stored in *r_in*.
       * To model buckets with different RF voltage :pycode:`unfold_beam=False`
         has to be used. The beam can be unfolded using the function
         :py:func:`.unfold_beam`. This function takes into account
         the true voltage in each bucket and distributes the particles in the
         bunches defined by :code:`ring.fillpattern` using a 6D orbit search.
    """

    use_mp = kwargs.pop("use_mp", False)
    use_mp, use_gpu = _mp_table[use_mp]

    trackdata = {}
    trackparam = {}
    part_kw = ["energy", "particle"]
    try:
        npart = numpy.shape(r_in)[1]
    except IndexError:
        npart = 1

    [trackparam.update((kw, kwargs.get(kw))) for kw in kwargs if kw in part_kw]
    trackparam.update({"npart": npart})

    if not in_place:
        r_in = r_in.copy()

    lattice = initialize_lpass(lattice, nturns, kwargs)
    ldtype = [
        ("islost", numpy.bool_),
        ("turn", numpy.uint32),
        ("elem", numpy.uint32),
        ("coord", numpy.float64, (6,)),
    ]
    loss_map = numpy.recarray((npart,), ldtype)
    lat_kw = ["turn"]
    [trackparam.update((kw, kwargs.get(kw))) for kw in kwargs if kw in lat_kw]
    trackparam.update({"refpts": get_uint32_index(lattice, refpts), "nturns": nturns})

    start_method = kwargs.pop("start_method", None)
    pool_size = kwargs.pop("pool_size", None)
    if use_mp:
        kwargs.update({"pool_size": pool_size, "start_method": start_method})
        rout = _plattice_pass(lattice, r_in, nturns=nturns, refpts=refpts, **kwargs)
    else:
        rout = _lattice_pass(
            lattice,
            r_in,
            nturns=nturns,
            refpts=refpts,
            use_gpu=use_gpu,
            no_varelem=False,
            **kwargs,
        )

    if kwargs.get("losses", False):
        rout, lm = rout
        lm["coord"] = lm["coord"].T
        for k, v in lm.items():
            loss_map[k] = v

    trackdata.update({"loss_map": loss_map})
    trackparam.update({"rout": r_in})

    return rout, trackparam, trackdata


def element_track(element: Element, r_in, in_place: bool = False, **kwargs):
    """
    :py:func:`element_track` tracks particles through one element of a
    calling the element-specific tracking function specified in the
    Element's *PassMethod* field

    Usage:
      >>> element_track(element, r_in)
      >>> element.track(r_in)

    Parameters:
        element: element to track through
        r_in: (6, N) array: input coordinates of N particles.
          For the best efficiency, *r_in*
          should be given as F_CONTIGUOUS numpy array.

    Keyword arguments:
        in_place (bool): If True *r_in* is modified in-place and
          reports the coordinates at the end of the element.
          (default: False)
        omp_num_threads (int):  Number of OpenMP threads
          (default: automatic)
        particle (Optional[Particle]): circulating particle.
          Default: :code:`lattice.particle` if existing,
          otherwise :code:`Particle('relativistic')`
        energy (Optiona[float]): lattice energy. Default 0.

    Returns:
        r_out: (6, N, R, T) array containing output coordinates of N particles
          at R reference points for T turns
    """
    if not in_place:
        r_in = r_in.copy()

    rout = _element_pass(element, r_in, **kwargs)
    return rout


def gpu_core_count(gpuPool: list[int] | None):
    """
    :py:func:`gpu_core_count` returns number of GPU core.

    Parameters:
        gpuPool: List of selected GPUs. If None, first GPU is selected.
    """
    gpus = gpu_info()
    nbgpu = len(gpus)
    nbcore = 0
    if gpuPool is None:
        gpuPool = [0]
    for g in gpuPool:
        if g < 0 or g >= nbgpu:
            raise AtError(
                "Invalid GPU id " + str(g) + ", must be in [0.." + str(nbgpu - 1) + "]"
            )
        nbcore += gpus[g][2]  # Compute unit number
    return nbcore


def gpu_info():
    """
    :py:func:`gpu_info` returns list of GPU present on the system and their
    corresponding information. If GPU support is not enabled or if no capable device
    are present on the system, an empty list is returned.

    Returns:
        gpu: [gpu name,hardware version (CUDA device),compute unit number,platform].
          The number of compute units is not the so-called number of "CUDA cores".
          The number of threads that can be executed simultaneously depends on the
          hardware and on the type of used instructions. For best performance, it is
          recommended to track a number of particles which is multiple of 64 (or 128)
          times the number of compute units.
    """
    if iscuda() or isopencl():
        return _gpuinfo()
    else:
        return []


internal_lpass = _lattice_pass
internal_epass = _element_pass
internal_plpass = _plattice_pass
Lattice.track = lattice_track
Element.track = element_track
