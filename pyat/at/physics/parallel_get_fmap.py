"""This script calculates the frequency and diffussion maps."""

# oblanco ALBA CELLS 2026jan06 use custom parallelization

from __future__ import annotations

from collections import defaultdict
from multiprocessing import Manager, Process, Queue, Value, cpu_count, current_process
from time import sleep

import numpy as np

# psutil is not a requirement but could be useful
try:
    # noinspection PyPackageRequirements
    import psutil

    PSUTIL_FOUND = True
except (ImportError, RuntimeError):
    MSG = (
        "psutil is unavailable => memory management is disabled.\n"
        + 'To enable memory management, run "pip install psutil"'
    )
    print(MSG)
    PSUTIL_FOUND = False

import at
from at.lattice import AtError, Lattice
from at.tracking import MPMode

# Jaime Coello de Portugal (JCdP) frequency analysis implementation
from .harmonic_analysis import get_tunes_harmonic

__all__ = ["get_fmap"]

# reference system coordinate names and indexes
_coord_index = np.array(
    [("x", 0), ("px", 1), ("y", 2), ("py", 3), ("dp", 4), ("ct", 5)]
)


def get_mem_info(pid: int) -> dict[str, int]:
    """Get memory info.

    Args:
        pid: process ID

    Returns:
        A dictionary with the pss memory from psutil.
    """
    res = defaultdict(int)
    mmap = psutil.Process(pid).memory_full_info()
    res["pss"] = mmap.pss
    return res


def get_ram_info() -> tuple:
    """Get ram memory info.

    Returns:
        Tuple of total ram, available ram, used ram, free ram,
        and percentage used obtained from psutil
    """
    ram = psutil.virtual_memory()
    total_ram = ram.total  # total installed RAM
    available_ram = ram.available  # RAM available for processes
    used_ram = ram.used  # RAM used by processes
    free_ram = ram.free  # RAM not being used
    percent_used = ram.percent  # percentage of RAM used
    return total_ram, available_ram, used_ram, free_ram, percent_used


def get_swap_info() -> tuple:
    """Get swap memory info.

    Returns:
        Tuple of total swap, used swap, free swap and percentage
        used obtained from psutil
    """
    swap = psutil.swap_memory()
    total_swap = swap.total  # total installed RAM
    used_swap = swap.used  # RAM used by processes
    free_swap = swap.free  # RAM not being used
    percent_used = swap.percent  # percentage of RAM used

    return total_swap, used_swap, free_swap, percent_used


def estimate_memory_resources(**kwargs: dict[str:any]) -> tuple:
    """Estimate the memory resources.

    Arguments:
        kwargs:
            verbose: default False
            mem_marging: minimum memory to keep free. Default np.nan.

    Returns:
        A tuple with the estimated free memory and the safe margin
    """
    # estimate memory resources
    verbose = kwargs.pop("verbose", False)
    mem_margin = kwargs.pop("mem_margin", np.nan) * 1024 * 1024
    ramtot, _, _, ramfree, _ = get_ram_info()
    swaptot, _, swapfree, _ = get_swap_info()
    onehundredmegabytes = 100 * 1024 * 1024
    safe_margin = min(onehundredmegabytes, 1 / 32 * ramtot, 1 / 4 * swaptot)
    if not np.isnan(mem_margin):
        safe_margin = mem_margin
    estimated_freemem = ramfree + swapfree
    if verbose:
        msg = (
            f"Estimated free memory : {estimated_freemem/1024/1024:.3f} MB, "
            + f"Memory safe margin : {safe_margin/1024/1024:.3f} MB"
        )
        print(msg)
    return estimated_freemem, safe_margin


def track_queue(
    ring: at.Lattice,
    zin: np.ndarray,
    nparticles: int,
    nturns: int,
    n_moving_slices: int,
    task_to_accomplish: list,
    nuaxes: tuple,
    worker_id: list,
    mem_limit: float,
    used_mem: list,
    mem_margin: float,
    turns_per_particle: list,
    elem_particle_lost: list,
    x1freq: list,
    x2freq: list,
    x1diff: list,
    x2diff: list,
    nudiff: list,
    verbose: bool,
    kwargs: dict,
) -> None:
    """Track particles with index inside queue.

    Arguments:
        ring: an pyat ring
        zin: input particles
        nparticles: number of particles
        nturns: number of turns.
        n_moving_slices: moving slices for frequency analysis
        task_to_accomplish: list of particles not yet tracked
        nuaxes: axes for frequency analysis
        worker_id: cpu ID processing a particle
        mem_limit: max memory allowed
        used_mem: estimation of the memor used
        mem_margin: minimun free memory at any time
        turns_per_particle: number of turns before the particle is lost
        elem_particle_lost: element where the particle is lost
        x1freq: tune on axis1
        x2freq: tune on axis2
        x1diff: tune std on axis1
        x2diff: tune std on axis2
        nudiff: 0.5*log10(dnuaxis1**2 + dnuaxis2**2)
        verbose: print info
        kwargs: keywords to pass to the tracking function
    """
    while not task_to_accomplish.empty():
        estimated_freemem = np.inf
        safe_margin = 0
        if PSUTIL_FOUND:
            estimated_freemem, safe_margin = estimate_memory_resources(
                verbose=verbose, mem_margin=mem_margin
            )
            if not np.isnan(mem_margin):
                safe_margin = mem_margin
        while estimated_freemem < safe_margin:
            if verbose:
                print("waiting for free memory ...")
            sleep(5)
        while used_mem.value > mem_limit:
            if verbose:
                print("using too much memory, waiting ...")
            sleep(5)
        task_id = task_to_accomplish.get(block=True, timeout=10)
        current = current_process()
        worker_id[task_id] = current.pid
        # we track twice the number of nturns in order to get two frequencies
        zout, _, loss_info = ring.track(
            zin[:, task_id], 2 * nturns, losses=True, **kwargs
        )
        msg_mem_info = ""
        if PSUTIL_FOUND:
            monitor_data = get_mem_info(current.pid)["pss"]
            with used_mem.get_lock():
                used_mem.value = used_mem.value + monitor_data
                msg_mem_info = f"Used mem: {used_mem.value/1024/1024:.3f} MB, "
        if verbose:
            msg = (
                msg_mem_info
                + f"cpu id {current.pid}, "
                + f"running particle id {task_id} of {nparticles}"
            )
            print(msg)
        # if the particle is lost
        if loss_info["loss_map"]["islost"][0]:
            turns_per_particle[task_id] = loss_info["loss_map"]["turn"][0]
            elem_particle_lost[task_id] = loss_info["loss_map"]["elem"][0]
        else:
            # get the tunes
            zout = zout.squeeze()
            tune_per_window = np.empty((2, n_moving_slices))
            tune_per_window[:] = np.nan
            shift_n = np.floor_divide(2 * nturns, n_moving_slices)
            for i in range(n_moving_slices):
                lim_inf = i * shift_n
                lim_sup = (i + 1) * shift_n
                z_1 = zout[nuaxes[0], lim_inf:lim_sup]
                z_2 = zout[nuaxes[1], lim_inf:lim_sup]
                z_1 = z_1 - z_1.mean()
                z_2 = z_2 - z_2.mean()
                tune1 = get_tunes_harmonic(z_1)
                tune2 = get_tunes_harmonic(z_2)
                if len(tune1) == 0 or len(tune2) == 0:
                    continue
                tune_per_window[0, i] = tune1[0]
                tune_per_window[1, i] = tune2[0]
            # metric
            dnu1, dnu2 = np.nanstd(tune_per_window, axis=1)
            dnu_norm = 0.5 * np.log10(dnu1 * dnu1 + dnu2 * dnu2)
            # min max diff
            if dnu_norm > -2:
                dnu_norm = -2
            if dnu_norm < -10:
                dnu_norm = -10
            # update
            turns_per_particle[task_id] = 2 * nturns
            x1freq[task_id] = np.nanmean(tune_per_window[0, :])
            x2freq[task_id] = np.nanmean(tune_per_window[1, :])
            x1diff[task_id] = dnu1
            x2diff[task_id] = dnu2
            nudiff[task_id] = dnu_norm
        if PSUTIL_FOUND:
            with used_mem.get_lock():
                used_mem.value = used_mem.value - monitor_data


def check_parallel_resources(
    use_mp: bool | at.MPMode,
    pool_size: int,
    kwargs: dict[str, any],
    verbose: bool,
) -> tuple:
    """Check parallel resources.

    Arguments:
        use_mp: Choose parallelization mode. True (CPU), at.MPMode.CPU, at.MPMode.GPU
        pool_size: integer number of workers to use
        kwargs: dictionary to include the gpu_pool parameters
        verbose: print info

    Returns:
        A tuple with the number of workers to upse and a dictionary with
        parallelization parameters

    Raises:
        AtError when MPMode.GPU is chosen
    """
    # For backward compatibility (use_mp can be a boolean)
    if use_mp is True:
        use_mp = MPMode.CPU
    nproc = cpu_count()
    if use_mp is MPMode.CPU:
        if np.isnan(pool_size):
            nprocu = nproc
        else:
            nprocu = pool_size
        if verbose:
            msg = (
                f"{nprocu} of {nproc} cpu used frequency map calculation\n"
                + "Multi-process acceptance calculation selected..."
            )
            print(msg)
        if nproc == 1 and verbose:
            print("Consider use_mp=False for single core computations")
    elif use_mp is MPMode.GPU:
        msg = "\nGPU acceptance calculation selected, but not yet implemented."
        raise AtError(msg)
    else:
        nprocu = 1
        if verbose:
            print("Single process acceptance calculation selected...")
        if nproc > 1 and verbose:
            print("Consider use_mp=True for parallelized computations")

    return nprocu, kwargs


def generate_grid(
    grid_size: list | tuple,
    window: list | tuple,
    offset: np.ndarray,
    shift_zero: float,
    axes: np.ndarray,
) -> tuple:
    """Generate a particle 2D grid.

    Arguments:
        grid_size: the number of points per axis
        window: [min_axis1,max_axis1,min_axis2,max_axis2]
        offset: 6D offset for the particle tracking
        shift_zero: small displacement to the 2D particles
        axes: defines the planes for the grid

    Returns:
        A tuple with the number of particles, the number of dimension,
        the 2D points, and an array of 6D particles

    Raises:
        AtError if the window is not well defined
    """
    if window[0] == window[1] or window[2] == window[3]:
        raise AtError("Window is too narrow")

    n_x1, n_x2 = grid_size
    min_x1 = min(window[0:2])
    max_x1 = max(window[0:2])
    min_x2 = min(window[2:4])
    max_x2 = max(window[2:4])
    x1vals = np.linspace(min_x1, max_x1, n_x1)
    x2vals = np.linspace(min_x2, max_x2, n_x2)

    nparticles = n_x1 * n_x2
    points = np.zeros((nparticles, 2))

    for _i in range(n_x1):
        for _j in range(n_x2):
            points[_i * n_x2 + _j, 0] = x1vals[_i]
            points[_i * n_x2 + _j, 1] = x2vals[_j]

    ndims = 6
    particles = np.zeros((ndims, nparticles))
    particles = particles + offset.reshape((6, 1))
    particles[axes[0], :] = particles[axes[0], :] + points[:, 0] + shift_zero
    particles[axes[1], :] = particles[axes[1], :] + points[:, 1] + shift_zero

    return nparticles, ndims, points, particles


def adapt_particles(user_particles: np.ndarray, axes: np.ndarray) -> tuple:
    """Adapt the particles given by the user.

    Arguments:
        user_particles: (6,n) particle array
        axes: planes to consider for the grid

    Returns:
        Tuple with number of particles, the number of dimensions, the points
        of reference, and the particles provided by the user
    """
    ndims, nparticles = user_particles.shape
    points = user_particles[axes, :]
    return nparticles, ndims, points, user_particles


def coord_to_indexes(useraxes: list | tuple) -> np.ndarray:
    """Transform user coordinate names to pyat indexes.

    Arguments:
        useraxes: user axes

    Returns:
        Numpy array with indexes corresponding to the user coordinate names
    """
    user_prefers_indexes = np.all([isinstance(i, int) for i in useraxes])
    if not user_prefers_indexes:
        axes = np.zeros((len(useraxes)), dtype=int)
        for _i, _ax in enumerate(useraxes):
            mask = _coord_index[:, 0] == _ax
            axes[_i] = int(np.squeeze(_coord_index[mask, 1][0]))
    else:
        axes = np.array(useraxes)
    return axes


def get_fmap(ring: Lattice, **kwargs: dict[str, any]) -> dict[str, any]:
    """Compute the frequency and diffussion maps of a ring.

    Arguments:
        ring: Lattice definition
        kwargs: any of the following.
            nturns: Number of turns for the tracking. See notes.
            axes: a list or tuple of coordinates; allowed values are:
            ``'x'``, ``'px'``, ``'y'``, ``'py'``, ``'dp'``, ``'ct'``
            or numbers from 0 to 5
            nuaxes: index of axes to extract the tunes. Default (0,2)
            n_moving_slices: number of moving windows to extract the tune.
            Default 2. See notes and ref. [2].
            particles: (6,n) particle array to track. If not given a rectangular
            grid of particles is created using 'window' and 'grid_size'
            window: list or tuple with (axis1_min,axis1_max,axis2_min,axis2_max)
            Note that window is only used if particles are not given.
            grid_size: (n_axis1,n_axis2) number of points per plane
            Note that grid_size is only used if particles are not given.
            offset: 6D offset for all particles. Default np.zeros((6))
            Note that offset is only used if particles are not given.
            shift_zero: small offset applied only on the axes. Default 10e-5
            Note that shift_zero is not used if particles are not given.
            verbose: print additional info, e.g. the cpu ID, the particle counter.
            use_mp: Default True. Only CPU parallelization is implemente, not GPU.
            pool_size: number of CPUs to use, otherwise it uses all available
            max_mem: memory usage limit in MB. If given, memory_margin is ignoredr,.
            otherwise, all the memory available is estimated and used as max.
            Requires psutil
            memory_margin: minimum memory to keep free in MB. If max_mem is not given
            it estimates the available and sets a margin. See notes.
            Requires psutil

    Returns:
        A dictionary containing
            name_of_axis1: the starting point on axis1 for all tracked particles
            name_of_axis2: the starting point on axis2 for all tracked particles
            "axes_indexes": array with the indexes of axis1 and axis2
            "axes_names": array with the axis1 and axis2 names
            "nu_axes_indexes": axes indexes used for frequency analysis
            "nu_axes_names": names of the axis for frequency analysis
            "nturns": number of tracked turns. The algorithm required 2*nturns.
            "nu_+name_of_axis1": tune on the axis1 per particle
            "nu_+name_of_axis2": tune on the axis2 per particle
            "dnu_+name_of_axis1": std of the tune on axis1
            "dnu_+name_of_axis2": std of the tune on axis2
            "dnulog10": log10 of the sqrt of the quadratic sum of dnu
            "turns_per_particle": number of turns before the particles was lost
            "elem_particle_lost": element where the particle was lost
            "worker_id": cpu ID who tracked the particle

    Examples:
        >>> fmapdata = get_fmap(
                ring,
                window=(-10e-3,10e-3,-10e-3,10e-3),
                offset=np.array([0,0,0,0,eoffset,0]),
                grid_size=(100,100),
                axes=('x','y'),
                n_moving_slices=10,
                shift_zero=1e-9,
                verbose=False,
                use_mp=at.MPMode.CPU,
                pool_size=pool_size,
                nturns=1024,
                max_mem=2048,
                mem_margin=0,
                )

        >>> # mask to filter only the surving particles
        >>> mask = ~np.isnan(fmapdata['dnulog10'])

    .. note::
       * The algorithm requires to track twice the number of turns to create
         at least two shifted windows where to evaluate the tune, and calculate
         the std.
       * The total number of turns tracked is 2*nturns. While the tune is extracted
         ``n_moving_slices`` times using a moving window. The tune variation is
         calculated as the standard deviation per plane. The variation per plane
         is added quadratically, then, log10(sqrt(dnu1**2+dnu2**2)) is returned
         as dnulog10.  More details on the moving window are in ref. [2]
       * ``use_mp=True`` or at.MPMode.CPUs are the only options available.
       * Memory arguments require psutil to have an effect
       * mem_margin is set as the minimum between 100 MB, 1/32 of total ram, and 1/4
         of swap memory
       * This routing also returns the number of turns the particle did in case
         it was lost.  This is complementary information to the frequency map as
         shown in ref [3].

       [1] J. Laskar, in Proceedings of the 20th Particle Accelerator Conference,
       Portland, OR, 2003 (IEEE, New York, 2003), p. 378.
       [2] D. Shatilov, E. Levichev, E. Simonov, M.. Zobov. Application of frequency map
       analysis to beam-beam effects study in crab waist collision scheme.
       Phys. Rev. ST Accel. Beams 14, 014001, 2011.
       [3] E. Serra-Carbonell, O Blanco, T. Guentzel. Application of fast algorithms to
       calculate dynamic and momentum aperture to the design of ALBA II. IPAC25, Taipei,
       Taiwan. MOPB065.
    """
    # initialize variables
    nturns = kwargs.pop("nturns", 1000)
    useraxes = kwargs.pop("axes", ("x", "px"))
    usernuaxes = kwargs.pop("nuaxes", (0, 2))
    n_moving_slices = kwargs.pop("n_moving_slices", 2)
    user_particles = kwargs.pop("particles", np.nan)
    window = kwargs.pop("window", (-10e-3, 10e-3, -10e-3, 10e-3))
    grid_size = kwargs.pop("grid_size", (10, 10))
    offset = kwargs.pop("offset", np.zeros((6)))
    shift_zero = kwargs.pop("shift_zero", 1e-5)
    verbose = kwargs.pop("verbose", False)
    use_mp = kwargs.pop("use_mp", True)
    pool_size = kwargs.pop("pool_size", np.nan)
    max_mem = kwargs.pop("max_mem", np.nan)
    mem_margin = kwargs.pop("mem_margin", np.nan * 1024 * 1024)

    verboseprint = print if verbose else lambda *a, **k: None

    nprocu, kwargs = check_parallel_resources(use_mp, pool_size, kwargs, verbose)

    axes = coord_to_indexes(useraxes)
    verboseprint(f"Using axes {_coord_index[axes,0]}")
    nuaxes = coord_to_indexes(usernuaxes)
    verboseprint(f"Using tune axes {_coord_index[nuaxes,0]}")

    if not np.any(np.isnan(user_particles)):
        nparticles, ndims, points, particles = adapt_particles(user_particles, axes)
    else:
        nparticles, ndims, points, particles = generate_grid(
            grid_size, window, offset, shift_zero, axes
        )

    verboseprint(f"Number of grid points: {nparticles}")

    task_to_accomplish = Queue()
    processes = []
    manager = Manager()
    for i in range(nparticles):
        task_to_accomplish.put(i)

    estimated_freemem = np.inf
    safe_margin = 0
    if PSUTIL_FOUND:
        estimated_freemem, safe_margin = estimate_memory_resources(
            verbose=verbose, mem_margin=mem_margin
        )
        if np.isnan(max_mem):
            if not np.isnan(mem_margin):
                print("User memory margin {mem_margin/1024/1024} MB")
                safe_margin = mem_margin
            max_mem = estimated_freemem - safe_margin
    mem_limit = max_mem * 1024 * 1024
    verboseprint(f"Memory usage limit set to {mem_limit/1024/1024} MB")

    # estimate memory usage
    floatsize = 8  # bytes in python
    n_in_parallel = min(nparticles, nprocu)
    needed_mem = (
        floatsize * ndims * n_in_parallel * nturns * 2
    )  # we track twice the nturns
    verboseprint(
        f"Estimated max. memory necessary for tracking: {needed_mem/1024/1024:.3} MB"
    )

    verboseprint("Tracking...")

    # creating processes
    worker_id = manager.list(nparticles * [np.nan])
    turns_per_particle = manager.list(nparticles * [np.nan])
    elem_particle_lost = manager.list(nparticles * [np.nan])
    x1freq = manager.list(nparticles * [np.nan])
    x2freq = manager.list(nparticles * [np.nan])
    x1diff = manager.list(nparticles * [np.nan])
    x2diff = manager.list(nparticles * [np.nan])
    nudiff = manager.list(nparticles * [np.nan])
    used_memory = Value("i", 0)
    for _worker in range(nprocu):
        prc = Process(
            target=track_queue,
            args=(
                ring,
                particles,
                nparticles,
                nturns,
                n_moving_slices,
                task_to_accomplish,
                nuaxes,
                worker_id,
                mem_limit,
                used_memory,
                mem_margin,
                turns_per_particle,
                elem_particle_lost,
                x1freq,
                x2freq,
                x1diff,
                x2diff,
                nudiff,
                verbose,
                kwargs,
            ),
        )
        processes.append(prc)
        prc.start()

    # completing process
    for prc in processes:
        prc.join()

    verboseprint("done")

    return {
        _coord_index[axes[0], 0]: points[:, 0],
        _coord_index[axes[1], 0]: points[:, 1],
        "axes_indexes": np.array(axes),
        "axes_names": np.array(_coord_index[axes, 0]),
        "nu_axes_indexes": np.array(nuaxes),
        "nu_axes_names": np.array(_coord_index[axes, 0]),
        "nturns": np.array(nturns),
        "nu_" + _coord_index[axes[0], 0]: np.array(x1freq),
        "nu_" + _coord_index[axes[1], 0]: np.array(x2freq),
        "dnu_" + _coord_index[axes[0], 0]: np.array(x1diff),
        "dnu_" + _coord_index[axes[1], 0]: np.array(x2diff),
        "dnulog10": np.array(nudiff),
        "turns_per_particle": np.array(turns_per_particle),
        "elem_particle_lost": np.array(elem_particle_lost),
        "worker_id": np.array(worker_id),
    }
