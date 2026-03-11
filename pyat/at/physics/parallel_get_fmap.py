"""This script calculates the frequency and diffussion maps."""

# oblanco 2026jan06

from __future__ import annotations
from collections import defaultdict

import numpy as np

try:
    # noinspection PyPackageRequirements
    import psutil
    psutil_found = True
except (ImportError, RuntimeError):
    msg = 'psutil is unavailable => memory management is disabled.\n' + \
          'To enable memory management, run "pip install psutil"'
    print(msg)
    psutil_found = False
import os
import at
from time import sleep
from at.lattice import Lattice, Refpts
from at.tracking import MPMode, gpu_core_count
from multiprocessing import Manager, Process, Queue, cpu_count, current_process, Value

__all__ = ["get_fmap"]

_coord_index = np.array([("x",0), ("px", 1), ("y", 2), ("py", 3), ("dp", 4), ("ct", 5)])

def get_mem_info(pid: int) -> dict[str, int]:
  res = defaultdict(int)
  mmap = psutil.Process(pid).memory_full_info()
  res['pss'] = mmap.pss
  return res


def get_ram_info():
    ram = psutil.virtual_memory()
    total_ram = ram.total # total installed RAM
    available_ram = ram.available # RAM available for processes
    used_ram = ram.used # RAM used by processes
    free_ram = ram.free # RAM not being used
    percent_used = ram.percent # percentage of RAM used
    return total_ram, available_ram, used_ram, free_ram, percent_used

def get_swap_info():
    swap = psutil.swap_memory()
    total_swap = swap.total # total installed RAM
    used_swap = swap.used # RAM used by processes
    free_swap = swap.free # RAM not being used
    percent_used = swap.percent # percentage of RAM used

    return total_swap, used_swap, free_swap, percent_used

def estimate_memory_resources(**kwargs):
    # estimate memory resources
    verbose = kwargs.pop("verbose", False)
    mem_margin = kwargs.pop("mem_margin", np.nan)*1024*1024
    ramtot,_,_,ramfree,_ = get_ram_info()
    swaptot,_,swapfree,_ = get_swap_info()
    onehundredmegabytes = 100*1024*1024
    safe_margin = min(onehundredmegabytes,1/32*ramtot,1/4*swaptot)
    if not np.isnan(mem_margin):
        safe_margin = mem_margin
    estimated_freemem = ramfree + swapfree
    if verbose:
        msg = f'Estimated free memory : {estimated_freemem/1024/1024:.3f} MB, ' + \
              f'Memory safe margin : {safe_margin/1024/1024:.3f} MB'
        print(msg)
    return estimated_freemem, safe_margin

def grid_configuration(planes,axes,window,shift_zero=1.0e-6):
    pass


def track_queue(
    ring: at.Lattice,
    zin: np.ndarray,
    nparticles: int,
    nturns: int,
    n_moving_slices: int,
    task_to_accomplish: list,
    nuaxes: tuple,
    worker_id: list,
    needed_mem: float,
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

    Parameters:
        ring: pyat ring
        zin: particles coordinates (6,nparticles)
        nturns: number of tracking turns
        task_to_accomplish: list with indexes of particles to track
    """
    while not task_to_accomplish.empty():

        estimated_freemem = np.inf
        safe_margin = 0
        if psutil_found:
            estimated_freemem, safe_margin = estimate_memory_resources(verbose=verbose,mem_margin=mem_margin)
            if not np.isnan(mem_margin):
                safe_margin = mem_margin
        while estimated_freemem < safe_margin:
            if verbose:
                print('waiting for free memory ...')
            sleep(5)
        while used_mem.value > mem_limit:
            if verbose:
                print('using too much memory, waiting ...')
            sleep(5)
        task_id = task_to_accomplish.get(block=True, timeout=10)
        current = current_process()
        worker_id[task_id] = current.pid
        # we track twice the number of nturns in order to get two frequencies
        zout, _, loss_info = ring.track(zin[:, task_id], 2*nturns, losses=True, **kwargs)
        msg_mem_info = ''
        if psutil_found:
            monitor_data = get_mem_info(current.pid)['pss']
            with used_mem.get_lock():
                used_mem.value = used_mem.value + monitor_data
                msg_mem_info = f'Used mem: {used_mem.value/1024/1024:.3f} MB, '
        if verbose:
           msg = msg_mem_info + f'cpu id {current.pid}, ' + \
                 f'running particle id {task_id} of {nparticles}'
           print(msg)
        # if the particle is lost
        if loss_info['loss_map']['islost'][0]:
            print('particle is lost')
            turns_per_particle[task_id] = loss_info['loss_map']['turn'][0]
            elem_particle_lost[task_id] = loss_info['loss_map']['elem'][0]
        else:
        # get the tunes
            zout = zout.squeeze()
            tune_per_window = np.empty((2,n_moving_slices))
            tune_per_window[:] = np.nan
            shift_n = np.floor_divide(2*nturns,n_moving_slices)
            for i in range(n_moving_slices):
                lim_inf = i*shift_n
                lim_sup = (i+1)*shift_n
                z1 = zout[nuaxes[0],lim_inf:lim_sup]
                z2 = zout[nuaxes[1],lim_inf:lim_sup]
                z1 = z1 - z1.mean()
                z2 = z2 - z2.mean()
                tune1 = at.harmonic_analysis.get_tunes_harmonic(z1)
                tune2 = at.harmonic_analysis.get_tunes_harmonic(z2)
                if len(tune1) == 0 or len(tune2) == 0:
                    continue
                tune_per_window[0,i] = tune1[0]
                tune_per_window[1,i] = tune2[0]
            # metric
            dnu1,dnu2 = np.nanstd(tune_per_window,axis=1)
            dnu_norm = 0.5*np.log10(dnu1*dnu1 + dnu2*dnu2)
            # min max diff
            if dnu_norm > -2:
                dnu_norm = -2
            if dnu_norm < -10:
                dnu_norm = -10
            #update
            turns_per_particle[task_id] = 2*nturns
            x1freq[task_id] = np.nanmean(tune_per_window[0,:])
            x2freq[task_id] = np.nanmean(tune_per_window[1,:])
            x2freq[task_id] = tune2
            x1diff[task_id] = dnu1
            x2diff[task_id] = dnu2
            nudiff[task_id] = dnu_norm
        if psutil_found:
            with used_mem.get_lock():
                used_mem.value = used_mem.value - monitor_data


def check_parallel_resources(use_mp,pool_size,gpu_pool,kwargs,verbose):
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
            msg = f"{nprocu} of {nproc} cpu used frequency map calculation\n" + \
                   "Multi-process acceptance calculation selected..."
            print(msg)
        if nproc == 1:
            if verbose:
                print("Consider use_mp=False for single core computations")
    elif use_mp is MPMode.GPU:
        print(type(gpu_pool))
        print(gpu_pool)
        nprocu = gpu_core_count(gpu_pool)
        kwargs["gpu_pool"] = gpu_pool if gpu_pool is not None else [0]
        if verbose:
            msg = f"\n{nprocu} GPU cores found" + \
                  "GPU acceptance calculation selected..."
            print(msg)
    else:
        nprocu = 1
        if verbose:
            print("Single process acceptance calculation selected...")
        if nproc > 1 and verbose:
            print("Consider use_mp=True for parallelized computations")

    return nprocu,kwargs


def generate_grid(grid_size,window,offset,shift_zero,axes):

    if window[0] == window[1] or window[2] == window[3]:
        verboseprint("Window is too narrow")

    n_x1, n_x2 = grid_size
    min_x1 = min(window[0:2])
    max_x1 = max(window[0:2])
    min_x2 = min(window[2:4])
    max_x2 = max(window[2:4])
    x1vals = np.linspace(min_x1, max_x1, n_x1)
    x2vals = np.linspace(min_x2, max_x2, n_x2)

    nparticles = n_x1 * n_x2
    points = np.zeros((nparticles,2))

    for i,x1val in enumerate(x1vals):
        for j,x2val in enumerate(x2vals):
            points[i*n_x2+j,0] = x1vals[i]
            points[i*n_x2+j,1] = x2vals[j]

    ndims = 6
    particles = np.zeros((ndims,nparticles))
    particles = particles + offset.reshape((6,1))
    particles[axes[0], :] = particles[axes[0], :] + points[:, 0] + shift_zero
    particles[axes[1], :] = particles[axes[1], :] + points[:, 1] + shift_zero

    return nparticles,ndims,points,particles

def adapt_particles(user_particles,axes):
    ndims, nparticles = user_particles.shape
    points = user_particles[axes,:]
    return nparticles,ndims,points,user_particles

def coord_to_indexes(useraxes):
    user_prefers_indexes = np.all([isinstance(i,int) for i in useraxes])
    if not user_prefers_indexes:
        axes = np.zeros((len(useraxes)),dtype=int)
        for i,ax in enumerate(useraxes):
            mask = _coord_index[:,0] == ax
            axes[i] = int(np.squeeze(_coord_index[mask,1][0]))
    else:
        axes = np.array(useraxes)
    return axes

def get_fmap(
        ring: Lattice,
        **kwargs):
    r"""Compute the frequency and diffussion maps of a ring.

    Parameters:
        ring:           Lattice definition
        planes:         max. dimension 2, Plane(s) to scan for the acceptance.
          Allowed values are: ``'x'``, ``'xp'``, ``'y'``,
          ``'yp'``, ``'dp'``, ``'ct'``
        npoints:        (len(planes),) array: number of points in each
          dimension
        amplitudes:     (len(planes),) array: set the search range:
        nturns:         Number of turns for the tracking
        refpts:         Observation points. Default: start of the machine
        dp:             static momentum offset
        offset:         initial orbit. Default: closed orbit
        use_mp:         Flag to activate CPU or GPU multiprocessing
          (default: False)
        gpu_pool:       List of GPU id to use when use_mp is
          :py:attr:`at.tracking.MPMode.GPU`. If None specified, if gets
          first GPU.
        verbose:        Print out some information
        shift_zero: Epsilon offset applied on all 6 coordinates
        start_method:   Python multiprocessing start method. The default
          ``None`` uses the python default that is considered safe.
          Available parameters: ``'fork'``, ``'spawn'``, ``'forkserver'``.
          The default for linux is ``'fork'``, the default for macOS and
          Windows is ``'spawn'``. ``'fork'`` may be used on macOS to speed up
          the calculation or to solve runtime errors, however it is
          considered unsafe.

    _pdict = {"x": 0, "xp": 1, "y": 2, "yp": 3, "dp": 4, "ct": 5}
    Returns:
        boundary:   (2,n) array: 2D acceptance
        survived:   (2,n) array: Coordinates of surviving particles
        tracked:    (2,n) array: Coordinates of tracked particles

    In case of multiple refpts, return values are lists of arrays, with one
    array per ref. point.

    Examples:

        >>> bf, sf, gf = ring.get_acceptance(planes, npoints, amplitudes)
        >>> plt.plot(*gf, ".")
        >>> plt.plot(*sf, ".")
        >>> plt.plot(*bf)
        >>> plt.show()

    .. note::

       * When``use_mp=True`` all the available CPUs will be used.
         This behavior can be changed by setting
         ``at.DConstant.patpass_poolsize`` to the desired value
       * When multiple ``refpts`` are provided particles are first
         projected to the beginning of the ring with tracking. Then,
         all particles are tracked up to ``nturns``. This allows to
         do most of the work in a single function call and allows for
         full parallelization.
    """
    #if start_method is not None:
    #    kwargs["start_method"] = start_method


    # initialize variables
    nturns = kwargs.pop("nturns", 1000)
    useraxes = kwargs.pop("axes", ('x','px'))
    usernuaxes = kwargs.pop("nuaxes", (0, 2))
    n_moving_slices = kwargs.pop('n_moving_slices',2)
    user_particles = kwargs.pop("particles", np.nan)
    window = kwargs.pop("window", (-10e-3, 10e-3, -10e-3, 10e-3))
    grid_size = kwargs.pop("grid_size", (10, 10))
    offset = kwargs.pop("offset", np.zeros((6)))
    shift_zero = kwargs.pop("shift_zero",1e-5)
    verbose = kwargs.pop("verbose", False)
    use_mp = kwargs.pop("use_mp", False)
    pool_size = kwargs.pop("pool_size", np.nan)
    gpu_pool = kwargs.pop('gpu_pool',[0])
    max_mem = kwargs.pop("max_mem", np.nan)
    mem_margin = kwargs.pop("mem_margin", np.nan*1024*1024)

    verboseprint = print if verbose else lambda *a, **k: None

    nprocu,kwargs = check_parallel_resources(use_mp,pool_size,gpu_pool,kwargs,verbose)

    axes = coord_to_indexes(useraxes)
    verboseprint(f'Using axes {_coord_index[axes,0]}')
    nuaxes = coord_to_indexes(usernuaxes)
    verboseprint(f'Using tune axes {_coord_index[nuaxes,0]}')

    if not np.any(np.isnan(user_particles)):
        nparticles,ndims,points,particles = adapt_particles(user_particles,axes)
    else:
        nparticles,ndims,points,particles = generate_grid(grid_size,window,offset,shift_zero,axes)


    verboseprint(f"Number of grid points: {nparticles}")

    task_to_accomplish = Queue()
    processes = []
    manager = Manager()
    for i in range(nparticles):
        task_to_accomplish.put(i)

    estimated_freemem = np.inf
    safe_margin = 0
    max_mem = estimated_freemem
    if psutil_found:
        estimated_freemem, safe_margin = estimate_memory_resources(verbose=verbose,mem_margin=mem_margin)
        if np.isnan(max_mem):
            if not np.isnan(mem_margin):
                print('User memory margin {mem_margin/1024/1024} MB')
                safe_margin = mem_margin
            max_mem = estimated_freemem - safe_margin
    mem_limit = max_mem*1024*1024
    verboseprint(f'Memory usage limit set to {mem_limit/1024/1024} MB')

    # estimate memory usage
    floatsize = 8 # bytes in python
    n_in_parallel = min(nparticles,nprocu)
    needed_mem = floatsize*ndims*n_in_parallel*nturns*2 # we track twice the nturns
    verboseprint(f"Estimated max. memory necessary for tracking: {needed_mem/1024/1024:.3} MB")

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
    used_memory = Value('i',0)
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
                needed_mem,
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

    return {_coord_index[axes[0],0]: points[:,0],
            _coord_index[axes[1],0]: points[:,1],
            'axes_indexes': np.array(axes),
            'axes_names': np.array(_coord_index[axes,0]),
            'nu_axes_indexes': np.array(nuaxes),
            'nu_axes_names': np.array(_coord_index[axes,0]),
            'nturns': np.array(nturns),
            'nu_'+_coord_index[axes[0],0]: np.array(x1freq),
            'nu_'+_coord_index[axes[1],0]: np.array(x2freq),
            'dnu_'+_coord_index[axes[0],0]: np.array(x1diff),
            'dnu_'+_coord_index[axes[1],0]: np.array(x2diff),
            'dnulog10': np.array(nudiff),
            'turns_per_particle': np.array(turns_per_particle),
            'elem_particle_lost': np.array(elem_particle_lost),
            'worker_id': np.array(worker_id),
            }

