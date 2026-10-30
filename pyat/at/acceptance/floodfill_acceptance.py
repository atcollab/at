"""This is an implementation of the Flood Fill algorithm for pyat."""

from multiprocessing import Manager, Process, Queue

import numpy as np
from at.lattice import Lattice, axisdef
from at.tracking.track import MPMode

# Author : E. Serra,  UAB and ALBA,  2025 original version in python
#                                    See. IPAC2025, MOPB065.
#                                    APPLICATION OF FAST ALGORITHMS TO
#                                    CALCULATE DYNAMIC AND
#                                    MOMENTUM APERTURE TO THE DESIGN OF
#                                    ALBA II.
#                                    JACoW-IPAC25-MOPB065
# Edited : O. Blanco, ALBA,          2025 pyat parallel version

__all__ = ["floodfill", "floodfill_1d"]


def floodfill(
    ring: Lattice,
    nturns: int = 1024,
    planes: list | tuple = ("x", "y"),
    amplitudes: list | tuple = (10e-3, 10e-3),
    bounds: list | tuple = ((-1, 1), (0, 1)),
    npoints: list | tuple = (10, 10),
    offset: list | np.ndarray | None = None,
    verbose: bool = False,
    use_mp: bool | MPMode = True,
    pool_size: int = 10,
) -> np.ndarray:
    """Find the 2D acceptance of the ring using Flood Fill [1]_.

    Flood fill [1]_ tracks particles from the exterior to the border of the
    acceptance.

    Parameters:
        ring: pyat lattice
        nturns: Number of turns for the tracking. Default: 1024
        planes: max. dimension 2, Plane(s) to scan for the acceptance.
          Allowed values are: ``'x'``, ``'px'``, ``'y'``,
          ``'py'``, ``'dp'``, ``'ct'``
        amplitudes: (2,) array, set the search range per plane.
          Default [10e-3,10e-3]
        bounds: (2,2) array, defines the tracked range: range=bounds*amplitude
            Default ((-1,1),(0,1))
        npoints: Number of steps per axis. Default (10,10).
        offset: Offset to be added. Default np.zeros((6)).
          This is useful to study off-axis acceptance on any plane,
          or off-momentum acceptance by adding dp to the 5th coord.,
          to track particles on the closed orbit, or to add
          a small deviation to the tracked coordinates,
          e.g. [10e-5 10e-5] in the transverse planes.
        verbose: Print extra info. Default 0.
        use_mp: Parallel tracking, default True. It uses at.MPMode.CPU,
          not at.MPMode.GPU.
        pool_size:  Number of cpus.

    Returns:
        A (4,n) array with the following rows:
          Position on axis 1 of the n tracked particles
          Position on axis 2 of the n tracked particles
          Number of turns the n particle did during tracking
          Flag set to 1 if the nth particle is lost, otherwise 0

    Example:
        >>> ff_data = floodfill(ring, nturns=500)
        >>> alive_mask = ff_data[3,:] == 0

    .. Note::

       This method is recommended for single or low number of CPUs, and does
       not scale well for parallel computing.

    References:
      .. [1] B. Riemann, M. Aiba, J. Kallestrup, and A. Streun, "Efficient
         algorithms for dynamic aperture and momentum acceptance
         calculation in synchrotron light sources",
         Phys. Rev. Accel. Beams, vol. 27, no. 9, p. 094 002, 2024.
         doi:10.1103/PhysRevAccelBeams.27.094002

    """

    def track_queue(
        ring: Lattice,
        zin: np.ndarray,
        nturns: int,
        n_x: int,
        n_y: int,
        task_to_accomplish: list,
        islost: list,
        final_turn: list,
        registered_for_tracking: list,
    ) -> None:
        """Track particles with index inside queue.

        Parameters:
            ring: pyat ring
            zin: particles coordinates (6,nparticles)
            nturns: number of tracking turns
            n_x: horizontal grid size
            n_y: vertical grid size
            task_to_accomplish: list with indexes of particles to track
            islost: list with particle loss result from tracking
            final_turn: turn when the particle is lost
            registered_for_tracking: particle tracked or to be tracked
        """
        while not task_to_accomplish.empty():
            task_id = task_to_accomplish.get(block=True, timeout=10)
            registered_for_tracking[task_id] = True
            *_, losses_data = ring.track(zin[:, task_id], nturns, losses=True)
            islost[task_id] = losses_data["loss_map"]["islost"][0]
            final_turn[task_id] = losses_data["loss_map"]["turn"][0]
            if islost[task_id]:
                # next particles
                id_move = np.array(
                    [task_id + 1, task_id - 1, task_id + n_y, task_id - n_y]
                )
                # use only valid index
                mask = np.logical_and(id_move >= 0, id_move < (n_x * n_y))
                newid = id_move[mask]
                for i in newid:
                    if not registered_for_tracking[i]:
                        registered_for_tracking[i] = True
                        task_to_accomplish.put(i)

    # Initialize variables
    if offset is None:
        offset = 6 * [0]
    offset = np.array(offset)
    planesi = np.atleast_1d(axisdef.axis_(tuple(planes), key="index"))

    # Initialize output in case we return earlier
    data_tracked = np.zeros((4, 0))

    verboseprint = print if verbose else lambda *a, **k: None

    verboseprint("Flood fill starts.")
    verboseprint(f"Maximum number of turns: {nturns}")

    window = np.ravel((np.array(bounds).reshape((2, 2)).T * np.array(amplitudes)).T)
    if window[0] == window[1] or window[2] == window[3]:
        verboseprint("Window is too narrow")
        return data_tracked

    # set the grid

    if np.any(np.asarray(npoints) < 1):
        print("Horizontal and vertical grid size should be more than 1")
        return data_tracked

    n_x, n_y = npoints
    min_x = min(window[0:2])
    max_x = max(window[0:2])
    min_y = min(window[2:4])
    max_y = max(window[2:4])
    xvals = np.linspace(min_x, max_x, n_x)
    yvals = np.linspace(min_y, max_y, n_y)

    nparticles = n_x * n_y
    points = np.zeros((2, nparticles))

    verboseprint(f"Number of grid points: {nparticles}")

    # set the order
    ii_ = 0
    for ix_ in range(n_x):
        for iy_ in range(n_y):
            points[0, ii_] = xvals[ix_]
            points[1, ii_] = yvals[iy_]
            ii_ = ii_ + 1

    ndims = 6
    particles = np.zeros((nparticles, ndims))
    particles = particles + offset
    particles[:, planesi[0]] = particles[:, planesi[0]] + points[0, :]
    particles[:, planesi[1]] = particles[:, planesi[1]] + points[1, :]

    # parallel parameters
    nproc = pool_size if use_mp else 1
    verboseprint(f"Number of processors: {nproc}")
    task_to_accomplish = Queue()
    processes = []

    leftside = np.arange(1, n_y - 1)
    rightside = np.arange((n_x - 1) * n_y + 1, n_x * n_y - 1)
    lowerline = np.asarray([*range(0, n_y * n_x, n_y)])
    upperline = lowerline + n_y - 1

    manager = Manager()
    registered_for_tracking = manager.list(nparticles * [False])
    for i in np.concatenate((upperline, leftside, lowerline, rightside)):
        task_to_accomplish.put(i)
        registered_for_tracking[i] = True

    verboseprint("Tracking...")

    # creating processes
    islost = manager.list(nparticles * [False])
    final_turn = manager.list(nparticles * [-1])
    for _worker in range(nproc):
        prc = Process(
            target=track_queue,
            args=(
                ring,
                particles.T,
                nturns,
                n_x,
                n_y,
                task_to_accomplish,
                islost,
                final_turn,
                registered_for_tracking,
            ),
        )
        processes.append(prc)
        prc.start()

    # completing process
    for prc in processes:
        prc.join()

    verboseprint("done")
    ft_ = np.array(final_turn)
    mask_tracked = ft_ != -1
    mask_islost = np.array(islost)
    turns_per_point = np.zeros(sum(mask_tracked))

    points_tracked = points[:, mask_tracked]
    lost_after_track = mask_islost[mask_tracked]
    turns_per_point[lost_after_track] = ft_[mask_islost]
    turns_per_point[~lost_after_track] = nturns
    return np.vstack((points_tracked, turns_per_point, lost_after_track))


def floodfill_1d(
    ring: Lattice,
    nturns: int = 1024,
    plane: list | tuple = ("dp"),
    amplitude: list | tuple = (10e-2,),
    bounds: list | tuple = ((-1, 1),),
    npoints: list | tuple = (10,),
    offset: list | np.ndarray | None = None,
    obspt: list | np.ndarray | None = None,
    verbose: bool = False,
    use_mp: bool | MPMode = True,
    pool_size: int = 10,
    **kwargs,
) -> np.ndarray:
    """Find the momentum acceptance of the ring using Flood Fill.

    Flood fill tracks particles from the exterior to the border of the
    acceptance.

    Parameters:
        ring: pyat lattice
        nturns: Number of turns for the tracking. Default: 1024
        plane: max. dimension 1, Plane to scan for the acceptance.
          Allowed values are: ``'x'``, ``'px'``, ``'y'``,
          ``'py'``, ``'dp'``, ``'ct'``
        amplitude: (2,) array, set the search range.
          Default [10e-3]
        bounds: (2,2) array, defines the tracked range: range=bounds*amplitude
            Default ((-1,1),(0,1))
        npoints: Number of steps per axis. Default (10,10).
        offset: Offset to be added. Default is the closed orbit per obspt.
          This is useful to study to add a small deviation in the transverse plane
          to the tracked coordinates
          e.g. 10e-5 in the transverse planes.
        obspt: reference points (default len(ring))
        verbose: Print extra info. Default 0.
        use_mp: Parallel tracking, default True. It uses at.MPMode.CPU,
          not at.MPMode.GPU.
        pool_size:  Number of cpus.

    Returns:
        A (5,n) array with the following rows:
          1) Observation points.
          2) energy of the tracked points
          3) `s` postions of the observation points.
          4) Number of turns the n particle did during tracking.
          5) Flag set to 1 if the nth particle is lost, otherwise 0.

    Example:
        >>> data_ff1d = floodfill_1d(ring, plane=["dp"], nturns=1024)
        >>> alive_mask = data_ff1d[4,:] == 0

        >>> particle_turns = data_ff1d[3,:]
        >>> particle_refpts = data_ff1d[0,:]
        >>> particle_spos = data_ff1d[2,:]
        >>> particle_offset = data_ff1d[1,:]
        >>> plt.scatter(particle_spos,particle_offset,c=particle_turns)
        >>> plt.show()


    .. Note::

       This method is recommended for single or low number of CPUs, and does
       not scale well for parallel computing.
    """

    def track_queue(
        ring: Lattice,
        points: np.ndarray,
        zin: np.ndarray,
        nturns: int,
        n_x: int,
        n_y: int,
        task_to_accomplish: list,
        islost: list,
        final_turn: list,
        registered_for_tracking: list,
    ) -> None:
        """Track particles with index inside queue.

        Parameters:
            ring: pyat ring
            points : array containing the pairs (obspt, energy offset) per
            input particle
            zin: particles coordinates (6,nparticles)
            nturns: number of tracking turns
            n_x: horizontal grid size
            n_y: vertical grid size
            task_to_accomplish: list with indexes of particles to track
            islost: list with particle loss result from tracking
            final_turn: turn when the particle is lost
            registered_for_tracking: particle tracked or to be tracked
        """
        while not task_to_accomplish.empty():
            task_id = task_to_accomplish.get(block=True, timeout=10)
            registered_for_tracking[task_id] = True
            newring = ring.rotate(int(points[0, task_id]))
            *_, losses_data = newring.track(zin[:, task_id], nturns, losses=True)
            islost[task_id] = losses_data["loss_map"]["islost"][0]
            final_turn[task_id] = losses_data["loss_map"]["turn"][0]
            if islost[task_id]:
                # next particles
                id_move = np.array(
                    [task_id + 1, task_id - 1, task_id + n_y, task_id - n_y]
                )
                # use only valid index
                mask = np.logical_and(id_move >= 0, id_move < (n_x * n_y))
                newid = id_move[mask]
                for i in newid:
                    if not registered_for_tracking[i]:
                        registered_for_tracking[i] = True
                        task_to_accomplish.put(i)

    # Initialize output in case we return earlier
    data_tracked = np.zeros((5, 0))

    # obspt
    if obspt is None:
        obspt = np.arange(len(ring))
    obspt = np.asarray(obspt)
    s_obspt = ring.get_s_pos(refpts=obspt)
    s_unique, idx_unique = np.unique(s_obspt, return_index=True)
    obspt_unique = obspt[idx_unique]

    verboseprint = print if verbose else lambda *a, **k: None

    verboseprint("Flood fill starts.")
    verboseprint(f"Maximum number of turns: {nturns}")

    window = np.ravel((np.array(bounds).reshape((1, 2)).T * np.array(amplitude)).T)
    if window[0] == window[1]:
        verboseprint("Window is too narrow")
        return data_tracked

    # set the grid
    n_y = npoints[0]
    min_y = min(window[0:2])
    max_y = max(window[0:2])
    xvals = obspt_unique
    n_x = len(obspt_unique)
    yvals = np.linspace(min_y, max_y, n_y)

    nparticles = n_x * n_y
    points = np.zeros((2, nparticles))

    verboseprint(f"Number of grid points: {nparticles}")

    # include option dp
    dp = kwargs.get("dp", 0.0)
    verboseprint(f"Using dp {dp:.6g}")

    # set the order
    ii_ = 0
    for ix_ in range(n_x):
        for iy_ in range(n_y):
            points[0, ii_] = int(xvals[ix_])
            points[1, ii_] = yvals[iy_]
            ii_ = ii_ + 1

    if offset is None:
        offset = [None]
    if None in offset:
        _, co = ring.find_orbit(obspt_unique)
        offset_repeat = np.repeat(co, n_y, axis=0)
    else:
        offset = offset[obspt_unique]
        offset_repeat = np.repeat(offset, n_y, axis=0)

    planesi = np.atleast_1d(axisdef.axis_(tuple(plane), key="index"))

    ndims = 6
    particles = offset_repeat
    particles[4, :] = particles[4, :] + dp
    particles[:, planesi[0]] = particles[:, planesi[0]] + points[1, :]

    # parallel parameters
    nproc = pool_size if use_mp else 1
    verboseprint(f"Number of processors: {nproc}")
    task_to_accomplish = Queue()
    processes = []

    leftside = np.arange(1, n_y - 1)
    rightside = np.arange((n_x - 1) * n_y + 1, n_x * n_y - 1)
    lowerline = np.asarray([*range(0, n_y * n_x, n_y)])
    upperline = lowerline + n_y - 1

    manager = Manager()
    registered_for_tracking = manager.list(nparticles * [False])
    for i in np.concatenate((upperline, leftside, lowerline, rightside)):
        task_to_accomplish.put(i)
        registered_for_tracking[i] = True

    verboseprint("Tracking...")

    # creating processes
    islost = manager.list(nparticles * [False])
    final_turn = manager.list(nparticles * [-1])
    for _worker in range(nproc):
        prc = Process(
            target=track_queue,
            args=(
                ring,
                points,
                particles.T,
                nturns,
                n_x,
                n_y,
                task_to_accomplish,
                islost,
                final_turn,
                registered_for_tracking,
            ),
        )
        processes.append(prc)
        prc.start()

    # completing process
    for prc in processes:
        prc.join()

    verboseprint("done")
    ft_ = np.array(final_turn)
    mask_tracked = ft_ != -1
    mask_islost = np.array(islost)
    turns_per_point = np.zeros(sum(mask_tracked))

    points_tracked = points[:, mask_tracked]
    s_unique_rep = np.repeat(s_unique, n_y, axis=0)
    s_tracked = s_unique_rep[mask_tracked]
    lost_after_track = mask_islost[mask_tracked]
    turns_per_point[lost_after_track] = ft_[mask_islost]
    turns_per_point[~lost_after_track] = nturns
    return np.vstack((points_tracked, s_tracked, turns_per_point, lost_after_track))
