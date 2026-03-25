"""This is an implementation of the Flood Fill algorithm for pyat."""

from multiprocessing import Manager, Process, Queue

import at
import numpy as np

# Author : E. Serra,  UAB and ALBA,  2025 original version in python
#                                    See. IPAC2025, MOPB065.
#                                    APPLICATION OF FAST ALGORITHMS TO
#                                    CALCULATE DYNAMIC AND
#                                    MOMENTUM APERTURE TO THE DESIGN OF
#                                    ALBA II.
#                                    JACoW-IPAC25-MOPB065
# Edited : O. Blanco, ALBA,          2025 pyat parallel version

__all__ = ["floodfill"]


def track_queue(
    ring: at.Lattice,
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
            id_move = np.array([task_id + 1, task_id - 1, task_id + n_y, task_id - n_y])
            # use only valid index
            mask = np.logical_and(id_move >= 0, id_move < (n_x * n_y))
            newid = id_move[mask]
            for i in newid:
                if not registered_for_tracking[i]:
                    registered_for_tracking[i] = True
                    task_to_accomplish.put(i)


def floodfill(ring: at.Lattice, **kwargs: dict[str, any]) -> tuple:
    """Find the 2D acceptance of the ring using Flood Fill.

    Flood fill tracks particles from the exterior to the border of the
    acceptance.
    The lost particles are returned in plost.
    The not lost particles are returned in pnotlost.

    Parameters:
        ring: pyat lattice.
        kwargs: see below.

    kwargs:
        nturns: Number of turns for the tracking. Default: 1000
        window: Min and max coordinate range, [Axis1min,Axis1max,
        Axis2min,Axis2max]. Default [-10e-3,10e-3,-10e-3,10e-3].
        Axis1 and Axis2 are defined by 'axes'.
        grid_size: Number of steps per axis. Default [10,10].
        axes: Indexes of axes to be scanned. Default is [0,2], i.e. x-y.
        six_doffset: Offset to be added. Default np.zeros((6,1)).
        This is useful to study off-axis acceptance on any plane,
        or off-momentum acceptance by adding dp to the 5th coord.,
        to track particles on the closed orbit, or to add
        a small deviation to the tracked coordinates,
        e.g. [10e-5 10e-5] in the transverse planes.
        verbose:    Print extra info. Default 0.
        use_mp:     Parallel tracking with queue. Only CPU, not GPU.
        pool_size:  Number of cpus.

    Returns:
        Tuple with (pnotlost, plost).
        pnotlost: array of size (2,n_not_lost) containing the 2D offsets of
        n_not_lost tracked particles that survived.
        plost: array of size (3,n_lost) containing the 2D offsets of
        n_lost trackedp articles that did not survive; and
        in the third row the turn on which each particle is lost.

    Example:
        >>> pnl, pl = floodfill(ring, nturns=500)

    Notes:
    This method is recomended for single or low number of CPUs, and,
    it does not scale well for parallel computing.

    Based on the article,
      B. Riemann, M. Aiba, J. Kallestrup, and A. Streun, "Efficient
      algorithms for dynamic aperture and momentum acceptance
      calculation in synchrotron light sources", Phys. Rev. Accel.
      Beams, vol. 27, no. 9, p. 094 002, 2024.
      doi:10.1103/PhysRevAccelBeams.27.094002
    """
    # initialize variables
    nturns = kwargs.pop("nturns", 1000)
    window = kwargs.pop("window", (-10e-3, 10e-3, -5e-3, 5e-3))
    grid_size = kwargs.pop("grid_size", (10, 10))
    axes = kwargs.pop("axes", (0, 2))
    sixd_offset = kwargs.pop("sixd_offset", np.zeros((6)))
    verbose = kwargs.pop("verbose", False)
    pool_size = kwargs.pop("pool_size", 10)
    use_mp = kwargs.pop("use_mp", True)

    # Initialize output in case we return earlier
    points_not_lost = np.zeros((2, 0))
    points_lost_turns = np.zeros((3, 0))

    verboseprint = print if verbose else lambda *a, **k: None

    verboseprint("Flood fill starts.")
    verboseprint(f"Maximum number of turns: {nturns}")

    if window[0] == window[1] or window[2] == window[3]:
        verboseprint("Window is too narrow")
        return points_not_lost, points_lost_turns

    # set the grid

    if np.any(np.asarray(grid_size) < 1):
        print("Horizontal and vertical grid size should be more than 1")
        return points_not_lost, points_lost_turns

    n_x, n_y = grid_size
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
    particles = particles + sixd_offset
    particles[:, axes[0]] = particles[:, axes[0]] + points[0, :]
    particles[:, axes[1]] = particles[:, axes[1]] + points[1, :]

    # parallel parameters
    if use_mp:
        nproc = pool_size
    else:
        nproc = 1
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
    mask_tracked = np.array(final_turn) != -1
    mask_islost = np.array(islost)

    points_not_lost = points[:, np.logical_and(~mask_islost, mask_tracked)]
    points_lost = points[:, np.logical_and(mask_islost, mask_tracked)]
    ft_ = np.array(final_turn)
    turns_lost = ft_[np.logical_and(mask_islost, mask_tracked)]
    points_lost_turns = np.vstack((points_lost, turns_lost))
    return points_not_lost, points_lost_turns
