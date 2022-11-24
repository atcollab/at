"""
Utility functions for tracking simulations
"""
import numpy
from at.lattice import Lattice, DConstant
from typing import Optional


__all__ = ['get_bunches', 'get_bunches_std_mean', 'unfold_beam']


def get_bunches(r_in: numpy.ndarray, nbunch: int,
                selected_bunches: Optional[numpy.ndarray] = None):
    """Function to get the bunches particles 6D coordinates

    Parameters:
      r_in: 6 x n_particles Fortran-ordered numpy array.
      nbunch: integer, total number of bunches
      selected_bunches: integer or array of integers, index
      of the selected bunches

    Returns:
      List of ndarray containing the 6 x n particles coordinates
      of the selected bunches
    """
    if selected_bunches is None:
        selected_bunches = numpy.arange(nbunch)
    else:
        selected_bunches = numpy.atleast_1d(selected_bunches)
    bunches = [r_in.T[ib::nbunch].T for ib in selected_bunches]
    return bunches


def get_bunches_std_mean(r_in: numpy.ndarray, nbunch: int,
                         selected_bunches: Optional[numpy.ndarray] = None):
    """Function to get the bunches standard deviation and center
    of mass

    Parameters:
      r_in: 6 x n_particles Fortran-ordered numpy array.
      nbunch: integer, total number of bunches
      selected_bunches: integer or array of integers, index
      of the selected bunches

    Returns:
      Lists of ndarray containing the 6D standard deviation
      and center of mass (std, mean)
    """
    bunches = get_bunches(r_in, nbunch, selected_bunches)
    std = [numpy.nanstd(b, axis=1) for b in bunches]
    mean = [numpy.nanmean(b, axis=1) for b in bunches]
    return std, mean


def unfold_beam(ring: Lattice, beam: numpy.ndarray, **kwargs):
    """Function to unfold the beam based on the ring fill pattern.
    The input particle distribution has to be in on bucket 0.
    Particle are distributed in bunches using ``i%ring.nbunch``
    where i is the particle index.
    For each bunches the absolute ``ct`` is computed using the 6D
    closed orbit search, this closed orbit is added to the input
    particles.
    The particle coordinates are modified in-place.

    Parameters:
        ring: Lattice description
        beam: array with shape(6, nparticles)

    Keyword Arguments:
        convergence (float):    Convergence criterion for 6D orbit
          Default: :py:data:`DConstant.OrbConvergence <.DConstant>`
        max_iterations (int):   Maximum number of iterations for
          6D orbit.
          Default: :py:data:`DConstant.OrbMaxIter <.DConstant>`
    """
    conv = kwargs.pop('convergence', DConstant.OrbConvergence)
    maxiter = kwargs.pop('max_iterations', DConstant.OrbMaxIter)
    o60, _ = ring.find_orbit(max_iterations=maxiter, convergence=conv)
    for i, spos in enumerate(ring.bunch_spos):
        guess = o60
        guess[5] += spos
        o6, _ = ring.find_orbit(guess=guess, max_iterations=maxiter,
                                convergence=conv)
        beam[:, i::ring.nbunch] = (beam[:, i::ring.nbunch].T + o6).T
