"""
Utility functions for tracking simulations.
"""

from __future__ import annotations
import numpy as np
import functools
from collections.abc import Sequence, Iterable
from ..lattice import Lattice, Element
from ..lattice import BeamMoments, Collective, QuantumDiffusion
from ..lattice import SimpleQuantDiff, VariableThinMultipole
from ..lattice import elements, refpts_iterator, set_value_refpts
from ..lattice import DConstant, checktype, checkattr, get_bool_index


__all__ = [
    "disable_varelem",
    "format_results",
    "fortran_align",
    "get_bunches",
    "get_bunches_std_mean",
    "has_collective",
    "initialize_lpass",
    "unfold_beam",
    "variable_refs",
]


DIMENSION_ERROR = "Input to lattice_pass() must be a 6xN array."
_COLLECTIVE_ELEMS = (BeamMoments, Collective)
_VAR_ELEMS = (QuantumDiffusion, SimpleQuantDiff, VariableThinMultipole)
_DISABLE_ELEMS = _COLLECTIVE_ELEMS + _VAR_ELEMS


def _set_beam_monitors(ring: Sequence[Element], nbunch: int, nturns: int):
    """Function to initialize the beam monitors."""
    monitors = list(refpts_iterator(ring, elements.BeamMoments))
    monitors += list(refpts_iterator(ring, elements.SliceMoments))
    for m in monitors:
        m.set_buffers(nturns, nbunch)
    return len(monitors) == 0


def variable_refs(ring):
    idref = get_bool_index(ring, checkattr("PassMethod", "IdentityPass"))
    varref = get_bool_index(ring, checktype(_DISABLE_ELEMS))
    varref = varref & ~idref
    return varref


def has_collective(ring) -> bool:
    """True if any element involves collective effects."""
    refs = get_bool_index(ring, checktype(_COLLECTIVE_ELEMS))
    return sum(refs) > 0


def disable_varelem(ring):
    """Function to disable collective effects elements."""
    refs = variable_refs(ring)
    if sum(refs) > 0:
        ring = set_value_refpts(ring, refs, "PassMethod", "IdentityPass", copy=True)
    return ring


def initialize_lpass(lattice: Iterable[Element], nturns: int, kwargs) -> list[Element]:
    """Function to initialize keyword arguments for lattice tracking."""
    if not isinstance(lattice, list):
        lattice = list(lattice)
    nbunch = getattr(lattice, "nbunch", 1)
    if not kwargs.pop("unfold_beam", True):
        kwargs["bunch_spos"] = np.zeros(nbunch)
    no_bm = _set_beam_monitors(lattice, nbunch, nturns)
    kwargs["keep_lattice"] = kwargs.get("keep_lattice", False) and no_bm
    if not kwargs.get("use_mp", False):
        kwargs.pop("pool_size", None)
        kwargs.pop("start_method", None)
    return lattice


def fortran_align(func):
    # noinspection PyShadowingNames
    """decorator to ensure that *r_in* is Fortran-aligned.

    :py:func:`fortran_align` ensures that the 2nd argument (usually *r_in*) of
    the decorated function is Fortran-aligned before calling the function

    Example:

        >>> @fortran_align
        ... def element_pass(element: Element, r_in, **kwargs):
        ... ...

        Ensure that *r_in* is Fortran-aligned
    """

    @functools.wraps(func)
    def wrapper(lattice, r_in, *args, **kwargs):
        assert r_in.shape[0] == 6 and r_in.ndim in (1, 2), DIMENSION_ERROR
        if r_in.flags.f_contiguous:
            return func(lattice, r_in, *args, **kwargs)
        else:
            r_fin = np.asfortranarray(r_in)
            r_out = func(lattice, r_fin, *args, **kwargs)
            r_in[:] = r_fin[:]
            return r_out

    return wrapper


def format_results(results, r_in, losses):
    """Function to format the ouput of parallelized tracking."""
    lin, lout = zip(*results, strict=True)
    # Update r_in with values at the end of tracking
    np.concatenate(lin, out=r_in, axis=1)
    if losses:
        lout, ldic = zip(*lout, strict=True)
        keys = ldic[0].keys()
        dicout = {k: np.hstack([li[k] for li in ldic]) for k in keys}
        return np.concatenate(lout, axis=1), dicout
    else:
        return np.concatenate(lout, axis=1)


def get_bunches(
    r_in: np.ndarray, nbunch: int, selected_bunches: np.ndarray | None = None
):
    """Function to get the bunches particles 6D coordinates.

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
        selected_bunches = np.arange(nbunch)
    else:
        selected_bunches = np.atleast_1d(selected_bunches)
    bunches = [r_in.T[ib::nbunch].T for ib in selected_bunches]
    return bunches


def get_bunches_std_mean(
    r_in: np.ndarray, nbunch: int, selected_bunches: np.ndarray | None = None
):
    """Function to get the bunches standard deviation and center of mass.

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
    std = [np.nanstd(b, axis=1) for b in bunches]
    mean = [np.nanmean(b, axis=1) for b in bunches]
    return std, mean


def unfold_beam(ring: Lattice, beam: np.ndarray, **kwargs) -> np.ndarray:
    """Function to unfold the beam based on the ring fill pattern.
    The input particle distribution has to be in on bucket 0.
    Particle are distributed in bunches using ``i%ring.nbunch``
    where i is the particle index.
    For each bunches the absolute ``ct`` is computed using the 6D
    closed orbit search, this closed orbit is added to the input
    particles.

    Parameters:
        ring: Lattice description
        beam: array with shape(6, nparticles)

    Keyword Arguments:
        convergence (float):    Convergence criterion for 6D orbit
          Default: :py:data:`DConstant.OrbConvergence <.DConstant>`
        max_iterations (int):   Maximum number of iterations for
          6D orbit.
          Default: :py:data:`DConstant.OrbMaxIter <.DConstant>`
    Return:
        beam (numpy.ndarray): unfolded beam
    """
    conv = kwargs.pop("convergence", DConstant.OrbConvergence)
    maxiter = kwargs.pop("max_iterations", DConstant.OrbMaxIter)
    o60, _ = ring.find_orbit(max_iterations=maxiter, convergence=conv)
    unfolded_beam = np.zeros(beam.shape)
    for i, spos in enumerate(ring.bunch_spos[::-1]):
        guess = o60.copy()
        guess[5] -= spos
        o6, _ = ring.find_orbit(guess=guess, max_iterations=maxiter, convergence=conv)
        unfolded_beam[:, i :: ring.nbunch] = (beam[:, i :: ring.nbunch].T + o6).T
    return unfolded_beam
