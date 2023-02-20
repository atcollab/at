from __future__ import annotations
from ..lattice import Lattice
from ..latticetools import SvdResponse, ResponseMatrix
from typing import Optional, Union
from collections.abc import Sequence
import matplotlib.pyplot as plt
from matplotlib.axes import Axes


def plot_norm(resp: SvdResponse,
              ax: Optional[tuple[Axes, Axes]] = None) -> None:
    r"""Plot the norm of :py:class:`.Observable`\ s and
    :py:class:`.Variable`\ s of a response matrix

    Plot the norm of the lines of the weighted response matrix
    (:py:class:`.Observable`\ s) and of its columns
    (:py:class:`.Variable`\ s)

    For a stable solution, the norms should have the same order of magnitude.
    If not, the weights of observables and variables should be adjusted.

    Args:
        resp:           Response matrix object
        ax:             tuple of :py:class:`~.matplotlib.axes.Axes`. If given,
          plots will be drawn in these axes.
    """
    obs, var = resp.check_norm()
    if ax is None:
        fig, (ax1, ax2) = plt.subplots(nrows=2)
    else:
        ax1, ax2 = ax[:2]
    ax1.bar(range(len(obs)), obs)
    ax1.set_title("Norm of weighted observables")
    ax1.set_xlabel("Observable #")
    ax2.bar(range(len(var)), var)
    ax2.set_title("Norm of weighted variables")
    ax2.set_xlabel("Variable #")


def plot_singular_values(resp: SvdResponse, ax: Axes = None,
                         logscale: bool = True) -> None:
    r"""Plot the singular values of a response matrix

    Args:
        resp:           Response matrix object
        logscale:       If :py:obj:`True`, use log scale
        ax:             If given, plots will be drawn in these axes.
    """
    singvals = resp.singular_values
    if ax is None:
        fig, ax = plt.subplots()
    ax.bar(range(len(singvals)), singvals)
    if logscale:
        ax.set_yscale('log')
    ax.set_title("Singular values")


def plot_contents(resp: ResponseMatrix, error: Union[Lattice, Sequence[float]],
                  ax: Axes = None, logscale: bool = True) -> None:
    """Plot the decomposition of an error vector on the basis of singular
    vectors

    Args:
        resp:           Response matrix object
        error:          Error to be analysed. *error* may be:

          * a :py:class:`.Lattice`: The response matrix observables will be
            evaluated for this :py:class:`.Lattice` and the deviation from
            target will be decomposed on the basis of singular vectors,
          * a Sequence of float: It will be decomposed on the basis of singular
            vectors
        logscale:       If :py:obj:`True`, use log scale
        ax:             If given, plots will be drawn in these axes.
    """
    if resp.singular_values is None:
        resp.solve()
    if isinstance(error, Lattice):
        obs = resp.observables
        obs.evaluate(error, r_in=resp.r_in)
        error = obs.flat_deviations
    corr = resp.uh @ error
    if ax is None:
        fig, ax = plt.subplots()
    ax.bar(range(len(corr)), corr)
    if logscale:
        ax.set_yscale('log')
    ax.set_title("SVD decomposition")
    ax.set_xlabel("Singular vector #")


SvdResponse.plot_norm = plot_norm
SvdResponse.plot_singular_values = plot_singular_values
ResponseMatrix.plot_contents = plot_contents
