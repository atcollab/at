from __future__ import annotations
from ..latticetools import SvdResponse
from typing import Optional
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


SvdResponse.plot_norm = plot_norm
SvdResponse.plot_singular_values = plot_singular_values
