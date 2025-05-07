from __future__ import annotations
from ..lattice import Lattice
from ..latticetools import ResponseMatrix
from typing import Optional
import matplotlib.pyplot as plt
from matplotlib.axes import Axes


def plot_norm(resp: ResponseMatrix, ax: Optional[tuple[Axes, Axes]] = None) -> None:
    r"""Plot the norm of the lines and columns of the weighted response matrix

    For a stable solution, the norms should have the same order of magnitude.
    If not, the weights of observables and variables should be adjusted.

    Args:
        resp:           Response matrix object
        ax:             tuple of :py:class:`~.matplotlib.axes.Axes`. If given,
          plots will be drawn in these axes.
    """
    obs, var = resp.check_norm()
    if ax is None:
        fig, (ax1, ax2) = plt.subplots(nrows=2, layout="constrained")
    else:
        ax1, ax2 = ax[:2]
    ax1.bar(range(len(obs)), obs)
    ax1.set_title("Norm of weighted observables")
    ax1.set_xlabel("Observable #")
    ax2.bar(range(len(var)), var)
    ax2.set_title("Norm of weighted variables")
    ax2.set_xlabel("Variable #")


def plot_singular_values(
    resp: ResponseMatrix, ax: Axes = None, logscale: bool = True
) -> None:
    r"""Plot the singular values of a response matrix

    Args:
        resp:           Response matrix object
        logscale:       If :py:obj:`True`, use log scale
        ax:             If given, plots will be drawn in these axes.
    """
    if resp.singular_values is None:
        resp.solve()
    singvals = resp.singular_values
    if ax is None:
        fig, ax = plt.subplots()
    ax.bar(range(len(singvals)), singvals)
    if logscale:
        ax.set_yscale("log")
    ax.set_title("Singular values")


def plot_obs_analysis(
    resp: ResponseMatrix, lattice: Lattice, ax: Axes = None, logscale: bool = True
) -> None:
    """Plot the decomposition of an error vector on the basis of singular
    vectors

    Args:
        resp:           Response matrix object
        lattice:        Lattice description. The response matrix observables
          will be evaluated for this :py:class:`.Lattice` and the deviation
          from   target will be decomposed on the basis of singular vectors,
        logscale:       If :py:obj:`True`, use log scale
        ax:             If given, plots will be drawn in these axes.
    """
    if resp.singular_values is None:
        resp.solve()
    obs = resp.observables
    # noinspection PyProtectedMember
    obs.evaluate(lattice, **resp._eval_args)
    corr = resp._uh @ obs.flat_deviations
    if ax is None:
        fig, ax = plt.subplots()
    ax.bar(range(len(corr)), corr)
    if logscale:
        ax.set_yscale("log")
    ax.set_title("SVD decomposition")
    ax.set_xlabel("Singular vector #")


def plot_var_analysis(
    resp: ResponseMatrix, lattice: Lattice, ax: Axes = None, logscale: bool = False
) -> None:
    """Plot the decomposition of a correction vector on the basis of singular
    vectors

    Args:
        resp:           Response matrix object
        lattice:        Lattice description. The variables will be evaluated
          for this :py:class:`.Lattice` and will be decomposed on the basis
          of singular vectors,
        logscale:       If :py:obj:`True`, use log scale
        ax:             If given, plots will be drawn in these axes.
    """
    if resp.singular_values is None:
        resp.solve()
    var = resp.variables
    if ax is None:
        fig, ax = plt.subplots()
    corr = (resp._v * resp.singular_values).T @ var.get(lattice)
    ax.bar(range(len(corr)), corr)
    if logscale:
        ax.set_yscale("log")
    ax.set_title("SVD decomposition")
    ax.set_xlabel("Singular vector #")


ResponseMatrix.plot_norm = plot_norm
ResponseMatrix.plot_singular_values = plot_singular_values
ResponseMatrix.plot_obs_analysis = plot_obs_analysis
ResponseMatrix.plot_var_analysis = plot_var_analysis
