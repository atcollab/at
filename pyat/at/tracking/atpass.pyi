"""Stub file for the 'atpass' extension."""

from __future__ import annotations

import numpy as np

from at.lattice import Element, Particle

_defref = np.array([], dtype=np.uint32)

def atpass(
    line: list[Element],
    r_in: np.ndarray,
    nturns: int,
    refpts: np.ndarray = _defref,
    turn: int | None = None,
    energy: float | None = None,
    particle: Particle | None = None,
    keep_counter: bool = False,
    reuse: bool = False,
    omp_num_thread: int = 0,
    losses: bool = False,
    bunch_spos=None,
    bunch_current=None,
): ...
def elempass(
    element: Element,
    r_in: np.ndarray,
    energy: float | None = None,
    particle: Particle | None = None,
): ...
def diffusion_matrix(
    element: Element,
    r_in: np.ndarray,
    energy: float | None = None,
    particle: Particle | None = None,
): ...
def reset_rng(*, rank: int = 0, seed: int | None = None) -> None: ...
def common_rng() -> float: ...
def thread_rng() -> float: ...
