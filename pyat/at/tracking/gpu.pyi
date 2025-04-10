"""Stub file for the 'gpu' extension."""

from __future__ import annotations

import numpy as np

from at.lattice import Element, Particle

def gpuinfo(): ...
def gpupass(
    line: list[Element],
    r_in: np.ndarray,
    nturns: int,
    refpts: np.ndarray = ...,
    turn: int | None = None,
    energy: float | None = None,
    particle: Particle | None = None,
    keep_counter: bool = False,
    reuse: bool = False,
    losses: bool = False,
    bunch_spos=None,
    bunch_current=None,
    gpu_pool: list[int] | None = None,
    tracking_starts=None,
    integrator=4,
    verbose: bool = False
): ...
