import numpy as np

from at.lattice import Element, Particle

_defref = np.array([], dtype=np.uint32)

def line_store(
    line: list[Element],
    energy: float | None = None,
    particle: Particle | None = None,
    bunch_spos=None,
): ...

def line_pass(
    line: list[Element],
    r_in: np.ndarray,
    nturns: int,
    refpts: np.ndarray = ...,
    turn: int | None = None,
    keep_counter: bool = False,
    reuse: bool = False,
    omp_num_thread: int = 0,
    losses: bool = False,
    start_elem: int | None = None,
    end_elem: int | None = None,
): ...

def elem_pass(
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
