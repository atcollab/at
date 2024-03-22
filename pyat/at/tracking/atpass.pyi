"""Stub file for the 'atpass' extension"""

import numpy as np
from typing import List, Optional
from at.lattice import Element, Particle

def atpass(line: List[Element], r_in: np.ndarray, nturns: int,
           refpts: np.ndarray,
           turn: Optional[int] = None,
           energy: Optional[float] = None,
           particle: Optional[Particle] = None,
           keep_counter: bool = False,
           reuse: bool = False,
           omp_num_thread: int = 0,
           losses: bool = False,
           bunch_spos = None, bunch_current = None): ...

def elempass(element: Element, r_in,
             energy: Optional[float] = None,
             particle: Optional[Particle] = None,
             ): ...

def reset_rng(*, rank: int = 0, seed: Optional[int] = None) -> None: ...
def common_rng() -> float: ...
def thread_rng() -> float: ...
