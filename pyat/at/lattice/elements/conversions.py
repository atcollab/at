"""Conversion functions for Element attributes"""

from __future__ import annotations

__all__ = ["_array", "_array66", "_int", "_float"]

import numpy as np


def _array(value, shape=(-1,), dtype=np.float64):
    # Ensure proper ordering(F) and alignment(A) for "C" access in integrators
    return np.require(value, dtype=dtype, requirements=["F", "A"]).reshape(
        shape, order="F"
    )


def _array66(value):
    return _array(value, shape=(6, 6))


def _float(value) -> float:
    return float(value)


def _int(value, vmin: int | None = None, vmax: int | None = None) -> int:
    intv = int(value)
    if vmin is not None and intv < vmin:
        raise ValueError(f"Value must be greater of equal to {vmin}")
    if vmax is not None and intv > vmax:
        raise ValueError(f"Value must be smaller of equal to {vmax}")
    return intv
