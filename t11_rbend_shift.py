"""Is AT's X0ref recoverable as Xsuite's rbend_shift?

xsuite.py documents X0ref as ignored, noting Xsuite has a "similar attribute
rbend_shift but there is no analytical conversion possible between both".
Scan rbend_shift on a single built line (a scalar change needs no kernel
rebuild) and see whether any value restores agreement, and if so where it sits
relative to X0ref.
"""
import numpy as np
import at
from xs_compare import (build_line, make_particles, track_at, track_xs,
                        rect_bend, ENERGY, NAMES)

e = rect_bend(k1=0.3)
print(f"AT X0ref = {e.X0ref:+.6e}   RefDZ = {e.RefDZ:+.6e}")

lat = at.Lattice([e.copy()], energy=ENERGY, periodicity=1)
line = build_line(lat)
b = line[0]
print("xsuite element:", type(b).__name__,
      "has rbend_shift:", hasattr(b, "rbend_shift"),
      "| current:", getattr(b, "rbend_shift", None))

parts = make_particles(7, 5e-3)
at_out = track_at(lat, parts)


def resid(shift):
    try:
        line[0].rbend_shift = shift
    except Exception as exc:                      # noqa: BLE001
        return None, f"cannot set: {exc}"
    xs = track_xs(line, parts, ENERGY)
    d = at_out - xs
    fin = np.isfinite(d).all(axis=0)
    if not fin.any():
        return None, "all lost"
    return np.nanmax(np.abs(d[:, fin]), axis=1), None


base, err = resid(0.0)
if err:
    print("rbend_shift unusable:", err)
else:
    print(f"\nrbend_shift = 0        -> {np.max(base):.4e}")
    x0 = e.X0ref
    for name, val in [("+X0ref", x0), ("-X0ref", -x0),
                      ("+X0ref/2", x0 / 2), ("-X0ref/2", -x0 / 2)]:
        r, err = resid(val)
        if r is not None:
            print(f"rbend_shift = {name:9s} ({val:+.5e}) -> {np.max(r):.4e}")

    print("\nfine scan around the best region:")
    best, bestv = None, None
    for val in np.linspace(-3 * abs(x0), 3 * abs(x0), 25):
        r, err = resid(val)
        if r is None:
            continue
        m = np.max(r)
        if best is None or m < best:
            best, bestv = m, val
    print(f"  minimum {best:.4e} at rbend_shift = {bestv:+.6e}")
    print(f"  X0ref = {x0:+.6e}   ratio = {bestv / x0 if x0 else float('nan'):.4f}")
