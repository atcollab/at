"""Per-order field error x misalignment, one order at a time.

This is the matrix that actually answers "does order N alone still agree once
the element is misaligned".  Every earlier per-order scan was misalignment-free,
and every earlier misalignment scan lumped all orders together.

Each order is run three ways: no misalignment, misaligned (exact), and the
difference tells whether misalignment introduces a cross-term beyond whatever
the order contributes on its own.
"""
import numpy as np
from xs_compare import (compare, sector_bend, multipole, rect_bend,
                        add_field_errors, apply_misalign)

MIS = dict(dx=2e-3, dy=-1e-3, dz=1e-3, tilt=1e-2, pitch=-5e-3, yaw=8e-3)
REL = 1e-4

BUILDERS = {
    "ExactSectorBend": lambda: sector_bend(),
    "ExactMultipole": lambda: multipole(),
    "ExactRectangularBend": lambda: rect_bend(),
}

for pmname, builder in BUILDERS.items():
    # A Dipole accepts MaxOrder up to len(Polynom)-1, and add_field_errors pads
    # to max(orders)+1, so orders 2-3 are testable on bends too (verified).
    maxord = 3
    print(f"\n===== {pmname}Pass =====", flush=True)
    print(f"{'error':<14}{'alone':>13}{'+misaligned':>14}{'misalign only':>15}",
          flush=True)
    print("-" * 56, flush=True)

    base = compare(builder(), label="", verbose=False)
    mis_only = compare(apply_misalign(builder(), exact=True, **MIS),
                       label="", verbose=False)
    print(f"{'(clean)':<14}{np.max(base):>13.3e}{'-':>14}"
          f"{np.max(mis_only):>15.3e}", flush=True)

    for o in range(0, maxord + 1):
        for kind, kw in ((f"B{o} normal", dict(skew=False)),
                         (f"A{o} skew  ", dict(normal=False))):
            alone = compare(add_field_errors(builder(), rel=REL, orders=(o,), **kw),
                            label="", verbose=False)
            e = add_field_errors(builder(), rel=REL, orders=(o,), **kw)
            both = compare(apply_misalign(e, exact=True, **MIS),
                           label="", verbose=False)
            a = np.max(alone) if alone is not None else np.nan
            b = np.max(both) if both is not None else np.nan
            print(f"{kind:<14}{a:>13.3e}{b:>14.3e}{'':>15}", flush=True)
