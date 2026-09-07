"""The main matrix: misalignments and field errors, separately and combined,
for ExactMultipolePass / ExactSectorBendPass / ExactRectangularBendPass.

For every case the element is tracked twice: once with AT's linearised R1/T1
misalignment (the current behaviour) and once with the exact transform.  Both
are compared against the Xsuite conversion of the same AT lattice.
"""
import sys
import numpy as np
from xs_compare import compare, sector_bend, multipole, rect_bend, apply_misalign

MIS_SMALL = dict(dx=1e-3, dy=1e-3, dz=1e-3, tilt=1e-3, pitch=1e-3, yaw=1e-3)
MIS_LARGE = dict(dx=2e-3, dy=-1e-3, dz=1e-3, tilt=1e-2, pitch=-5e-3, yaw=8e-3)


def add_field_errors(elem, rel=1e-4, orders=(0, 1, 2, 3)):
    """Add skew+normal errors of relative size `rel` to the given orders."""
    scale = max(abs(elem.PolynomB[1]), 1.0) if len(elem.PolynomB) > 1 else 1.0
    a = np.array(elem.PolynomA, dtype=float)
    b = np.array(elem.PolynomB, dtype=float)
    if len(a) < 4:
        a = np.pad(a, (0, 4 - len(a)))
        b = np.pad(b, (0, 4 - len(b)))
    for n in orders:
        a[n] += rel * scale
        b[n] += rel * scale
    elem.PolynomA = a
    elem.PolynomB = b
    elem.MaxOrder = max(elem.MaxOrder, max(orders))
    return elem


BUILDERS = {
    "ExactMultipole": multipole,
    "ExactSectorBend": sector_bend,
    "ExactRectangularBend": rect_bend,
}

CASES = [
    ("clean",                    None,      None),
    ("field errors 1e-4",        None,      dict(rel=1e-4)),
    ("field errors 1e-2",        None,      dict(rel=1e-2)),
    ("misalign small",           MIS_SMALL, None),
    ("misalign large",           MIS_LARGE, None),
    ("misalign small + fields",  MIS_SMALL, dict(rel=1e-4)),
    ("misalign large + fields",  MIS_LARGE, dict(rel=1e-4)),
    ("misalign large + fields2", MIS_LARGE, dict(rel=1e-2)),
]

for pmname, builder in BUILDERS.items():
    print(f"\n===== {pmname}Pass =====", flush=True)
    print(f"{'case':<26}{'linear (R1/T1)':>17}{'exact (new)':>15}{'gain':>9}",
          flush=True)
    print("-" * 67, flush=True)
    for label, mis, fld in CASES:
        res = {}
        for mode in ("linear", "exact"):
            e = builder()
            if fld:
                add_field_errors(e, **fld)
            if mis:
                apply_misalign(e, exact=(mode == "exact"), **mis)
            res[mode] = compare(e, label=label, verbose=False)
        lin, exa = res["linear"], res["exact"]
        if lin is None or exa is None:
            print(f"{label:<26}{'LOST':>17}{'LOST':>15}", flush=True)
            continue
        l, x = np.max(lin), np.max(exa)
        gain = f"{l / x:>8.0f}x" if x > 0 and mis else "       -"
        print(f"{label:<26}{l:>17.3e}{x:>15.3e}{gain}", flush=True)
        sys.stdout.flush()
