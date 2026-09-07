"""ExactRectangularBendPass, now that rbendtune() works on numpy>=2.

Also checks the anchor question: Xsuite's RBend anchors misalignment on
length_straight (the chord) while the AT side defaults to Length/2 (the arc).
If those differ, MisalignAnchor should show it.
"""
import numpy as np
from xs_compare import compare, rect_bend, add_field_errors, apply_misalign

MIS = dict(dx=2e-3, dy=-1e-3, dz=1e-3, tilt=1e-2, pitch=-5e-3, yaw=8e-3)

print("=== baseline ===")
compare(rect_bend(k1=0.0, tune=False), label="k1=0 (no tune needed)", show_xs=True)
compare(rect_bend(), label="k1=0.3 + rbendtune   ", show_xs=True)

print("\n=== per-order field errors (rel 1e-4), alone vs misaligned ===")
print(f"{'error':<14}{'alone':>13}{'+misaligned':>14}")
print("-" * 41)
mo = compare(apply_misalign(rect_bend(), exact=True, **MIS), verbose=False)
print(f"{'(misalign)':<14}{'-':>13}{np.max(mo):>14.3e}")
for o in range(4):
    for kind, kw in ((f"B{o} normal", dict(skew=False)), (f"A{o} skew  ", dict(normal=False))):
        a = compare(add_field_errors(rect_bend(), rel=1e-4, orders=(o,), **kw), verbose=False)
        e = add_field_errors(rect_bend(), rel=1e-4, orders=(o,), **kw)
        b = compare(apply_misalign(e, exact=True, **MIS), verbose=False)
        av = np.max(a) if a is not None else np.nan
        bv = np.max(b) if b is not None else np.nan
        print(f"{kind:<14}{av:>13.3e}{bv:>14.3e}")

print("\n=== anchor sensitivity (dx only, which exposed the s-shift before) ===")
for anchor in (None, 2.0 / 2, None):
    e = rect_bend()
    apply_misalign(e, exact=True, dx=1e-3)
    if anchor is not None:
        e.MisalignAnchor = anchor
    compare(e, label=f"dx, anchor={anchor}", verbose=True)
