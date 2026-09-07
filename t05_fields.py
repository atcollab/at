"""Isolate the sector-bend field-error disagreement, one order at a time,
and check the rectangular-bend baseline.

MaxOrder is now raised before the polynomials are assigned, so orders 2-3
really are applied (t03 silently truncated them).
"""
import numpy as np
from xs_compare import compare, sector_bend, rect_bend, multipole, add_field_errors

print("=== ExactSectorBendPass: one order at a time, rel = 1e-4 ===")
compare(sector_bend(), label="no errors      ", verbose=True)
for o in (0, 1, 2, 3):
    for kind, kw in (("normal B%d" % o, dict(skew=False)),
                     ("skew   A%d" % o, dict(normal=False))):
        e = add_field_errors(sector_bend(), rel=1e-4, orders=(o,), **kw)
        compare(e, label=f"{kind}      ", verbose=True)

print("\n=== ExactSectorBendPass: B0 magnitude scan ===")
for rel in (1e-5, 1e-4, 1e-3, 1e-2):
    e = add_field_errors(sector_bend(), rel=rel, orders=(0,), skew=False)
    compare(e, label=f"B0 rel={rel:.0e}", verbose=True)

print("\n=== ExactSectorBendPass: does the face angle matter? (B0 rel=1e-4) ===")
for face in (0.0, 0.5):
    e = add_field_errors(sector_bend(face=face), rel=1e-4, orders=(0,), skew=False)
    compare(e, label=f"face={face}    ", verbose=True)

print("\n=== ExactRectangularBendPass baseline (no errors at all) ===")
compare(rect_bend(), label="rect clean     ", verbose=True, show_xs=True)
compare(rect_bend(k1=0.0), label="rect k1=0      ", verbose=True, show_xs=True)
compare(rect_bend(k1=0.0, angle=0.0), label="rect angle=0   ", verbose=True, show_xs=True)
