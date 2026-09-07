"""Why does a pure dx regress on a curved element?

The residual is exactly dx*2*sin(angle/2), so this prints the per-coordinate
breakdown and scans the bend angle to identify which term carries it.
"""
import numpy as np
from xs_compare import compare, sector_bend, multipole, apply_misalign, NAMES

print("=== straight control: ExactMultipolePass, dx = 1 mm ===")
compare(apply_misalign(multipole(), exact=False, dx=1e-3), label="linear")
compare(apply_misalign(multipole(), exact=True, dx=1e-3), label="exact ")

for angle in (1.0, 0.4, 0.1):
    print(f"\n=== ExactSectorBendPass, angle = {angle} rad, dx = 1 mm ===")
    print(f"    predicted dx*2*sin(angle/2) = {1e-3 * 2 * np.sin(angle / 2):.4e}")
    compare(apply_misalign(sector_bend(angle=angle), exact=False, dx=1e-3),
            label="linear")
    compare(apply_misalign(sector_bend(angle=angle), exact=True, dx=1e-3),
            label="exact ")

print("\n=== dz on a bend (works) for contrast, angle = 1 rad ===")
compare(apply_misalign(sector_bend(), exact=True, dz=1e-3), label="exact dz")
