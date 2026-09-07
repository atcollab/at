"""Misalignment DOF scan: AT's linear R1/T1 vs the new exact transform,
both measured against Xsuite."""
import numpy as np
from xs_compare import compare, sector_bend, apply_misalign

CASES = [
    ("dx   = 1 mm",      dict(dx=1e-3)),
    ("dy   = 1 mm",      dict(dy=1e-3)),
    ("dz   = 1 mm",      dict(dz=1e-3)),
    ("tilt = 1 mrad",    dict(tilt=1e-3)),
    ("pitch= 1 mrad",    dict(pitch=1e-3)),
    ("yaw  = 1 mrad",    dict(yaw=1e-3)),
    ("tilt = 10 mrad",   dict(tilt=1e-2)),
    ("pitch= 10 mrad",   dict(pitch=1e-2)),
    ("yaw  = 10 mrad",   dict(yaw=1e-2)),
    ("all six",          dict(dx=1e-3, dy=1e-3, dz=1e-3,
                              tilt=1e-3, pitch=1e-3, yaw=1e-3)),
    ("all six, large",   dict(dx=2e-3, dy=-1e-3, dz=1e-3,
                              tilt=1e-2, pitch=-5e-3, yaw=8e-3)),
]

print(f"{'case':<18}{'linear (R1/T1)':>18}{'exact (new)':>18}{'gain':>10}")
print("-" * 64)
for label, kw in CASES:
    lin = compare(apply_misalign(sector_bend(), exact=False, **kw),
                  label=label, verbose=False)
    exa = compare(apply_misalign(sector_bend(), exact=True, **kw),
                  label=label, verbose=False)
    if lin is None or exa is None:
        print(f"{label:<18}{'LOST':>18}{'LOST':>18}")
        continue
    l, e = np.max(lin), np.max(exa)
    gain = l / e if e > 0 else np.inf
    print(f"{label:<18}{l:>18.3e}{e:>18.3e}{gain:>9.0f}x")
