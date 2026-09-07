"""Root-cause the order-0 (B0/A0) sector-bend disagreement, finish the order
scan, and re-check the rectangular bend now that rbendtune() is called.

The B0 residual sits in y/py, which a normal dipole error cannot produce
through the body kick alone -- so the suspicion is the edge/fringe, not the
K0h term. Toggle them independently.
"""
from xs_compare import compare, sector_bend, rect_bend, add_field_errors

def sb(**kw):
    return sector_bend(**kw)

print("=== B0 (rel 1e-4): isolate body vs face vs fringes ===")
for label, kw in [
    ("face=.5 fb=1 fq=1 (default)", dict(face=0.5, fringe_bend=1, fringe_quad=1)),
    ("face=.5 fb=0 fq=1          ", dict(face=0.5, fringe_bend=0, fringe_quad=1)),
    ("face=.5 fb=1 fq=0          ", dict(face=0.5, fringe_bend=1, fringe_quad=0)),
    ("face=.5 fb=0 fq=0          ", dict(face=0.5, fringe_bend=0, fringe_quad=0)),
    ("face=0  fb=1 fq=1          ", dict(face=0.0, fringe_bend=1, fringe_quad=1)),
    ("face=0  fb=0 fq=0          ", dict(face=0.0, fringe_bend=0, fringe_quad=0)),
]:
    e = add_field_errors(sb(**kw), rel=1e-4, orders=(0,), skew=False)
    compare(e, label=label, verbose=True)

print("\n=== same toggles, CLEAN (no B0) -- baseline for each config ===")
for label, kw in [
    ("face=.5 fb=0 fq=0          ", dict(face=0.5, fringe_bend=0, fringe_quad=0)),
    ("face=0  fb=0 fq=0          ", dict(face=0.0, fringe_bend=0, fringe_quad=0)),
]:
    compare(sb(**kw), label=label, verbose=True)

print("\n=== B0 magnitude scan (default config) ===")
for rel in (1e-5, 1e-4, 1e-3):
    e = add_field_errors(sb(), rel=rel, orders=(0,), skew=False)
    compare(e, label=f"B0 rel={rel:.0e}   ", verbose=True)

print("\n=== remaining orders on the sector bend (rel 1e-4) ===")
for o in (2, 3):
    for kind, kw in ((f"normal B{o}", dict(skew=False)), (f"skew   A{o}", dict(normal=False))):
        e = add_field_errors(sb(), rel=1e-4, orders=(o,), **kw)
        compare(e, label=f"{kind}      ", verbose=True)

print("\n=== ExactRectangularBendPass, now with rbendtune() ===")
compare(rect_bend(), label="rect clean     ", verbose=True)
compare(rect_bend(k1=0.0), label="rect k1=0      ", verbose=True)
compare(add_field_errors(rect_bend(), rel=1e-4, orders=(0, 1)),
        label="rect + fields  ", verbose=True)
