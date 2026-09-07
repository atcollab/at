"""Root-cause the A0 (skew dipole) disagreement the same way B0 was done.

For B0 the residual vanished at face=0 with fringes off, proving the body kick
was fine and the edge was at fault.  A0 has never been through this isolation:
its origin is currently only a hypothesis (that Xsuite's Bend never computes a
skew-dipole fringe at all).
"""
from xs_compare import compare, sector_bend, add_field_errors


def a0(rel=1e-4, **kw):
    return add_field_errors(sector_bend(**kw), rel=rel, orders=(0,), normal=False)


print("=== A0 (rel 1e-4): isolate body vs face vs fringes ===")
for label, kw in [
    ("face=.5 fb=1 fq=1 (default)", dict(face=0.5, fringe_bend=1, fringe_quad=1)),
    ("face=.5 fb=1 fq=0          ", dict(face=0.5, fringe_bend=1, fringe_quad=0)),
    ("face=.5 fb=0 fq=1          ", dict(face=0.5, fringe_bend=0, fringe_quad=1)),
    ("face=0  fb=1 fq=1          ", dict(face=0.0, fringe_bend=1, fringe_quad=1)),
    ("face=0  fb=1 fq=0          ", dict(face=0.0, fringe_bend=1, fringe_quad=0)),
    ("face=0  fb=0 fq=1          ", dict(face=0.0, fringe_bend=0, fringe_quad=1)),
    ("face=0  fb=0 fq=0          ", dict(face=0.0, fringe_bend=0, fringe_quad=0)),
]:
    compare(a0(**kw), label=label, verbose=True)

print("\n=== A0 magnitude scan (default config) ===")
for rel in (1e-5, 1e-4, 1e-3):
    compare(a0(rel=rel), label=f"A0 rel={rel:.0e}   ", verbose=True)

print("\n=== control: same configs with NO A0, to get each baseline ===")
for label, kw in [
    ("face=0  fb=1 fq=1          ", dict(face=0.0, fringe_bend=1, fringe_quad=1)),
    ("face=0  fb=0 fq=1          ", dict(face=0.0, fringe_bend=0, fringe_quad=1)),
    ("face=0  fb=1 fq=0          ", dict(face=0.0, fringe_bend=1, fringe_quad=0)),
]:
    compare(sector_bend(**kw), label=label, verbose=True)
