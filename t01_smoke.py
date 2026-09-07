"""Smoke test: does the harness reproduce agreement for the clean element?"""
from xs_compare import compare, sector_bend, apply_misalign

print("=== ExactSectorBendPass, no errors ===")
compare(sector_bend(), label="clean", show_xs=True)

print()
print("=== dy = 1 mm, LINEAR (R1/T1, AT default) ===")
compare(apply_misalign(sector_bend(), exact=False, dy=1e-3),
        label="linear dy", show_xs=True)

print()
print("=== dy = 1 mm, EXACT (new) ===")
compare(apply_misalign(sector_bend(), exact=True, dy=1e-3),
        label="exact dy", show_xs=True)
